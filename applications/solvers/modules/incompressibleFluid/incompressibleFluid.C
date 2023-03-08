/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "incompressibleFluid.H"
#include "localEulerDdtScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleFluid, 0);
    addToRunTimeSelectionTable(solver, incompressibleFluid, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::incompressibleFluid::correctCoNum()
{
    fluidSolver::correctCoNum(phi);
}


void Foam::solvers::incompressibleFluid::continuityErrors()
{
    fluidSolver::continuityErrors(phi);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleFluid::incompressibleFluid(fvMesh& mesh)
:
    fluidSolver(mesh),

    p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    pressureReference(p, pimple.dict()),

    U
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    phi
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U) & mesh.Sf()
    ),

    viscosity(viscosityModel::New(mesh)),

    momentumTransport
    (
        incompressible::momentumTransportModel::New
        (
            U,
            phi,
            viscosity
        )
    ),

    MRF(mesh)
{
    mesh.schemes().setFluxRequired(p.name());

    momentumTransport->validate();

    if (mesh.dynamic())
    {
        Info<< "Constructing face momentum Uf" << endl;

        Uf = new surfaceVectorField
        (
            IOobject
            (
                "Uf",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::interpolate(U)
        );
    }

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        Info<< "Using LTS" << endl;

        trDeltaT = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    runTime.name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleFluid::~incompressibleFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleFluid::preSolve()
{
    // Read the controls
    read();

    fvModels().preUpdateMesh();

    if (transient())
    {
        correctCoNum();
    }
    else if (LTS)
    {
        setRDeltaT();
    }

    // Update the mesh for topology change, mesh to mesh mapping
    mesh.update();
}


void Foam::solvers::incompressibleFluid::prePredictor()
{
    fvModels().correct();

    if (pimple.predictTransport())
    {
        momentumTransport->predict();
    }
}


void Foam::solvers::incompressibleFluid::thermophysicalPredictor()
{}


void Foam::solvers::incompressibleFluid::pressureCorrector()
{
    while (pimple.correct())
    {
        correctPressure();
    }

    tUEqn.clear();
}


void Foam::solvers::incompressibleFluid::postCorrector()
{
    if (pimple.correctTransport())
    {
        viscosity->correct();
        momentumTransport->correct();
    }
}

#if 0
void Foam::solvers::momentumTransport_t::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*G*epsilon_()/k_()
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      + epsilonSource()
      + fvModels.source(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

    correctNut();
}
#endif

void Foam::solvers::incompressibleFluid::postSolve()
{}


// ************************************************************************* //

#if 0
#include <cxxabi.h>

Foam::string get_typename(const std::type_info &t_info)
{
    char *name = abi::__cxa_demangle(t_info.name(), 0, 0, 0);
    Foam::string res(name);
    free(name); // 省事
    return res;
}


record_t record; // record 的定义


#include "error.H"
#include "OStringStream.H"
#include "OSspecific.H"
#include "IFstream.H"

#include <inttypes.h>
#include <cxxabi.h>
#include <execinfo.h>
#include <dlfcn.h>

namespace Foam {
word demangleSymbol(const char* sn);
fileName absolutePath(const char* fn);
void printSourceFileAndLine
(
    Ostream& os,
    const fileName& filename,
    Dl_info *info,
    void *addr
);

void Foam::Foam_error::printStack(Ostream& os)
{
    // Get raw stack symbols
    const size_t CALLSTACK_SIZE = 128;

    void *callstack[CALLSTACK_SIZE];
    size_t size = backtrace(callstack, CALLSTACK_SIZE);

    Dl_info *info = new Dl_info;

    fileName fname = "???";
    word address;

    for(size_t i=0; i<size; i++)
    {
        int st = dladdr(callstack[i], info);

        os << '#' << label(i) << "  ";
        if (st != 0 && info->dli_fname != nullptr && info->dli_fname[0] != '\0')
        {
            fname = absolutePath(info->dli_fname);

            os <<
            (
                (info->dli_sname != nullptr)
              ? demangleSymbol(info->dli_sname)
              : "?"
            );
        }
        else
        {
            os << "?";
        }

        printSourceFileAndLine(os, fname, info, callstack[i]);
        os << nl;
    }

    delete info;
}

} // namespace Foam
#endif
