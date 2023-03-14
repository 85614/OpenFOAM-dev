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
#include "epsilonWallFunctionFvPatchScalarField.H"

#include "record.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleFluid, 0);
    addToRunTimeSelectionTable(solver, incompressibleFluid, fvMesh);
}

struct Foam_error
    : public Foam::error
{
    static void printStack(Ostream &os, size_t i);
};

} // namespace Foam 

void printCaller() 
{
    Foam::Foam_error::printStack(Foam::Info, 1);
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

struct {
    const Foam::Time *ptr = nullptr;
} __use_assert;

bool use_assert() 
{
    return __use_assert.ptr->timeIndex() <= 5;
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

    viscosity(dynamic_cast<viscosity_t*>(viscosityModel::New(mesh).ptr())),

    momentumTransport
    (
        static_cast<momentumTransport_t*>(
            dynamic_cast<momentumTransport_t::_MyBase*>(
        incompressible::momentumTransportModel::New
        (
            U,
            phi,
            viscosity
        )
        .ptr()
        ))
    ),

    MRF(mesh)
{
    record.insert_clone(PAIR(U));
    record.insert_clone(PAIR(p));
    record.insert_clone(PAIR(phi));
    record.insert_clone("k", (momentumTransport->k()()), _WHERE);
    record.insert_clone("nut", (momentumTransport->nut()()), _WHERE);
    record.insert_clone("nu", (momentumTransport->nu()()), _WHERE);
    record.insert_clone("epsilon", (momentumTransport->epsilon()()), _WHERE);
    record.insert_clone("sigma", (momentumTransport->sigma()()), _WHERE);
    assert(typeid(momentumTransport_t::_MyBase) == typeid(momentumTransport()));
    assert(typeid(viscosity_t) == typeid(viscosity()));

    __use_assert.ptr = &runTime;

    if (0)
    {
        // in foamRun 
        auto &&solver = *this;
        solver.preSolve();
        {
            solver.moveMesh(); // do nothing
            solver.prePredictor(); // do nothing
            solver.momentumPredictor();
            solver.thermophysicalPredictor(); // do nothing
            solver.pressureCorrector(); // 压力修正 内循环
            solver.postCorrector(); // 
        }

        solver.postSolve();
    }
    mesh.schemes().setFluxRequired(p.name());
    assert(typeid(Foam::fvSchemes) == typeid(mesh.schemes()));
    momentumTransport->validate();
    {
        // momentumTransport->correctNut();
        auto &&Cmu_ = momentumTransport->Cmu_;
        auto &&k_ = momentumTransport->k();
        auto &&epsilon_ = momentumTransport->epsilon();
        auto &&nut = momentumTransport->nut();
        assert(nut().internalField() == Cmu_*sqr(k_)/epsilon_);
        record.assign("nut", (momentumTransport->nut()()), _WHERE);
    }

    assert(!mesh.dynamic());
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
    assert(transient());

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
    return; // 什么也没做

    fvModels().correct();

    {
        PtrListDictionary<fvModel>& modelList(fvModels());
        assert(modelList.size() == 0);
    }
    assert(pimple.predictTransport());
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
    assert(pimple.correctTransport());
    if (pimple.correctTransport())
    {
        viscosity->correct(); // do nothing
        momentumTransport->momentumTransport_t::correct();
    }
}


void Foam::solvers::momentumTransport_t::correct()
{
    if (use_assert()) {
        assert(epsilon_.boundaryFieldRef()[0].manipulatedMatrix() ==  record.get<volScalarField>("epsilon").boundaryFieldRef()[0].manipulatedMatrix());
    }
    // 和 _MyBase::correct 一样
    assert(this->turbulence_);
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
    {
        // BasicMomentumTransportModel::correct();
        // do nothing
    }

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
    {
        assert(tgradU().v() == tgradU().internalField());
    }
    tgradU.clear();
    if (use_assert()) 
    {
        assert(record.equal("epsilon", epsilon_));
    }
    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
    if (use_assert())
    {
        // internalField 也改了
        assert(!(record.get<volScalarField>("epsilon").internalField() == epsilon_.internalField()));
        assert(!record._equal(record.get<volScalarField>("epsilon").boundaryField(), epsilon_.boundaryField()));

        record.assign("epsilon", epsilon_, _WHERE);
        assert(record.equal("epsilon", epsilon_));
    }
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
    if (use_assert()) 
    {
        static_assert(std::is_same<Foam::geometricOneField, std::decay<decltype(alpha)>::type>::value);
        static_assert(std::is_same<Foam::geometricOneField, std::decay<decltype(rho)>::type>::value);
        
        static_assert(std::is_same<Foam::geometricOneField::Internal, std::decay<decltype(alpha())>::type>::value);
        static_assert(std::is_same<Foam::geometricOneField::Internal, std::decay<decltype(rho())>::type>::value);
        
        assert(record.equal("phi", alphaRhoPhi));
        assert(record.equal(PAIR(U)));
        assert(record.equal(PAIR(nut)));
        // assert(record.equal("epsilon", epsilon_, _WHERE));
        assert(record._equal(DepsilonEff(), this->nut_/sigmaEps_ + this->nu()));
        static const auto sigmaEps_const = sigmaEps_;

        static const auto const_record = std::make_tuple(Cmu_.value(), C1_.value(), C2_.value(), C3_.value(), sigmak_.value(), sigmaEps_.value());

        assert(const_record == std::make_tuple(Cmu_.value(), C1_.value(), C2_.value(), C3_.value(), sigmak_.value(), sigmaEps_.value()));

        fvScalarMatrix zeroMatrix(epsilon_, epsEqn->dimensions());
        assert(record._equal(epsilonSource()(), zeroMatrix));
        assert(record._equal(fvModels.source(alpha, rho, epsilon_)(), zeroMatrix));


        auto &epsilon_ = record.get<volScalarField>("epsilon");
        tmp<fvScalarMatrix> epsEqn
        (
            fvm::ddt(alpha, rho, epsilon_)
        + fvm::div(alphaRhoPhi, epsilon_)
        - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
        ==
            C1_*alpha()*rho()*G*epsilon_()/k_()
        - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilon_)
        - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
        //   + epsilonSource()
        //   + fvModels.source(alpha, rho, epsilon_)
        );

        record.insert_clone_or_assign(PAIR(epsEqn()));
    }

    epsEqn.ref().relax();
    assert(0 == epsEqn->relaxationFactor());
    if (use_assert()) 
    {
        assert(record.equal(PAIR(epsEqn())));
    }
    fvConstraints.constrain(epsEqn.ref());
    if (use_assert()) 
    {
        assert(record.equal(PAIR(epsEqn())));
        assert(record.equal("epsilon", epsilon_));
        assert(epsilon_.boundaryFieldRef()[0].manipulatedMatrix() ==  record.get<volScalarField>("epsilon").boundaryFieldRef()[0].manipulatedMatrix());

    }
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());

    if (use_assert()) 
    {
        assert(record.equal("epsilon", epsilon_));
        assert(!record.equal(PAIR(epsEqn())));
        auto &epsEqn_record = record.get(PAIR(epsEqn()));

        {
            auto &epsilon_record = record.get<decltype(epsilon_)>("epsilon");
            assert(record._equal(epsilon_record, epsilon_));

            // 内部实现 epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
            auto &bFields = epsilon_record.boundaryFieldRef();
            forAll(bFields, patchi)
            {
                // PRINT_EXPR(get_typename(typeid(bFields[patchi]))); // Foam::epsilonWallFunctionFvPatchScalarField
                using patch_t = Foam::epsilonWallFunctionFvPatchScalarField;
                if (patchi + 1 < bFields.size())
                {
                    assert(typeid(bFields[patchi]) == typeid(patch_t));
                    auto &patch_ref = dynamic_cast<patch_t &>(bFields[patchi]);
                    patch_ref.manipulateMatrix(epsEqn_record);
                }
            }
        }
        // record.assign(PAIR(epsEqn()));
        assert(record.equal(PAIR(epsEqn())));
    }
    assert(epsilon_.boundaryFieldRef()[0].manipulatedMatrix());
    solve(epsEqn);
    assert(!epsilon_.boundaryFieldRef()[0].manipulatedMatrix());
    fvConstraints.constrain(epsilon_);
    if (use_assert()) {
        auto &epsEqn_record = record.get<fvScalarMatrix>("epsEqn()");
        solve(epsEqn_record);
        assert(record.equal("epsilon", epsilon_));
    }
    bound(epsilon_, this->epsilonMin_);
    record.assign("epsilon", epsilon_, _WHERE);

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
    
    if (use_assert()) {
        record.assign("k", k_, _WHERE);
        record.assign(PAIR(nut));
        assert(record.equal("epsilon", epsilon_));
        assert(epsilon_.boundaryFieldRef()[0].manipulatedMatrix() ==  record.get<volScalarField>("epsilon").boundaryFieldRef()[0].manipulatedMatrix());
    }
}

void Foam::solvers::momentumTransport_t::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}

void Foam::solvers::incompressibleFluid::postSolve()
{}


// ************************************************************************* //

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

void Foam::Foam_error::printStack(Ostream& os, size_t i)
{
    // Get raw stack symbols
    const size_t CALLSTACK_SIZE = 8;

    void *callstack[CALLSTACK_SIZE];
    size_t size = backtrace(callstack, CALLSTACK_SIZE);

    Dl_info *info = new Dl_info;

    fileName fname = "???";
    word address;

    // for(size_t i=0; i<size; i++)
    if (i < size)
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
