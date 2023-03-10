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

Class
    Foam::solvers::incompressibleFluid

Description
    Solver module for steady or transient turbulent flow of incompressible
    isothermal fluids with optional mesh motion and change.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, constraining or limiting
    the solution.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

SourceFiles
    incompressibleFluid.C

See also
    Foam::solvers::fluidSolver
    Foam::solvers::isothermalFluid

\*---------------------------------------------------------------------------*/

#ifndef incompressibleFluid_H
#define incompressibleFluid_H

#include "fluidSolver.H"
#include "viscosityModel.H"
#include "incompressibleMomentumTransportModels.H"
#include "pressureReference.H"
#include "IOMRFZoneList.H"
#include "kEpsilon.H"
#include "constantViscosityModel.H"
#include <cassert>
#include <unordered_map>

#define PRINT_EXPR(x) (void)(Info << #x " = " << (x) << endl)
#define PRINT(x) (void)(Info << #x " begin:\n" << (x) << "\n" #x " end"<< endl)
Foam::string get_typename(const std::type_info &t_info);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

struct momentumTransport_t
    : Foam::RASModels::kEpsilon<Foam::incompressibleMomentumTransportModel>
{
    using _MyBase = Foam::RASModels::kEpsilon<Foam::incompressibleMomentumTransportModel>;
    using _MyBase::kEpsilon;
    using _MyBase::operator=;
    using BasicMomentumTransportModel = Foam::incompressibleMomentumTransportModel;
    using _MyBase::Cmu_; // 改为 public
    void correct(); // 搬过来，方便修改
    void correctNut(); // 同上


};

struct fvModels_t
    : Foam::fvModels
{
    using _MyBase = Foam::fvModels;
    using _MyBase::fvModels;
    using _MyBase::source; // 公开所有source
};

static_assert(std::has_virtual_destructor<momentumTransport_t>::value);

/*---------------------------------------------------------------------------*\
                          Class incompressibleFluid Declaration
\*---------------------------------------------------------------------------*/

class incompressibleFluid
:
    public fluidSolver
{
    using _MyBase = fluidSolver;
protected:

    // Pressure

        //- Pressure field
        volScalarField p;

        //- Pressure reference
        Foam::pressureReference pressureReference;


    // Kinematic properties

        //- Velocity field
        volVectorField U;

        //- Mass-flux field
        surfaceScalarField phi;


    // Momentum transport

        //- Kinematic viscosity model
        // autoPtr<viscosityModel> viscosity;
        using viscosity_t = Foam::viscosityModels::constant;
        autoPtr<viscosity_t> viscosity;

        //- Pointer to the momentum transport model
        // autoPtr<incompressible::momentumTransportModel> momentumTransport;
        autoPtr<momentumTransport_t> momentumTransport;


    // Continuity properties

        //- Current cumulative continuity error
        scalar cumulativeContErr;


    // Optional models

        //- MRF zone list
        IOMRFZoneList MRF;


    // Cached temporary fields

        //- Pointer to the surface momentum field
        //  used to recreate the flux after mesh-change
        autoPtr<surfaceVectorField> Uf;

        //- Optional LTS reciprocal time-step field
        tmp<volScalarField> trDeltaT;

        //- Cached momentum matrix
        //  shared between the momentum predictor and pressure corrector
        tmp<fvVectorMatrix> tUEqn;


private:

    // Private Member Functions

        //- Set rDeltaT for LTS
        void setRDeltaT();

        //- Correct the cached Courant numbers
        void correctCoNum();

        //- Calculate and print the continuity errors
        void continuityErrors();

        //- Construct the pressure equation
        //  and correct the pressure and velocity
        void correctPressure();


public:

    //- Runtime type information
    TypeName("incompressibleFluid");


    // Constructors

        //- Construct from region mesh
        incompressibleFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        incompressibleFluid(const incompressibleFluid&) = delete;


    //- Destructor
    virtual ~incompressibleFluid();


    // Member Functions

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual bool moveMesh();

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum and thermophysical transport modelling
        virtual void postCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();

        fvModels_t &fvModels() const
        {
            return static_cast<fvModels_t &>(_MyBase::fvModels());
        }

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const incompressibleFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
