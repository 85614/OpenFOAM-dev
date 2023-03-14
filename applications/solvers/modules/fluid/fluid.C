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

#include "fluid.H"
#include "addToRunTimeSelectionTable.H"
#include "tool.H"
#include "eddyDiffusivity.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(fluid, 0);
    addToRunTimeSelectionTable(solver, fluid, fvMesh);
}

template<class _Base>
struct with_base: public _Base 
{
    using _MyBase = _Base;
    using _Base::_Base;
};

#define MAKE_SPECIAL_T(var, base_t) struct var##_t: public with_base<base_t>

struct thermophysicalTransport_t
    :public with_base<Foam::turbulenceThermophysicalTransportModels::eddyDiffusivity<Foam::RASThermophysicalTransportModel<Foam::ThermophysicalTransportModel<Foam::compressibleMomentumTransportModel, Foam::fluidThermo>>>>
{

};
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::fluid::fluid(fvMesh& mesh)
:
    isothermalFluid(mesh),

    thermophysicalTransport
    (
        fluidThermoThermophysicalTransportModel::New
        (
            momentumTransport(),
            thermo
        )
    )
{
    thermo.validate(type(), "h", "e");
    // PRINT_EXPR(get_typename(typeid(thermophysicalTransport())));
    assert(typeid(thermophysicalTransport()) == typeid(thermophysicalTransport_t::_MyBase));
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::fluid::~fluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::fluid::prePredictor()
{
    isothermalFluid::prePredictor();

    if (pimple.predictTransport())
    {
        thermophysicalTransport->predict();
    }
}


void Foam::solvers::fluid::postCorrector()
{
    isothermalFluid::postCorrector();

    if (pimple.correctTransport())
    {
        thermophysicalTransport->correct();
    }
}

#include "tool.C"
// ************************************************************************* //
