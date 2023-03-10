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
#include "fvmDiv.H"
#include "record.H"

bool use_assert();
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleFluid::momentumPredictor()
{
    tUEqn =
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + momentumTransport->divDevSigma(U)
     ==
        fvModels().source(U)
    );
    if (use_assert()) {
        // 确定没有在意料之外的地方被修改
        assert(record.equal(PAIR(U)));
        assert(record.equal(PAIR(phi)));
        assert(record.equal("nu", (momentumTransport->nu()())));
        assert(record.equal("nut", (momentumTransport->nut()())));

        assert(record._equal(momentumTransport->nuEff(), momentumTransport->nu() + momentumTransport->nut()));

        static_assert(std::is_same<Foam::geometricOneField, std::decay<decltype(momentumTransport->alpha())>::type>::value);
        static_assert(std::is_same<Foam::geometricOneField, std::decay<decltype(momentumTransport->rho())>::type>::value);

        auto mytUEqn = fvm::ddt(U) + fvm::div(phi, U);
        // momentumTransport->divDevSigma(U)
        mytUEqn = mytUEqn +
            (
                - fvc::div((momentumTransport->alpha()*momentumTransport->rho()*momentumTransport->nuEff())*dev2(T(fvc::grad(U))))
                - fvm::laplacian(momentumTransport->alpha()*momentumTransport->rho()*momentumTransport->nuEff(), U)
            ); 
        // MRF.DDt(U) 和 fvModels().source(U) 为 0
        
        assert(record._equal(tUEqn, mytUEqn));
        record.insert_clone_or_assign(PAIR(tUEqn()));
    }

    if (0) 
    {
        auto T_grad_U = T(fvc::grad(U));
        const Foam::tensor &T_grad_U_ref_0 = T_grad_U.ref().primitiveField()[0];
        Foam::dev2(T_grad_U_ref_0);
        {
            const Foam::tensor &t = T_grad_U_ref_0;
            using Cmpt = Foam::scalar;
            t - SphericalTensor<Cmpt>::twoThirdsI*tr(t);
        }
        volVectorField U_copy = U;
        record._equal(U, U);
    }
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();
    if (use_assert()) {
        // UEqn.relax(UEqn.relaxationFactor());
        assert(UEqn.relaxationFactor() == 0); // 没设置，默认是0
        // UEqn.relax(Foam::scalar(0)); // <= 0 直接 return 了 ...
        assert(record.equal(PAIR(tUEqn())));
    }

    fvConstraints().constrain(UEqn);
    
    if (use_assert()) {
        const PtrListDictionary<fvConstraint>& constraintList(fvConstraints());
        assert(constraintList.empty());
        // 什么也没做
        assert(record.equal(PAIR(tUEqn())));
    }

    assert(pimple.momentumPredictor());
    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
        record.assign(PAIR(U));

        fvConstraints().constrain(U);
        // assert(!fvConstraints().constrain(U));
    }
    if (use_assert()) {
        assert(record.equal(PAIR(p)));
        assert(record.equal(PAIR(U)));
    }
}


// ************************************************************************* //
