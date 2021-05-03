/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Application
    ajointShapeOptimizationFoam

Description
    Steady-state solver for incompressible, turbulent flow of non-Newtonian
    fluids with optimisation of duct shape by applying "blockage" in regions
    causing pressure loss as estimated using an adjoint formulation.

    References:
    \verbatim
        "Implementation of a continuous adjoint for topology optimization of
         ducted flows"
        C. Othmer,
        E. de Villiers,
        H.G. Weller
        AIAA-2007-3947
        http://pdf.aiaa.org/preview/CDReadyMCFD07_1379/PV2007_3947.pdf
    \endverbatim

    Note that this solver optimises for total pressure loss whereas the
    above paper describes the method for optimising power-loss.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "fvIOoptionList.H"

template<class Type>
void zeroCells
(
    GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells
)
{
    forAll(cells, i)
    {
        vf[cells[i]] = pTraits<Type>::zero;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "initAdjointContinuityErrs.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        laminarTransport.lookup("lambda") >> lambda;

        #include "updateTauR.C"

        //zeroCells(Tau, inletCells);
        //zeroCells(Tau, outletCells);

        // Pressure-velocity SIMPLE corrector
        int currentSimpleIter = 0;
        // Adjoint Pressure-velocity SIMPLE corrector
        while (currentSimpleIter < nSimpleLoop)
        {
            // Momentum predictor

            tmp<fvVectorMatrix> UEqn
            (
                fvm::div(phi, U)
              - fvm::laplacian(nu,U)
              + fvc::div(dev(tauR))
              ==
                fvOptions(U)
              + explicitGradient
            );

            UEqn().relax();

            fvOptions.constrain(UEqn());

            solve(UEqn() == -fvc::grad(p));

            fvOptions.correct(U);

            volScalarField rAU(1.0/UEqn().A());
            volVectorField HbyA("HbyA", U);
            HbyA = rAU*UEqn().H();
            UEqn.clear();
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::interpolate(HbyA) & mesh.Sf()
            );

            fvOptions.makeRelative(phiHbyA);

            adjustPhi(phiHbyA, U, p);

            // Non-orthogonal pressure corrector loop
            while (simple.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (simple.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            currentSimpleIter ++;
        }

        int currentAdjointIter = 0;
        forAll(Ua, cellI)
        {
            Ua[cellI] = 0*Ua[cellI];
            pa[cellI] = 0*pa[cellI];
        }

        // Adjoint Pressure-velocity SIMPLE corrector
        while (currentAdjointIter < nAdjointLoop)
        {
            // Adjoint Momentum predictor
            if (adjointType == "QiqiSimple")
            {
                // Qiqi simple form
                solve(-fvm::laplacian(nu, Ua) + 2*(U - UDNS)/unitTime);
            }
            else if (adjointType == "QiqiComplex")
            {
                // Qiqi complex form
                tmp<fvVectorMatrix> UaEqn
                (
                   (fvc::grad(U) & Ua)
                  -(U & fvc::grad(Ua))
                  - fvm::laplacian(nu, Ua)
                  - fvc::div(nu*dev(T(fvc::grad(Ua))))
                  // dimension needs test
                  + 2 * (U - UDNS) / unitTime      // derivative of J (objective function)
                );

                #include "UaEqn.H"
            }
            else if (adjointType == "Othmer")
            {
                // Othmer form
                tmp<fvVectorMatrix> UaEqn
                (
                   fvm::div(-phi,Ua)
                  -(fvc::grad(Ua) & U)
                  - fvm::laplacian(nu, Ua)
                  - fvc::div(nu*dev(T(fvc::grad(Ua))))
                  // dimension needs test
                  + 2*(U - UDNS)/unitTime      // derivative of J (objective function)
                );

                #include "UaEqn.H"
            }

            currentAdjointIter ++;
        }

        // turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
