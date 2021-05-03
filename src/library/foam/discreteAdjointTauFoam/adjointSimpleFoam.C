/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "createMRF.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    #include "adjointSettings.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dco::ga1s<double>::global_tape = dco::ga1s<double>::tape_t::create();
    forAll(tauR,i){
        for(int j = 0; j < 5; j++){
            dco::ga1s<double>::global_tape->register_variable(tauR[i].component(j));
        }
    }

    for(int i = 0; i < tauR.size(); i++){
        if(mesh.V()[i] == 0)
            Info << "zero V in mesh!" << endl;
    }

    Info<< "\nStarting time loop\n" << endl;

    // run until end time reached / converged
    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << (scalar)runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " <<   (scalar)runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "Calculating Cost Function" << endl;
    // calc cost function
    scalar J = 0;
    forAll(U,cellI)
    {
        vector Udiff = U[cellI] - UDNS[cellI];
        J += (Udiff & Udiff)*mesh.V()[cellI];
    }
    // Set adjoint of J
    if(Pstream::master())
        dco::derivative(J)=1;

    Info<< "Interpreting tape" << endl;
    dco::ga1s<double>::global_tape->interpret_adjoint();

    // get adjoint sensitivitys, scale with cell volume, write to sens
    forAll(tauR,i){
        double tmp;
        for(int j = 0; j < 5; j++){
            sens[i].component(j) = dco::derivative(tauR[i].component(j))/mesh.V()[i];
        }
    }

    runTime.setTime(0,(label)0);
    sens.write();
    Info << "Sens sum: " << gSum(sens) << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
