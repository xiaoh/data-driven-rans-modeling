    // Create Exact Solutions

    forAll(mesh.C(), cellI)
    {
        scalar x = mesh.C()[cellI].component(0);
        scalar y = mesh.C()[cellI].component(1);
        scalar z = mesh.C()[cellI].component(2);
        U[cellI].x() = Foam::sin(x);
        U[cellI].y() = Foam::cos(x);
        U[cellI].z() = Foam::cos(2*x);
        p[cellI] = Foam::exp(y);
    }

    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& currPatch = mesh.boundary()[patchI];
        forAll(currPatch, faceI)
        {
            scalar x = mesh.C().boundaryField()[patchI][faceI].component(0);
            scalar y = mesh.C().boundaryField()[patchI][faceI].component(1);
            scalar z = mesh.C().boundaryField()[patchI][faceI].component(2);
            U.boundaryField()[patchI][faceI].x() = Foam::sin(x);
            U.boundaryField()[patchI][faceI].y() = Foam::cos(x);
            U.boundaryField()[patchI][faceI].z() = Foam::cos(2*x);
            p.boundaryField()[patchI][faceI] = Foam::exp(y);
        }
    }

    //Write fields out
    //if (runTime.timeIndex() == 0)
    //PressureGradient = (U & fvc::grad(p)) / (sqrt((fvc::grad(p) & fvc::grad(p)) * (U & U)) + mag(U & fvc::grad(p)));
