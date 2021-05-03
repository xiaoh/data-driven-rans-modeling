        sensitiveJ = fvc::grad(Ua);

        volTensorField laplacianTauR("laplacianTauR", fvc::laplacian(tauR));
        volTensorField deltaTauR
        (
            "deltaTauR",
            mesh.fieldRelaxationFactor("tauR")
           *(
                lambda*sensitiveJ
              + 2*regulationCoeff*laplacianTauR
            )
        );

        tauR += deltaTauR;
        volVectorField divTauR("divTauR", fvc::div(tauR));

        if (runTime.outputTime())
        {
            deltaTauR.write();
            divTauR.write();
        }

        volScalarField J("J", (U - UDNS) & (U - UDNS));

        scalar JTotal(0.0);
        forAll(J, cellI)
        {
            JTotal += J.internalField()[cellI]*mesh.V()[cellI];
        }

        Info<< "J is: " << JTotal << nl << endl;
