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
        pGradAlongStream_Exact[cellI] = Foam::cos(x) * Foam::exp(y)
                        / (Foam::sqrt( Foam::exp(y) * Foam::exp(y) * ( 1 + Foam::sqr(Foam::cos(2*x))))
                        + Foam::mag(Foam::cos(x) * Foam::exp(y)));
        scalar UG_numerator = Foam::mag(Foam::sin(x)*Foam::cos(x)*Foam::sin(x)
                            - Foam::sin(x)*Foam::sin(x)*Foam::cos(x)
                            - Foam::sin(x)*2*Foam::sin(2*x)*Foam::cos(2*x));
        scalar UG_denominator = (Foam::sqr(Foam::sin(x)) 
                              + Foam::sqr(Foam::cos(x))
                              + Foam::sqr(Foam::cos(2*x)))
                              * Foam::sqr(Foam::sin(x)*Foam::cos(x)
                                         -Foam::cos(x)*Foam::sin(x)
                                         -Foam::cos(2*x)*2*Foam::sin(2*x));
        UGradMisalignment_Exact[cellI] = UG_numerator
                                      / (Foam::sqrt(UG_denominator) + Foam::mag(UG_numerator));
        // UGradMisalignment_Exact[cellI] = (Foam::sqrt(UG_denominator) + Foam::mag(UG_numerator));
        // UGradMisalignment_Exact[cellI] = UG_numerator;
        // UGradMisalignment_Exact[cellI] = Foam::sqrt(UG_denominator);
        scalar magSqrS = Foam::sqr(Foam::cos(x))
                       + 2 * Foam::sqr(0.5*Foam::sin(x))
                       + 2 * Foam::sqr(Foam::sin(2*x));
        scalar magSqrOmega = 2 * Foam::sqr(0.5*Foam::sin(x))
                           + 2 * Foam::sqr(Foam::sin(2*x));
        QCriterion_Exact[cellI] = (magSqrOmega - magSqrS)
                                / (magSqrOmega + magSqrS);
        ReT_Exact[cellI] = Foam::min((Foam::sqrt(k[cellI]) * dist[cellI]) / (50.0 * nu.value()),
                                     maxReT.value());
        TurbulenceIntensity_Exact[cellI] = k[cellI] /
                                           (
                                           0.5 * (Foam::sqr(Foam::sin(x)) 
                                           + Foam::sqr(Foam::cos(x))
                                           + Foam::sqr(Foam::cos(2*x)))
                                           + k[cellI]
                                           );
        TurbulenceTime_Exact[cellI] = Foam::sqrt(magSqrS) * k[cellI]
                                    / (Foam::sqrt(magSqrS) * k[cellI]
                                      + epsilon[cellI]);
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
            pGradAlongStream_Exact.boundaryField()[patchI][faceI] = Foam::cos(x) * Foam::exp(y)
                                                    / (Foam::sqrt( Foam::exp(y) * Foam::exp(y) 
                                                    * ( 1 + Foam::sqr(Foam::cos(2*x))) )
                                                    + Foam::mag(Foam::cos(x) * Foam::exp(y)));
            scalar UG_numerator = Foam::mag(Foam::sin(x)*Foam::cos(x)*Foam::sin(x)
                                - Foam::sin(x)*Foam::sin(x)*Foam::cos(x)
                                - Foam::sin(x)*2*Foam::sin(2*x)*Foam::cos(2*x));
            scalar UG_denominator = (Foam::sqr(Foam::sin(x)) 
                                  + Foam::sqr(Foam::cos(x))
                                  + Foam::sqr(Foam::cos(2*x)))
                                  * Foam::sqr(Foam::sin(x)*Foam::cos(x)
                                             -Foam::cos(x)*Foam::sin(x)
                                             -Foam::cos(2*x)*2*Foam::sin(2*x));
            // UGradMisalignment_Exact.boundaryField()[patchI][faceI] = UG_numerator
            //                               / (Foam::sqrt(UG_denominator) + Foam::mag(UG_numerator));
            UGradMisalignment_Exact.boundaryField()[patchI][faceI] = (Foam::sqrt(UG_denominator) + Foam::mag(UG_numerator));
            // UGradMisalignment_Exact.boundaryField()[patchI][faceI] = UG_numerator;
            // UGradMisalignment_Exact.boundaryField()[patchI][faceI] = Foam::sqrt(UG_denominator);
            scalar magSqrS = Foam::magSqr(Foam::sin(x))
                           + 2 * Foam::magSqr(0.5*Foam::cos(x))
                           + 2 * Foam::magSqr(Foam::cos(2*x));
            scalar magSqrOmega = 2 * Foam::magSqr(0.5*Foam::cos(x))
                               + 2 * Foam::magSqr(Foam::cos(2*x));
            QCriterion_Exact.boundaryField()[patchI][faceI] = (magSqrOmega - magSqrS)
                                    / (magSqrOmega + magSqrS);
        }
    }

    //Write fields out
    //if (runTime.timeIndex() == 0)
    //PressureGradient = (U & fvc::grad(p)) / (sqrt((fvc::grad(p) & fvc::grad(p)) * (U & U)) + mag(U & fvc::grad(p)));
