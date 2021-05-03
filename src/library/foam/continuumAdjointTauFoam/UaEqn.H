                UaEqn().relax();

                solve(UaEqn() == - fvc::grad(pa));

                volScalarField rAUa(1.0/UaEqn().A());
                volVectorField HbyAa("HbyAa", Ua);
                HbyAa = rAUa*UaEqn().H();
                UaEqn.clear();
                surfaceScalarField phiHbyAa
                (
                    "phiHbyAa",
                    fvc::interpolate(HbyAa) & mesh.Sf()
                );
                adjustPhi(phiHbyAa, Ua, pa);

                // Non-orthogonal pressure corrector loop
                while (simple.correctNonOrthogonal())
                {
                    fvScalarMatrix paEqn
                    (
                        fvm::laplacian(rAUa, pa) == fvc::div(phiHbyAa)
                    );

                    paEqn.setReference(paRefCell, paRefValue);
                    paEqn.solve();

                    if (simple.finalNonOrthogonalIter())
                    {
                        phia = phiHbyAa - paEqn.flux();
                    }
                }

                #include "adjointContinuityErrs.H"

                // Explicitly relax pressure for adjoint momentum corrector
                pa.relax();

                // Adjoint momentum corrector
                Ua = HbyAa - rAUa*fvc::grad(pa);
                Ua.correctBoundaryConditions();
