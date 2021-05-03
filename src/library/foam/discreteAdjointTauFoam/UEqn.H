    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      // + MRF.DDt(U)
      - fvm::laplacian(turbulence->nu(),U)
      + fvc::div(dev(tauR))
      ==
        fvOptions(U)
    );

    UEqn().relax();

    fvOptions.constrain(UEqn());

    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
