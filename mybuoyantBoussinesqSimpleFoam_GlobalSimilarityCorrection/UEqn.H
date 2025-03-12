    Info<< "UEqn.H\n" << endl;

    // Solve the momentum equation

    nut.correctBoundaryConditions();

    volScalarField nuEff("nuEff", nu + nut);
    

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      - fvm::laplacian(nuEff, U)
      - fvc::div(nuEff*dev(fvc::grad(U)().T()))
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
