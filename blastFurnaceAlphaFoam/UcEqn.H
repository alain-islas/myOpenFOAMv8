    MRF.correctBoundaryVelocity(Uc);

    fvVectorMatrix UcEqn
    (
        fvm::ddt(alphac, rhoc, Uc) + fvm::div(alphaRhoPhic, Uc)
      + MRF.DDt(alphac*rhoc, Uc)
      + turbulence->divDevTau(Uc)
     ==
        alphac*rhoc()*g
      + biomassMPPICParcels.SU(Uc)
      + ironOreMPPICParcels.SU(Uc)
      + limestoneMPPICParcels.SU(Uc)
      + fvOptions(rhoc, Uc)
    );

    UcEqn.relax();

    fvOptions.constrain(UcEqn);

    if (pimple.momentumPredictor())
    {
        solve(UcEqn == -alphac*fvc::grad(p));

        fvOptions.correct(Uc);
        Kc = 0.5*magSqr(Uc);
    }
