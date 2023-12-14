{
    volScalarField& he = thermo.he();
    // (1) Added
    const volScalarField alphaEff("alphaEff", thermophysicalTransport->alphaEff());

    fvScalarMatrix EEqn
    (
        fvm::ddt(alphac, rhoc, he) + mvConvection->fvmDiv(alphaRhoPhic, he)
      + fvc::ddt(alphac, rhoc, Kc) + fvc::div(alphaRhoPhic, Kc)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(alphaRhoPhic/fvc::interpolate(rhoc), Uc),
                p,
                "div(phiv,p)"
            )
          : -dpdt*alphac
        )
	// (2) Note the diffusion term has been re-designed since OpenFOAM 8
	// https://cfd.direct/openfoam/free-software/redesign-for-the-future/
       //+ thermophysicalTransport->divq(he)
       -fvm::laplacian
	(
	   alphac*alphaEff,
           he
	)
     ==
        alphac*rhoc*(Uc&g)
      + combustion->Qdot()
      + biomassMPPICParcels.Sh(he)
      + ironOreMPPICParcels.Sh(he)
      + limestoneMPPICParcels.Sh(he)
      + radiation->Sh(thermo, he)
      + fvOptions(rhoc, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();

    radiativeSourceTerm = radiation->Sh(thermo, he) & he;

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
