{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      + thermophysicalTransport.divq(he)
     ==
        rho*(U&g)
      + reaction.Qdot()
      + coalParcels.Sh(he)
      + limestoneParcels.Sh(he)
      + rad.Sh(thermo, he)
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    rad.correct();

    Info<< "T gas min/max   = " << min(thermo.T()).value() << ", "
        << max(thermo.T()).value() << endl;
}
