tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        alphaRhoPhic,
        mesh.divScheme("div(alphaRhoPhi,Yi_h)")
    )
);


{
    combustion->correct();
    volScalarField Yt(0.0*Y[0]);
    // (1) Added
    const volScalarField muEff("muEff", turbulence->muEff());

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(alphac, rhoc, Yi)
              + mvConvection->fvmDiv(alphaRhoPhic, Yi)
	      // (2) Note the diffusion term has been re-designed since OpenFOAM 8
	      // https://cfd.direct/openfoam/free-software/redesign-for-the-future/
              //+ thermophysicalTransport->divj(Yi)
	      -fvm::laplacian
	      (
	          alphac*muEff,
	          Yi
	      )
              ==
                biomassMPPICParcels.SYi(i, Yi)
	      + ironOreMPPICParcels.SYi(i, Yi)
	      + limestoneMPPICParcels.SYi(i, Yi)
              + combustion->R(Yi)
              + fvOptions(rhoc, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
