Info<< "\nConstructing coal cloud" << endl;
coalCloud coalParcels
(
    "coalCloud1",
/*    rhoFluid[i],
    UFluid[i],
    gFluid[i],
    slgThermoFluid[i]
*/
    rho,
    U,
    g,
    slgThermo
);

Info<< "\nConstructing limestone cloud" << endl;
basicThermoCloud limestoneParcels
(
    "limestoneCloud1",
/*    rhoFluid[i],
    UFluid[i],
    gFluid[i],
    slgThermoFluid[i]
*/
    rho,
    U,
    g,
    slgThermo
);
