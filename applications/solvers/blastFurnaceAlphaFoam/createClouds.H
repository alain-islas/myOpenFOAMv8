Info<< "\nConstructing biomassMPPIC cloud" << endl;
coalMPPICCloud biomassMPPICParcels
(
    "biomassCloud1",
    rhoc,
    Uc,
    g,
    slgThermo
);

Info<< "\nConstructing ironOreMPPIC cloud" << endl;
coalMPPICCloud ironOreMPPICParcels
(
    "ironOreCloud1",
    rhoc,
    Uc,
    g,
    slgThermo
);

Info<< "\nConstructing limestoneMPPIC cloud" << endl;
coalMPPICCloud limestoneMPPICParcels
(
    "limestoneCloud1",
    rhoc,
    Uc,
    g,
    slgThermo
);
