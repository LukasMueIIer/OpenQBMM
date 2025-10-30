//Field creation for the transport of water vapour
volScalarField massFractionWater
(
    IOobject
    (
        "Y_waterVapour",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& T = thermo.T();

dimensionedScalar PrtVapourH2O("PrtVapourH2O",dimless,0.7);
dimensionedScalar RVapour("RVapour",dimGasConstant,461.6);
dimensionedScalar fracT("fractTvapour",dimTemperature,373.16);  //Reference Temperature for saturation
dimensionedScalar capPressure("saturationPressureCap",dimPressure,1e5); //Used to cap the saturation pressure


Info << "Creating Saturation Field wrt. Water" << nl << endl;
volScalarField saturationWater
(
    IOobject
    (
        "S_water",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    0,
    dimless
);

Info << "Creating Vapour Pressure Field" << nl << endl;
volScalarField pVapour
(
    IOobject
    (
        "pVapour",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho * massFractionWater * RVapour * T
);

Info << "Creating Saturation Pressure wrt. Water Field" << nl << endl;
volScalarField pSaturationWater
(
    IOobject
    (
        "pSaturationWater",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho * massFractionWater * RVapour * T
);

#include "updateWaterSaturation.H"

volScalarField waterAbsorption
(
    IOobject
    (
        "waterAbsorption",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimDensity/dimTime,0)
);


