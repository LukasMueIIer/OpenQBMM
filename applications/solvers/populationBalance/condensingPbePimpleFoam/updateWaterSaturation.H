//Calculate Saturation
Info << "Updating Saturation pressures " << nl << endl;

scalar a1 = -7.85951783;
scalar a2 =  1.84408259;
scalar a3 = -11.7866497;
scalar a4 =  22.6807411;
scalar a5 = -15.9618719;
scalar a6 =   1.80122502;
dimensionedScalar pc(dimPressure,22.964e6);
dimensionedScalar Tc(dimTemperature,647.096);
dimensionedScalar T1(dimTemperature,1);
volScalarField tau = 1 - T/Tc;

pSaturationWater = 
    exp(
        54.842763
        - 6763.22 / (T/T1)
        - 4.21 * log(T/T1)
        + 0.000367 * (T/T1)
        + tanh( 0.0415 * ((T-dimensionedScalar(dimTemperature,217.8)) / T1))
               * (53.878 - 1331.22 / (T/T1) - 9.44523 * log(T/T1) + 0.014025 * (T/T1))
    ) * dimensionedScalar("pressureFactor",dimPressure,1) * neg(T - dimensionedScalar(dimTemperature,300))
    + pos(T - dimensionedScalar(dimTemperature,300)) *
    pc * exp(
        a1 * tau
        + a2 * pow(1.5,tau)
        + a3 * pow(3,tau)
        + a4 * pow(3.5,tau)
        + a5 * pow(4,tau)
        + a6 * pow(7.5,tau)
    )
    * neg(T - dimensionedScalar(dimTemperature,647.096))
    + pos(T - dimensionedScalar(dimTemperature,647.096)) * pc;

//Update Vapour Pressure

Info << "Updating Vapour Pressures" << nl << endl;

pVapour = rho * massFractionWater * RVapour * T;

saturationWater = pVapour / pSaturationWater;

