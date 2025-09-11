//Solving the transport Equation for the water vapour mass fraction and calculating the saturation

Info << "Calculating vapour transport " << nl << endl;


volScalarField DEff
(
    (turbulence->nut() / PrtVapourH2O)//Hand Calculated diffusivity
);


fvScalarMatrix waterMassFractionTransport
(
    fvm::div(phi, massFractionWater) - fvm::laplacian(DEff * rho, massFractionWater)
);

waterMassFractionTransport.relax();

waterMassFractionTransport.solve();

fvOptions.correct(massFractionWater);


//Bound negative Values to zero
forAll(mesh.C(),cellID)
{
    if(neg(massFractionWater[cellID]))
    {
        massFractionWater[cellID] = scalar(0.);
    }
}


//Calculate Saturation
Info << "Updating Saturation pressures " << nl << endl;


pSaturationWater = 
pow(10, -7.90298*(fracT/T-1) 
+  5.02808 * log10(fracT/T)
- 1.3816e-7 * (pow(10,11.344*(1-T/fracT)) - 1)
+ 8.1328e-3 * (pow(10,-3.49149*(fracT/T-1) ) - 1) 
+ log10(dimensionedScalar(dimless,1013.246))) * dimensionedScalar("pressureFactor",dimPressure,100) 
* neg(T - fracT)
+ capPressure * pos(T - fracT);

//Update Vapour Pressure

Info << "Updating Vapour Pressures" << nl << endl;

pVapour = rho * massFractionWater * RVapour * T;

saturationWater = pVapour / pSaturationWater;

