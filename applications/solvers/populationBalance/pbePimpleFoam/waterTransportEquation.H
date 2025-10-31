//Solving the transport Equation for the water vapour mass fraction and calculating the saturation

Info << "Calculating vapour transport " << nl << endl;


volScalarField DEff
(
    (turbulence->nut() / PrtVapourH2O)//Hand Calculated diffusivity
);


fvScalarMatrix waterMassFractionTransport
(
    fvm::ddt(rho, massFractionWater) + fvm::div(phi, massFractionWater) - fvm::laplacian(DEff * rho, massFractionWater) + (waterAbsorption)
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
#include "updateWaterSaturation.H"


