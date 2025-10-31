//Solving the transport Equation for the water vapour mass fraction and calculating the saturation

Info << "Calculating vapour transport " << nl << endl;


volScalarField DEff
(
    (turbulence->nut() / PrtVapourH2O)//Hand Calculated diffusivity
);

//Store old value for relaxation
volScalarField oldMassFractionWater(massFractionWater);

scalar relaxationFactor(0.3);

fvScalarMatrix waterMassFractionTransport
(
    fvm::ddt(rho, massFractionWater) + fvm::div(phi, massFractionWater) - fvm::laplacian(DEff * rho, massFractionWater) + (waterAbsorption)
);

waterMassFractionTransport.relax();

waterMassFractionTransport.solve();

//Relax solution

massFractionWater = relaxationFactor * massFractionWater + (1 - relaxationFactor) * oldMassFractionWater;

//fvOptions.correct(massFractionWater);


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


