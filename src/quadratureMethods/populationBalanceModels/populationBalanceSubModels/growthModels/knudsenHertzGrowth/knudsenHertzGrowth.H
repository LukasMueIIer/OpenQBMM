/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015 by Matteo Icardi and 2017 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 Alberto Passalacqua
-------------------------------------------------------------------------------
2017-03-28 Alberto Passalacqua: Adapted to single scalar calculation.
2020-03-12 Alberto Passalacqua: Added access to fvMesh.
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::populationBalanceSubModels::growthModels::knudsenHertzGrowth

Description
    Constant growth model, which is only activated when local supersaturation is present. 
    The value of the constant is read from dictionary
    or assumed equal to 1. The minimum and maximum value of the abscissa to
    which the model is applied can be defined by the user or it is assumed to
    be 1.0e-10 and 1.0e10.

SourceFiles
    knudsenHertzGrowth.C

\*---------------------------------------------------------------------------*/

#ifndef knudsenHertzGrowth_H
#define knudsenHertzGrowth_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
                    Class knudsenHertzGrowth Declaration
\*---------------------------------------------------------------------------*/

class knudsenHertzGrowth
:
    public growthModel
{
private:

    // Private data

        //- Minimum value of the abscissa to apply the model
        scalar minAbscissa_;

        //- Maximum value of the abscissa to apply the model
        scalar maxAbscissa_;

        //- condensation coefficient
        scalar alphaCondensation_;

        //- gas constant of water vapour
        scalar RVapour_;

        //- Access to the saturation field, moved to Kg;
        mutable refPtr<volScalarField> saturationWater_;
        mutable refPtr<volScalarField> pVapour_; //Vaper pressure of water
        mutable refPtr<volScalarField> pSaturationWater_; //Saturation Vapour Pressure wrt. liquid Water
        mutable refPtr<volScalarField> T_; //Temperature Field
        mutable label lookedUpSaturation;        //indicator if we allready performed the lookup

        void lookUpSaturation() const; //function to look up saturation
        
        //Access to the sink term field of the Water Transport
        mutable refPtr<volScalarField> waterAbsorption_;
        //will be looked up at the same time as the Saturation, so use same flag

        void calculateSinkTerm
        (
            const labelList& momentOrder,
            const label celli,
            const scalarQuadratureApproximation& quadrature
        );

        void calculateSinkTerm
        (
            const labelList& momentOrder,
            const label celli,
            const velocityQuadratureApproximation& quadrature
        ); 
       

public:

        //- Runtime type information
        TypeName("knudsenHertzGrowth");


    // Constructors

        //- Construct from components
        knudsenHertzGrowth
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~knudsenHertzGrowth();


    // Member Functions

        //- Growth model
        scalar Kg
        (
            const scalar& abscissa,
            const bool lengthBased,
            const label environment = 0
        ) const;

        //Overwrite the Base Class Functions
        scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const scalarQuadratureApproximation& quadrature
        );

        //- Return phase space convection
        scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const velocityQuadratureApproximation& quadrature
        );


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace growthModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
