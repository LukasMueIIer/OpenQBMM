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

\*---------------------------------------------------------------------------*/

#include "saturationActivatedGrowth.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(saturationActivatedGrowth, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        saturationActivatedGrowth,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::saturationActivatedGrowth
::saturationActivatedGrowth
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    growthModel(dict, mesh),
    minAbscissa_(dict.lookupOrDefault("minAbscissa", scalar(0))),
    maxAbscissa_(dict.lookupOrDefault("maxAbscissa", GREAT)),
    lookedUpSaturation(0),
    saturationWater_(nullptr),
    waterAbsorption_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::saturationActivatedGrowth
::~saturationActivatedGrowth()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::populationBalanceSubModels::growthModels::saturationActivatedGrowth::lookUpSaturation
(

) const
{
    saturationWater_.reset(
        mesh_.getObjectPtr<volScalarField>(
            word("S_water")
        )
    );
    
    waterAbsorption_.reset(
        mesh_.getObjectPtr<volScalarField>(
            word("waterAbsoption")
        )
    );

    lookedUpSaturation = 1;
}

Foam::scalar
Foam::populationBalanceSubModels::growthModels::saturationActivatedGrowth::Kg
(
    const scalar& abscissa,
    const bool lengthBased,
    const label environment
) const
{   
    if(lookedUpSaturation == 0)
    {
       lookUpSaturation(); 
    }
    
    return Cg_.value()
          *pos0(abscissa - minAbscissa_)
          *neg0(abscissa - maxAbscissa_)
          *pos0((*saturationWater_)[environment] - 1); //only activate if sat>1
}

Foam::scalar
Foam::populationBalanceSubModels::growthModels::saturationActivatedGrowth::phaseSpaceConvection
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)
{
    if(lookedUpSaturation == 0)
    {
       lookUpSaturation(); 
    }

    const PtrList<volScalarNode>& nodes = quadrature.nodes();
    label sizeIndex = nodes[0].sizeIndex();

    if(sizeIndex != -1)
    {

        //Calculate the order of the calculation
        label sizeOrder = momentOrder[sizeIndex];
        bool lengthBased = nodes[0].lengthBased();
        bool volumeFraction = nodes[0].useVolumeFraction();



        if (volumeFraction)
        {
            if (lengthBased)
            {
                sizeOrder += 3;
            }
            else
            {
                sizeOrder += 1;
            }
        }

        if(
            sizeOrder == 2
        )
        {
        this->calculateSinkTerm
        (
            momentOrder,
            celli,
            quadrature
        ); 
        }
    }

    return growthModel::phaseSpaceConvection
    (
        momentOrder,
        celli,
        quadrature
    );
}


Foam::scalar
Foam::populationBalanceSubModels::growthModels::saturationActivatedGrowth::phaseSpaceConvection
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation& quadrature
)
{
    if(lookedUpSaturation == 0)
    {
       lookUpSaturation(); 
    }

    const PtrList<volVelocityNode>& nodes = quadrature.nodes();
    label sizeIndex = nodes[0].sizeIndex();

    if(sizeIndex != -1)
    {

        //Calculate the order of the calculation
        label sizeOrder = momentOrder[sizeIndex];
        bool lengthBased = nodes[0].lengthBased();
        bool volumeFraction = nodes[0].useVolumeFraction();



        if (volumeFraction)
        {
            if (lengthBased)
            {
                sizeOrder += 3;
            }
            else
            {
                sizeOrder += 1;
            }
        }

        if(
            sizeOrder == 2
        )
        {
        calculateSinkTerm
        (
            momentOrder,
            celli,
            quadrature
        ); 
        }
    }

    return growthModel::phaseSpaceConvection
    (
        momentOrder,
        celli,
        quadrature
    );
}

void
Foam::populationBalanceSubModels::growthModels::saturationActivatedGrowth::calculateSinkTerm
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)
{
    if(lookedUpSaturation == 0)
    {
       lookUpSaturation(); 
    }

    scalar gSource(0);

    const PtrList<volScalarNode>& nodes = quadrature.nodes();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
        (*waterAbsorption_)[celli] = gSource;
    }

    label sizeOrder = momentOrder[sizeIndex];
    bool lengthBased = nodes[0].lengthBased();
    bool volumeFraction = nodes[0].useVolumeFraction();

    if (volumeFraction)
    {
        if (lengthBased)
        {
            sizeOrder += 3;
        }
        else
        {
            sizeOrder += 1;
        }
    }

    if (sizeOrder < 1)
    {
        (*waterAbsorption_)[celli] = gSource;
    }

    const labelList& scalarIndexes = nodes[0].scalarIndexes();

    if (!nodes[0].extended())
    {
        forAll(nodes, pNodeI)
        {
            const volScalarNode& node = nodes[pNodeI];

            scalar bAbscissa =
                max(node.primaryAbscissae()[sizeIndex][celli], scalar(0));

            scalar d = node.d(celli, bAbscissa);

            scalar n =
                node.n(celli, node.primaryWeight()[celli], bAbscissa);

            scalar gSourcei =
                n
               *Kg(d, lengthBased, celli)
               *pow(bAbscissa, 2);

            forAll(scalarIndexes, nodei)
            {
                if (scalarIndexes[nodei] != sizeIndex)
                {
                    gSourcei *=
                        pow
                        (
                            node.primaryAbscissae()[nodei][celli],
                            momentOrder[scalarIndexes[nodei]]
                        );
                }
            }

            gSource += gSourcei;
        }

        (*waterAbsorption_)[celli] = gSource;
    }

    forAll(nodes, pNodeI)
    {
        const volScalarNode& node = nodes[pNodeI];

        forAll(node.secondaryWeights()[sizeIndex], sNodei)
        {
            scalar bAbscissa =
                max
                (
                    node.secondaryAbscissae()[sizeIndex][sNodei][celli],
                    scalar(0)
                );

            scalar d = node.d(celli, bAbscissa);

            scalar n =
                node.n(celli, node.primaryWeight()[celli], bAbscissa)
               *node.secondaryWeights()[sizeIndex][sNodei][celli];

            scalar gSourcei =
                n
               *Kg(d, lengthBased, celli)
               *pow(bAbscissa, 2);

            forAll(scalarIndexes, cmpt)
            {
                if (scalarIndexes[cmpt] != sizeIndex)
                {
                    gSourcei *=
                        node.secondaryWeights()[cmpt][sNodei][celli]
                       *pow
                        (
                            node.secondaryAbscissae()[cmpt][sNodei][celli],
                            momentOrder[scalarIndexes[cmpt]]
                        );
                }
            }

            gSource += gSourcei;
        }
    }

    (*waterAbsorption_)[celli] = gSource;
}


void
Foam::populationBalanceSubModels::growthModels::saturationActivatedGrowth::calculateSinkTerm
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation& quadrature
)
{
    if(lookedUpSaturation == 0)
    {
       lookUpSaturation(); 
    }

    scalar gSource(0);

    const PtrList<volVelocityNode>& nodes = quadrature.nodes();
    label sizeIndex = nodes[0].sizeIndex();

    if (sizeIndex == -1)
    {
      (*waterAbsorption_)[celli] = gSource;
    }

    label sizeOrder = momentOrder[sizeIndex];
    bool lengthBased = nodes[0].lengthBased();
    bool volumeFraction = nodes[0].useVolumeFraction();

    if (volumeFraction)
    {
        if (lengthBased)
        {
            sizeOrder += 3;
        }
        else
        {
            sizeOrder += 1;
        }
    }

    if (sizeOrder < 1)
    {
        (*waterAbsorption_)[celli] = gSource;
    }

    const labelList& scalarIndexes = nodes[0].scalarIndexes();
    const labelList& velocityIndexes = nodes[0].velocityIndexes();

    forAll(nodes, pNodeI)
    {
        const volVelocityNode& node = nodes[pNodeI];

        scalar bAbscissa =
            max(node.primaryAbscissae()[sizeIndex][celli], scalar(0));

        scalar d = node.d(celli, bAbscissa);

        scalar n =
            node.n(celli, node.primaryWeight()[celli], bAbscissa);

        scalar gSourcei =
            n*Kg(d, lengthBased, celli)*pow(bAbscissa, 2);

        forAll(scalarIndexes, nodei)
        {
            if (scalarIndexes[nodei] != sizeIndex)
            {
                gSourcei *=
                    pow
                    (
                        node.primaryAbscissae()[nodei][celli],
                        momentOrder[scalarIndexes[nodei]]
                    );
            }
        }

        forAll(velocityIndexes, cmpt)
        {
            gSourcei *=
                pow
                (
                    node.velocityAbscissae()[celli][cmpt],
                    momentOrder[velocityIndexes[cmpt]]
                );
        }

        gSource += gSourcei;
    }

    (*waterAbsorption_)[celli] = gSource;
}




// ************************************************************************* //
