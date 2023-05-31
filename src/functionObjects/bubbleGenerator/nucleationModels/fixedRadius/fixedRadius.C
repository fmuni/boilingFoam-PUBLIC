/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "fixedRadius.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace nucleationModels
{
    defineTypeNameAndDebug(fixedRadius, 0);
    addToRunTimeSelectionTable(nucleationModel, fixedRadius, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nucleationModels::fixedRadius::fixedRadius
(
    const dictionary dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),
    rad_(readScalar(dict.lookup("radius"))),
    imposeT_(dict.lookupOrDefault<bool>("imposeT",false)),
    radT_(0),
    alphaTol_(dict.lookupOrDefault<scalar>("alphaTol",1e-3)),
    phaseName_(dict.lookupOrDefault<word>("phaseName","liquid")),
    Tint_(0),
    displ_(dict.lookupOrDefault("displacement", vector::zero))

{
    if(imposeT_)
    {
        Tint_ = readScalar(dict.lookup("internalT"));
        radT_ = readScalar(dict.lookup("radiusT"));
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nucleationModels::fixedRadius::~fixedRadius()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::nucleationModels::fixedRadius::
nucleate(const vector& location) const
{
    label vapourFound(0);

    vector bLoc = location + displ_;

    volScalarField& alpha
    (
        const_cast<volScalarField&>
        (
            mesh_.lookupObject<volScalarField>("alpha." + phaseName_)
        )
    );

    volScalarField& T
    (
        const_cast<volScalarField&>
        (
            mesh_.lookupObject<volScalarField>("T")
        )
    );

    List<label>     cellIds;

    forAll(alpha,cellI)
    {
        scalar dist(mag(mesh_.C()[cellI] - bLoc));
        if(dist < max(radT_,rad_))
        {
            if(alpha[cellI] < 1.0 - alphaTol_)
            {
                vapourFound = 1;
                break;
            }

            cellIds.append(cellI);
        }
    }

    reduce(vapourFound,maxOp<label>());

    if(vapourFound > 0)
    {
        return false;
    }

    Info<< "Nucleating bubble at site : "
        << location << endl;

    forAll(cellIds,id)
    {
        label cellI = cellIds[id];

        if(imposeT_)
        {
            T[cellI] = Tint_;
        }

        scalar dist(mag(mesh_.C()[cellI] - bLoc));
        if(dist < rad_)
        {
            alpha[cellI] = 0;
        }

    }

    return true;
}


// ************************************************************************* //
