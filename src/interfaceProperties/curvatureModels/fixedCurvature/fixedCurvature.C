/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

#include "fixedCurvature.H"
#include "addToRunTimeSelectionTable.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "surfaceFields.H"
#include "surfaceMesh.H"
#include "unitConversion.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace curvatureModels
{
    defineTypeNameAndDebug(fixed, 0);
    addToRunTimeSelectionTable(curvatureModel, fixed, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::curvatureModels::fixed::fixed
(
    const dictionary& dict,
    const fvMesh& mesh,
    reconstructionSchemes& surf
)
:
    curvatureModel(dict,mesh,surf),
    value_("value",dimK,readScalar(dict.lookup("value")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::curvatureModels::fixed::~fixed()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::curvatureModels::fixed::calculateK
(
    const volScalarField& alpha1,
    const volVectorField& U
)
{

    K_ = value_;

}

bool Foam::curvatureModels::fixed::readDict(const dictionary& dict)
{
    value_ = dimensionedScalar("value",dimK,readScalar(dict.lookup("value")));
    return true;
}


bool Foam::curvatureModels::fixed::writeData(Ostream& os) const
{

    return false;
}


// ************************************************************************* //
