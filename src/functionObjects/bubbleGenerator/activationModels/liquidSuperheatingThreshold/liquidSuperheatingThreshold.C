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

#include "liquidSuperheatingThreshold.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace activationModels
{
    defineTypeNameAndDebug(liquidSuperheatingThreshold, 0);
    addToRunTimeSelectionTable(activationModel, liquidSuperheatingThreshold, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activationModels::liquidSuperheatingThreshold::liquidSuperheatingThreshold
(
    const dictionary dict,
    const fvMesh& mesh
)
:
    activationModel(dict, mesh),
    Tact_(readScalar(dict.lookup("Tactivation"))),
    alphaTol_(dict.lookupOrDefault<scalar>("alphaTol",1e-3)),
    phaseName_(dict.lookupOrDefault<word>("phaseName","liquid"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModels::liquidSuperheatingThreshold::~liquidSuperheatingThreshold()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::activationModels::liquidSuperheatingThreshold::
isActive(const vector& location) const
{
    const volScalarField& alpha
    (
        mesh_.lookupObject<volScalarField>("alpha." + phaseName_)
    );

    const volScalarField& T
    (
        mesh_.lookupObject<volScalarField>("T")
    );

    label cellId(-1);

    cellId = mesh_.findCell(location);

    scalar alphaValue(0);
    scalar Tvalue(0);

    autoPtr<interpolation<scalar> > interpT =
        interpolation<scalar>::New("cellPointFace", T);

    if(cellId > -1)
    {
        alphaValue = alpha[cellId];
        Tvalue = interpT->interpolate(location,cellId);
    }

    reduce(alphaValue,maxOp<scalar>());
    reduce(Tvalue,maxOp<scalar>());

    if(alphaValue > 1.0 - alphaTol_ && Tvalue > Tact_)
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
