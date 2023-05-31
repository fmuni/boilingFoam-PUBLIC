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

#include "gradAlphaCurvature.H"
#include "addToRunTimeSelectionTable.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "surfaceFields.H"
#include "surfaceMesh.H"
#include "unitConversion.H"
#include "fvcAverage.H"
#include "fvcReconstruct.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace curvatureModels
{
    defineTypeNameAndDebug(gradAlphaCurvature, 0);
    addToRunTimeSelectionTable(curvatureModel, gradAlphaCurvature, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::curvatureModels::gradAlphaCurvature::gradAlphaCurvature
(
    const dictionary& dict,
    const fvMesh& mesh,
    reconstructionSchemes& surf
)
:
    curvatureModel(dict,mesh,surf),
    nSmoothIter_(dict.lookupOrDefault("KSmoothingIter",0))
{
    if(nSmoothIter_ > 0)
    {
        Info<<"gradAlphaCurvature: Using alpha-based curvature smoothing with "
            << nSmoothIter_ <<" iterations."
            <<endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::curvatureModels::gradAlphaCurvature::~gradAlphaCurvature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::curvatureModels::gradAlphaCurvature::calculateK
(
    const volScalarField& alpha1,
    const volVectorField& U
)
{

    const surfaceVectorField& Sf = mesh_.Sf();

    // Cell gradient of alpha
    volVectorField gradAlpha(fvc::grad(alpha1, "nHat"));

    //- Use smoothed alpha if required
    if (nSmoothIter_>0)
    {
        volScalarField alphaSmooth(alpha1);
        for(label iter=0;iter<nSmoothIter_;iter++)
        {
            alphaSmooth = fvc::average(alphaSmooth);
        }

        gradAlpha = fvc::grad(alphaSmooth, "nHat");
    }

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    //gradAlphaf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
    // surfaceVectorField nHatfv
    // (
    //     (gradAlphaf + deltaN_*vector(0, 0, 1)
    //    *sign(gradAlphaf.component(vector::Z)))/(mag(gradAlphaf) + deltaN_)
    // );
    correctContactAngle
    (
        nHatfv.boundaryFieldRef(),
        gradAlphaf.boundaryField(),
        U,
        alpha1
    );


    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatfv & Sf);

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    /*
    volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryFieldRef()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */
}

bool Foam::curvatureModels::gradAlphaCurvature::readDict(const dictionary& dict)
{
    return true;
}


bool Foam::curvatureModels::gradAlphaCurvature::writeData(Ostream& os) const
{

    return false;
}


// ************************************************************************* //
