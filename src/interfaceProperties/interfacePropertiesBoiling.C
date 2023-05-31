/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "interfacePropertiesBoiling.H"

#include "mathematicalConstants.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Boiling::interfaceProperties::interfaceProperties
(
    volScalarField& alpha1,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        alpha1.mesh().solverDict(alpha1.name()).get<scalar>("cAlpha")
    ),

    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),
    surf_(reconstructionSchemes::New(alpha1, phi, U, dict)),
    KPtr_(curvatureModel::New(dict, alpha1.mesh(), surf_())),

    alpha1_(alpha1),
    U_(U),



    rhoCorrection_
    (
        transportPropertiesDict_.lookupOrDefault("sigmaRhoCorrection",false)
    )

{
    if(rhoCorrection_)
    {
        Info<<"interfaceProperties: Using density correction when computing"
            <<" the surface tension force."
            <<endl;
    }

    calculateK();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::Boiling::interfaceProperties::sigmaK() const
{
    return cFact()*(sigmaPtr_->sigma()*KPtr_->K());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::Boiling::interfaceProperties::surfaceTensionForce() const
{

    return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_);
}


Foam::tmp<Foam::volScalarField>
Foam::Boiling::interfaceProperties::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}

void Foam::Boiling::interfaceProperties::calculateK()
{
    KPtr_->calculateK(alpha1_,U_);
}

void Foam::Boiling::interfaceProperties::correct()
{
    //- Force reconstruction of the interface
    surf_->reconstruct();
    calculateK();
}

Foam::tmp<Foam::volScalarField>
Foam::Boiling::interfaceProperties::cFact() const
{
    tmp<volScalarField> tCorr
    (
        new volScalarField
        (
            IOobject
            (
                "rhoCorr",
                alpha1_.mesh().time().timeName(),
                alpha1_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            alpha1_.mesh(),
            dimensionedScalar("corr",dimless,1.)
        )
    );

    volScalarField& corr = tCorr.ref();

    if(rhoCorrection_)
    {
        wordList phases(transportPropertiesDict_.lookup("phases"));

        dimensionedScalar rhol
        (
            "rho",
            dimDensity,
            transportPropertiesDict_.subDict(phases[0])
        );

        dimensionedScalar rhov
        (
            "rho",
            dimDensity,
            transportPropertiesDict_.subDict(phases[1])
        );

        const  volScalarField& rho
        (
            alpha1_.mesh().lookupObject<volScalarField>("rho")
        );

        dimensionedScalar rhoMean((rhol+rhov)/2.0);
        corr = rho/rhoMean;
    }

    return tCorr;
}


bool Foam::Boiling::interfaceProperties::read()
{
    alpha1_.mesh().solverDict(alpha1_.name()).readEntry("cAlpha", cAlpha_);
    sigmaPtr_->readDict(transportPropertiesDict_);
    KPtr_->readDict(transportPropertiesDict_);

    return true;
}


// ************************************************************************* //
