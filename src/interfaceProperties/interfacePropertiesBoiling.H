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

Class
    Foam::interfaceProperties

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    interfacePropertiesBoiling.C

\*---------------------------------------------------------------------------*/

#ifndef interfacePropertiesBoiling_H
#define interfacePropertiesBoiling_H

#include "IOdictionary.H"
#include "surfaceTensionModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "curvatureModel.H"
#include "reconstructionSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Boiling
{
/*---------------------------------------------------------------------------*\
                     Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class interfaceProperties
{
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Compression coefficient
        scalar cAlpha_;

        //- Surface tension
        autoPtr<surfaceTensionModel> sigmaPtr_;

        //- Reconstruction scheme
        autoPtr<reconstructionSchemes>   surf_;
        //- Curvature model
        autoPtr<curvatureModel> KPtr_;

        volScalarField& alpha1_;
        const volVectorField& U_;



        //- Switch for density correction
        Switch rhoCorrection_;


    // Private Member Functions

        //- No copy construct
        interfaceProperties(const interfaceProperties&) = delete;

        //- No copy assignment
        void operator=(const interfaceProperties&) = delete;

        //- Re-calculate the interface curvature
        void calculateK();

        //- Compute the correction factor to surface tension force
        tmp<volScalarField> cFact() const;

public:

    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfaceProperties
        (
            volScalarField& alpha1,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const IOdictionary&
        );


    // Member Functions

        scalar cAlpha() const
        {
            return cAlpha_;
        }

        const dimensionedScalar& deltaN() const
        {
            return KPtr_->deltaN();
        }

        const surfaceScalarField& nHatf() const
        {
            return KPtr_->nHatf();
        }

        reconstructionSchemes& surf()
        {
            return surf_();
        }

        tmp<volScalarField> sigmaK() const;

        tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        void correct();

        //- Read transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
