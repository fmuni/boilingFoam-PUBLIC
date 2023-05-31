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

Class
    Foam::curvatureModels::gradAlphaCurvature

Description
    gradAlphaCurvature curvature model.

Usage
    Example of the surface tension specification:
    \verbatim
        sigma
        {
            type                gradAlphaCurvature;
            sigma               0.07;
        }
    \endverbatim

See also
    Foam::curvatureModel

SourceFiles
    gradAlphaCurvaturecurvature.C

\*---------------------------------------------------------------------------*/

#ifndef gradAlphaCurvatur_H
#define gradAlphaCurvature_H

#include "curvatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace curvatureModels
{

/*---------------------------------------------------------------------------*\
                           Class gradAlphaCurvature Declaration
\*---------------------------------------------------------------------------*/

class gradAlphaCurvature
:
    public curvatureModel
{
    // Private data


        label nSmoothIter_;

public:

    //- Runtime type information
    TypeName("gradAlpha");


    // Constructors

        //- Construct from dictionary and mesh
        gradAlphaCurvature
        (
            const dictionary& dict,
            const fvMesh& mesh,
            reconstructionSchemes& surf

        );


    //- Destructor
    virtual ~gradAlphaCurvature();


    // Member Functions

        /// Calculate K
        virtual void calculateK
        (
            const volScalarField& alpha1,
            const volVectorField& U
        );

        //- Update surface tension coefficient from given dictionary
        virtual bool readDict(const dictionary& dict);

        //- Write in dictionary format
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace curvatureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
