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
    Foam::curvatureModel

Description
    Abstract base-class for curvature models which return the curvature field.

Usage
    Example of the surface tension specification:
    \verbatim
        sigma
        {
            type                <surface tension model type>;
            <coefficient name>  <coefficient value>;
            .
            .
            .
        }
    \endverbatim
    For simplicity and backward-compatibility the constant value format is
    also supported, e.g.
    \verbatim
        sigma           0.07;
    \endverbatim

SourceFiles
    curvatureModel.C
    newcurvatureModel.C

\*---------------------------------------------------------------------------*/

#ifndef curvatureModel_H
#define curvatureModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "regIOobject.H"
#include "dimensionedTypes.H"
#include "surfaceFields.H"
#include "surfaceInterpolate.H"
#include "reconstructionSchemes.H"


#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class curvatureModel Declaration
\*---------------------------------------------------------------------------*/

class curvatureModel
:
    public regIOobject
{

    const reconstructionSchemes& surf_;

protected:

    // Protected member data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Curvature field
        mutable volScalarField K_;

        surfaceScalarField nHatf_;

         //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;


        // //- Smoothing iterations
        // label nSmoothIter_;


    // Protected member functions

        static const dictionary& KDict(const dictionary& dict)
        {
            return dict.subDict("curvatureModel");
        }

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::Boundary& nHat,
            const surfaceVectorField::Boundary& gradAlphaf,
            const volVectorField& U,
            const volScalarField& alpha1
        ) const;


public:

    //- Runtime type information
    TypeName("curvatureModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            curvatureModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                reconstructionSchemes& surf
            ),
            (dict, mesh, surf)
        );


    // Static data members

        //- Surface tension coefficient dimensions
        static const dimensionSet dimK;


    // Constructors

        // Construct from mesh
        curvatureModel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            reconstructionSchemes& surf
        );


    //- Destructor
    virtual ~curvatureModel();


    // Selectors

        static autoPtr<curvatureModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            reconstructionSchemes& surf
        );


    // Member Functions

        //- Return urvature
        virtual const volScalarField& K() const
        {
            // smoothK();
            // return Ksmt_;
            return K_;
        };

         //- Return urvature
        virtual const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        };

        //- Return urvature
       virtual const surfaceScalarField& nHatf() const
       {
           return nHatf_;
       };

        /// Calculate K
        virtual void calculateK
        (
            const volScalarField& alpha1,
            const volVectorField& U
        ) = 0;

        //- Update surface tension coefficient from given dictionary
        virtual bool readDict(const dictionary& dict) = 0;

        //- Write in dictionary format
        virtual bool writeData(Ostream& os) const = 0;

        const reconstructionSchemes& surf()
        {
            return surf_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
