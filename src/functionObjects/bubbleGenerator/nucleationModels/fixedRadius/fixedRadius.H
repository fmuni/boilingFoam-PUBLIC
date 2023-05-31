/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::nucleationModels::fixedRadius

Description
    Use alpha field and temperature as thresholds

SourceFiles
    fixedRadius.C

\*---------------------------------------------------------------------------*/

#ifndef fixedRadius_H
#define fixedRadius_H

#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                           Class fixedRadius Declaration
\*---------------------------------------------------------------------------*/

class fixedRadius
:
    public nucleationModel
{
    const scalar    rad_;
    const bool      imposeT_;
    mutable scalar    radT_;
    const scalar    alphaTol_;
    const word      phaseName_;
    mutable scalar    Tint_;


    // Displacement from the site locus
    const vector displ_;

public:

    //- Runtime type information
    TypeName("fixedRadius");


    // Constructors

        fixedRadius
        (
            const dictionary dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~fixedRadius();


    // Member Functions

        virtual bool nucleate(const vector& location) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
