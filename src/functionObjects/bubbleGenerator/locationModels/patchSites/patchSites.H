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
    Foam::locationModels::patchSites

Description
    Sites are the internal cells of a user-defined patch.

SourceFiles
    patchSites.C

\*---------------------------------------------------------------------------*/

#ifndef patchSites_H
#define patchSites_H

#include "locationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace locationModels
{

/*---------------------------------------------------------------------------*\
                           Class patchSites Declaration
\*---------------------------------------------------------------------------*/

class patchSites
:
    public locationModel
{
    vectorField sites_;
    const word patchName_;

public:

    //- Runtime type information
    TypeName("patchSites");


    // Constructors

        patchSites
        (
            const dictionary dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~patchSites();


    // Member Functions

        virtual const vectorField& sites() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace locationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
