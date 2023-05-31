/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::locationModel

Description

SourceFiles
    locationModel.C

\*---------------------------------------------------------------------------*/

#ifndef locationModel_H
#define locationModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class locationModel Declaration
\*---------------------------------------------------------------------------*/

class locationModel
{
protected:

    // Protected Data
    const dictionary    dict_;
    const fvMesh&       mesh_;


public:

    //- Runtime type information
    TypeName("locationModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            locationModel,
            dictionary,
            (
                const dictionary dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Static data members

    // Constructors

        locationModel
        (
            const dictionary dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~locationModel();


    // Selectors

        static autoPtr<locationModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        //- location coefficient
        virtual const vectorField& sites() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
