/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR, AFFILIATION
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
    Foam::functionObjects::postBoiling


Description
    Compute quantities of interest for cases with phase change.


Usage
    \verbatim
    postData
    {
        // Mandatory entries (unmodifiable)
        type           postBoiling;
        libs           (boilingFunctionObject);

        alphaName   alpha.liquid;

        //- Patches were relevant quantities are computed
        patches ( fluid_to_solid );

        //- Saturation temperature
        Tsat    373.15;

        //- Reference length;
        Lr  229e-6;
    }
    \endverbatim

SourceFiles
    postBoiling.C

\*---------------------------------------------------------------------------*/

#ifndef postBoiling_H
#define postBoiling_H

#include "fvMeshFunctionObject.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class postBoiling Declaration
\*---------------------------------------------------------------------------*/

class postBoiling
:
    public fvMeshFunctionObject
{

    //- Name of the liquid phase
    word alphaName_;

    //- List of patches
    wordList patches_;

    //- Saturation temperature
    scalar Tsat_;

    //- Reference length
    scalar Lr_;

    //- Wall nusselt (zero inside the domain)
    volScalarField Nu_;

    //- Optional entries for splitting patch into two (for multiregion)
    Switch splitPatch_;
    vector splitVec_;



public:

    //- Runtime type information
    TypeName("postBoiling");


    // Constructors

        //- Construct from Time and dictionary
        postBoiling
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        postBoiling(const postBoiling&) = delete;

        //- No copy assignment
        void operator=(const postBoiling&) = delete;


    //- Destructor
    virtual ~postBoiling() = default;


    // Member Functions

        //- Read the postBoiling data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the postBoiling
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
