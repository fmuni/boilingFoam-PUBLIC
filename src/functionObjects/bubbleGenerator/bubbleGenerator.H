/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR, AFFILIATION
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
    Foam::functionObjects::bubbleGenerator

Group
    grpFieldFunctionObjects

Description
    Function object used to generate bubbles in multiphase simulations.

Usage

    \verbatim
    bubbleGenerator1
    {
        // Mandatory entries (unmodifiable)
        type           bubbleGenerator;
        libs           (boilingFunctionObject);

        // Mandatory entries (runtime modifiable)

        // Submodels
        locationModel   fixedLocations;

        frequencyModel  fixedFrequency;

        initializationModel fixedRadii;

        fixedLocationsDict
        {
            locations List<vector>
            2
            (
                (1 0 0)
                (1 1 0)
            );

        }

        fixedFrequencyDict
        {
            frequency   1e4;
        }

        fixedRadiiDict
        {
            radius  2.29e-4;
            radiusT 2.5e-4;
            Tsat    373.14;
        }

    }
    \endverbatim


Note
    - This should be used in conjunction with icoBoilingFoam

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    bubbleGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef bubbleGenerator_H
#define bubbleGenerator_H

#include "fvMeshFunctionObject.H"
#include "locationModel.H"
#include "activationModel.H"
#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class bubbleGenerator Declaration
\*---------------------------------------------------------------------------*/

class bubbleGenerator
:
    public fvMeshFunctionObject
{
    // Private Data
    const Time&                      runTime_;

    autoPtr<locationModel>          location_;

    autoPtr<activationModel>      activation_;

    autoPtr<nucleationModel>      nucleation_;

    const bool      fixDT_;
    mutable scalar  newDT_;

    void adjustDT(bool adjust) const;


public:

    //- Runtime type information
    TypeName("bubbleGenerator");


    // Constructors

        //- Construct from Time and dictionary
        bubbleGenerator
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        bubbleGenerator(const bubbleGenerator&) = delete;

        //- No copy assignment
        void operator=(const bubbleGenerator&) = delete;


    //- Destructor
    virtual ~bubbleGenerator() = default;


    // Member Functions

        //- Read the bubbleGenerator data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the bubbleGenerator
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
