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
    Foam::functionObjects::bubbleTracker

Group
    grpFieldFunctionObjects

Description
    Function object used to identify and track bubbles in multiphase simulations.

Usage

    \verbatim
    bubbleTracker1
    {
        // Mandatory entries (unmodifiable)
        type           bubbleTracker;
        libs           (boilingFunctionObject);

        // Mandatory entries (runtime modifiable)

        alphaName       alpha.liquid;
        alphaTreshold   1e-3;

    }
    \endverbatim


Note
    - This should be used in conjunction with icoBoilingFoam

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    bubbleTracker.C

\*---------------------------------------------------------------------------*/

#ifndef bubbleTracker_H
#define bubbleTracker_H

#include "fvMeshFunctionObject.H"
#include "locationModel.H"
#include "activationModel.H"
#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

    struct bubble
    {
        scalar  volume;
        vector  centreOfMass;
        vector  velocityAve;
        vector  velocityDev;
        scalar  temperatureAve;
        scalar  temperatureDev;
        Switch  active;
    };
/*---------------------------------------------------------------------------*\
                   Class bubbleTracker Declaration
\*---------------------------------------------------------------------------*/

class bubbleTracker
:
    public fvMeshFunctionObject
{
    // Private Data
    const Time&     runTime_;

    const word      alphaName_;
    const scalar    alphaTol_;
    label           maxNofBubbles_;

    volScalarField  bubbleIDs_;

    List<bubble>    bubbleList_;

    void updateColorField();

    inline label getBubbleLocId
    (
        label& globalId,
        labelList& bubbleIDsConversion
//        label& minId,
//        label& maxNofBubbles
    )
    {
        // return (globalId > maxNofBubbles-1) ?
        //     globalId - minId + maxNofBubbles
        //     :
        //     globalId;
        forAll(bubbleIDsConversion,id)
        {
            if(bubbleIDsConversion[id] == globalId)
            {
                return id;
            }
        }

        return -1;
    };

    inline label getBubbleLocId
    (
        label& globalId,
        label& minId,
        label& maxNofBubbles
    )
    {
        if(globalId > maxNofBubbles-1)
        {
            return globalId - minId + maxNofBubbles;
        }

        return globalId;
    };

    inline label getBubbleGlobalId
    (
        label& localId,
        label& minId,
        label& maxNofBubbles
    )
    {

        if(localId > maxNofBubbles-1)
        {
            return minId + localId;
        }

        return localId;
    };


public:

    //- Runtime type information
    TypeName("bubbleTracker");


    // Constructors

        //- Construct from Time and dictionary
        bubbleTracker
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        bubbleTracker(const bubbleTracker&) = delete;

        //- No copy assignment
        void operator=(const bubbleTracker&) = delete;


    //- Destructor
    virtual ~bubbleTracker() = default;


    // Member Functions

        //- Read the bubbleTracker data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the bubbleTracker
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
