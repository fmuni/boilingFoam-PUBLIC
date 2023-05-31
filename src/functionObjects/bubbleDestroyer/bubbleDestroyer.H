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
    Foam::functionObjects::bubbleDestroyer

Group
    grpFieldFunctionObjects

Description
    Destroys bubbles in the domain when they reach a location.
    The alpha field is interpolated at several probes along the domain and
    (if the value of alpha at the last probe is <1) is reset to 1 in the
    region corresponding to the bubble.

Usage
\c system/controlDict.functions:
    \verbatim
    bubbleDestroyer1
    {
        // Mandatory entries (unmodifiable)
        type           bubbleDestroyer;
        libs           (bubbleDestroyerFunctionObject);

        // Mandatory entries (runtime modifiable)
        numberOfProbes  20;
        alphaName       alpha.liquid;
        direction       (1 0 0);

        // Optional entries (runtime modifiable)
        tolerance 1e-3;

    }
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::bubbleDestroyer
    ...

SourceFiles
    bubbleDestroyer.C
    bubbleDestroyerTEMPLATES.C
    ...

\*---------------------------------------------------------------------------*/

#ifndef bubbleDestroyer_H
#define bubbleDestroyer_H

#include "fvMeshFunctionObject.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class bubbleDestroyer Declaration
\*---------------------------------------------------------------------------*/

class bubbleDestroyer
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Number of probes to be used
        label numberOfProbes_;

        //- Name of the alpha field
        word alphaName_;

        //- Direction of the aligned probes
        vector direction_;

        //- Offset for the last bubble. Used when the domain is symmetric
        //  or when the bubble must be destroyed before the end of the domain.
        vector offset_;

        //- Tolerance on alpha
        scalar tol_;

        //- List of probes
        pointField probeList_;

        //- Remove momentum when the bubble is destroyed
        Switch removeMomentum_;

        //- Switch for resetting alpha everywhere
        Switch resetAlpha_;

        //- Switch for resetting U everywhere
        Switch resetU_;

        //- Initial time for reset
        word initialT_;

        //- Delay before reinitializing alpha
        scalar alphaReinitDelay_;

        //- Counter for time passed from bubble destruction
        scalar counterT_;

        //- Switch for resetting U delayed
        Switch alphaMustBeRestored_;

        //- Reset delayed alpha
        void resetAlphaDelayed();
public:

    //- Runtime type information
    TypeName("bubbleDestroyer");


    // Constructors

        //- Construct from Time and dictionary
        bubbleDestroyer
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        bubbleDestroyer(const bubbleDestroyer&) = delete;

        //- No copy assignment
        void operator=(const bubbleDestroyer&) = delete;


    //- Destructor
    virtual ~bubbleDestroyer() = default;


    // Member Functions

        //- Read the bubbleDestroyer data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the bubbleDestroyer
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
