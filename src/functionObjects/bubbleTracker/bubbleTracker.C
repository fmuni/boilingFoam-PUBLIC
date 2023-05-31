/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 Federico Municchi. All rights reserved.
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

\*---------------------------------------------------------------------------*/

#include "bubbleTracker.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "coupledFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(bubbleTracker, 0);
    addToRunTimeSelectionTable(functionObject, bubbleTracker, dictionary);
}
}

// * * * * * * * * * * * * * *  Private Functions  * * * * * * * * * * * * * //
void Foam::functionObjects::bubbleTracker::updateColorField()
{
    //- Get alpha field
    const volScalarField& alpha
    (
        mesh_.lookupObject<volScalarField>(alphaName_)
    );

    //- Create objects for the tracking of the local ids

    //- This is the minimum ID.
    //  Values between 0 and (maxNofBubbles_ - 1)  are reserved for the global
    //  (true) indexing of the bubbles.
    label minId((1 + Pstream::myProcNo())*maxNofBubbles_);

    //- This list is used to ensure consistency in the bubble IDs.
    labelList bubbleIDsConversion(2*maxNofBubbles_,-1);
    forAll(bubbleList_,id)
    {
        bubbleIDsConversion[id] = id;
    }

    //- Keep track of the max id assigned
    label maxId(minId);

    //- Now build the new colour field.
    //  At this stage, it doesn't matter if a bubble has multiple indexes.
    forAll(alpha,cellI)
    {

        //- Check if there is some vapour
        if(alpha[cellI] < (1.0 - alphaTol_))
        {
            //- Check if belongs to a bubble. In that case do nothing.
            if(bubbleIDs_[cellI] > SMALL)
            {
                continue;
            }

            //- Check the neighbours to understand if it belongs to a bubble.
            forAll(mesh_.cellCells()[cellI],nId)
            {
                label nI(mesh_.cellCells()[cellI][nId]);

                if(bubbleIDs_[cellI] < SMALL)
                {
                    //- If the neighbour belongs to a bubble
                    //  add it to that bubble.
                    if(bubbleIDs_[nI] > SMALL)
                    {
                        bubbleIDs_[cellI] = bubbleIDs_[nI];
                    }
                }
                else
                {
                    if(bubbleIDs_[nI] > SMALL)
                    {
                        //- Assign the minimum id
                        bubbleIDs_[cellI] =
                            min(bubbleIDs_[cellI],bubbleIDs_[nI]);

                        //- Update conversion table
                        label maxofTwos(max(bubbleIDs_[cellI],bubbleIDs_[nI]));
                        label maxofTwosId
                        (
                            getBubbleLocId
                            (
                                maxofTwos,
                                minId,
                                maxNofBubbles_
                            )
                        );

                        bubbleIDsConversion[maxofTwosId] = bubbleIDs_[cellI];

                        // label absVal(bubbleIDs_[cellI]);
                        // bubbleIDsConversion[maxofTwosId] =
                        // getBubbleLocId
                        // (
                        //     absVal,
                        //     minId,
                        //     maxNofBubbles_
                        // );

                        bubbleIDs_[nI] = bubbleIDs_[cellI];
                    }
                }
            }

            //- If the bubble ID is still 0, then create a "new id" increasing
            //  over the maximum.
            if(bubbleIDs_[cellI] < SMALL)
            {
                bubbleIDs_[cellI] = maxId;
                bubbleIDsConversion[maxId] = maxId;
                maxId++;
            }
        }
        else
        {
            bubbleIDs_[cellI] = 0;
        }



    }

    //- Update colour field with the new map
    forAll(bubbleIDs_,cellI)
    {
        if(bubbleIDs_[cellI] > SMALL)
        {
            label currId(bubbleIDs_[cellI]);
            bubbleIDs_[cellI] =
                bubbleIDsConversion
                [
                    getBubbleLocId
                    (
                        currId,
                        minId,
                        maxNofBubbles_
                    )
                ];
        }
    }

    bubbleIDs_.correctBoundaryConditions();

    //- Now check processor boundaries to understand if the bubble ID should be
    //  changed.
    label iter(0);
    label res(1);
    while(res>0 && iter < 100)
    {
        iter++;
        res = 0;
        labelList oldTable(bubbleIDsConversion);
        forAll(bubbleIDs_.boundaryField(), patchI)
        {
            //- Check if it is a coupled patch
            if
            (
                isA<coupledFvPatchScalarField>
                (
                    bubbleIDs_.boundaryField()[patchI]
                )
            )
            {
                //- Get patch as coupled
                coupledFvPatchScalarField& bf
                (
                    refCast<coupledFvPatchScalarField>
                    (
                        bubbleIDs_.boundaryFieldRef()[patchI]
                    )
                );

                //- Get interface fields
                scalarField neighField(bf.patchNeighbourField());
                scalarField intField(bf.patchInternalField());

                //- Loop to correct bubble ids
                forAll(intField,faceId)
                {
                    //- Continue loop if no bubble is there
                    if(neighField[faceId] < SMALL || intField[faceId] < SMALL)
                    {
                        continue;
                    }

                    //- Correct IDs mapping if necessary
                    label valueId(intField[faceId]);
                    if(valueId > neighField[faceId])
                    {
                        label maxofTwosId
                        (
                            getBubbleLocId
                            (
                                valueId,
                                minId,
                                maxNofBubbles_
                            )
                        );
                        bubbleIDsConversion[maxofTwosId] = neighField[faceId];
                    }

                }

            }

        }

        //- Compute residual based on the new conversion table
        forAll(oldTable,tabId)
        {
            res += mag(oldTable[tabId] - bubbleIDsConversion[tabId]);
        }

        reduce(res,maxOp<label>());

        if(res == 0)
        {
            Info<<"\nConverged with " << iter << " iterations" << endl;
        }
        else if(iter == 100)
        {
            Info<<"\nReached 100 iterations!\n";
        }

        //- Update colour field with the new map
        forAll(bubbleIDs_,cellI)
        {
            if(bubbleIDs_[cellI] > SMALL)
            {
                label currId(bubbleIDs_[cellI]);
                bubbleIDs_[cellI] =
                    bubbleIDsConversion
                    [
                        getBubbleLocId
                        (
                            currId,
                            bubbleIDsConversion
                        )
                    ];
            }
        }

        bubbleIDs_.correctBoundaryConditions();
    }

    //- Reorganise bubbles to get meaningful ids

    //- First, create a list of new bubbles belonging to this processor
    labelList   newBubbleList;
    for(label bId(0);bId<maxNofBubbles_;bId++)
    {
        if
        (
            bubbleIDsConversion[maxNofBubbles_ + bId]
            ==
            minId + bId
        )
        {
            newBubbleList.append(minId +bId);
        }
    }

    //- Now communicate the list
    labelList displ(UPstream::nProcs(),0);
    labelList numfrags(UPstream::nProcs(),0);
    label   remoteSize = 0;

    numfrags[UPstream::myProcNo()] = newBubbleList.size();

    reduce(numfrags,sumOp<labelList>());

    forAll(displ,proc)
    {
        displ[proc] = remoteSize;
        remoteSize += numfrags[proc];
    }

    label remoteDispl_ = displ[UPstream::myProcNo()];

    labelList globalNewBubbleList(remoteSize,0);

    forAll(newBubbleList,nbId)
    {
        globalNewBubbleList[remoteDispl_+nbId] = newBubbleList[nbId];
    }

    reduce(globalNewBubbleList,maxOp<labelList >());

    //- Now everyone has a list of new bubbles

    //- Update the colour field with the new ids
    forAll(bubbleIDs_,cellI)
    {
        if(bubbleIDs_[cellI] > maxNofBubbles_ - 1)
        {
            label currId(-1);
            forAll(globalNewBubbleList,bId)
            {
                if(globalNewBubbleList[bId] == bubbleIDs_[cellI])
                {
                    currId = bId;
                }
            }
            bubbleIDs_[cellI] = bubbleList_.size() + currId;
        }
    }

    //-  Add new bubbles in the same order.
    forAll(globalNewBubbleList,bId)
    {
        bubbleList_.append(bubble());
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::bubbleTracker::bubbleTracker
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    runTime_(runTime),
    alphaName_(dict.lookup("alphaName")),
    alphaTol_(dict.lookupOrDefault("alphaTol",1e-1)),
    maxNofBubbles_(dict.lookupOrDefault("maxNofBubblesPerProc",20)),
    bubbleIDs_
    (
        IOobject
        (
            "bubbleIDs",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("bubbleIDs",dimless,0)
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::bubbleTracker::read(const dictionary& dict)
{

    return true;
}


bool Foam::functionObjects::bubbleTracker::execute()
{

    updateColorField();

    return true;
}


bool Foam::functionObjects::bubbleTracker::end()
{
    return true;
}


bool Foam::functionObjects::bubbleTracker::write()
{
    return true;
}

// ************************************************************************* //
