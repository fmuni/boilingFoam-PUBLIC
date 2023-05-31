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

#include "bubbleGenerator.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(bubbleGenerator, 0);
    addToRunTimeSelectionTable(functionObject, bubbleGenerator, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::bubbleGenerator::bubbleGenerator
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    runTime_(runTime),
    fixDT_(dict.lookupOrDefault<bool>("fixTimeStepAfterNucleation",false))
{
    read(dict);

    location_.reset
    (
        locationModel::New
        (
            dict,
            mesh_
        )
    );

    activation_.reset
    (
        activationModel::New
        (
            dict,
            mesh_
        )
    );

    nucleation_.reset
    (
        nucleationModel::New
        (
            dict,
            mesh_
        )
    );

    if(fixDT_)
    {
        newDT_ = readScalar(dict.lookup("deltaT"));
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::bubbleGenerator::read(const dictionary& dict)
{

    return true;
}


bool Foam::functionObjects::bubbleGenerator::execute()
{

    bool nucleated(false);

    const vectorField& sites = location_->sites();

    forAll(sites,siteId)
    {
        if(activation_->isActive(sites[siteId]))
        {
            nucleated =
            (
                nucleated
                ||
                nucleation_->nucleate(sites[siteId])
            );
        }
    }

    adjustDT(nucleated & fixDT_);

    return true;
}


bool Foam::functionObjects::bubbleGenerator::end()
{
    return true;
}


bool Foam::functionObjects::bubbleGenerator::write()
{
    return true;
}

void Foam::functionObjects::bubbleGenerator::adjustDT(bool adjust) const
{
    Time& runTime = const_cast<Time&>(runTime_);

    if(adjust)
    {
        Info<< "BubbleGenerator : Setting deltaT to " << newDT_
            << " s" << endl;
        runTime.setDeltaT(newDT_);
    }
}

// ************************************************************************* //
