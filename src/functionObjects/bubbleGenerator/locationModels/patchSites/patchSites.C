/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "patchSites.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace locationModels
{
    defineTypeNameAndDebug(patchSites, 0);
    addToRunTimeSelectionTable(locationModel, patchSites, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::locationModels::patchSites::patchSites
(
    const dictionary dict,
    const fvMesh& mesh
)
:
    locationModel(dict, mesh),
    sites_(),
    patchName_(dict.lookup("patchName"))
{
    label patchId = mesh.boundaryMesh().findPatchID(patchName_);
    
    if(patchId != -1)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchId];

        const labelUList& pcells = pp.faceCells();

        sites_.resize(pcells.size());

        forAll(pcells,celli)
        {
            sites_[celli] = mesh.C()[pcells[celli]];
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::locationModels::patchSites::~patchSites()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField&
Foam::locationModels::patchSites::sites() const
{

    //- NOTE: sites will always correspond to the initial mesh!

    return sites_;
}


// ************************************************************************* //
