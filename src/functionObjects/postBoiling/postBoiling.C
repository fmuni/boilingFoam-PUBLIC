/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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

#include "fvCFD.H"
#include "postBoiling.H"
#include "error.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include <iomanip>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(postBoiling, 0);
    addToRunTimeSelectionTable(functionObject, postBoiling, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::postBoiling::postBoiling
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    alphaName_(dict.get<word>("alphaName")),
    patches_(dict.get<wordList>("patches")),
    Tsat_(dict.get<scalar>("Tsat")),
    Lr_(dict.get<scalar>("Lr")),
    Nu_
    (
        IOobject
        (
            "Nu",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Nu",dimless,0)
    ),
    splitPatch_(dict.getOrDefault<Switch>("splitPatch",false)),
    splitVec_(vector::zero)
{
    read(dict);

    if(splitPatch_)
    {
        splitVec_ = vector(dict.get<vector>("splitVector"));
    }

    //- Create output file for csv
    if (Pstream::master())
    {
        std::ofstream file;
        file.open ("boilingPost.csv");
        if(file.is_open())
        {
            //- Write header
            file<< "Time, "
                << "volume of vapour";

            forAll(patches_,patchi)
            {
                word patchName(patches_[patchi]);
                file<< ", " << patchName << " Area"
                    << ", " << patchName << " Nu"
                    << ", " << patchName << " qw"
                    << ", " << patchName << " T"
                    << ", " << patchName << " dry to total";

                if(splitPatch_)
                {
                    file<< ", " << patchName << "2 Area"
                        << ", " << patchName << "2 Nu"
                        << ", " << patchName << "2 qw"
                        << ", " << patchName << "2 T"
                        << ", " << patchName << "2 dry to total";
                }
            }

            file<< "\n";
        }
        else
        {
            FatalErrorInFunction<< "\n Unable to open file boilingPost.csv"
                                << abort(FatalError);
        }

        file.close();
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::postBoiling::read(const dictionary& dict)
{

    return true;
}


bool Foam::functionObjects::postBoiling::execute()
{
    //- Prepare file for master
    std::ofstream file;
    if (Pstream::master())
    {

        file.open
        (
            "boilingPost.csv",
            std::ofstream::out | std::ofstream::app
        );

        file.setf(std::ios::scientific);
        file.precision(8);

        if(file.is_open())
        {
            file<< mesh_.time().value();
        }
        else
        {
            FatalErrorInFunction<< "\n Unable to open file boilingPost.csv"
                                << abort(FatalError);
        }

    }

    //- Compute vapour volume
    const volScalarField& alphal =
        mesh_.lookupObject<volScalarField>(alphaName_);

    dimensionedScalar bubbleV
    (
        fvc::domainIntegrate(1.0 - alphal)
    );

    if (Pstream::master())
    {
        file<< ", " << bubbleV.value();
    }

    //- Compute wall heat flux
    const volScalarField& kappa =
        mesh_.lookupObject<volScalarField>("kappa");

    const volScalarField& T =
        mesh_.lookupObject<volScalarField>("T");


    forAll(patches_,patchi)
    {
        scalar Tave(0), Nuave(0), dry(0), qw(0), totS(SMALL);
        //- More for splitPatch
        scalar Tave2(0), Nuave2(0), dry2(0), qw2(0), totS2(SMALL);
        label patchId = mesh_.boundaryMesh().findPatchID(patches_[patchi]);

        //- Check that at least one processor has the patch
        label patchIdglobal(patchId);

        reduce(patchIdglobal,maxOp<label>());
        if(patchIdglobal == -1)
        {
            FatalErrorInFunction<< "\n Unable to open file patch "
                                << patches_[patchi] << endl
                                << abort(FatalError);
        }
        //- Assume someone does not have this patch.
        //  Only processors with this patch should perform the computations.
        if(patchId > -1)
        {
            scalarField deltaT(Tsat_ - T.boundaryField()[patchId]);

            //- Correct deltaT to be always != 0
            forAll(deltaT,fc)
            {
                if(mag(deltaT[fc]) < SMALL)
                {
                    deltaT[fc] = SMALL;
                }
            }

            //- Correct local patch Nusselt
            Nu_.boundaryFieldRef()[patchId] =
            (
                -T.boundaryField()[patchId].snGrad()*Lr_/deltaT
            );

            //- Now compute the integral quantities
            scalarField magS(mag(mesh_.boundaryMesh()[patchId].faceAreas()));
            vectorField Sc(mesh_.boundaryMesh()[patchId].faceCentres());

            forAll(deltaT,fc)
            {
                if
                (
                    splitPatch_
                    &&
                    ((Sc[fc]&splitVec_) > magSqr(splitVec_))
                )
                {
                    totS2 += magS[fc];
                    Tave2 += (Tsat_ - deltaT[fc])*magS[fc];
                    Nuave2 += Nu_.boundaryField()[patchId][fc]*magS[fc];
                    qw2 +=
                        Nu_.boundaryField()[patchId][fc]*magS[fc]*deltaT[fc]/Lr_
                        *kappa.boundaryField()[patchId][fc];
                    dry2 +=
                        (1.0 - alphal.boundaryField()[patchId][fc])*magS[fc];

                }
                else
                {
                    totS += magS[fc];
                    Tave += (Tsat_ - deltaT[fc])*magS[fc];
                    Nuave += Nu_.boundaryField()[patchId][fc]*magS[fc];
                    qw +=
                        Nu_.boundaryField()[patchId][fc]*magS[fc]*deltaT[fc]/Lr_
                        *kappa.boundaryField()[patchId][fc];
                    dry +=
                        (1.0 - alphal.boundaryField()[patchId][fc])*magS[fc];

                }

            }
        }

        reduce(totS,sumOp<scalar>());
        reduce(Tave,sumOp<scalar>());
        reduce(qw,sumOp<scalar>());
        reduce(dry,sumOp<scalar>());
        reduce(Nuave,sumOp<scalar>());

        Tave /= totS;
        qw /= totS;
        dry /= totS;
        Nuave /= totS;

        if(splitPatch_)
        {
            reduce(totS2,sumOp<scalar>());
            reduce(Tave2,sumOp<scalar>());
            reduce(qw2,sumOp<scalar>());
            reduce(dry2,sumOp<scalar>());
            reduce(Nuave2,sumOp<scalar>());

            Tave2 /= totS2;
            qw2 /= totS2;
            dry2 /= totS2;
            Nuave2 /= totS2;
        }

        //- Finally write to file
        if (Pstream::master())
        {

            file<<", " << totS
                <<", " << Nuave
                <<", " << qw
                <<", " << Tave
                <<", " << dry;

            if(splitPatch_)
            {
                file<<", " << totS2
                    <<", " << Nuave2
                    <<", " << qw2
                    <<", " << Tave2
                    <<", " << dry2;
            }
        }

    }

    //- New line and close file
    if (Pstream::master())
    {

        file<<"\n";
        file.close();
    }

    return true;
}


bool Foam::functionObjects::postBoiling::end()
{
    return true;
}


bool Foam::functionObjects::postBoiling::write()
{
    return true;
}


// ************************************************************************* //
