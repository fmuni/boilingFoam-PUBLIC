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

#include "bubbleDestroyer.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"
#include "interpolation.H"
#include "surfaceFields.H"




#define DEBUG false

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(bubbleDestroyer, 0);
    addToRunTimeSelectionTable(functionObject, bubbleDestroyer, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::bubbleDestroyer::bubbleDestroyer
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    numberOfProbes_(dict.get<label>("numberOfProbes")),
    alphaName_(dict.get<word>("alphaName")),
    direction_(dict.get<vector>("direction")),
    offset_(dict.getOrDefault<vector>("offset",vector::zero)),
    tol_(dict.getOrDefault<scalar>("tolerance", 1e-3)),
    probeList_(numberOfProbes_,vector::zero),
    removeMomentum_(dict.getOrDefault<Switch>("removeMomentum",false)),
    resetAlpha_(dict.getOrDefault<Switch>("resetAlpha",false)),
    resetU_(dict.getOrDefault<Switch>("resetU",false)),
    initialT_(dict.getOrDefault<word>("initialTime","0")),
    alphaReinitDelay_(dict.getOrDefault<scalar>("alphaReinitDelay",-1)),
    counterT_(0),
    alphaMustBeRestored_(false)
{
    Info<<"\nUsing functionObject: bubbleDestroyer " << endl;

    read(dict);

    //- Normalize the direction
    direction_ /= mag(direction_);

    //- Get global mesh bounding box
    const boundBox bounds(mesh_.bounds());
    vector bMax(bounds.max() + offset_);
    vector bMin(bounds.min() + offset_);

    //- Get the maximum and minimum bounds
    scalar maxB(bMax&direction_);
    scalar minB(bMin&direction_);

    //- Get the channel centre
    vector maxCC(bMax - (maxB*direction_));
    vector minCC(bMin - (minB*direction_));

    scalar maxCCx(maxCC.x());
    scalar maxCCy(maxCC.y());
    scalar maxCCz(maxCC.z());
    scalar minCCx(minCC.x());
    scalar minCCy(minCC.y());
    scalar minCCz(minCC.z());

    reduce(maxB,maxOp<scalar>());
    reduce(maxCCx,maxOp<scalar>());
    reduce(maxCCy,maxOp<scalar>());
    reduce(maxCCz,maxOp<scalar>());

    reduce(minB,minOp<scalar>());
    reduce(minCCx,minOp<scalar>());
    reduce(minCCy,minOp<scalar>());
    reduce(minCCz,minOp<scalar>());

    vector channelCC
    (
        (maxCCx + minCCx)/2.0,
        (maxCCy + minCCy)/2.0,
        (maxCCz + minCCz)/2.0
    );

    if(DEBUG)
    {
        Info<<"\nmaxB = " << maxB << "  "
            <<"minB = " << minB << " "
            <<"channelCC = " << channelCC
            << endl;
    }

    scalar delta((maxB - minB)/numberOfProbes_);

    //- Create list of probes
    forAll(probeList_,prb)
    {
        probeList_[prb] = (minB + delta*prb)*direction_ + channelCC;
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::bubbleDestroyer::read(const dictionary& dict)
{
    return true;
}

void Foam::functionObjects::bubbleDestroyer::resetAlphaDelayed()
{
    if(!resetAlpha_ || !alphaMustBeRestored_)
    {
        return;
    }

    counterT_ += mesh_.time().deltaT().value();

    if(counterT_ > alphaReinitDelay_)
    {
        Info << "\nRestoring alpha" <<endl;
        volScalarField& alpha
        (
            const_cast<volScalarField&>
            (
                mesh_.lookupObject<volScalarField>(alphaName_)
            )
        );

        volScalarField alpha0
        (
            IOobject
            (
                alpha.name(),
                initialT_,
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );

        alpha = alpha0;
        alpha.correctBoundaryConditions();

        alphaMustBeRestored_ = false;
    }
}


bool Foam::functionObjects::bubbleDestroyer::execute()
{

    resetAlphaDelayed();

    //- Create a list of cells corresponding to the interpolation points
    labelList cellList(numberOfProbes_,-1);

    forAll(probeList_,prb)
    {
        cellList[prb] = mesh_.findCell(probeList_[prb]);
    }

    labelList cellListGlobal(cellList);

    reduce(cellListGlobal,maxOp<labelList>());

    if(DEBUG)
    {
        forAll(cellListGlobal,cellI)
        {
            if(cellListGlobal[cellI]==-1)
            {
                Info<<"bubbleDestroyer: cannot find probe "
                    << cellI << " corresponding to point "
                    << probeList_[cellI]
                    << endl;
            }
        }
    }

    //- Now interpolate the alpha field.
    //  Const_cast is required because we might need to modify the field later.
    volScalarField& alpha
    (
        const_cast<volScalarField&>
        (
            mesh_.lookupObject<volScalarField>(alphaName_)
        )
    );

    //- Get reference to velocity in case needs to be reset
    volVectorField& U
    (
        const_cast<volVectorField&>
        (
            mesh_.lookupObject<volVectorField>("U")
        )
    );


    autoPtr<interpolation<scalar>> interpolator
    (
        interpolation<scalar>::New("cellPointFace", alpha)
    );

    //- Interpolate at the last probe to check if the bubble should be removed
    label lastValidPrb(-1);
    for(label cellI=numberOfProbes_-1;cellI>0;cellI--)
    {
        if(cellListGlobal[cellI] != -1)
        {
            lastValidPrb = cellI;
            break;
        }
    }

    if(lastValidPrb == -1)
    {
        FatalErrorInFunction<<"No valid last probe in bubbleDestroyer!"
                            << abort(FatalError);
    }

    scalar alphaEnd(2.0);
    if(cellList[lastValidPrb] > -1)
    {
        alphaEnd = interpolator->interpolate
        (
            probeList_[lastValidPrb],
            cellList[lastValidPrb]
        );
    }

    reduce(alphaEnd,minOp<scalar>());

    //- End if the last probe has only liquid
    if(alphaEnd > 1.0-tol_)
    {
        return true;
    }

    Info << "Found bubble at the end of the domain" << endl;

    if(alphaReinitDelay_ > SMALL)
    {
        alphaMustBeRestored_ = true;
    }

    if(resetU_)
    {
        Info << "\nRestoring U" <<endl;
        volVectorField U0
        (
            IOobject
            (
                U.name(),
                initialT_,
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );

        U = U0;
        U.correctBoundaryConditions();
    }

    if(resetAlpha_ && alphaReinitDelay_ < SMALL)
    {
        Info << "\nRestoring alpha" <<endl;
        volScalarField alpha0
        (
            IOobject
            (
                alpha.name(),
                initialT_,
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );

        alpha = alpha0;
        alpha.correctBoundaryConditions();

        return true;
    }

    if(resetAlpha_)
    {
        return true;
    }

    //- Interpolate at the other probes
    scalarField intrpAlpha(numberOfProbes_-1,2.0);

    forAll(cellList,prb)
    {
        if(prb == numberOfProbes_ -1)
        {
            break;
        }

        if(cellList[prb] > -1)
        {
            intrpAlpha[prb] = interpolator->interpolate
            (
                probeList_[prb],
                cellList[prb]
            );
        }
    }

    reduce(intrpAlpha,minOp<scalarField>());

    //- Now set alpha to 1 for all cells after the last liquid
    scalar lastLiq(0);

    forAll(intrpAlpha,prb)
    {
        if(intrpAlpha[prb]<1.0-tol_)
        {
            break;
        }

        lastLiq = probeList_[prb]&direction_;
    }


    Info << "\nRemoving bubbles after L = " << lastLiq << endl;

    forAll(mesh_.C(),cellI)
    {
        if((mesh_.C()[cellI]&direction_) > lastLiq)
        {
            if(removeMomentum_)
            {
                if(alpha[cellI]<1.0-tol_)
                {
                    U[cellI] *= 0.;
                }
            }

            alpha[cellI] = 1.0;
        }
    }

    //- Also correct patches
    volScalarField::Boundary& aBf = alpha.boundaryFieldRef();
    const surfaceVectorField::Boundary& cfBf = mesh_.Cf().boundaryField();

    forAll(aBf,patchi)
    {
        forAll(aBf[patchi],facei)
        {
            if((cfBf[patchi][facei]&direction_)>lastLiq)
            {
                aBf[patchi][facei] = 1.0;
            }
        }
    }

    return true;
}


bool Foam::functionObjects::bubbleDestroyer::end()
{
    return true;
}


bool Foam::functionObjects::bubbleDestroyer::write()
{
    return true;
}


// ************************************************************************* //
