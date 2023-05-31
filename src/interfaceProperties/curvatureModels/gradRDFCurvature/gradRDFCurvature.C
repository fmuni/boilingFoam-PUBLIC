/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

#include "gradRDFCurvature.H"
#include "addToRunTimeSelectionTable.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "surfaceFields.H"
#include "surfaceMesh.H"
#include "unitConversion.H"
#include "interpolationCellPoint.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "Pstream.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace curvatureModels
{
    defineTypeNameAndDebug(gradRDFCurvature, 0);
    addToRunTimeSelectionTable(curvatureModel, gradRDFCurvature, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::curvatureModels::gradRDFCurvature::gradRDFCurvature
(
    const dictionary& dict,
    const fvMesh& mesh,
    reconstructionSchemes& surf
)
:
    curvatureModel(dict,mesh,surf)
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::curvatureModels::gradRDFCurvature::~gradRDFCurvature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::curvatureModels::gradRDFCurvature::calculateK
(
    const volScalarField& alpha1,
    const volVectorField& U
)
{
    const volScalarField& RDF(mesh_.lookupObject<volScalarField>("RDF"));
    const surfaceVectorField& Sf = mesh_.Sf();


    //- Get centres and normals
    const volVectorField& centres(surf().centre());

    const volVectorField& normals(surf().normal());


    volVectorField gradRDF(fvc::grad(RDF));

    surfaceVectorField gradRDFf(fvc::interpolate(gradRDF));

    // Face unit interface normal
    surfaceVectorField nHatfv
    (
        gradRDFf
        /
        ( mag(gradRDFf) + deltaNa() )
    );

    correctContactAngle
    (
        nHatfv.boundaryFieldRef(),
        gradRDFf.boundaryField(),
        U,
        alpha1
    );


    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;
    // correctContactAngle
    // (
    //     nHatfv.boundaryFieldRef(),
    //     gradRDFf.boundaryField(),
    //     U,
    //     RDF
    // );

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);

    //- Impose max/min based on the values at the interface
    scalar maxK(-GREAT);
    scalar minK(GREAT);

    forAll(K_,cellI)
    {
        //- Cells at the interface
        if (mag(normals[cellI]) > SMALL)
        {
            maxK = max(maxK,K_[cellI]);
            minK = min(minK,K_[cellI]);
        }
    }

    reduce(maxK,maxOp<scalar>());
    reduce(minK,minOp<scalar>());

    K_ = min(K_,dimensionedScalar("maxK",K_.dimensions(),maxK));
    K_ = max(K_,dimensionedScalar("minK",K_.dimensions(),minK));

    //- Now that a first curvature is computed, interpolate at the
    //  interface.
/*
    //- Interpolate at the interface
    interpolationCellPoint<scalar> interpolK(K_);

    scalarField KatInterFace;
    vectorField centresList;

    forAll(K_,cellI)
    {
           //- Cells at the interface
           if (mag(normals[cellI]) > SMALL)
           {
               //- Substitute the interpolated value
               K_[cellI]=  interpolK.interpolate(centres[cellI],cellI,-1);

               //- Store the interpolated value
               KatInterFace.append(K_[cellI]);

               //- Also store the centre
               centresList.append(centres[cellI]);
           }
           else
           {
                //- Set to zero if not at the interface
               K_[cellI]= 0;
           }
    }

    //- Create global list of centres and curvature values
    //  Then, communicate to all processors as global data.
    labelList nInterPoints(Pstream::nProcs(),0);
    nInterPoints[Pstream::myProcNo()] = centresList.size();

    reduce(nInterPoints,sumOp<labelList>());

    labelList stride(Pstream::nProcs(),0);
    label totInterPoints(0);

    forAll(nInterPoints,proc)
    {
        stride[proc] = totInterPoints;
        totInterPoints += nInterPoints[proc];
    }

    vectorField globalCentres(totInterPoints,vector::zero);
    scalarField globalCurvs(totInterPoints,0);

    label myStride(stride[Pstream::myProcNo()]);
    forAll(centresList,intId)
    {
        globalCentres[myStride + intId] = centresList[intId];
        globalCurvs[myStride + intId] = KatInterFace[intId];
    }

    reduce(globalCentres,sumOp<vectorField>());
    reduce(globalCurvs,sumOp<scalarField>());

    //- This is for checking the close neighbors

    Random rndGen(17301893);

    // Slightly boiling bb. Slightly off-centred just so on symmetric
    // geometry there are less face/edge aligned items.
    treeBoundBox bb
    (
        treeBoundBox(globalCentres).extend(rndGen, 1e-4)
    );

    // bb.min() -= point(1e-8, 1e-8, 1e-8);
    // bb.max() += point(1e-8, 1e-8, 1e-8);

    indexedOctree<treeDataPoint> surfaceTree
    (
        treeDataPoint
        (
            globalCentres
        ),
        bb,     // bb
        8,      // maxLevel
        10,     // leafsize
        3.0     // duplicity
    );

    forAll(K_,cellI)
    {
        //- Cells NOT at the interface but in the region where the RDF is
        //  defined
        if (mag(normals[cellI]) < SMALL && mag(RDF[cellI]) > SMALL)
        {
            pointIndexHit pHit =
                surfaceTree.findNearest(mesh_.C()[cellI], GREAT);

            //- Get the id of the corresponding hit node
            const label idx = pHit.index();

            //- If cell is found
            if (idx != -1)
            {
                K_[cellI]=  globalCurvs[idx];
            }
        }
    }
*/
    K_.correctBoundaryConditions();

}

bool Foam::curvatureModels::gradRDFCurvature::readDict(const dictionary& dict)
{
    return true;
}


bool Foam::curvatureModels::gradRDFCurvature::writeData(Ostream& os) const
{

    return false;
}


// ************************************************************************* //
