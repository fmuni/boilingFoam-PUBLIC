/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 DLR
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

#include "plicRDF.H"
#include "interpolationCellPoint.H"
#include "fvc.H"
#include "leastSquareGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace reconstruction
{
    defineTypeNameAndDebug(plicRDF, 0);
    addToRunTimeSelectionTable(reconstructionSchemes,plicRDF, components);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::reconstruction::plicRDF::interpolateNormal()
{
    scalar dt = mesh_.time().deltaTValue();

    leastSquareGrad<scalar> lsGrad("polyDegree1",mesh_.geometricD());

    exchangeFields_.setUpCommforZone(interfaceCell_,false);

    Map<vector> mapCentre
    (
        exchangeFields_.getDatafromOtherProc(interfaceCell_, centre_)
    );
    Map<vector> mapNormal
    (
        exchangeFields_.getDatafromOtherProc(interfaceCell_, normal_)
    );

    Map<vector> mapCC
    (
        exchangeFields_.getDatafromOtherProc(interfaceCell_, mesh_.C())
    );
    Map<scalar> mapAlpha
    (
        exchangeFields_.getDatafromOtherProc(interfaceCell_, alpha1_)
    );

    DynamicField<vector > cellCentre(100);
    DynamicField<scalar > alphaValues(100);

    DynamicList<vector> foundNormals(30);

    const labelListList& stencil = exchangeFields_.getStencil();

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        vector estimatedNormal{Zero};
        scalar weight{0};
        foundNormals.clear();
        forAll(stencil[celli], i)
        {
            const label gblIdx = stencil[celli][i];
            vector n =
                exchangeFields_.getValue(normal_, mapNormal, gblIdx);
            point p = mesh_.C()[celli]-U_[celli]*dt;
            if (mag(n) != 0)
            {
                n /= mag(n);
                vector centre =
                    exchangeFields_.getValue(centre_, mapCentre, gblIdx);
                vector distanceToIntSeg = (tensor::I- n*n) & (p - centre);
                estimatedNormal += n /max(mag(distanceToIntSeg), SMALL);
                weight += 1/max(mag(distanceToIntSeg), SMALL);
                foundNormals.append(n);
            }
        }

        if (weight != 0 && mag(estimatedNormal) != 0)
        {
            estimatedNormal /= weight;
            estimatedNormal /= mag(estimatedNormal);
        }

        bool tooCoarse = false;

        if (foundNormals.size() > 1 && mag(estimatedNormal) != 0)
        {
            forAll(foundNormals, i)
            {
                // all have the length of 1
                // to coarse if normal angle is bigger than 10 deg
                if ((estimatedNormal & foundNormals[i]) <= 0.98)
                {
                    tooCoarse = true;
                }
            }
        }
        else
        {
            tooCoarse = true;
        }

        // if a normal was found and the interface is fine enough
        // smallDist is always smallDist
        if (mag(estimatedNormal) != 0 && !tooCoarse)
        {
            interfaceNormal_[i] = estimatedNormal;
        }
        else
        {
            cellCentre.clear();
            alphaValues.clear();

            forAll(stencil[celli],i)
            {
                const label gblIdx = stencil[celli][i];
                cellCentre.append
                (
                    exchangeFields_.getValue(mesh_.C(), mapCC, gblIdx)
                );
                alphaValues.append
                (
                    exchangeFields_.getValue(alpha1_, mapAlpha, gblIdx)
                );
            }
            cellCentre -= mesh_.C()[celli];
            interfaceNormal_[i] = lsGrad.grad(cellCentre, alphaValues);
        }

    }
}

void Foam::reconstruction::plicRDF::gradSurf(const volScalarField& phi)
{
    leastSquareGrad<scalar> lsGrad("polyDegree1", mesh_.geometricD());

    exchangeFields_.setUpCommforZone(interfaceCell_, false);

    Map<vector> mapCC
    (
        exchangeFields_.getDatafromOtherProc(interfaceCell_, mesh_.C())
    );
    Map<scalar> mapPhi
    (
        exchangeFields_.getDatafromOtherProc(interfaceCell_, phi)
    );

    DynamicField<vector> cellCentre(100);
    DynamicField<scalar> phiValues(100);

    const labelListList& stencil = exchangeFields_.getStencil();

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];

        cellCentre.clear();
        phiValues.clear();

        for (const label gblIdx : stencil[celli])
        {
            cellCentre.append
            (
                exchangeFields_.getValue(mesh_.C(), mapCC, gblIdx)
            );
            phiValues.append
            (
                exchangeFields_.getValue(phi, mapPhi, gblIdx)
            );
        }

        cellCentre -= mesh_.C()[celli];
        interfaceNormal_[i] = lsGrad.grad(cellCentre, phiValues);
    }
}


void Foam::reconstruction::plicRDF::setInitNormals(bool interpolate)
{
    interfaceLabels_.clear();

    forAll(alpha1_, celli)
    {
        if (sIterPLIC_.isASurfaceCell(alpha1_[celli]))
        {
            interfaceCell_[celli] = true; // is set to false earlier
            interfaceLabels_.append(celli);
        }
    }
    interfaceNormal_.setSize(interfaceLabels_.size());

    RDF_.markCellsNearSurf(interfaceCell_, RDFStencil_);
    const boolList& nextToInterface_ = RDF_.nextToInterface();
    exchangeFields_.updateStencil(nextToInterface_);

    if (interpolate)
    {
        interpolateNormal();
    }
    else
    {
        gradSurf(alpha1_);
    }
}


void Foam::reconstruction::plicRDF::calcResidual
(
    Map<scalar>& normalResidual,
    Map<scalar>& avgAngle
)
{
    exchangeFields_.setUpCommforZone(interfaceCell_,false);

    Map<vector> mapNormal
    (
        exchangeFields_.getDatafromOtherProc(interfaceCell_, normal_)
    );

    const labelListList& stencil = exchangeFields_.getStencil();

    normalResidual.clear();

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        if (mag(normal_[celli]) == 0 || mag(interfaceNormal_[i]) == 0)
        {
            continue;
        }

        scalar avgDiffNormal = 0;
        scalar maxDiffNormal = GREAT;
        scalar weight= 0;
        const vector cellNormal = normal_[celli]/mag(normal_[celli]);

        for (const label gblIdx : stencil[celli])
        {
            vector normal =
                exchangeFields_.getValue(normal_, mapNormal, gblIdx);

            if (mag(normal) != 0 && i != 0)
            {
                vector n = normal/mag(normal);
                scalar cosAngle = max(min((cellNormal & n), 1), -1);
                avgDiffNormal += acos(cosAngle) * mag(normal);
                weight += mag(normal);
                if (cosAngle < maxDiffNormal)
                {
                    maxDiffNormal = cosAngle;
                }
            }
        }

        if (weight != 0)
        {
            avgDiffNormal /= weight;
        }
        else
        {
            avgDiffNormal = 0;
        }

        vector newCellNormal = normalised(interfaceNormal_[i]);

        scalar normalRes = (1 - (cellNormal & newCellNormal));
        avgAngle.insert(celli, avgDiffNormal);
        normalResidual.insert(celli, normalRes);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstruction::plicRDF::plicRDF
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    reconstructionSchemes
    (
        typeName,
        alpha1,
        phi,
        U,
        dict
    ),
    mesh_(alpha1.mesh()),

    interfaceNormal_(0.2*mesh_.nCells()),

    isoFaceTol_(modelDict().getOrDefault<scalar>("isoFaceTol", 1e-8)),
    surfCellTol_(modelDict().getOrDefault<scalar>("surfCellTol", 1e-8)),
    tol_(modelDict().getOrDefault<scalar>("tol", 1e-6)),
    relTol_(modelDict().getOrDefault<scalar>("relTol", 0.1)),
    iteration_(modelDict().getOrDefault<label>("iterations", 5)),
    interpolateNormal_(modelDict().getOrDefault("interpolateNormal", true)),
    RDF_(mesh_),
    exchangeFields_(zoneDistribute::New(mesh_)),
    sIterPLIC_(mesh_,surfCellTol_),
    RDFStencil_(readLabel(modelDict().lookup("RDFStencil")))
{
    setInitNormals(false);

    centre_ = dimensionedVector("centre", dimLength, Zero);
    normal_ = dimensionedVector("normal", dimArea, Zero);
    
    Info<<"\nBoiling RDF stencil cells = " << RDFStencil_<< endl;

    forAll(interfaceLabels_, i)
    {
        const label celli = interfaceLabels_[i];
        if (mag(interfaceNormal_[i]) == 0)
        {
            continue;
        }
        sIterPLIC_.vofCutCell
        (
            celli,
            alpha1_[celli],
            isoFaceTol_,
            100,
            interfaceNormal_[i]
        );

        if (sIterPLIC_.cellStatus() == 0)
        {
            normal_[celli] = sIterPLIC_.surfaceArea();
            centre_[celli] = sIterPLIC_.surfaceCentre();
            if (mag(normal_[celli]) == 0)
            {
                normal_[celli] = Zero;
                centre_[celli] = Zero;
            }
        }
        else
        {
            normal_[celli] = Zero;
            centre_[celli] = Zero;
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstruction::plicRDF::reconstruct(bool forceUpdate)
{
    const bool uptodate = alreadyReconstructed(forceUpdate);

    if (uptodate && !forceUpdate)
    {
        return;
    }

    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if (interfaceCell_.size() != mesh_.nCells())
        {
            interfaceCell_.resize(mesh_.nCells());
        }
    }
    interfaceCell_ = false;

    // Sets interfaceCell_ and interfaceNormal
    setInitNormals(interpolateNormal_);

    centre_ = dimensionedVector("centre", dimLength, Zero);
    normal_ = dimensionedVector("normal", dimArea, Zero);

    // nextToInterface is update on setInitNormals
    const boolList& nextToInterface_ = RDF_.nextToInterface();

    labelHashSet tooCoarse;

    for (int iter=0; iter<iteration_; ++iter)
    {
        forAll(interfaceLabels_, i)
        {
            const label celli = interfaceLabels_[i];
            if (mag(interfaceNormal_[i]) == 0 || tooCoarse.found(celli))
            {
                continue;
            }
            sIterPLIC_.vofCutCell
            (
                celli,
                alpha1_[celli],
                isoFaceTol_,
                100,
                interfaceNormal_[i]
            );

            if (sIterPLIC_.cellStatus() == 0)
            {

                normal_[celli] = sIterPLIC_.surfaceArea();
                centre_[celli] = sIterPLIC_.surfaceCentre();
                if (mag(normal_[celli]) == 0)
                {
                    normal_[celli] = Zero;
                    centre_[celli] = Zero;
                }
            }
            else
            {
                normal_[celli] = Zero;
                centre_[celli] = Zero;
            }
        }

        normal_.correctBoundaryConditions();
        centre_.correctBoundaryConditions();
        Map<scalar> residual;
        Map<scalar> avgAngle;

        surfaceVectorField::Boundary nHatb(mesh_.Sf().boundaryField());
        nHatb *= 1/(mesh_.magSf().boundaryField());

        {
            RDF_.constructRDF
            (
                nextToInterface_,
                centre_,
                normal_,
                exchangeFields_,
                false
            );
            RDF_.updateContactAngle(alpha1_, U_, nHatb);
            gradSurf(RDF_);
            calcResidual(residual, avgAngle);
        }


        label resCounter = 0;
        scalar avgRes = 0;
        scalar avgNormRes = 0;

        Map<scalar>::iterator resIter = residual.begin();
        Map<scalar>::iterator avgAngleIter = avgAngle.begin();

        while (resIter.found())
        {
            if (avgAngleIter() > 0.26 && iter > 0) // 15 deg
            {
                tooCoarse.set(resIter.key());
            }
            else
            {
                avgRes += resIter();
                scalar normRes = 0;
                scalar discreteError = 0.01*sqr(avgAngleIter());
                if (discreteError != 0)
                {
                    normRes= resIter()/max(discreteError, tol_);
                }
                else
                {
                    normRes= resIter()/tol_;
                }
                avgNormRes += normRes;
                resCounter++;

            }

            ++resIter;
            ++avgAngleIter;
        }

        reduce(avgRes,sumOp<scalar>());
        reduce(avgNormRes,sumOp<scalar>());
        reduce(resCounter,sumOp<label>());

        if (resCounter == 0) // avoid division  by zero and leave loop
        {
            resCounter = 1;
            avgRes = 0;
            avgNormRes = 0;
        }


        if (iter == 0)
        {
            DebugInfo
                << "initial residual absolute = "
                << avgRes/resCounter << nl
                << "initial residual normalized = "
                << avgNormRes/resCounter << nl;
        }

        if
        (
            (
                (avgNormRes/resCounter < relTol_ || avgRes/resCounter < tol_)
             && (iter >= 1 )
            )
         || iter + 1  == iteration_
        )
        {
            DebugInfo
                << "iterations = " << iter << nl
                << "final residual absolute = "
                << avgRes/resCounter << nl
                << "final residual normalized = " << avgNormRes/resCounter
                << endl;

            break;
        }
    }
}


void Foam::reconstruction::plicRDF::mapAlphaField() const
{
    // without it we seem to get a race condition
    mesh_.C();

    cutCellPLIC cutCell(mesh_);

    forAll(normal_, celli)
    {
        if (mag(normal_[celli]) != 0)
        {
            vector n = normal_[celli]/mag(normal_[celli]);
            scalar cutValue = (centre_[celli] - mesh_.C()[celli]) & (n);
            cutCell.calcSubCell
            (
                celli,
                cutValue,
                n
            );
            alpha1_[celli] = cutCell.VolumeOfFluid();

        }
    }
    alpha1_.correctBoundaryConditions();
    alpha1_.oldTime () = alpha1_;
    alpha1_.oldTime().correctBoundaryConditions();
}


// ************************************************************************* //
