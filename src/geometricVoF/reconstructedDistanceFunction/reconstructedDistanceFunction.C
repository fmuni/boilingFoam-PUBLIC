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

#include "emptyPolyPatch.H"
#include "reconstructedDistanceFunction.H"
#include "processorPolyPatch.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "wedgePolyPatch.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::reconstructedDistanceFunction::coupledFacesPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (isA<coupledPolyPatch>(pp))
        {
            nCoupled += pp.size();
        }
    }
    labelList nCoupledFaces(nCoupled);
    nCoupled = 0;

    for (const polyPatch& pp : patches)
    {
        if (isA<coupledPolyPatch>(pp))
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                nCoupledFaces[nCoupled++] = facei++;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>
        (
            mesh_.faces(),
            nCoupledFaces
        ),
        mesh_.points()
    );
}


void Foam::reconstructedDistanceFunction::markCellsNearSurf
(
    const boolList& interfaceCells,
    const label neiRingLevel
)
{
    neiRingLevel_ = neiRingLevel;

    // performance might be improved by increasing the saving last iterations
    // cells in a Map and loop over the map
    if (mesh_.topoChanging())
    {
        // Introduced resizing to cope with changing meshes
        if (nextToInterface_.size() != mesh_.nCells())
        {
            nextToInterface_.resize(mesh_.nCells());
        }
        coupledBoundaryPoints_ = coupledFacesPatch()().meshPoints();
    }

    const labelListList& pCells = mesh_.cellPoints();
    const labelListList& cPoints = mesh_.pointCells();

    boolList alreadyMarkedPoint(mesh_.nPoints(), false);
    nextToInterface_ = false;

    // do coupled face first
    Map<bool> syncMap;

    for (int level=0;level<=neiRingLevel;level++)
    {
        // parallel
        if (level > 0)
        {
            forAll(coupledBoundaryPoints_, i)
            {
                const label pi = coupledBoundaryPoints_[i];
                forAll(mesh_.pointCells()[pi], j)
                {
                    const label celli = cPoints[pi][j];
                    if (cellDistLevel_[celli] == level-1)
                    {
                        syncMap.insert(pi, true);
                        break;
                    }
                }
            }

            syncTools::syncPointMap(mesh_, syncMap, orEqOp<bool>());

            // mark parallel points first
            forAllConstIters(syncMap, iter)
            {
                const label pi = iter.key();

                if (!alreadyMarkedPoint[pi])
                {
                    // loop over all cells attached to the point
                    forAll(cPoints[pi], j)
                    {
                        const label pCelli = cPoints[pi][j];
                        if (cellDistLevel_[pCelli] == -1)
                        {
                            cellDistLevel_[pCelli] = level;
                            nextToInterface_[pCelli] = true;
                        }
                    }
                }
                alreadyMarkedPoint[pi] = true;
            }
        }


        forAll(cellDistLevel_, celli)
        {
            if (level == 0)
            {
                if (interfaceCells[celli])
                {
                    cellDistLevel_[celli] = 0;
                    nextToInterface_[celli] = true;
                }
                else
                {
                    cellDistLevel_[celli] = -1;
                }
            }
            else
            {
                if (cellDistLevel_[celli] == level-1)
                {
                    forAll(pCells[celli], i)
                    {
                        const label pI = pCells[celli][i];

                        if (!alreadyMarkedPoint[pI])
                        {
                            forAll(cPoints[pI], j)
                            {
                                const label pCelli = cPoints[pI][j];
                                if (cellDistLevel_[pCelli] == -1)
                                {
                                    cellDistLevel_[pCelli] = level;
                                    nextToInterface_[pCelli] = true;
                                }
                            }
                        }
                        alreadyMarkedPoint[pI] = true;
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstructedDistanceFunction::reconstructedDistanceFunction
(
    const fvMesh& mesh
)
:
    volScalarField
    (
        IOobject
        (
            "RDF",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimLength, Zero)
    ),
    mesh_(mesh),
    coupledBoundaryPoints_(coupledFacesPatch()().meshPoints()),
    cellDistLevel_
    (
        IOobject
        (
            "cellDistLevel",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("cellDistLevel", dimless, -1)
    ),
    nextToInterface_(mesh.nCells(), false),
    neiRingLevel_(1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField&  Foam::reconstructedDistanceFunction::constructRDF
(
    const boolList& nextToInterface,
    const volVectorField& centre,
    const volVectorField& normal,
    zoneDistribute& distribute,
    bool updateStencil
)
{
    volScalarField& reconDistFunc = *this;

    if (nextToInterface.size() != centre.size())
    {
        FatalErrorInFunction
            << "size of nextToInterface: " << nextToInterface.size()
            << "size of centre:" <<  centre.size()
            << "do not match. Did the mesh change?"
            << exit(FatalError);
        return reconDistFunc;
    }


    distribute.setUpCommforZone(nextToInterface, updateStencil);

    Map<vector> mapCentres =
        distribute.getDatafromOtherProc(nextToInterface, centre);
    Map<vector> mapNormal =
        distribute.getDatafromOtherProc(nextToInterface, normal);

    const labelListList& stencil = distribute.getStencil();

//- Old algorithm
/*
   forAll(nextToInterface,celli)
    {
        if (nextToInterface[celli])
        {
            if (mag(normal[celli]) != 0) // interface cell
            {
                vector n = -normal[celli]/mag(normal[celli]);
                scalar dist = (centre[celli] - mesh_.C()[celli]) & n;
                reconDistFunc[celli] = dist;
            }
            else // nextToInterfaceCell or level == 1 cell
            {
                scalar averageDist = 0;
                scalar avgWeight = 0;
                const point p = mesh_.C()[celli];

                forAll(stencil[celli], i)
                {
                    const label gblIdx = stencil[celli][i];
                    vector n = -distribute.getValue(normal, mapNormal, gblIdx);
                    if (mag(n) != 0)
                    {
                        n /= mag(n);
                        vector c = distribute.getValue(centre,mapCentres,gblIdx);
                        vector distanceToIntSeg = (c - p);
                        scalar distToSurf = distanceToIntSeg & (n);
                        scalar weight = 0;

                        if (mag(distanceToIntSeg) != 0)
                        {
                            distanceToIntSeg /= mag(distanceToIntSeg);
                            weight = sqr(mag(distanceToIntSeg & n));
                        }
                        else // exactly on the center
                        {
                            weight = 1;
                        }
                        averageDist += distToSurf * weight;
                        avgWeight += weight;
                    }

                }

                if (avgWeight != 0)
                {
                    reconDistFunc[celli] = averageDist / avgWeight;
                }
            }
        }
        else
        {
            reconDistFunc[celli] = 0;
        }
    }

    forAll(reconDistFunc.boundaryField(), patchI)
    {
        fvPatchScalarField& pRDF = reconDistFunc.boundaryFieldRef()[patchI];
        if (isA<calculatedFvPatchScalarField>(pRDF))
        {
            const polyPatch& pp = pRDF.patch().patch();
            forAll(pRDF, i)
            {
                const label pCellI = pp.faceCells()[i];

                if (nextToInterface_[pCellI])
                {
                    scalar averageDist = 0;
                    scalar avgWeight = 0;
                    const point p = mesh_.C().boundaryField()[patchI][i];

                    forAll(stencil[pCellI], j)
                    {
                        const label gblIdx = stencil[pCellI][j];
                        vector n = -distribute.getValue(normal, mapNormal, gblIdx);
                        if (mag(n) != 0)
                        {
                            n /= mag(n);
                            vector c =
                                distribute.getValue(centre, mapCentres, gblIdx);
                            vector distanceToIntSeg = (c - p);
                            scalar distToSurf = distanceToIntSeg & (n);
                            scalar weight = 0;

                            if (mag(distanceToIntSeg) != 0)
                            {
                                distanceToIntSeg /= mag(distanceToIntSeg);
                                weight = sqr(mag(distanceToIntSeg & n));
                            }
                            else // exactly on the center
                            {
                                weight = 1;
                            }
                            averageDist += distToSurf * weight;
                            avgWeight += weight;
                        }
                    }

                    if (avgWeight != 0)
                    {
                        pRDF[i] = averageDist / avgWeight;
                    }
                    else
                    {
                        pRDF[i] = 0;
                    }
                }
                else
                {
                    pRDF[i] = 0;
                }
            }
        }
    } */

    //- FM: new iterative algorithm
    //  The RDF is computed in successive steps where only node neighbours
    //  communicate.

    //- Store normal direction field (this will be non zero also far from)
    //  the interface.
    volVectorField  localNormal(normal);

    //- Set the RDF to zero
    reconDistFunc *= 0.;

    //- Cell centers should be accessible within the stencil
    Map<vector> mapCellC =
        distribute.getDatafromOtherProc(nextToInterface, mesh_.C());

    //- Create new maps for exchange of updated fields
    Map<vector> mapLocalNormal =
        distribute.getDatafromOtherProc(nextToInterface, localNormal);

    Map<scalar> mapRDF =
        distribute.getDatafromOtherProc(nextToInterface, reconDistFunc);

    for(label level(0); level < neiRingLevel_ + 1; level++ )
    {

        forAll(nextToInterface,celli)
        {
            //- Select one level at the time to save computations
            if
            (
                cellDistLevel_[celli] > level-1
                &&
                cellDistLevel_[celli] < level+1
            )
            {
                if (mag(normal[celli]) > SMALL)
                {
                    //- This only works for cells at the interface
                    vector n = -normal[celli]/mag(normal[celli]);
                    scalar dist = (centre[celli] - mesh_.C()[celli]) & n;
                    reconDistFunc[celli] = dist;
                    localNormal[celli] = n;
                }
                else if
                (
                    level == 1
                    ||
                    level == 0
                )// nextToInterfaceCell or level == 1 cell
                {
                    scalar averageDist = 0;
                    scalar avgWeight = 0;
                    const point p = mesh_.C()[celli];

                    forAll(stencil[celli], i)
                    {
                        const label gblIdx = stencil[celli][i];
                        vector n = -distribute.getValue
                        (
                            normal,
                            mapNormal,
                            gblIdx
                        );

                        if (mag(n) > SMALL)
                        {
                            n /= mag(n);
                            vector c = distribute.getValue
                            (
                                centre,
                                mapCentres,
                                gblIdx
                            );
                            vector distanceToIntSeg = (c - p);
                            scalar distToSurf = distanceToIntSeg & (n);
                            scalar weight(0.);

                            if (mag(distanceToIntSeg) > SMALL)
                            {
                                distanceToIntSeg /= mag(distanceToIntSeg);
                                weight = sqr(mag(distanceToIntSeg & n));
                            }
                            else // exactly on the center
                            {
                                weight = 1.0;
                            }
                            averageDist += distToSurf * weight;
                            avgWeight += weight;

                            //- The local normal should be the weighted average
                            //  of the normals
                            localNormal[celli] += n;
                        }

                    }

                    if (avgWeight > SMALL)
                    {
                        reconDistFunc[celli] = averageDist / avgWeight;

                        //- Also store the local normal
                        localNormal[celli] /= avgWeight;
                    }
                    else
                    {
                        reconDistFunc[celli] = 0;
                    }

                }
                else
                {
                    //- For furhter cells the RDF is built based on the RDF
                    //  of the neighbours and their local normals.
                    scalar averageDist = 0;
                    scalar avgWeight = 0;
                    const point p = mesh_.C()[celli];

                    forAll(stencil[celli], i)
                    {
                        const label gblIdx = stencil[celli][i];
                        vector n =
                            distribute.getValue
                            (
                                localNormal,
                                mapLocalNormal,
                                gblIdx
                            );

                        if (mag(n) > SMALL)
                        {
                            n /= mag(n);

                            //- Get neighbour RDF
                            scalar rdfN =
                                distribute.getValue
                                (
                                    reconDistFunc,
                                    mapRDF,
                                    gblIdx
                                );

                            //- Get neighbour cell centre
                            vector c =
                                distribute.getValue
                                (
                                    mesh_.C(),
                                    mapCellC,
                                    gblIdx
                                );

                            vector distanceToIntSeg = (c - p);
                            scalar distToSurf =
                                rdfN + (distanceToIntSeg & (n));
                            scalar weight(0.);

                            if (mag(distanceToIntSeg) > SMALL)
                            {
                                distanceToIntSeg /= mag(distanceToIntSeg);
                                weight = sqr(mag(distanceToIntSeg & n));
                            }
                            else // exactly on the center
                            {
                                weight = 1.;
                            }
                            averageDist += distToSurf * weight;
                            avgWeight += weight;

                            //- The local normal should be the weighted average
                            //  of the normals
                            localNormal[celli] += n;
                        }

                    }

                    if (avgWeight > SMALL)
                    {
                        reconDistFunc[celli] = averageDist / avgWeight;

                        //- Also store the local normal
                        localNormal[celli] /= avgWeight;
                    }
                }
            }

        }

        //- Update new maps for exchange of updated fields
        mapLocalNormal =
            distribute.getDatafromOtherProc(nextToInterface, localNormal);

        mapRDF =
            distribute.getDatafromOtherProc(nextToInterface, reconDistFunc);

    }


    //- Now the same algorithm is applied to the boundary field
    //  There is no need to loop over the levels because the maps are filled
    //  already
    forAll(reconDistFunc.boundaryField(), patchI)
    {
        fvPatchScalarField& pRDF = reconDistFunc.boundaryFieldRef()[patchI];
        if (isA<calculatedFvPatchScalarField>(pRDF))
        {
            const polyPatch& pp = pRDF.patch().patch();
            forAll(pRDF, i)
            {
                const label pCellI = pp.faceCells()[i];

                if (mag(normal[pCellI]) > SMALL)
                {
                    //- To be consistent with cells at the interface
                    vector n = -normal[pCellI]/mag(normal[pCellI]);
                    scalar dist =
                    (
                        centre[pCellI] - mesh_.C().boundaryField()[patchI][i]
                    ) & n;

                    pRDF[i] = dist;
                }
                else
                {

                    //- On the next neighbour (standard algorithm)
                    scalar averageDist = 0;
                    scalar avgWeight = 0;
                    const point p = mesh_.C().boundaryField()[patchI][i];

                    forAll(stencil[pCellI], j)
                    {
                        const label gblIdx = stencil[pCellI][j];
                        vector n = distribute.getValue
                        (
                            localNormal,
                            mapLocalNormal,
                            gblIdx
                        );

                        if (mag(n) > SMALL)
                        {
                            n /= mag(n);
                            //- Get neighbour RDF
                            scalar rdfN =
                                distribute.getValue
                                (
                                    reconDistFunc,
                                    mapRDF,
                                    gblIdx
                                );

                            //- Get neighbour cell centre
                            vector c =
                                distribute.getValue
                                (
                                    mesh_.C(),
                                    mapCellC,
                                    gblIdx
                                );
                            vector distanceToIntSeg = (c - p);

                            scalar distToSurf =
                                rdfN + (distanceToIntSeg & (n));
                            scalar weight = 0;

                            if (mag(distanceToIntSeg) > SMALL)
                            {
                                distanceToIntSeg /= mag(distanceToIntSeg);
                                weight = sqr(mag(distanceToIntSeg & n));
                            }
                            else // exactly on the center
                            {
                                weight = 1.;
                            }
                            averageDist += distToSurf * weight;
                            avgWeight += weight;
                        }
                    }

                    if (avgWeight > SMALL)
                    {
                        pRDF[i] = averageDist / avgWeight;
                    }
                    else
                    {
                        pRDF[i] = 0;
                    }

                }
            }
        }
    }

    reconDistFunc.correctBoundaryConditions();

    return reconDistFunc;
}


void Foam::reconstructedDistanceFunction::updateContactAngle
(
    const volScalarField& alpha,
    const volVectorField& U,
    surfaceVectorField::Boundary& nHatb
)
{
    const fvMesh& mesh = alpha.mesh();
    const volScalarField::Boundary& abf = alpha.boundaryField();
    volScalarField::Boundary& RDFbf = this->boundaryFieldRef();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<alphaContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                degToRad()*acap.theta(U.boundaryField()[patchi], nHatp)
            );

            RDFbf[patchi] =
                1.0/acap.patch().deltaCoeffs()*cos(theta)
              + RDFbf[patchi].patchInternalField();
        }
    }
}


// ************************************************************************* //
