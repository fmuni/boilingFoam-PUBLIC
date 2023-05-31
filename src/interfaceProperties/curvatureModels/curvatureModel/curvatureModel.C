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

#include "curvatureModel.H"
#include "fvMesh.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "surfaceFields.H"
#include "surfaceMesh.H"
#include "unitConversion.H"
#include "fvcAverage.H"
#include "fvcReconstruct.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(curvatureModel, 0);
    defineRunTimeSelectionTable(curvatureModel, dictionary);
}

const Foam::dimensionSet Foam::curvatureModel::dimK(0, -1, 0, 0, 0);

void Foam::curvatureModel::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf,
    const volVectorField& U,
    const volScalarField& alpha1
) const
{
    const fvMesh& mesh = alpha1.mesh();
    const volScalarField::Boundary& abf = alpha1.boundaryField();

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
                degToRad() * acap.theta(U.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::curvatureModel::curvatureModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    reconstructionSchemes& surf
)
:
    regIOobject
    (
        IOobject
        (
            typeName, mesh.name(),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    surf_(surf),
    mesh_(mesh),
    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimK, Zero)
    ),
    nHatf_
    (
        IOobject
        (
            "nHatf",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimArea, Zero)
    ),


    deltaN_
    (
        "deltaN",
        1e-8/cbrt(average(mesh_.V()))
    )
    // nSmoothIter_
    // (
    //     dict.lookupOrDefault("nSmoothIter",0)
    // )
{
    // if(nSmoothIter_ > 0)
    // {
    //     Info<<"\nUsing curvature smoothing with " << nSmoothIter_
    //         <<" iterations."
    //         <<endl;
    // }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::curvatureModel::~curvatureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::curvatureModel::writeData(Ostream& os) const
{
    return os.good();
}

// void Foam::curvatureModel::smoothK() const
// {
//     Ksmt_ = K_;
//     for(label iter(0);iter<nSmoothIter_;iter++)
//     {
//         //surfaceScalarField Ksmtf(fvc::interpolate(Ksmt_)*K_.mesh().magSf());
//         Ksmt_ = fvc::average(Ksmt_); //fvc::reconstructMag(Ksmtf);
//         Ksmt_.correctBoundaryConditions();
//     }
// }


// ************************************************************************* //
