/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "turbulentTemperatureBoilingMicrolayerFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureBoilingMicrolayerFvPatchScalarField::
turbulentTemperatureBoilingMicrolayerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase
    (
        patch(),
        "undefined",
        "undefined",
        "undefined-K",
        "undefined-alpha"
    ),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this
    ),
    TnbrName_("undefined-Tnbr"),
    coeffs_(0,0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureBoilingMicrolayerFvPatchScalarField::
turbulentTemperatureBoilingMicrolayerFvPatchScalarField
(
    const turbulentTemperatureBoilingMicrolayerFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        ptf
    ),
    TnbrName_(ptf.TnbrName_),
    coeffs_(ptf.coeffs_)
{}


turbulentTemperatureBoilingMicrolayerFvPatchScalarField::
turbulentTemperatureBoilingMicrolayerFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        dict
    ),
    TnbrName_(dict.get<word>("Tnbr")),
    coeffs_(dict.get<scalarField>("coefficients"))
{

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }

    // Store patch value as initial guess when running in database mode
    mappedPatchFieldBase<scalar>::initRetrieveField
    (
        this->internalField().name(),
        *this
    );
    mappedPatchFieldBase<scalar>::initRetrieveField
    (
        this->internalField().name() + "_weights",
        this->patch().deltaCoeffs()
    );
}


turbulentTemperatureBoilingMicrolayerFvPatchScalarField::
turbulentTemperatureBoilingMicrolayerFvPatchScalarField
(
    const turbulentTemperatureBoilingMicrolayerFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    temperatureCoupledBase(patch(), wtcsf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), iF),
        *this,
        wtcsf
    ),
    TnbrName_(wtcsf.TnbrName_),
    coeffs_(wtcsf.coeffs_)
{}


turbulentTemperatureBoilingMicrolayerFvPatchScalarField::
turbulentTemperatureBoilingMicrolayerFvPatchScalarField
(
    const turbulentTemperatureBoilingMicrolayerFvPatchScalarField& wtcsf
)
:
    mixedFvPatchScalarField(wtcsf),
    temperatureCoupledBase(patch(), wtcsf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), wtcsf.internalField()),
        *this,
        wtcsf
    ),
    TnbrName_(wtcsf.TnbrName_),
    coeffs_(wtcsf.coeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureBoilingMicrolayerFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        mappedPatchFieldBase<scalar>::mapper
        (
            patch(),
            this->internalField()
        );

    const tmp<scalarField> myKDelta = kappa(*this)*patch().deltaCoeffs();
    tmp<scalarField> nbrIntFld;
    tmp<scalarField> nbrKDelta;
    scalarField nbrCoeffs;

    //Pout<< "updateCoeffs() : mpp.sameWorld():" << mpp.sameWorld() << endl;


    if (mpp.sameWorld())
    {
        // Same world so lookup
        const auto& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
        const label nbrPatchID = mpp.samplePolyPatch().index();
        const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];

        const turbulentTemperatureBoilingMicrolayerFvPatchScalarField&
        nbrField =
        refCast
        <
        const turbulentTemperatureBoilingMicrolayerFvPatchScalarField
        >
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>
            (
                TnbrName_
            )
        );


        // Get neighbour internal data in local order. Does all
        // comms/reordering already
        nbrIntFld = this->mappedInternalField();
        // Calculate neighbour weighting (in neighbouring ordering)
        nbrKDelta = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();

        nbrCoeffs = nbrFld.layerCoeffs();
    }
    else
    {
        // Different world so use my region,patch. Distribution below will
        // do the reordering
        nbrIntFld = this->mappedInternalField();
        nbrKDelta = new scalarField(myKDelta());
    }


    //Pout<< "updateCoeffs() : nbrIntFld:" << flatOutput(nbrIntFld()) << endl;
    //Pout<< "updateCoeffs() : nbrKDelta BEFORE:" << flatOutput(nbrKDelta())
    //    << endl;

    distribute(this->internalField().name() + "_weights", nbrKDelta.ref());
    //Pout<< "updateCoeffs() : nbrKDelta AFTER:" << flatOutput(nbrKDelta())
    //    << endl;
    //Pout<< "updateCoeffs() : myKDelta:" << flatOutput(myKDelta())
    //    << endl;


    // Both sides agree on
    // - temperature : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
    // - gradient    : (temperature-fld)*delta
    // We've got a degree of freedom in how to implement this in a mixed bc.
    // (what gradient, what fixedValue and mixing coefficient)
    // Two reasonable choices:
    // 1. specify above temperature on one side (preferentially the high side)
    //    and above gradient on the other. So this will switch between pure
    //    fixedvalue and pure fixedgradient
    // 2. specify gradient and temperature such that the equations are the
    //    same on both sides. This leads to the choice of
    //    - refGradient = zero gradient
    //    - refValue = neighbour value
    //    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    this->refValue() = nbrIntFld();
    this->refGrad() = 0.0;
    this->valueFraction() = nbrKDelta()/(nbrKDelta() + myKDelta());

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << mpp.sampleRegion() << ':'
            << mpp.samplePatch() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


void turbulentTemperatureBoilingMicrolayerFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntry("Tnbr", TnbrName_);

    temperatureCoupledBase::write(os);
    mappedPatchFieldBase<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureBoilingMicrolayerFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
