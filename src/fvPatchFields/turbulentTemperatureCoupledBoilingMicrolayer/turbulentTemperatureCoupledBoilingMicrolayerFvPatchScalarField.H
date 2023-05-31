/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::compressible::
        turbulentTemperatureBoilingMicrolayerFvPatchScalarField

Description
    Boundary condition for the temperature in the presence of a microlayer.
    The source term from the microlayer is given as a polynomialfunction of the
    temperature at the boundary so that it can be treated semi-implicitly.

    Q_ml = Q' T_p + Q''

    This condition is fundamentally a fixed value, where:

    T_p = w/(1+q') T_n + (1-w)/(1+q') T_c + Q''/(\Delta*(1+q'))

    And the factor q' is given by:

    q' = Q/(\Delta (k_f + k_s))

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureBoilingMicrolayer;
        Tnbr            T;
        kappaMethod     lookup;
        kappa           kappa;
        Qcoeffs         3(1 2 3);
        value           uniform 300;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

See also
    Foam::temperatureCoupledBase

SourceFiles
    turbulentTemperatureBoilingMicrolayerFvPatchScalarField.C

Author
    Federico Municchi. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureBoilingMicrolayerFvPatchScalarField_H
#define turbulentTemperatureBoilingMicrolayerFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "mappedPatchFieldBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class turbulentTemperatureBoilingMicrolayerFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureBoilingMicrolayerFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase,
    public mappedPatchFieldBase<scalar>
{
    // Private data

        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Coefficients of the polynomium
        scalarField coeffs_;


public:

    //- Runtime type information
    TypeName("turbulentTemperatureBoilingMicrolayer");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureBoilingMicrolayerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureBoilingMicrolayerFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureBoilingMicrolayerFvPatchScalarField onto a
        //  new patch
        turbulentTemperatureBoilingMicrolayerFvPatchScalarField
        (
            const turbulentTemperatureBoilingMicrolayerFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        turbulentTemperatureBoilingMicrolayerFvPatchScalarField
        (
            const turbulentTemperatureBoilingMicrolayerFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentTemperatureBoilingMicrolayerFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentTemperatureBoilingMicrolayerFvPatchScalarField
        (
            const turbulentTemperatureBoilingMicrolayerFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentTemperatureBoilingMicrolayerFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get polynomium coefficients
        virtual scalarField layerCoeffs()
        {
            return coeffs_;
        }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
