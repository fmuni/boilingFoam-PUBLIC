/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::activationModels::liquidSuperheatingThreshold

Description
    Use alpha field and temperature as thresholds

SourceFiles
    liquidSuperheatingThreshold.C

\*---------------------------------------------------------------------------*/

#ifndef liquidSuperheatingThreshold_H
#define liquidSuperheatingThreshold_H

#include "activationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace activationModels
{

/*---------------------------------------------------------------------------*\
                           Class liquidSuperheatingThreshold Declaration
\*---------------------------------------------------------------------------*/

class liquidSuperheatingThreshold
:
    public activationModel
{
    const scalar    Tact_;
    const scalar    alphaTol_;
    const word      phaseName_;

public:

    //- Runtime type information
    TypeName("liquidSuperheatingThreshold");


    // Constructors

        liquidSuperheatingThreshold
        (
            const dictionary dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~liquidSuperheatingThreshold();


    // Member Functions

        virtual bool isActive(const vector& location) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace activationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
