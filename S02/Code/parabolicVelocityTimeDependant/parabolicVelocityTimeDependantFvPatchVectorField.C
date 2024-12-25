/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "parabolicVelocityTimeDependantFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parabolicVelocityTimeDependantFvPatchVectorField::
parabolicVelocityTimeDependantFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    u0_(dict.lookup<scalar>("u0")),
	A_(dict.lookup<scalar>("A")),
    omega_(dict.lookup<scalar>("omega")),
    radius_(dict.lookup<scalar>("radius"))
    
{}


Foam::parabolicVelocityTimeDependantFvPatchVectorField::
parabolicVelocityTimeDependantFvPatchVectorField
(
    const parabolicVelocityTimeDependantFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    u0_(ptf.u0_),
    A_(ptf.A_),
    omega_(ptf.omega_),
    radius_(ptf.radius_)
{}


Foam::parabolicVelocityTimeDependantFvPatchVectorField::
parabolicVelocityTimeDependantFvPatchVectorField
(
    const parabolicVelocityTimeDependantFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    u0_(ptf.u0_),
    A_(ptf.A_),
    omega_(ptf.omega_),
    radius_(ptf.radius_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::parabolicVelocityTimeDependantFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // const scalar t = this->db().time().value();
 /*   const scalarField ts(size(), t);

    // Compute geometry
    const vector axisHat = normalised(axis_);
    const vectorField d(patch().Cf() - origin_);
    const vectorField r(d - (axisHat & d)*axisHat);
    const scalarField magR(mag(r));
    const vectorField rHat(normalised(r));

    // Evaluate individual velocity components
    const scalarField axialVelocity(axialVelocity_->value(ts, magR));
    const scalarField radialVelocity(radialVelocity_->value(ts, magR));
    tmp<scalarField> tangentialVelocity;
    if (omega_.valid())
    {
        tangentialVelocity = omega_->value(t)*magR;
    }
    else
    {
        tangentialVelocity = tangentialVelocity_->value(ts, magR);
    }

    // Combine components the complete vector velocity
    operator==
    (
        axialVelocity*axisHat
      + radialVelocity*rHat
      + tangentialVelocity*(axisHat ^ rHat)
    );

    fixedValueFvPatchField<vector>::updateCoeffs();
    */
    const scalar t = this->db().time().value();
    
    vectorField r(patch().Cf() - origin_);
    scalarField magr(mag(r));
    
    scalar umax  = u0_ + A_ * sin(omega_ * t);
    
    scalarField parabolicVelocityTimeDependant = umax * (1 - sqr(magr/radius_));
    
    vector normalizedAxis = axis_ / mag(axis_);
    
    operator == (parabolicVelocityTimeDependant * normalizedAxis);
    
    fixedValueFvPatchField<vector>::updateCoeffs();
    
}


void Foam::parabolicVelocityTimeDependantFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry(os, "origin", origin_);
    writeEntry(os, "axis", axis_);
    writeEntry(os, "u0", u0_);
    writeEntry(os, "A", A_);
    writeEntry(os, "omega", omega_);
    writeEntry(os, "radius", radius_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       parabolicVelocityTimeDependantFvPatchVectorField
   );
}


// ************************************************************************* //
