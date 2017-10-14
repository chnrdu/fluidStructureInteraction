/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "directionMixedDisplacementFvPatchVectorField.H"
#include "symmTransformField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

directionMixedDisplacementFvPatchVectorField
::directionMixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF)
{}


directionMixedDisplacementFvPatchVectorField
::directionMixedDisplacementFvPatchVectorField
(
    const directionMixedDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper)
{}


directionMixedDisplacementFvPatchVectorField
::directionMixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF, dict)
{
    Info << "Direction mixed boundary condition with non-orthogonal correction"
        << endl;
    directionMixedFvPatchVectorField::evaluate();
}


directionMixedDisplacementFvPatchVectorField
::directionMixedDisplacementFvPatchVectorField
(
    const directionMixedDisplacementFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void directionMixedDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


void directionMixedDisplacementFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);

//     const directionMixedDisplacementFvPatchVectorField& dmptf =
//         refCast<const directionMixedDisplacementFvPatchVectorField >(ptf);
}


tmp<Field<vector> > directionMixedDisplacementFvPatchVectorField
::snGrad() const
{
    Field<vector> pif = this->patchInternalField();

    Field<vector> normalValue = 
        transform(this->valueFraction(), this->refValue());

    vectorField n = this->patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);

    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + this->dimensionedInternalField().name() + ")"
        );

    Field<vector> gradValue = 
        pif + (k&gradD.patchInternalField())
      + this->refGrad()/this->patch().deltaCoeffs();

    Field<vector> transformGradValue =
        transform(I - this->valueFraction(), gradValue);

//     vectorField sng =
//     (
//         normalValue + transformGradValue 
//       - (pif + (k&gradD.patchInternalField()))
//     )
//    *this->patch().deltaCoeffs();

//     Info << sng << endl;

    return
    (
        normalValue + transformGradValue 
      - (pif + (k&gradD.patchInternalField()))
    )
   *this->patch().deltaCoeffs();
}


void directionMixedDisplacementFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField n = this->patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);

    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + this->dimensionedInternalField().name() + ")"
        );

    Field<vector> normalValue = 
        transform(this->valueFraction(), this->refValue());

    Field<vector> gradValue =
        this->patchInternalField()
      + (k&gradD.patchInternalField())
      + this->refGrad()/this->patch().deltaCoeffs();

    Field<vector> transformGradValue =
        transform(I - this->valueFraction(), gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);

    fvPatchField<vector>::evaluate();
}


void directionMixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField, 
    directionMixedDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
