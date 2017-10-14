/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "tractionDisplacementIncrementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "constitutiveModel.H"
#include "stressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionDisplacementIncrementFvPatchVectorField::
tractionDisplacementIncrementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionDisplacementIncrementFvPatchVectorField::
tractionDisplacementIncrementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size())
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;

    Info << "Creating traction displacement incr boundary conditions" << endl;
}


tractionDisplacementIncrementFvPatchVectorField::
tractionDisplacementIncrementFvPatchVectorField
(
    const tractionDisplacementIncrementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper)
{}


tractionDisplacementIncrementFvPatchVectorField::
tractionDisplacementIncrementFvPatchVectorField
(
    const tractionDisplacementIncrementFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


tractionDisplacementIncrementFvPatchVectorField::
tractionDisplacementIncrementFvPatchVectorField
(
    const tractionDisplacementIncrementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionDisplacementIncrementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void tractionDisplacementIncrementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionDisplacementIncrementFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementIncrementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void tractionDisplacementIncrementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Looking up stress model
    const stressModel& stress =
        this->db().objectRegistry::lookupObject<stressModel>
        (
            "stressProperties"
        );

    Switch nonLinear
    (
        stress.stressProperties().lookup("nonLinear")
    );

    Switch enforceLinear
    (
        stress.stressProperties().lookup("enforceLinear")
    );

    word DDName = this->dimensionedInternalField().name();

    const fvsPatchField<tensor>& gradDD =
        patch().lookupPatchField<surfaceTensorField, tensor>
        (
            "grad" + DDName + "f"
        );

    const fvsPatchField<tensor>& gradD =
        patch().lookupPatchField<surfaceTensorField, tensor>
        (
            "gradDf"
        );

    const fvsPatchField<symmTensor>& sigma =
        patch().lookupPatchField<surfaceSymmTensorField, symmTensor>
        (
            "sigmaf"
        );

    const fvsPatchField<scalar>& mu =
        patch().lookupPatchField<surfaceScalarField, scalar>
        (
            "muf"
        );

    const fvsPatchField<scalar>& lambda =
        patch().lookupPatchField<surfaceScalarField, scalar>
        (
            "lambdaf"
        );

    vectorField n = patch().nf();

    {
        vectorField t = traction_;

        if (nonLinear && !enforceLinear)
        {
            tensorField F = I + gradD + gradDD;

            scalarField J = det(F);

            tensorField invF = hinv(F);

            scalarField SoS0 = mag(J*(invF & n));

            vectorField nCurrent = (invF & n);
            nCurrent /= mag(nCurrent);

            // Cauchy traction
            t -= pressure_*nCurrent;

            // 2nd Piola-Kirchhoff traction
            t = (t & invF)*SoS0;
        }
        else
        {
            t -= pressure_*n;
        }

        // 2nd Piola-Kirchhoff traction increment
        vectorField DTraction = t - (n&sigma);

        gradient() =
            DTraction
          - (n & (mu*gradDD.T() - (mu + lambda)*gradDD))
          - n*lambda*tr(gradDD);

        if (nonLinear && !enforceLinear)
        {
            gradient() -=
                (n & (mu*(gradDD & gradDD.T())))
              + (n & (mu*(gradDD & gradD.T())))
              + (n & (mu*(gradD & gradDD.T())))
              + 0.5*n*lambda*tr(gradDD & gradDD.T())
              + 0.5*n*lambda*tr(gradDD & gradD.T())
              + 0.5*n*lambda*tr(gradD & gradDD.T());
        }

        if (stress.rheology().plasticityActive())
        {
            gradient() +=
                2*mu
               *(
                    n 
                  & stress.rheology()
                   .DEpsilonP().boundaryField()[patch().index()]
                );
        }

        gradient() /= (2.0*mu + lambda);
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}

void tractionDisplacementIncrementFvPatchVectorField
::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField n = patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);

    word DDName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradDD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + DDName + ")"
        );

    Field<vector>::operator=
    (
        this->patchInternalField() 
      + (k&gradDD.patchInternalField())
      + gradient()/this->patch().deltaCoeffs()
    );

//     vectorField dUP = (k&gradU.patchInternalField());

//     vectorField nGradUP = (n&gradU.patchInternalField());

//     Field<vector>::operator=
//     (
//         this->patchInternalField() + dUP
//       + 0.5*(gradient() + nGradUP)
//        /this->patch().deltaCoeffs()
//     );

    fvPatchField<vector>::evaluate();
}


// Write
void tractionDisplacementIncrementFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField, 
    tractionDisplacementIncrementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
