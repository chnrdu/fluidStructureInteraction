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

#include "tractionDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "constitutiveModel.H"
#include "stressModel.H"
#include "pRveUnsTotalLagrangianStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
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


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
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

    Info << "Creating traction displacement boundary conditions" << endl;
}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper)
{}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void tractionDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionDisplacementFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void tractionDisplacementFvPatchVectorField::updateCoeffs()
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

    bool thermalStress = stress.thermalStress();

    word DName = this->dimensionedInternalField().name();

    const fvsPatchField<tensor>& gradDf =
        patch().lookupPatchField<surfaceTensorField, tensor>
        (
            "grad" + DName + "f"
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

    vectorField t = traction_;

    {
        if (nonLinear && !enforceLinear)
        {
            tensorField F = I + gradDf;

            scalarField J = det(F);

            tensorField invF = hinv(F);

            scalarField SoS0 = mag(J*(invF & n));

            vectorField nCurrent = (invF & n);
            nCurrent /= mag(nCurrent);

            t -= pressure_*nCurrent;

            t = (t & invF)*SoS0;
        }
        else
        {
            t -= pressure_*n;
        }

        gradient() =
            t
          - (n & (mu*gradDf.T() - (mu + lambda)*gradDf))
          - n*lambda*tr(gradDf);

        if (nonLinear && !enforceLinear)
        {
            gradient() -=
                (n & (mu*(gradDf & gradDf.T())))
              + 0.5*n*lambda*tr(gradDf & gradDf.T());
        }
        
        if (thermalStress)
        {
            const fvPatchField<scalar>& DT =
                patch().lookupPatchField<volScalarField, scalar>("DT");

            const fvsPatchField<scalar>& threeK =
                patch().lookupPatchField<surfaceScalarField, scalar>
                (
                    "threeKf"
                );

            const fvsPatchField<scalar>& alpha =
                patch().lookupPatchField<surfaceScalarField, scalar>
                (
                    "alphaf"
                );

            gradient() += n*threeK*alpha*DT;
        }

        gradient() /= (2.0*mu + lambda);
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}

void tractionDisplacementFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField n = patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);

    word UName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + UName + ")"
        );

    Field<vector>::operator=
    (
        this->patchInternalField() 
      + (k&gradU.patchInternalField())
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
void tractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    writeEntry("value", os);
}


// void tractionDisplacementFvPatchVectorField::operator=
// (
//     const fvPatchField<vector>& ptf
// )
// {
//     fvPatchField<vector>::operator=(ptf);

//     gradient() = fvPatchField::snGrad();

//     if (ptf.type() == tractionDisplacementFvPatchVectorField::typeName)
//     {
//         const tractionDisplacementFvPatchVectorField& dmptf =
//             refCast<const tractionDisplacementFvPatchVectorField>(ptf);

//         gradient() = dmptf.gradient();
//     }
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, tractionDisplacementFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
