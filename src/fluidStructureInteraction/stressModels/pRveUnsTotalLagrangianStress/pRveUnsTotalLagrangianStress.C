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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "pRveUnsTotalLagrangianStress.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stressModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pRveUnsTotalLagrangianStress, 0);
addToRunTimeSelectionTable
(
    stressModel, 
    pRveUnsTotalLagrangianStress, 
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


pRveUnsTotalLagrangianStress::pRveUnsTotalLagrangianStress(const fvMesh& mesh)
:
    unsTotalLagrangianStress(mesh),
    avgDeformationGradient_
    (
        stressProperties().lookup("avgDeformationGradient")
    ),
    totPointD_
    (
        IOobject
        (
            "totPointD",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    totSigma_
    (
        IOobject
        (
            "totSigma",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    ),
    totEpsilon_
    (
        IOobject
        (
            "totEpsilon",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pRveUnsTotalLagrangianStress::updateTotalFields()
{
    Info << "Update total fields" << endl;

    totPointD_.internalField() = 
        pointD().internalField()
      + ((avgDeformationGradient_ - I) & mesh().points());

    Switch nonLinear(stressProperties().lookup("nonLinear"));
    Switch enforceLinear(stressProperties().lookup("enforceLinear"));

    // Calc average strain
    symmTensor avgE = 
        0.5*symm(avgDeformationGradient_ + avgDeformationGradient_.T()) - I;
    if (nonLinear && !enforceLinear)
    {
        avgE = (symm(avgDeformationGradient_.T()&avgDeformationGradient_) - I);
    }

    totEpsilon_ = epsilon() + avgE;

    // Calc average stress
    volSymmTensorField avgSigma = 
        2*mu()*avgE + lambda()*tr(avgE)*I;

    totSigma_ = sigma() + avgSigma;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stressModels
} // End namespace Foam

// ************************************************************************* //
