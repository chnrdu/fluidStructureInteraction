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

#include "pRveUnsIncrTotalLagrangianStress.H"

#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"

#include "addToRunTimeSelectionTable.H"

// #include "fvcGradf.H"
// #include "tractionDisplacementIncrementFvPatchVectorField.H"
// #include "skewCorrectionVectors.H"
// #include "multiMaterial.H"
// #include "twoDPointCorrector.H"

// #include "componentReferenceList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stressModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pRveUnsIncrTotalLagrangianStress, 0);
addToRunTimeSelectionTable
(
    stressModel,
    pRveUnsIncrTotalLagrangianStress,
    dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pRveUnsIncrTotalLagrangianStress::pRveUnsIncrTotalLagrangianStress
(
    const fvMesh& mesh
)
:
    unsIncrTotalLagrangianStress(mesh),
    avgDEpsilon_
    (
        stressProperties().lookup("avgDEpsilon")
    ),
    avgEpsilon_
    (
        stressProperties().lookup("avgEpsilon")
    ),
//     avgDeformationGradientIncrement_
//     (
//         stressProperties().lookup("avgDeformationGradientIncrement")
//     ),
//     avgDeformationGradient_
//     (
//         stressProperties().lookup("avgDeformationGradient")
//     ),
    totD_
    (
        IOobject
        (
            "totD",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
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
    totSigmaf_
    (
        IOobject
        (
            "totSigmaf",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// symmTensor pRveUnsIncrTotalLagrangianStress::avgDEpsilon() const
// {
//     symmTensor avgDEpsilon = symmTensor::zero;

//     Switch nonLinear
//     (
//         stressProperties().lookup("nonLinear")
//     );

//     Switch enforceLinear
//     (
//         stressProperties().lookup("enforceLinear")
//     );
    
//     tensor avgDeformationGradientNew = 
//         avgDeformationGradient()
//       + avgDeformationGradientIncrement();

//     symmTensor avgEnew = 
//         0.5
//        *symm
//         (
//             avgDeformationGradientNew 
//           + avgDeformationGradientNew.T()
//         ) - I;
        
//     symmTensor avgEold = 
//         0.5
//        *symm
//         (
//             avgDeformationGradient() 
//           + avgDeformationGradient().T()
//         ) - I;

//     if (nonLinear && !enforceLinear)
//     {
//         avgEnew = 
//             (
//                 symm
//                 (
//                     avgDeformationGradientNew.T()
//                   & avgDeformationGradientNew
//                 ) - I
//             );

//         avgEold = 
//             (
//                 symm
//                 (
//                     avgDeformationGradient_.T()
//                   & avgDeformationGradient_
//                 ) - I
//             );
//     }

//     avgDEpsilon = avgEnew - avgEold;
    
//     return avgDEpsilon;    
// }


// symmTensor pRveUnsIncrTotalLagrangianStress::avgEpsilon() const
// {
//     Switch nonLinear
//     (
//         stressProperties().lookup("nonLinear")
//     );

//     Switch enforceLinear
//     (
//         stressProperties().lookup("enforceLinear")
//     );
    
//     tensor avgDeformationGradientNew = 
//         avgDeformationGradient()
//       + avgDeformationGradientIncrement();

//     symmTensor avgEnew = 
//         0.5
//        *symm
//         (
//             avgDeformationGradientNew 
//           + avgDeformationGradientNew.T()
//         ) - I;
        
//     if (nonLinear && !enforceLinear)
//     {
//         avgEnew = 
//             (
//                 symm
//                 (
//                     avgDeformationGradientNew.T()
//                   & avgDeformationGradientNew
//                 ) - I
//             );
//     }
    
//     return avgEnew;    
// }


// bool pRveUnsIncrTotalLagrangianStress::evolve()
// {
//     Info << "Evolving stress model" << endl;

//     int nCorr
//     (
//         readInt(stressProperties().lookup("nCorrectors"))
//     );

//     scalar convergenceTolerance
//     (
//         readScalar(stressProperties().lookup("convergenceTolerance"))
//     );

//     Switch nonLinear(stressProperties().lookup("nonLinear"));
//     Switch debug(stressProperties().lookup("debug"));

//     componentReferenceList cr
//     (
//         stressProperties().lookup("componentReference"),
//         componentReference::iNew(mesh())
//     );

//     int iCorr = 0;
//     scalar initialResidual = 0;
//     lduMatrix::solverPerformance solverPerf;
//     scalar res = GREAT;

//     lduMatrix::debug = debug;

//     bool enforceLinear = false;
//     stressProperties().set("enforceLinear", enforceLinear);

//     do
//     {
//         if (lduMatrix::debug)
//         {
//             Info<< "Time: " << runTime().timeName() 
//                 << ", outer iteration: " << iCorr << endl;
//         }

//         DD().storePrevIter();

//         fvVectorMatrix DDEqn
//         (
//             rho()*fvm::d2dt2(DD())
//           - fvm::laplacian(2*muf() + lambdaf(), DD(), "laplacian(DDD,DD)")
//          == fvc::div
//             (
//                 mesh().Sf()
//               & (
//                   - (muf() + lambdaf())*gradDDf()
//                   + muf()*gradDDf().T() + lambdaf()*(I*tr(gradDDf()))
//                 )
//             )
//         );

//         // Update strain increment
//         DEpsilonf() = symm(gradDDf());
//         if (nonLinear && !enforceLinear)
//         {
//             DEpsilonf() += 0.5*symm(gradDDf() & gradDDf().T());
//             DEpsilonf() += 0.5*symm(gradDDf() & gradDf().T());
//             DEpsilonf() += 0.5*symm(gradDf() & gradDDf().T());
//         }

//         if (nonLinear && !enforceLinear)
//         {
//             DSigmaf() = 2*muf()*DEpsilonf() + I*(lambdaf()*tr(DEpsilonf()));

//             if (rheology().plasticityActive())
//             {
//                 DSigmaf() -= 2*muf()*rheology().DEpsilonPf();
//             }

//             DDEqn -=
//                 fvc::div
//                 (
//                     muf()*(mesh().Sf() & (gradDDf() & gradDf().T()))
//                   + muf()*(mesh().Sf() & (gradDf() & gradDDf().T()))
//                   + muf()*(mesh().Sf() & (gradDDf() & gradDDf().T()))
//                   + 0.5*lambdaf()*tr(gradDDf() & gradDf().T())*mesh().Sf()
//                   + 0.5*lambdaf()*tr(gradDf() & gradDDf().T())*mesh().Sf()
//                   + 0.5*lambdaf()*tr(gradDDf() & gradDDf().T())*mesh().Sf()
//                 )
//               + fvc::div(mesh().Sf() & (DSigmaf() & gradDf()))
//               + fvc::div(mesh().Sf() & ((sigmaf() + DSigmaf()) & gradDDf()));
//         }

//         if (rheology().plasticityActive())
//         {
//             DDEqn += fvc::div(2*muf()*(mesh().Sf() & rheology().DEpsilonPf()));
//         }

//         if (interface().valid())
//         {
//             interface()->correct(DDEqn);
//         }

//         forAll (cr, crI)
//         {
//             DDEqn.setComponentReference
//             (
//                 cr[crI].patchIndex(),
//                 cr[crI].faceIndex(),
//                 cr[crI].dir(),
//                 cr[crI].value()
//             );
//         }

//         solverPerf = DDEqn.solve();

//         if(iCorr == 0)
//         {
//             initialResidual = solverPerf.initialResidual();
//         }

//         DD().relax();

//         if (interface().valid())
//         {
//             interface()->updateDisplacementIncrement(pointDD());
//             interface()->updateDisplacementIncrementGradient
//             (
//                 gradDD(), 
//                 gradDDf()
//             );
//         }
//         else
//         {
//             volToPoint().interpolate(DD(), pointDD());
//             gradDD() = fvc::grad(DD(), pointDD());
//             gradDDf() = fvc::fGrad(DD(), pointDD());
//         }

//         if (nonLinear && !enforceLinear)
//         {
//             surfaceScalarField Det = det(I+gradDf()+gradDDf());

//             scalar minDetFf = min(Det).value();
//             reduce(minDetFf, minOp<scalar>());

//             scalar maxDetFf = max(Det).value();
//             reduce(maxDetFf, maxOp<scalar>());

//             if ( (minDetFf<0.01) || (maxDetFf>100) )
//             {
//                 Info << minDetFf << ", " << maxDetFf << endl;
//                 enforceLinear = true;
//                 stressProperties().set("enforceLinear", enforceLinear);
//             }
//         }

//         // Calculate momentu residual
//         res = residual();

//         if (lduMatrix::debug)
//         {
//             Info << "Relative residual = " << res << endl;
//         }

//         // Calculate strain increment
//         {
//             DEpsilon() = symm(gradDD());

//             if(nonLinear && !enforceLinear)
//             {
//                 DEpsilon() += 0.5*symm(gradDD() & gradDD().T());
//                 DEpsilon() += 0.5*symm(gradDD() & gradD().T());
//                 DEpsilon() += 0.5*symm(gradD() & gradDD().T());
//             }
//         }

//         // Correct plasticity term
//         rheology().correct();
//     }
//     while
//     (
//         (res > convergenceTolerance) 
//      && (++iCorr < nCorr)
//     );

//     DU() = fvc::ddt(DD());

//     // Calculate second Piola-Kirchhoff stress increment
//     {
//         DSigma() = 2*mu()*DEpsilon() + I*(lambda()*tr(DEpsilon()));

//         if (rheology().plasticityActive())
//         {
//             DSigma() -= 2*mu()*rheology().DEpsilonP();
//         }
//     }

//     Info << solverPerf.solverName() << ": Solving for " << DD().name()
//         << ", Initial residula = " << initialResidual
//         << ", Final residual = " << solverPerf.initialResidual()
//         << ", No outer iterations = " << iCorr 
//         << ", Relative momentum residual = " << res 
//         << ", enforceLinear = " << enforceLinear << endl;

//     lduMatrix::debug = 1;

//     if (nonLinear && enforceLinear)
//     {
//         return false;
//     }

//     return true;
// }


void pRveUnsIncrTotalLagrangianStress::updateTotalFields()
{
    unsIncrTotalLagrangianStress::updateTotalFields();

    avgEpsilon_ += avgDEpsilon_;

    totPointD_.internalField() +=
        pointDD().internalField()
      + (avgDEpsilon_ & mesh().points());

    totD_ += DD() + (avgDEpsilon_ & mesh().C());

//     D() += (avgDeformationGradientIncrement_&mesh().C());

//     pointD().internalField() += 
//         (avgDeformationGradientIncrement_&mesh().points());

//     Switch nonLinear(stressProperties().lookup("nonLinear"));
//     Switch enforceLinear(stressProperties().lookup("enforceLinear"));

    // Calc average strain increment

//     epsilon() += avgDEpsilon();
//     epsilonf() += avgDEpsilon();

    // Calc average stress increment

//     symmTensor avgEnew = avgEpsilon();

    volSymmTensorField avgDSigma = 
        2*mu()*avgDEpsilon_ + lambda()*tr(avgDEpsilon_)*I;

    totSigma_ += DSigma() + avgDSigma;

    surfaceSymmTensorField avgDSigmaf = 
        2*muf()*avgDEpsilon_ + lambdaf()*tr(avgDEpsilon_)*I;

    totSigmaf_ += DSigmaf() + avgDSigmaf;

//     avgDeformationGradient_ += avgDeformationGradientIncrement_;
}


bool pRveUnsIncrTotalLagrangianStress::writeObject
(
    IOstream::streamFormat,
    IOstream::versionNumber,
    IOstream::compressionType
) const
{
    volScalarField totSigmaEq
    (
        IOobject
        (
            "totSigmaEq",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(totSigma_)))
    );
    totSigmaEq.write();

    Info<< "Max sigmaEq = " << max(totSigmaEq).value() << endl;

    Info<< "totSigmaEq, max: " << gMax(totSigmaEq.internalField()) 
        << ", avg: " << gAverage(totSigmaEq.internalField()) 
        << ", min: " << gMin(totSigmaEq.internalField()) << endl;

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stressModels
} // End namespace Foam

// ************************************************************************* //
