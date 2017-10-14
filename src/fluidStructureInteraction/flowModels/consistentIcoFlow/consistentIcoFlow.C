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

#include "consistentIcoFlow.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "adjustPhi.H"

#include "findRefCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flowModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(consistentIcoFlow, 0);
addToRunTimeSelectionTable(flowModel, consistentIcoFlow, dictionary);
// addToRunTimeSelectionTable(icoFlow, consistentIcoFlow, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void consistentIcoFlow::makeSf() const
{
    // Find global face zones
    if (SfPtr_)
    {
        FatalErrorIn
        (
            "void flowModel::makeSf() const"
        )
            << "Face surface vectors alrady created"
                << abort(FatalError);
    }

//     const fvMesh& mesh = flowModel::mesh();

    IOobject SfHeader
    (
        "Sf",
        runTime().timeName(),
        mesh(),
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );

    SfPtr_ = 
        new surfaceVectorField
        (
            SfHeader,
            mesh(),
            dimensionedVector("0", dimArea, vector::zero)
        );
    surfaceVectorField& Sf = *SfPtr_;

    if(!SfHeader.headerOk())
    {
        const vectorField& allFaceAreas = mesh().faceAreas();

        Sf.internalField() =
            vectorField::subField(allFaceAreas, mesh().nInternalFaces());

        const fvPatchList& patches = mesh().boundary();

        forAll (patches, patchI)
        {
            Sf.boundaryField()[patchI] = 
                patches[patchI].patchSlice(allFaceAreas);
        }
    }
}


void consistentIcoFlow::updateSf()
{
    Sf().oldTime();

    const vectorField& allFaceAreas = mesh().faceAreas();

    Sf().internalField() =
        vectorField::subField(allFaceAreas, mesh().nInternalFaces());

    const fvPatchList& patches = mesh().boundary();

    forAll (patches, patchI)
    {
        Sf().boundaryField()[patchI] = 
            patches[patchI].patchSlice(allFaceAreas);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

consistentIcoFlow::consistentIcoFlow(const fvMesh& mesh)
:
    icoFlow(mesh),
    SfPtr_(NULL)
{
    phi().oldTime();
    updateSf();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void consistentIcoFlow::evolve()
{
    Info << "Evolving flow model: " << this->type() << endl;

    const fvMesh& mesh = flowModel::mesh();

    updateSf();

    int nCorr(readInt(flowProperties().lookup("nCorrectors")));

    int nNonOrthCorr =
        readInt(flowProperties().lookup("nNonOrthogonalCorrectors"));

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p(), flowProperties(), pRefCell, pRefValue);

    if(mesh.moving())
    {
        // Make the fluxes relative
        phi() -= fvc::meshPhi(U());
    }

    // CourantNo
    {
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;
        scalar velMag = 0.0;

        if (mesh.nInternalFaces())
        {
            surfaceScalarField SfUfbyDelta =
                mesh.surfaceInterpolation::deltaCoeffs()*mag(phi());

            CoNum = max(SfUfbyDelta/mesh.magSf())
                .value()*runTime().deltaT().value();

            meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf()))
                .value()*runTime().deltaT().value();

            velMag = max(mag(phi())/mesh.magSf()).value();
        }

        Info<< "Courant Number mean: " << meanCoNum
            << " max: " << CoNum
            << " velocity magnitude: " << velMag << endl;
    }

    fvVectorMatrix UEqn
    (
        fvm::ddt(U())
      + fvm::div(phi(), U())
      - fvm::laplacian(nu(), U())
    );

    solve(UEqn == -gradp());

    // --- PISO loop

//     volScalarField rUA = 1.0/UEqn.A();

    for (int corr=0; corr<nCorr; corr++)
    {
        volScalarField AU = UEqn.A();
        volVectorField HU = UEqn.H();

        U() = HU/AU;

//         U() = rUA*UEqn.H();
//         phi() = (fvc::interpolate(U())&mesh.Sf());

#       include "calcPhi.H"

        adjustPhi(phi(), U(), p());

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian
                (
                    1.0/fvc::interpolate(AU), p(), 
                    "laplacian((1|A(U)),p)"
                )
             == fvc::div(phi())
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phi() -= pEqn.flux();
            }
        }

        // Continuity error
        {
            volScalarField contErr = fvc::div(phi());

            scalar sumLocalContErr = runTime().deltaT().value()*
                mag(contErr)().weightedAverage(mesh.V()).value();

            scalar globalContErr = runTime().deltaT().value()*
                contErr.weightedAverage(mesh.V()).value();

            Info<< "time step continuity errors : sum local = "
                << sumLocalContErr << ", global = " << globalContErr << endl;
        }

        gradp() = fvc::grad(p());

        U() -= gradp()/AU;
        U().correctBoundaryConditions();
    }
}

const surfaceVectorField& consistentIcoFlow::Sf() const
{
    if (!SfPtr_)
    {
        makeSf();
    }

    return *SfPtr_;
}


surfaceVectorField& consistentIcoFlow::Sf()
{
    if (!SfPtr_)
    {
        makeSf();
    }

    return *SfPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowModels
} // End namespace Foam

// ************************************************************************* //
