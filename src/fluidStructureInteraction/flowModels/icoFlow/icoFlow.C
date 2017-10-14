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

#include "icoFlow.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "findRefCell.H"
#include "adjustPhi.H"
#include "fluidStructureInterface.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flowModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(icoFlow, 0);
addToRunTimeSelectionTable(flowModel, icoFlow, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

icoFlow::icoFlow(const fvMesh& mesh)
:
    flowModel(this->typeName, mesh),
    U_
    (
        IOobject
        (
            "U",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    p_
    (
        IOobject
        (
            "p",
            runTime().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    gradp_(fvc::grad(p_)),
    phi_
    (
        IOobject
        (
            "phi",
            runTime().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U_) & mesh.Sf()
    ),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            runTime().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nu_(transportProperties_.lookup("nu")),
    rho_(transportProperties_.lookup("rho"))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const volVectorField& icoFlow::U() const
{
    return U_;
}


const volScalarField& icoFlow::p() const
{
    return p_;
}


//- Patch viscous force (N/m2)
tmp<vectorField> icoFlow::patchViscousForce(const label patchID) const
{
    tmp<vectorField> tvF
    (
        new vectorField(mesh().boundary()[patchID].size(), vector::zero)
    );

    tvF() = rho_.value()*nu_.value()*U().boundaryField()[patchID].snGrad();

//     vectorField n = mesh().boundary()[patchID].nf();
//     tvF() -= n*(n&tvF());

    return tvF;
}

//- Patch pressure force (N/m2)
tmp<scalarField> icoFlow::patchPressureForce(const label patchID) const
{
    tmp<scalarField> tpF
    (
        new scalarField(mesh().boundary()[patchID].size(), 0)
    );

    tpF() = rho_.value()*p().boundaryField()[patchID];

    return tpF;
}

//- Patch viscous force (N/m2)
tmp<vectorField> icoFlow::faceZoneViscousForce
(
    const label zoneID,
    const label patchID
) const
{
    vectorField pVF = patchViscousForce(patchID);

    tmp<vectorField> tvF
    (
        new vectorField(mesh().faceZones()[zoneID].size(), vector::zero)
    );
    vectorField& vF = tvF();

    const label patchStart = 
        mesh().boundaryMesh()[patchID].start();

    forAll(pVF, i)
    {
        vF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] = 
            pVF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(vF, sumOp<vectorField>());
    

    return tvF;
}

//- Patch pressure force (N/m2)
tmp<scalarField> icoFlow::faceZonePressureForce
(
    const label zoneID, 
    const label patchID
) const
{
    scalarField pPF = patchPressureForce(patchID);

    tmp<scalarField> tpF
    (
        new scalarField(mesh().faceZones()[zoneID].size(), 0)
    );
    scalarField& pF = tpF();

    const label patchStart = 
        mesh().boundaryMesh()[patchID].start();

    forAll(pPF, i)
    {
        pF[mesh().faceZones()[zoneID].whichFace(patchStart + i)] = 
            pPF[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(pF, sumOp<scalarField>());

    return tpF;
}

tmp<scalarField> icoFlow::faceZoneMuEff
(
    const label zoneID,
    const label patchID                
) const
{
    tmp<scalarField> tMuEff
    (
        new scalarField
        (
            mesh().faceZones()[zoneID].size(), 
            rho_.value()*nu_.value()
        )
    );

    return tMuEff;
}

void icoFlow::evolve()
{
    Info << "Evolving flow model" << endl;

    const fvMesh& mesh = flowModel::mesh();

//     // Looking up fluid structure interface
//     const fluidStructureInterface& fsi =
//         mesh.lookupObject<fluidStructureInterface>
//         (
//             "fsiProperties"
//         );

//     if (nGradUn_.size() != U_.boundaryField()[fsi.fluidPatchIndex()].size())
//     {
//         nGradUn_ = 
//             scalarField(U_.boundaryField()[fsi.fluidPatchIndex()].size(), 0);
//     }

    int nCorr(readInt(flowProperties().lookup("nCorrectors")));

    int nNonOrthCorr =
        readInt(flowProperties().lookup("nNonOrthogonalCorrectors"));

    // Prepare for the pressure solution
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p_, flowProperties(), pRefCell, pRefValue);

    if(mesh.moving())
    {
        // Make the fluxes relative
        phi_ -= fvc::meshPhi(U_);
    }

    // CourantNo
    {
        scalar CoNum = 0.0;
        scalar meanCoNum = 0.0;
        scalar velMag = 0.0;

        if (mesh.nInternalFaces())
        {
            surfaceScalarField SfUfbyDelta =
                mesh.surfaceInterpolation::deltaCoeffs()*mag(phi_);

            CoNum = max(SfUfbyDelta/mesh.magSf())
                .value()*runTime().deltaT().value();

            meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf()))
                .value()*runTime().deltaT().value();

            velMag = max(mag(phi_)/mesh.magSf()).value();
        }

        Info<< "Courant Number mean: " << meanCoNum
            << " max: " << CoNum
            << " velocity magnitude: " << velMag << endl;
    }

    fvVectorMatrix UEqn
    (
        fvm::ddt(U_)
      + fvm::div(phi_, U_)
      - fvm::laplacian(nu_, U_)
    );

    solve(UEqn == -gradp_);

    // --- PISO loop

    volScalarField rUA = 1.0/UEqn.A();

    for (int corr=0; corr<nCorr; corr++)
    {
        U_ = rUA*UEqn.H();
        phi_ = (fvc::interpolate(U_) & mesh.Sf());

//             + fvc::ddtPhiCorr(rUA, U_, phi_);

// #       include "scaleInterfacePhi.H"

// #       include "setInterfacePressureGradient.H"

        adjustPhi(phi_, U_, p_);

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rUA, p_) == fvc::div(phi_)
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phi_ -= pEqn.flux();
            }
        }

        // Continuity error
        {
            volScalarField contErr = fvc::div(phi_);

            scalar sumLocalContErr = runTime().deltaT().value()*
                mag(contErr)().weightedAverage(mesh.V()).value();

            scalar globalContErr = runTime().deltaT().value()*
                contErr.weightedAverage(mesh.V()).value();

            Info<< "time step continuity errors : sum local = "
                << sumLocalContErr << ", global = " << globalContErr << endl;
        }

        gradp_ = fvc::grad(p_);

        U_ -= rUA*gradp_;
        U_.correctBoundaryConditions();

// #       include "updateInterfaceNGradUn.H"
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowModels
} // End namespace Foam

// ************************************************************************* //
