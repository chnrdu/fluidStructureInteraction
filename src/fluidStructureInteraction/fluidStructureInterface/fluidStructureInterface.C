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

#include "fluidStructureInterface.H"
#include "volFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"
#include "RectangularMatrix.H"
#include "primitivePatchInterpolation.H"
#include "twoDPointCorrector.H"

#include "tetPointFields.H"
#include "fixedValueTetPolyPatchFields.H"
#include "tetPolyPatchInterpolation.H"
#include "tetFemMatrices.H"

#include "fixedValuePointPatchFields.H"
#include "RBFMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidStructureInterface, 0);
//     defineRunTimeSelectionTable(fluidStructureInterface, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fluidStructureInterface::calcCurrentSolidZonePoints() const
{
    // Find global face zones
    if (currentSolidZonePointsPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcCurrentSolidZonePoints() const"
        )
            << "Current solid zone points alarady exist"
                << abort(FatalError);
    }

    currentSolidZonePointsPtr_ =
        new vectorField(stress().currentFaceZonePoints(solidZoneIndex()));
}

void Foam::fluidStructureInterface::calcCurrentSolidZone2Points() const //zone2
{
    // Find global face zones
    if (currentSolidZone2PointsPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcCurrentSolidZone2Points() const"
        )
            << "Current solid zone2 points alarady exist"
                << abort(FatalError);
    }

    currentSolidZone2PointsPtr_ =
        new vectorField(stress().currentFaceZonePoints(solidZone2Index()));
}

void Foam::fluidStructureInterface::calcCurrentSolidZonePatch() const
{
    // Find global face zones
    if (currentSolidZonePatchPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcCurrentSolidZonePatch() const"
        )
            << "Current solid zone patch alarady exists"
                << abort(FatalError);
    }

    currentSolidZonePatchPtr_ =
        new PrimitivePatch<face, List, const pointField&>
        (
            solidMesh().faceZones()[solidZoneIndex_]().localFaces(),
            currentSolidZonePoints()
        );
}

void Foam::fluidStructureInterface::calcCurrentSolidZone2Patch() const //zones
{
    // Find global face zones
    if (currentSolidZone2PatchPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcCurrentSolidZone2Patch() const"
        )
            << "Current solid zone2 patch alarady exists"
                << abort(FatalError);
    }

    currentSolidZone2PatchPtr_ =
        new PrimitivePatch<face, List, const pointField&>
        (
            solidMesh().faceZones()[solidZone2Index_]().localFaces(),
            currentSolidZone2Points()
        );
}
void Foam::fluidStructureInterface::calcFluidToSolidInterpolator() const
{
    // Find global face zones
    if (fluidToSolidPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcFluidToSolidInterpolator() const"
        )
            << "Fluid to solid interpolator already exists"
                << abort(FatalError);
    }

    fluidToSolidPtr_ =
        new zoneToZoneInterpolation
        (
            fluidMesh().faceZones()[fluidZoneIndex_](),
            solidMesh().faceZones()[solidZoneIndex_](),
            intersection::VISIBLE
        );


    Info << "Checking fluid-to-solid interpolator" << endl;
    {
        vectorField fluidPatchFaceCentres =
            vectorField
            (
                fluidMesh().boundaryMesh()[fluidPatchIndex_].faceCentres()
            );

        vectorField fluidZoneFaceCentres
        (
            fluidMesh().faceZones()[fluidZoneIndex_].size(),
            vector::zero
        );

        const label fluidPatchStart =
            fluidMesh().boundaryMesh()[fluidPatchIndex_].start();

        forAll (fluidPatchFaceCentres, i)
        {
            fluidZoneFaceCentres
            [
                fluidMesh().faceZones()[fluidZoneIndex_].whichFace
                (
                    fluidPatchStart + i
                )
            ] =
                fluidPatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(fluidZoneFaceCentres, sumOp<vectorField>());

        vectorField solidZoneFaceCentres =
            fluidToSolidPtr_->faceInterpolate
            (
                fluidZoneFaceCentres
            );

        vectorField solidPatchFaceCentres
        (
            solidMesh().boundaryMesh()[solidPatchIndex_].size(),
            vector::zero
        );

        const label solidPatchStart =
            solidMesh().boundaryMesh()[solidPatchIndex_].start();

        forAll(solidPatchFaceCentres, i)
        {
            solidPatchFaceCentres[i] =
                solidZoneFaceCentres
                [
                    solidMesh().faceZones()[solidZoneIndex_]
                   .whichFace(solidPatchStart + i)
                ];
        }

        scalar maxDist = gMax
        (
            mag
            (
                solidPatchFaceCentres
              - solidMesh().boundaryMesh()[solidPatchIndex_].faceCentres()
            )
        );

        Info << "Fluid-to-solid face interpolation error: " << maxDist
            << endl;
    }
}


void Foam::fluidStructureInterface::calcFluid2ToSolidInterpolator() const
{
    // Find global face zones
    if (fluid2ToSolidPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcFluid2ToSolidInterpolator() const"
        )
            << "Fluid2 to solid interpolator already exists"
                << abort(FatalError);
    }

    fluid2ToSolidPtr_ =
        new zoneToZoneInterpolation
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_](),
            solidMesh().faceZones()[solidZone2Index_](),
            intersection::VISIBLE
        );


    Info << "Checking fluid2-to-solid interpolator" << endl;
    {
        vectorField fluid2PatchFaceCentres =
            vectorField
            (
                fluidMesh2().boundaryMesh()[fluid2PatchIndex_].faceCentres()
            );

        vectorField fluid2ZoneFaceCentres
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_].size(),
            vector::zero
        );

        const label fluid2PatchStart =
            fluidMesh2().boundaryMesh()[fluid2PatchIndex_].start();

        forAll (fluid2PatchFaceCentres, i)
        {
            fluid2ZoneFaceCentres
            [
                fluidMesh2().faceZones()[fluid2ZoneIndex_].whichFace
                (
                    fluid2PatchStart + i
                )
            ] =
                fluid2PatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(fluid2ZoneFaceCentres, sumOp<vectorField>());

        vectorField solidZone2FaceCentres =
            fluid2ToSolidPtr_->faceInterpolate
            (
                fluid2ZoneFaceCentres
            );

        vectorField solidPatch2FaceCentres
        (
            solidMesh().boundaryMesh()[solidPatch2Index_].size(),
            vector::zero
        );

        const label solidPatch2Start =
            solidMesh().boundaryMesh()[solidPatch2Index_].start();

        forAll(solidPatch2FaceCentres, i)
        {
            solidPatch2FaceCentres[i] =
                solidZone2FaceCentres
                [
                    solidMesh().faceZones()[solidZone2Index_]
                   .whichFace(solidPatch2Start + i)
                ];
        }

        scalar maxDist = gMax
        (
            mag
            (
                solidPatch2FaceCentres
              - solidMesh().boundaryMesh()[solidPatch2Index_].faceCentres()
            )
        );

        Info << "Fluid2-to-solid face interpolation error: " << maxDist
            << endl;
    }
}

void Foam::fluidStructureInterface::calcGgiFluidToSolidInterpolator() const
{
    // Find global face zones
    if (ggiFluidToSolidPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcGgiFluidToSolidInterpolator() const"
        )
            << "Ggi fluid to solid interpolator already exists"
                << abort(FatalError);
    }

    ggiFluidToSolidPtr_ =
        new ggiZoneInterpolation
        (
            fluidMesh().faceZones()[fluidZoneIndex_](),
            solidMesh().faceZones()[solidZoneIndex_](),
            tensorField(0),
            tensorField(0),
            vectorField(0), // Slave-to-master separation. Bug fix
            0,              // Non-overlapping face tolerances
            0,              // HJ, 24/Oct/2008
            true,           // Rescale weighting factors.  Bug fix, MB.
            ggiInterpolation::AABB //BB_OCTREE  // Octree search, MB.
        );


    Info << "Checking fluid-to-solid interpolator" << endl;
    {
        vectorField fluidPatchFaceCentres =
            vectorField
            (
                fluidMesh().boundaryMesh()[fluidPatchIndex_].faceCentres()
            );

        vectorField fluidZoneFaceCentres
        (
            fluidMesh().faceZones()[fluidZoneIndex_].size(),
            vector::zero
        );

        const label fluidPatchStart =
            fluidMesh().boundaryMesh()[fluidPatchIndex_].start();

        forAll (fluidPatchFaceCentres, i)
        {
            fluidZoneFaceCentres
            [
                fluidMesh().faceZones()[fluidZoneIndex_].whichFace
                (
                    fluidPatchStart + i
                )
            ] =
                fluidPatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(fluidZoneFaceCentres, sumOp<vectorField>());

        vectorField solidZoneFaceCentres =
            ggiFluidToSolidPtr_->masterToSlave
            (
                fluidZoneFaceCentres
            );

        vectorField solidPatchFaceCentres
        (
            solidMesh().boundaryMesh()[solidPatchIndex_].size(),
            vector::zero
        );

        const label solidPatchStart =
            solidMesh().boundaryMesh()[solidPatchIndex_].start();

        forAll(solidPatchFaceCentres, i)
        {
            solidPatchFaceCentres[i] =
                solidZoneFaceCentres
                [
                    solidMesh().faceZones()[solidZoneIndex_]
                   .whichFace(solidPatchStart + i)
                ];
        }

        scalar maxDist = gMax
        (
            mag
            (
                solidPatchFaceCentres
              - solidMesh().boundaryMesh()[solidPatchIndex_].faceCentres()
            )
        );

        Info << "Fluid-to-solid face interpolation error: " << maxDist
            << endl;
    }

    Info << "Number of uncovered master faces: "
        << ggiFluidToSolidPtr_->uncoveredMasterFaces().size() << endl;

    Info << "Number of uncovered slave faces: "
        << ggiFluidToSolidPtr_->uncoveredSlaveFaces().size() << endl;
}

void Foam::fluidStructureInterface::calcGgiFluid2ToSolidInterpolator() const
{
    // Find global face zones
    if (ggiFluid2ToSolidPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcGgiFluid2ToSolidInterpolator() const"
        )
            << "Ggi fluid2 to solid interpolator already exists"
                << abort(FatalError);
    }

    ggiFluid2ToSolidPtr_ =
        new ggiZoneInterpolation
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_](),
            solidMesh().faceZones()[solidZone2Index_](),
            tensorField(0),
            tensorField(0),
            vectorField(0), // Slave-to-master separation. Bug fix
            0,              // Non-overlapping face tolerances
            0,              // HJ, 24/Oct/2008
            true,           // Rescale weighting factors.  Bug fix, MB.
            ggiInterpolation::AABB //BB_OCTREE  // Octree search, MB.
        );


    Info << "Checking fluid2-to-solid interpolator" << endl;
    {
        vectorField fluid2PatchFaceCentres =
            vectorField
            (
                fluidMesh2().boundaryMesh()[fluid2PatchIndex_].faceCentres()
            );

        vectorField fluid2ZoneFaceCentres
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_].size(),
            vector::zero
        );

        const label fluid2PatchStart =
            fluidMesh2().boundaryMesh()[fluid2PatchIndex_].start();

        forAll (fluid2PatchFaceCentres, i)
        {
            fluid2ZoneFaceCentres
            [
                fluidMesh2().faceZones()[fluid2ZoneIndex_].whichFace
                (
                    fluid2PatchStart + i
                )
            ] =
                fluid2PatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(fluid2ZoneFaceCentres, sumOp<vectorField>());

        vectorField solidZone2FaceCentres =
            ggiFluid2ToSolidPtr_->masterToSlave
            (
                fluid2ZoneFaceCentres
            );

        vectorField solidPatch2FaceCentres
        (
            solidMesh().boundaryMesh()[solidPatch2Index_].size(),
            vector::zero
        );

        const label solidPatch2Start =
            solidMesh().boundaryMesh()[solidPatch2Index_].start();

        forAll(solidPatch2FaceCentres, i)
        {
            solidPatch2FaceCentres[i] =
                solidZone2FaceCentres
                [
                    solidMesh().faceZones()[solidZone2Index_]
                   .whichFace(solidPatch2Start + i)
                ];
        }

        scalar maxDist = gMax
        (
            mag
            (
                solidPatch2FaceCentres
              - solidMesh().boundaryMesh()[solidPatch2Index_].faceCentres()
            )
        );

        Info << "Fluid2-to-solid face interpolation error: " << maxDist
            << endl;
    }

    Info << "Number of uncovered master faces: "
        << ggiFluid2ToSolidPtr_->uncoveredMasterFaces().size() << endl;

    Info << "Number of uncovered slave faces: "
        << ggiFluid2ToSolidPtr_->uncoveredSlaveFaces().size() << endl;
}

void Foam::fluidStructureInterface::calcGgiInterpolator() const
{
    // Create extended ggi interpolation
    if (ggiInterpolatorPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcGgiInterpolator() const"
        )
            << "Ggi interpolator already exists"
                << abort(FatalError);
    }

    // Create copy of solid face zone primitive patch in current configuration

    deleteDemandDrivenData(currentSolidZonePatchPtr_);
    deleteDemandDrivenData(currentSolidZonePointsPtr_);

//     currentSolidZonePatch().movePoints(currentSolidZonePoints());

    Info << "Create extended GGI zone-to-zone interpolator" << endl;

    currentSolidZonePatch(),

    ggiInterpolatorPtr_ =
        new extendedGgiZoneInterpolation
        (
            fluidMesh().faceZones()[fluidZoneIndex_](),
            currentSolidZonePatch(),
//             solidMesh().faceZones()[solidZoneIndex_](),
            tensorField(0),
            tensorField(0),
            vectorField(0), // Slave-to-master separation. Bug fix
            0,              // Non-overlapping face tolerances
            0,              // HJ, 24/Oct/2008
            true,           // Rescale weighting factors.  Bug fix, MB.
            ggiInterpolation::AABB
            // N_SQUARED BB_OCTREE AABB THREE_D_DISTANCE
            // Octree search, MB.
        );


    Info << "Checking fluid-to-solid face interpolator" << endl;
    {
        vectorField fluidPatchFaceCentres =
            vectorField
            (
                fluidMesh().boundaryMesh()[fluidPatchIndex_].faceCentres()
            );

        vectorField fluidZoneFaceCentres
        (
            fluidMesh().faceZones()[fluidZoneIndex_].size(),
            vector::zero
        );

        const label fluidPatchStart =
            fluidMesh().boundaryMesh()[fluidPatchIndex_].start();

        forAll (fluidPatchFaceCentres, i)
        {
            fluidZoneFaceCentres
            [
                fluidMesh().faceZones()[fluidZoneIndex_].whichFace
                (
                    fluidPatchStart + i
                )
            ] =
                fluidPatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(fluidZoneFaceCentres, sumOp<vectorField>());

        vectorField solidZoneFaceCentres =
            ggiInterpolatorPtr_->masterToSlave
            (
                fluidZoneFaceCentres
            );

        vectorField solidPatchFaceCentres
        (
            solidMesh().boundaryMesh()[solidPatchIndex_].size(),
            vector::zero
        );

        const label solidPatchStart =
            solidMesh().boundaryMesh()[solidPatchIndex_].start();

        forAll(solidPatchFaceCentres, i)
        {
            solidPatchFaceCentres[i] =
                solidZoneFaceCentres
                [
                    solidMesh().faceZones()[solidZoneIndex_]
                   .whichFace(solidPatchStart + i)
                ];
        }

        scalar maxDist = gMax
        (
            mag
            (
                solidPatchFaceCentres
              - solidMesh().boundaryMesh()[solidPatchIndex_].faceCentres()
            )
        );

        Info << "Fluid-to-solid face interpolation error: " << maxDist
            << endl;
    }

    Info << "Checking solid-to-fluid point interpolator (GGI)" << endl;
    {
        vectorField solidZonePoints_ =
            solidMesh().faceZones()[solidZoneIndex_]().localPoints();

        vectorField solidZonePoints =
            ggiInterpolatorPtr_->slaveToMasterPointInterpolate
            (
                solidZonePoints_
            );

        vectorField fluidZonePoints =
            fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

        scalar maxDist = gMax
        (
            mag
            (
                fluidZonePoints
              - solidZonePoints
            )
        );

        Info << "Solid-to-fluid point interpolation error (GGI): " << maxDist
            << endl;
    }

    Info << "Number of uncovered master faces: "
        << ggiInterpolatorPtr_->uncoveredMasterFaces().size() << endl;

    Info << "Number of uncovered slave faces: "
        << ggiInterpolatorPtr_ ->uncoveredSlaveFaces().size() << endl;
}

void Foam::fluidStructureInterface::calcGgiInterpolator2() const
{
    // Create extended ggi interpolation
    if (ggiInterpolator2Ptr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcGgiInterpolator2() const"
        )
            << "Ggi interpolator already exists"
                << abort(FatalError);
    }

    // Create copy of solid face zone primitive patch in current configuration

    deleteDemandDrivenData(currentSolidZonePatchPtr_);
    deleteDemandDrivenData(currentSolidZonePointsPtr_);
    deleteDemandDrivenData(currentSolidZone2PatchPtr_);
    deleteDemandDrivenData(currentSolidZone2PointsPtr_);

//     currentSolidZonePatch().movePoints(currentSolidZonePoints());

    Info << "Create extended GGI zone-to-zone interpolator" << endl;

    currentSolidZone2Patch(),

    ggiInterpolator2Ptr_ =
        new extendedGgiZoneInterpolation
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_](),
            currentSolidZone2Patch(),
//             solidMesh().faceZones()[solidZoneIndex_](),
            tensorField(0),
            tensorField(0),
            vectorField(0), // Slave-to-master separation. Bug fix
            0,              // Non-overlapping face tolerances
            0,              // HJ, 24/Oct/2008
            true,           // Rescale weighting factors.  Bug fix, MB.
            ggiInterpolation::AABB
            // N_SQUARED BB_OCTREE AABB THREE_D_DISTANCE
            // Octree search, MB.
        );


    Info << "Checking fluid2-to-solid face interpolator" << endl;
    {
        vectorField fluid2PatchFaceCentres =
            vectorField
            (
                fluidMesh2().boundaryMesh()[fluid2PatchIndex_].faceCentres()
            );

        vectorField fluid2ZoneFaceCentres
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_].size(),
            vector::zero
        );

        const label fluid2PatchStart =
            fluidMesh2().boundaryMesh()[fluid2PatchIndex_].start();

        forAll (fluid2PatchFaceCentres, i)
        {
            fluid2ZoneFaceCentres
            [
                fluidMesh2().faceZones()[fluid2ZoneIndex_].whichFace
                (
                    fluid2PatchStart + i
                )
            ] =
                fluid2PatchFaceCentres[i];
        }

        // Parallel data exchange: collect faceCentres field on all processors
        reduce(fluid2ZoneFaceCentres, sumOp<vectorField>());

        vectorField solidZone2FaceCentres =
            ggiInterpolator2Ptr_->masterToSlave
            (
                fluid2ZoneFaceCentres
            );

        vectorField solidPatch2FaceCentres
        (
            solidMesh().boundaryMesh()[solidPatch2Index_].size(),
            vector::zero
        );

        const label solidPatch2Start =
            solidMesh().boundaryMesh()[solidPatch2Index_].start();

        forAll(solidPatch2FaceCentres, i)
        {
            solidPatch2FaceCentres[i] =
                solidZone2FaceCentres
                [
                    solidMesh().faceZones()[solidZone2Index_]
                   .whichFace(solidPatch2Start + i)
                ];
        }

        scalar maxDist = gMax
        (
            mag
            (
                solidPatch2FaceCentres
              - solidMesh().boundaryMesh()[solidPatch2Index_].faceCentres()
            )
        );

        Info << "Fluid2-to-solid face interpolation error: " << maxDist
            << endl;
    }

    Info << "Checking solid-to-fluid2 point interpolator (GGI)" << endl;
    {
        vectorField solidZone2Points_ =
            solidMesh().faceZones()[solidZone2Index_]().localPoints();

        vectorField solidZone2Points =
            ggiInterpolator2Ptr_->slaveToMasterPointInterpolate
            (
                solidZone2Points_
            );

        vectorField fluid2ZonePoints =
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().localPoints();

        scalar maxDist = gMax
        (
            mag
            (
                fluid2ZonePoints
              - solidZone2Points
            )
        );

        Info << "Solid-to-fluid2 point interpolation error (GGI): " << maxDist
            << endl;
    }

    Info << "Number of uncovered master faces: "
        << ggiInterpolator2Ptr_->uncoveredMasterFaces().size() << endl;

    Info << "Number of uncovered slave faces: "
        << ggiInterpolator2Ptr_ ->uncoveredSlaveFaces().size() << endl;
}

void Foam::fluidStructureInterface::calcSolidToFluidInterpolator() const
{
    // Find global face zones
    if (solidToFluidPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcSolidToFluidInterpolator() const"
        )
            << "Solid to fluid interpolator already exists"
                << abort(FatalError);
    }

    solidToFluidPtr_ =
        new zoneToZoneInterpolation
        (
            solidMesh().faceZones()[solidZoneIndex_](),
            fluidMesh().faceZones()[fluidZoneIndex_](),
            intersection::VISIBLE
        );

    Info << "Checking solid-to-fluid interpolator" << endl;
    {
        vectorField solidZonePoints_ =
            solidMesh().faceZones()[solidZoneIndex_]().localPoints();

        vectorField solidZonePoints =
            solidToFluidPtr_->pointInterpolate
            (
                solidZonePoints_
            );

        vectorField fluidZonePoints =
            fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

        scalar maxDist = gMax
        (
            mag
            (
                fluidZonePoints
              - solidZonePoints
            )
        );

        Info << "Solid-to-fluid point interpolation error: " << maxDist
            << endl;
    }
}

void Foam::fluidStructureInterface::calcSolidToFluid2Interpolator() const
{
    // Find global face zones
    if (solidToFluid2Ptr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcSolidToFluid2Interpolator() const"
        )
            << "Solid to fluid2 interpolator already exists"
                << abort(FatalError);
    }

    solidToFluid2Ptr_ =
        new zoneToZoneInterpolation
        (
            solidMesh().faceZones()[solidZone2Index_](),
            fluidMesh2().faceZones()[fluid2ZoneIndex_](),
            intersection::VISIBLE
        );

    Info << "Checking solid-to-fluid2 interpolator" << endl;
    {
        vectorField solidZone2Points_ =
            solidMesh().faceZones()[solidZone2Index_]().localPoints();

        vectorField solidZone2Points =
            solidToFluidPtr_->pointInterpolate
            (
                solidZone2Points_
            );

        vectorField fluid2ZonePoints =
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().localPoints();

        scalar maxDist = gMax
        (
            mag
            (
                fluid2ZonePoints
              - solidZone2Points
            )
        );

        Info << "Solid-to-fluid2 point interpolation error: " << maxDist
            << endl;
    }
}


void Foam::fluidStructureInterface::
calcAccumulatedFluidInterfaceDisplacement() const
{
    // Read accumulated displacement
    if (accumulatedFluidInterfaceDisplacementPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcAccumulatedFluidInterfaceDisplacement() const"
        )
            << "Accumulated displacement field already exists"
                << abort(FatalError);
    }

    // Accumulated fluid interface displacement
    IOobject accumulatedFluidInterfaceDisplacementHeader
    (
        "accumulatedFluidInterfaceDisplacement",
        flow().runTime().timeName(),
        fluidMesh(),
        IOobject::MUST_READ
    );

    if(accumulatedFluidInterfaceDisplacementHeader.headerOk())
    {
        Pout << "Reading accumulated fluid interface displacement" << endl;

        accumulatedFluidInterfaceDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "accumulatedFluidInterfaceDisplacement",
                    flow().runTime().timeName(),
                    fluidMesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
    else
    {
        Pout << "Creating accumulated fluid interface displacement" << endl;

        accumulatedFluidInterfaceDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "accumulatedFluidInterfaceDisplacement",
                    flow().runTime().timeName(),
                    fluidMesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                vectorField
                (
                    fluidMesh().boundaryMesh()[fluidPatchIndex()].nPoints(),
                    vector::zero
                )
            );
    }
}


void Foam::fluidStructureInterface::
calcAccumulatedFluid2InterfaceDisplacement() const
{
    // Read accumulated displacement
    if (accumulatedFluid2InterfaceDisplacementPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcAccumulatedFluid2InterfaceDisplacement() const"
        )
            << "Accumulated displacement field already exists"
                << abort(FatalError);
    }

    // Accumulated fluid interface displacement
    IOobject accumulatedFluid2InterfaceDisplacementHeader
    (
        "accumulatedFluid2InterfaceDisplacement",
        flow2().runTime().timeName(),
        fluidMesh2(),
        IOobject::MUST_READ
    );

    if(accumulatedFluid2InterfaceDisplacementHeader.headerOk())
    {
        Pout << "Reading accumulated fluid2 interface displacement" << endl;

        accumulatedFluid2InterfaceDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "accumulatedFluid2InterfaceDisplacement",
                    flow2().runTime().timeName(),
                    fluidMesh2(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
    else
    {
        Pout << "Creating accumulated fluid2 interface displacement" << endl;

        accumulatedFluid2InterfaceDisplacementPtr_ =
            new vectorIOField
            (
                IOobject
                (
                    "accumulatedFluid2InterfaceDisplacement",
                    flow2().runTime().timeName(),
                    fluidMesh2(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                vectorField
                (
                    fluidMesh2().boundaryMesh()[fluid2PatchIndex()].nPoints(),
                    vector::zero
                )
            );
    }
}

void Foam::fluidStructureInterface::calcMinEdgeLength() const
{
    // Read accumulated displacement
    if (minEdgeLengthPtr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcMinEdgeLength() const"
        )
            << "Minimal edge lengths already exist"
                << abort(FatalError);
    }

    minEdgeLengthPtr_ =
        new scalarField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            0
        );
    scalarField& minEdgeLength = *minEdgeLengthPtr_;


    const edgeList& edges =
        fluidMesh().faceZones()[fluidZoneIndex_]().edges();

    const vectorField& points =
        fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

    const labelListList& pointEdges =
        fluidMesh().faceZones()[fluidZoneIndex_]().pointEdges();

    forAll(points, pointI)
    {
        const labelList& curPointEdges = pointEdges[pointI];

        scalar minLength = GREAT;

        forAll(curPointEdges, edgeI)
        {
            const edge& curEdge = edges[curPointEdges[edgeI]];

            scalar Le = curEdge.mag(points);

            if (Le < minLength)
            {
                minLength = Le;
            }
        }

        minEdgeLength[pointI] = minLength;
    }
//     Pout << "Min edge length: " << min(minEdgeLength) << endl;
//     Pout << "gMin edge length: " << gMin(minEdgeLength) << endl;
}

void Foam::fluidStructureInterface::calcMinEdgeLength2() const
{
    // Read accumulated displacement
    if (minEdgeLength2Ptr_)
    {
        FatalErrorIn
        (
            "void fluidStructureInterface::"
            "calcMinEdgeLength2() const"
        )
            << "Minimal edge lengths already exist"
                << abort(FatalError);
    }

    minEdgeLength2Ptr_ =
        new scalarField
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().nPoints(),
            0
        );
    scalarField& minEdgeLength2 = *minEdgeLength2Ptr_;


    const edgeList& edges =
        fluidMesh2().faceZones()[fluid2ZoneIndex_]().edges();

    const vectorField& points =
        fluidMesh2().faceZones()[fluid2ZoneIndex_]().localPoints();

    const labelListList& pointEdges =
        fluidMesh2().faceZones()[fluid2ZoneIndex_]().pointEdges();

    forAll(points, pointI)
    {
        const labelList& curPointEdges = pointEdges[pointI];

        scalar minLength = GREAT;

        forAll(curPointEdges, edgeI)
        {
            const edge& curEdge = edges[curPointEdges[edgeI]];

            scalar Le = curEdge.mag(points);

            if (Le < minLength)
            {
                minLength = Le;
            }
        }

        minEdgeLength2[pointI] = minLength;
    }
//     Pout << "Min edge length: " << min(minEdgeLength) << endl;
//     Pout << "gMin edge length: " << gMin(minEdgeLength) << endl;
}

Foam::vectorIOField&
Foam::fluidStructureInterface::accumulatedFluidInterfaceDisplacement()
{
    if (!accumulatedFluidInterfaceDisplacementPtr_)
    {
        calcAccumulatedFluidInterfaceDisplacement();
    }

    return *accumulatedFluidInterfaceDisplacementPtr_;
}

Foam::vectorIOField&
Foam::fluidStructureInterface::accumulatedFluid2InterfaceDisplacement()
{
    if (!accumulatedFluid2InterfaceDisplacementPtr_)
    {
        calcAccumulatedFluid2InterfaceDisplacement();
    }

    return *accumulatedFluid2InterfaceDisplacementPtr_;
}

const Foam::scalarField& Foam::fluidStructureInterface::minEdgeLength() const
{
    if (!minEdgeLengthPtr_)
    {
        calcMinEdgeLength();
    }

    return *minEdgeLengthPtr_;
}
const Foam::scalarField& Foam::fluidStructureInterface::minEdgeLength2() const
{
    if (!minEdgeLength2Ptr_)
    {
        calcMinEdgeLength2();
    }

    return *minEdgeLength2Ptr_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidStructureInterface::fluidStructureInterface
(
    dynamicFvMesh& fMesh,
    dynamicFvMesh& fMesh2, //add mesh2
    fvMesh& sMesh
)
:
    IOdictionary
    (
        IOobject
        (
            "fsiProperties",
            fMesh.time().constant(),
            fMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )

    ),
    fluidMesh_(fMesh),
    flow_(flowModel::New(fluidMesh_)),
    fluidMesh2_(fMesh2),
    flow2_(flowModel::New(fluidMesh2_)), // initial fluid2
    solidMesh_(sMesh),
    stress_(stressModel::New(solidMesh_)),
    solidPatchIndex_(-1),
    solidZoneIndex_(-1),
    solidPatch2Index_(-1),
    solidZone2Index_(-1),
    fluidPatchIndex_(-1),
    fluidZoneIndex_(-1),
    fluid2PatchIndex_(-1),
    fluid2ZoneIndex_(-1), //fluid2 index
    currentSolidZonePointsPtr_(NULL),
    currentSolidZonePatchPtr_(NULL),
    currentSolidZone2PointsPtr_(NULL), //zone2
    currentSolidZone2PatchPtr_(NULL),
    fluidToSolidPtr_(NULL),
    fluid2ToSolidPtr_(NULL), //fluid2 ptr
    ggiFluidToSolidPtr_(NULL),
    ggiFluid2ToSolidPtr_(NULL), //ggifluid2 ptr
    ggiInterpolatorPtr_(NULL),
    ggiInterpolator2Ptr_(NULL),
    solidToFluidPtr_(NULL),
    solidToFluid2Ptr_(NULL), //solid fluid2 ptr
    couplingScheme_(lookup("couplingScheme")),
    relaxationFactor_(readScalar(lookup("relaxationFactor"))),
    aitkenRelaxationFactor_(relaxationFactor_),
    outerCorrTolerance_(readScalar(lookup("outerCorrTolerance"))),
    nOuterCorr_(readInt(lookup("nOuterCorr"))),
    coupled_(lookup("coupled")),
    predictor_(false), //(lookup("predictor")),
    couplingReuse_(readInt(lookup("couplingReuse"))),
    interfaceDeformationLimit_
    (
        readScalar(lookup("interfaceDeformationLimit"))
    ),
    fluidZonePointsDispl_(),
    fluidZonePointsDisplRef_(),
    fluidZonePointsDisplPrev_(),
    fluid2ZonePointsDispl_(),
    fluid2ZonePointsDisplRef_(),
    fluid2ZonePointsDisplPrev_(), //fluid2 displace
    solidZonePointsDispl_(),
    solidZonePointsDisplRef_(),
    solidZonePressure_(),
    solidZoneTraction_(),
    solidZoneTractionPrev_(),
    predictedSolidZoneTraction_(),
    solidZone2PointsDispl_(),  //zone2
    solidZone2PointsDisplRef_(),
    solidZone2Pressure_(),
    solidZone2Traction_(),
    solidZone2TractionPrev_(),
    predictedSolidZone2Traction_(),
    residual_(),
    residual2_(),
    residualPrev_(),
    residual2Prev_(),
    maxResidualNorm_(0),
    outerCorr_(0),

//     closedFluidDomain_(lookup("closedFluidDomain")),
//     refPressure_(0),
//     refPressureIncrement_(0),
//     compressibility_(readScalar(lookup("compressibility"))),
    interpolatorUpdateFrequency_
    (
        readInt(lookup("interpolatorUpdateFrequency"))
    ),
    fluidPatchPointsV_(),
    fluidPatchPointsW_(),
    fluidPatchPointsT_(),
    fluid2PatchPointsV_(), //add fluid2 patch
    fluid2PatchPointsW_(),
    fluid2PatchPointsT_(),
    accumulatedFluidInterfaceDisplacementPtr_(NULL),
    accumulatedFluid2InterfaceDisplacementPtr_(NULL),
    minEdgeLengthPtr_(NULL),
    minEdgeLength2Ptr_(NULL)
{
    // Solid patch index

    word solidPatchName(lookup("solidPatch"));

    polyPatchID solidPatch(solidPatchName, solidMesh().boundaryMesh());

    if (!solidPatch.active())
    {
        FatalErrorIn("fluidStructureInterface::fluidStructureInterface(...)")
            << "Solid patch name " << solidPatchName << " not found."
                << abort(FatalError);
    }

    solidPatchIndex_ = solidPatch.index();


    // Solid face zone index

    word solidZoneName(lookup("solidZone"));

    faceZoneID solidZone
    (
        solidZoneName,
        solidMesh().faceZones()
    );

    if (!solidZone.active())
    {
        FatalErrorIn("")
            << "Solid face zone name " << solidZoneName
                << " not found.  Please check your face zone definition."
                << abort(FatalError);
    }

    solidZoneIndex_ = solidZone.index();

    // Solid patch2 index

    word solidPatch2Name(lookup("solidPatch2")); //zone2

    polyPatchID solidPatch2(solidPatch2Name, solidMesh().boundaryMesh());

    if (!solidPatch2.active())
    {
        FatalErrorIn("fluidStructureInterface::fluidStructureInterface(...)")
            << "Solid patch2 name " << solidPatch2Name << " not found."
                << abort(FatalError);
    }

    solidPatch2Index_ = solidPatch2.index();


    // Solid face zone2 index

    word solidZone2Name(lookup("solidZone2"));

    faceZoneID solidZone2
    (
        solidZone2Name,
        solidMesh().faceZones()
    );

    if (!solidZone2.active())
    {
        FatalErrorIn("")
            << "Solid face zone2 name " << solidZone2Name
                << " not found.  Please check your face zone definition."
                << abort(FatalError);
    }

    solidZone2Index_ = solidZone2.index();
    // Fluid patch index

    word fluidPatchName(lookup("fluidPatch"));

    polyPatchID fluidPatch(fluidPatchName, fluidMesh().boundaryMesh());

    if (!fluidPatch.active())
    {
        FatalErrorIn("fluidStructureInterface::fluidStructureInterface(...)")
            << "Fluid patch name " << fluidPatchName << " not found."
                << abort(FatalError);
    }

    fluidPatchIndex_ = fluidPatch.index();


    // Fluid face zone index

    word fluidZoneName(lookup("fluidZone"));

    faceZoneID fluidZone
    (
        fluidZoneName,
        fluidMesh().faceZones()
    );

    if (!fluidZone.active())
    {
        FatalErrorIn("fluidStructureInterface::fluidStructureInterface(...)")
            << "Fluid face zone name " << fluidZoneName
                << " not found.  Please check your face zone definition."
                << abort(FatalError);
    }

    fluidZoneIndex_ = fluidZone.index();

    // Fluid2 patch index

    word fluid2PatchName(lookup("fluid2Patch"));

    polyPatchID fluid2Patch(fluid2PatchName, fluidMesh2().boundaryMesh());

    if (!fluid2Patch.active())
    {
        FatalErrorIn("fluidStructureInterface::fluidStructureInterface(...)")
            << "Fluid2 patch name " << fluid2PatchName << " not found."
                << abort(FatalError);
    }

    fluid2PatchIndex_ = fluid2Patch.index();


    // Fluid2 face zone index

    word fluid2ZoneName(lookup("fluid2Zone"));

    faceZoneID fluid2Zone
    (
        fluid2ZoneName,
        fluidMesh2().faceZones()
    );

    if (!fluid2Zone.active())
    {
        FatalErrorIn("fluidStructureInterface::fluidStructureInterface(...)")
            << "Fluid2 face zone name " << fluid2ZoneName
                << " not found.  Please check your face zone definition."
                << abort(FatalError);
    }

    fluid2ZoneIndex_ = fluid2Zone.index();


    // Check coupling scheme
    if
    (
        (couplingScheme_ == "IQN-ILS")
     || (couplingScheme_ == "Aitken")
     || (couplingScheme_ == "FixedRelaxation")
    )
    {
        Info<< "Selecting coupling scheme " << couplingScheme_ << endl;
    }
    else
    {
        FatalErrorIn
        (
            "fluidStructureInterface::fluidStructureInterface(...)"
        )   << "couplingScheme: " << couplingScheme_
            << " is not a valid choice. "
            << "Options are: IQN-ILS, Aitken, FixedRelaxation"
            << abort(FatalError);
    }

    // Initialize solid zone pressure
    solidZonePressure_ =
        scalarField(solidMesh().faceZones()[solidZoneIndex()].size(), 0.0);

    solidZoneTraction_ =
        vectorField
        (
            solidMesh().faceZones()[solidZoneIndex_]().size(),
            vector::zero
        );

    solidZoneTractionPrev_ =
        vectorField
        (
            solidMesh().faceZones()[solidZoneIndex_]().size(),
            vector::zero
        );

    predictedSolidZoneTraction_ =
        vectorField
        (
            solidMesh().faceZones()[solidZoneIndex_]().size(),
            vector::zero
        );
    solidZone2Pressure_ =
        scalarField(solidMesh().faceZones()[solidZone2Index()].size(), 0.0); //zone2

    solidZone2Traction_ =
        vectorField
        (
            solidMesh().faceZones()[solidZone2Index_]().size(),
            vector::zero
        );

    solidZone2TractionPrev_ =
        vectorField
        (
            solidMesh().faceZones()[solidZone2Index_]().size(),
            vector::zero
        );

    predictedSolidZone2Traction_ =
        vectorField
        (
            solidMesh().faceZones()[solidZone2Index_]().size(),
            vector::zero
        );
    residual_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    residual2_ =
        vectorField
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().nPoints(),
            vector::zero
        );
//     Pout << "fluid mesh: " << fluidMesh().allPoints().size() << ", "
//         << fluidMesh().points().size() << endl;

//     Pout << "solid mesh: " << solidMesh().allPoints().size() << ", "
//         << solidMesh().points().size() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidStructureInterface::~fluidStructureInterface()
{
    deleteDemandDrivenData(currentSolidZonePointsPtr_);
    deleteDemandDrivenData(currentSolidZonePatchPtr_);
    deleteDemandDrivenData(fluidToSolidPtr_);
    deleteDemandDrivenData(ggiFluidToSolidPtr_);
    deleteDemandDrivenData(ggiInterpolatorPtr_);
    deleteDemandDrivenData(solidToFluidPtr_);
    deleteDemandDrivenData(accumulatedFluidInterfaceDisplacementPtr_);
    deleteDemandDrivenData(minEdgeLengthPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField&
Foam::fluidStructureInterface::currentSolidZonePoints() const
{
    if (!currentSolidZonePointsPtr_)
    {
        calcCurrentSolidZonePoints();
    }

    return *currentSolidZonePointsPtr_;
}

const Foam::vectorField& //zone2
Foam::fluidStructureInterface::currentSolidZone2Points() const
{
    if (!currentSolidZone2PointsPtr_)
    {
        calcCurrentSolidZone2Points();
    }

    return *currentSolidZone2PointsPtr_;
}

const Foam::PrimitivePatch<Foam::face, Foam::List, const Foam::pointField&>&
Foam::fluidStructureInterface::currentSolidZonePatch() const
{
    if (!currentSolidZonePatchPtr_)
    {
        calcCurrentSolidZonePatch();
    }

    return *currentSolidZonePatchPtr_;
}


const Foam::PrimitivePatch<Foam::face, Foam::List, const Foam::pointField&>& //zone2
Foam::fluidStructureInterface::currentSolidZone2Patch() const
{
    if (!currentSolidZone2PatchPtr_)
    {
        calcCurrentSolidZone2Patch();
    }

    return *currentSolidZone2PatchPtr_;
}

const Foam::zoneToZoneInterpolation&
Foam::fluidStructureInterface::fluidToSolid() const
{
    if (!fluidToSolidPtr_)
    {
        calcFluidToSolidInterpolator();
    }

    return *fluidToSolidPtr_;
}

const Foam::zoneToZoneInterpolation&
Foam::fluidStructureInterface::fluid2ToSolid() const
{
    if (!fluid2ToSolidPtr_)
    {
        calcFluid2ToSolidInterpolator();
    }

    return *fluid2ToSolidPtr_;
}

const Foam::ggiZoneInterpolation&
Foam::fluidStructureInterface::ggiFluidToSolid() const
{
    if (!ggiFluidToSolidPtr_)
    {
        calcGgiFluidToSolidInterpolator();
    }

    return *ggiFluidToSolidPtr_;
}

const Foam::ggiZoneInterpolation&
Foam::fluidStructureInterface::ggiFluid2ToSolid() const
{
    if (!ggiFluid2ToSolidPtr_)
    {
        calcGgiFluid2ToSolidInterpolator();
    }

    return *ggiFluid2ToSolidPtr_;
}

const Foam::extendedGgiZoneInterpolation&
Foam::fluidStructureInterface::ggiInterpolator() const
{
    if (!ggiInterpolatorPtr_)
    {
        calcGgiInterpolator();
    }

    return *ggiInterpolatorPtr_;
}

const Foam::extendedGgiZoneInterpolation&
Foam::fluidStructureInterface::ggiInterpolator2() const
{
    if (!ggiInterpolator2Ptr_)
    {
        calcGgiInterpolator2();
    }

    return *ggiInterpolator2Ptr_;
}

const Foam::zoneToZoneInterpolation&
Foam::fluidStructureInterface::solidToFluid() const
{
    if (!solidToFluidPtr_)
    {
        calcSolidToFluidInterpolator();
    }

    return *solidToFluidPtr_;
}

const Foam::zoneToZoneInterpolation&
Foam::fluidStructureInterface::solidToFluid2() const
{
    if (!solidToFluid2Ptr_)
    {
        calcSolidToFluid2Interpolator();
    }

    return *solidToFluid2Ptr_;
}

void Foam::fluidStructureInterface::initializeFields()
{
    fluidZonePointsDispl_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    fluidZonePointsDisplRef_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    fluidZonePointsDisplPrev_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    solidZonePointsDispl_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    solidZonePointsDisplRef_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );

    solidZone2PointsDispl_ = //zone2
        vectorField
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().nPoints(),
            vector::zero
        );

    solidZone2PointsDisplRef_ =
        vectorField
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().nPoints(),
            vector::zero
        );

   fluid2ZonePointsDispl_ =
        vectorField
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().nPoints(),
            vector::zero
        );

    fluid2ZonePointsDisplRef_ =
        vectorField
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().nPoints(),
            vector::zero
        );

    fluid2ZonePointsDisplPrev_ =
        vectorField
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().nPoints(),
            vector::zero
        );

    solidZonePointsDispl_ =   //?bian liang ming bu bian
        vectorField
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().nPoints(),
            vector::zero
        );

    solidZonePointsDisplRef_ =
        vectorField
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().nPoints(),
            vector::zero
        );

        solidZone2PointsDispl_ =   //?bian liang ming bu bian // zone2
        vectorField
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().nPoints(),
            vector::zero
        );

    solidZone2PointsDisplRef_ =
        vectorField
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().nPoints(),
            vector::zero
        );

    residualPrev_ = residual_;
    residual2Prev_ = residual2_;
//         vectorField
//         (
//             fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
//             vector::zero
//         );

    residual_ =
        vectorField
        (
            fluidMesh().faceZones()[fluidZoneIndex_]().nPoints(),
            vector::zero
        );


    residual2_ =
        vectorField
        (
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().nPoints(),
            vector::zero
        );

    maxResidualNorm_ = 0;

    outerCorr_ = 0;

    nOuterCorr_ = readInt(lookup("nOuterCorr"));

    outerCorrTolerance_ = readScalar(lookup("outerCorrTolerance"));

    coupled_ = Switch(lookup("coupled"));

    couplingReuse_ = readInt(lookup("couplingReuse"));

    relaxationFactor_ = readScalar(lookup("relaxationFactor"));

//     refPressure_ += refPressureIncrement_;

//     if (timeIndex_ < runTime().timeIndex())
//     {
//         timeIndex_ = runTime().timeIndex();

//     }
}


void Foam::fluidStructureInterface::updateInterpolator()
{
//     label interpolatorUpdateFrequency_ = 2;

    if (interpolatorUpdateFrequency_ != 0)
    {
        if (((runTime().timeIndex()-1)%interpolatorUpdateFrequency_) == 0)
        {
            deleteDemandDrivenData(ggiInterpolatorPtr_);
            ggiInterpolator();
        }
    }
    else
    {
        if ((runTime().timeIndex()-1) == 0)
        {
            deleteDemandDrivenData(ggiInterpolatorPtr_);
            ggiInterpolator();
        }
    };
    if (interpolatorUpdateFrequency_ != 0)
    {
        if (((runTime().timeIndex()-1)%interpolatorUpdateFrequency_) == 0)
        {
            deleteDemandDrivenData(ggiInterpolator2Ptr_);
            ggiInterpolator2();
        }
    }
    else
    {
        if ((runTime().timeIndex()-1) == 0)
        {
            deleteDemandDrivenData(ggiInterpolator2Ptr_);
            ggiInterpolator2();
        }
    }
}


void Foam::fluidStructureInterface::updateDisplacement()
{
    Info << "\nTime = " << flow().runTime().timeName()
        << ", iteration: " << outerCorr() << endl;

//     if (outerCorr_ == 1)
//     {
//         // Cancel residual from previous time step
//         fluidZonePointsDisplPrev() = fluidZonePointsDispl();
//         fluidZonePointsDispl() += residualPrev();
//     }


    if (couplingScheme() == "FixedRelaxation")
    {
        Info << "Current fsi under-relaxation factor: "
            << relaxationFactor() << endl;

        fluidZonePointsDisplPrev() = fluidZonePointsDispl();

        // if (outerCorr_ == 1)
        // {
        //     // Cancel residual from previous time step
        //     fluidZonePointsDispl() += residualPrev();
        // }

        fluidZonePointsDispl() += relaxationFactor()*residual();

    }
    else if (couplingScheme() == "Aitken")
    {
        if (outerCorr_ < 3)
        {
            Info << "Current fsi under-relaxation factor: "
                << relaxationFactor() << endl;

            fluidZonePointsDisplPrev() = fluidZonePointsDispl();

            // if (outerCorr_ == 1)
            // {
            //     // Cancel residual from previous time step
            //     fluidZonePointsDispl() += residualPrev();
            // }

            fluidZonePointsDispl() += relaxationFactor()*residual();
        }
        else
        {
            aitkenRelaxationFactor() =
               -aitkenRelaxationFactor()
               *(
                    sum
                    (
                        residualPrev()
                      & (residual() - residualPrev())
                    )
                   /(
                        sum
                        (
                            (residual() - residualPrev())
                          & (residual() - residualPrev())
                        )
                    )
                );

            if (Pstream::parRun())
            {
                if(!Pstream::master())
                {
                    aitkenRelaxationFactor() = 0.0;
                }

                //- pass to all procs
                reduce(aitkenRelaxationFactor(), sumOp<scalar>());
            }

            aitkenRelaxationFactor() = mag(aitkenRelaxationFactor());

            if (aitkenRelaxationFactor()>1)
            {
                aitkenRelaxationFactor() = relaxationFactor();
            }

            Info << "Current fsi under-relaxation factor (Aitken): "
                << aitkenRelaxationFactor() << endl;

            fluidZonePointsDisplPrev() = fluidZonePointsDispl();

            fluidZonePointsDispl() +=
                aitkenRelaxationFactor()*residual();
        }
    }
    else if (couplingScheme() == "IQN-ILS")
    {
//      A fluid structure interaction solver with IQN-ILS
//      coupling algorithm (J. Degroote, K.-J. Bathe and J. Vierendeels.
//      Performance of a new partitioned procedure versus a monolithic
//      procedure in fluid-structure interaction. Computers & Structures

        // IQN-ILS
        if (outerCorr_ == 1)
        {
            // Clean up data from old time steps

            Info << "Modes before clean-up : " << fluidPatchPointsT_.size();

            while (true)
            {
                if (fluidPatchPointsT_.size())
                {
                    if
                    (
                        flow().runTime().timeIndex()-couplingReuse()
                      > fluidPatchPointsT_[0]
                    )
                    {
                        for (label i = 0; i < fluidPatchPointsT_.size()-1; i++)
                        {
                            fluidPatchPointsT_[i] = fluidPatchPointsT_[i+1];
                            fluidPatchPointsV_[i] = fluidPatchPointsV_[i+1];
                            fluidPatchPointsW_[i] = fluidPatchPointsW_[i+1];
                        }
                        fluidPatchPointsT_.remove();
                        fluidPatchPointsV_.remove();
                        fluidPatchPointsW_.remove();
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }

            Info << ", modes after clean-up : "
                << fluidPatchPointsT_.size() << endl;
        }
        else if (outerCorr_ == 2)
        {
            // Set reference in the first coupling iteration
            solidZonePointsDisplRef() = solidZonePointsDispl();
            fluidZonePointsDisplRef() = fluidZonePointsDispl();
        }
        else
        {
            // Reference has been set in the first coupling iteration
            fluidPatchPointsV_.append
            (
                (
                    solidZonePointsDispl()
                  - fluidZonePointsDispl()
                )
              - (
                    solidZonePointsDisplRef()
                  - fluidZonePointsDisplRef()
                )
            );

            fluidPatchPointsW_.append
            (
                solidZonePointsDispl()
              - solidZonePointsDisplRef()
            );

            fluidPatchPointsT_.append
            (
                flow().runTime().timeIndex()
            );
        }

        if (fluidPatchPointsT_.size() > 1)
        {
            updateDisplacementUsingIQNILS();
        }
        else
        {
            // Relax the interface displacement
            Info << "Current fsi under-relaxation factor: "
                << relaxationFactor() << endl;

            fluidZonePointsDisplPrev() = fluidZonePointsDispl();

            // if (outerCorr_ == 1)
            // {
            //     // Cancel residual from previous time step
            //     fluidZonePointsDispl() += residualPrev();
            // }

            fluidZonePointsDispl() += relaxationFactor()*residual();
        }
    }


    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    if (Pstream::parRun())
    {
        if(!Pstream::master())
        {
            fluidZonePointsDispl() *= 0.0;
        }

        //- pass to all procs
        reduce(fluidZonePointsDispl(), sumOp<vectorField>());

        label globalFluidZoneIndex =
            findIndex(flow().globalFaceZones(), fluidZoneIndex());

        if (globalFluidZoneIndex == -1)
        {
            FatalErrorIn
            (
                "fluidStructureInterface::updateDisplacement()"
            )   << "global zone point map is not availabel"
                << abort(FatalError);
        }

        const labelList& map =
            flow().globalToLocalFaceZonePointMap()[globalFluidZoneIndex];

        if(!Pstream::master())
        {
            vectorField fluidZonePointsDisplGlobal =
                fluidZonePointsDispl();

            forAll(fluidZonePointsDisplGlobal, globalPointI)
            {
                label localPoint = map[globalPointI];

                fluidZonePointsDispl()[localPoint] =
                    fluidZonePointsDisplGlobal[globalPointI];
            }
        }
    }
}

void Foam::fluidStructureInterface::updateDisplacement2()
{
    Info << "\nTime = " << flow2().runTime().timeName()
        << ", iteration: " << outerCorr() << endl;

//     if (outerCorr_ == 1)
//     {
//         // Cancel residual from previous time step
//         fluidZonePointsDisplPrev() = fluidZonePointsDispl();
//         fluidZonePointsDispl() += residualPrev();
//     }


    if (couplingScheme() == "FixedRelaxation")
    {
        Info << "Current fsi under-relaxation factor: "
            << relaxationFactor() << endl;

        fluid2ZonePointsDisplPrev() = fluid2ZonePointsDispl();

        // if (outerCorr_ == 1)
        // {
        //     // Cancel residual from previous time step
        //     fluidZonePointsDispl() += residualPrev();
        // }

        fluid2ZonePointsDispl() += relaxationFactor()*residual2();

    }
    else if (couplingScheme() == "Aitken")
    {
        if (outerCorr_ < 3)
        {
            Info << "Current fsi under-relaxation factor: "
                << relaxationFactor() << endl;

            fluid2ZonePointsDisplPrev() = fluid2ZonePointsDispl();

            // if (outerCorr_ == 1)
            // {
            //     // Cancel residual from previous time step
            //     fluidZonePointsDispl() += residualPrev();
            // }

            fluid2ZonePointsDispl() += relaxationFactor()*residual2();
        }
        else
        {
            aitkenRelaxationFactor() =
               -aitkenRelaxationFactor()
               *(
                    sum
                    (
                        residual2Prev()
                      & (residual2() - residual2Prev())
                    )
                   /(
                        sum
                        (
                            (residual2() - residual2Prev())
                          & (residual2() - residual2Prev())
                        )
                    )
                );

            if (Pstream::parRun())
            {
                if(!Pstream::master())
                {
                    aitkenRelaxationFactor() = 0.0;
                }

                //- pass to all procs
                reduce(aitkenRelaxationFactor(), sumOp<scalar>());
            }

            aitkenRelaxationFactor() = mag(aitkenRelaxationFactor());

            if (aitkenRelaxationFactor()>1)
            {
                aitkenRelaxationFactor() = relaxationFactor();
            }

            Info << "Current fsi under-relaxation factor (Aitken): "
                << aitkenRelaxationFactor() << endl;

            fluid2ZonePointsDisplPrev() = fluid2ZonePointsDispl();

            fluid2ZonePointsDispl() +=
                aitkenRelaxationFactor()*residual2();
        }
    }
    else if (couplingScheme() == "IQN-ILS")
    {
//      A fluid structure interaction solver with IQN-ILS
//      coupling algorithm (J. Degroote, K.-J. Bathe and J. Vierendeels.
//      Performance of a new partitioned procedure versus a monolithic
//      procedure in fluid-structure interaction. Computers & Structures

        // IQN-ILS
        if (outerCorr_ == 1)
        {
            // Clean up data from old time steps

            Info << "Modes before clean-up : " << fluid2PatchPointsT_.size();

            while (true)
            {
                if (fluid2PatchPointsT_.size())
                {
                    if
                    (
                        flow2().runTime().timeIndex()-couplingReuse()
                      > fluid2PatchPointsT_[0]
                    )
                    {
                        for (label i = 0; i < fluid2PatchPointsT_.size()-1; i++)
                        {
                            fluid2PatchPointsT_[i] = fluid2PatchPointsT_[i+1];
                            fluid2PatchPointsV_[i] = fluid2PatchPointsV_[i+1];
                            fluid2PatchPointsW_[i] = fluid2PatchPointsW_[i+1];
                        }
                        fluid2PatchPointsT_.remove();
                        fluid2PatchPointsV_.remove();
                        fluid2PatchPointsW_.remove();
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    break;
                }
            }

            Info << ", modes after clean-up : "
                << fluid2PatchPointsT_.size() << endl;
        }
        else if (outerCorr_ == 2)
        {
            // Set reference in the first coupling iteration
            solidZone2PointsDisplRef() = solidZone2PointsDispl();
            fluid2ZonePointsDisplRef() = fluid2ZonePointsDispl();
        }
        else
        {
            // Reference has been set in the first coupling iteration
            fluid2PatchPointsV_.append
            (
                (
                    solidZone2PointsDispl()
                  - fluid2ZonePointsDispl()
                )
              - (
                    solidZone2PointsDisplRef()
                  - fluid2ZonePointsDisplRef()
                )
            );

            fluid2PatchPointsW_.append
            (
                solidZone2PointsDispl()
              - solidZone2PointsDisplRef()
            );

            fluid2PatchPointsT_.append
            (
                flow2().runTime().timeIndex()
            );
        }

        if (fluid2PatchPointsT_.size() > 1)
        {
            updateDisplacement2UsingIQNILS();
        }
        else
        {
            // Relax the interface displacement
            Info << "Current fsi under-relaxation factor: "
                << relaxationFactor() << endl;

            fluid2ZonePointsDisplPrev() = fluid2ZonePointsDispl();

            // if (outerCorr_ == 1)
            // {
            //     // Cancel residual from previous time step
            //     fluidZonePointsDispl() += residualPrev();
            // }

            fluid2ZonePointsDispl() += relaxationFactor()*residual2();
        }
    }


    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    if (Pstream::parRun())
    {
        if(!Pstream::master())
        {
            fluid2ZonePointsDispl() *= 0.0;
        }

        //- pass to all procs
        reduce(fluid2ZonePointsDispl(), sumOp<vectorField>());

        label globalFluid2ZoneIndex =
            findIndex(flow2().globalFaceZones(), fluid2ZoneIndex());

        if (globalFluid2ZoneIndex == -1)
        {
            FatalErrorIn
            (
                "fluidStructureInterface::updateDisplacement2()"
            )   << "global zone point map is not availabel"
                << abort(FatalError);
        }

        const labelList& map =
            flow2().globalToLocalFaceZonePointMap()[globalFluid2ZoneIndex];

        if(!Pstream::master())
        {
            vectorField fluid2ZonePointsDisplGlobal =
                fluid2ZonePointsDispl();

            forAll(fluid2ZonePointsDisplGlobal, globalPointI)
            {
                label localPoint = map[globalPointI];

                fluid2ZonePointsDispl()[localPoint] =
                    fluid2ZonePointsDisplGlobal[globalPointI];
            }
        }
    }
}

void Foam::fluidStructureInterface::updateWeakDisplacement()
{
    vectorField solidZonePointsDisplAtSolid =
        stress().faceZonePointDisplacementIncrement(solidZoneIndex());

    solidZonePointsDispl() =
        ggiInterpolator().slaveToMasterPointInterpolate
        (
            solidZonePointsDisplAtSolid
        );

    residualPrev() = residual();

    residual() = solidZonePointsDispl() - fluidZonePointsDispl();

    fluidZonePointsDisplPrev() = fluidZonePointsDispl();

    fluidZonePointsDispl() += residual();

    // Make sure that displacement on all processors is equal to one
    // calculated on master processor
    if (Pstream::parRun())
    {
        if(!Pstream::master())
        {
            fluidZonePointsDispl() *= 0.0;
        }

        //- pass to all procs
        reduce(fluidZonePointsDispl(), sumOp<vectorField>());

        label globalFluidZoneIndex =
            findIndex(flow().globalFaceZones(), fluidZoneIndex());

        if (globalFluidZoneIndex == -1)
        {
            FatalErrorIn
            (
                "fluidStructureInterface::updateDisplacement()"
            )   << "global zone point map is not availabel"
                << abort(FatalError);
        }

        const labelList& map =
            flow().globalToLocalFaceZonePointMap()[globalFluidZoneIndex];

        if(!Pstream::master())
        {
            vectorField fluidZonePointsDisplGlobal =
                fluidZonePointsDispl();

            forAll(fluidZonePointsDisplGlobal, globalPointI)
            {
                label localPoint = map[globalPointI];

                fluidZonePointsDispl()[localPoint] =
                    fluidZonePointsDisplGlobal[globalPointI];
            }
        }
    }
}


void Foam::fluidStructureInterface::updateDisplacementUsingIQNILS()
{
    // Consider fluidPatchPointsV as a matrix V
    // with as columns the items
    // in the DynamicList and calculate the QR-decomposition of V
    // with modified Gram-Schmidt
    label cols = fluidPatchPointsV_.size();
    RectangularMatrix<scalar> R(cols, cols, 0.0);
    RectangularMatrix<scalar> C(cols, 1);
    RectangularMatrix<scalar> Rcolsum(1, cols);
    DynamicList<vectorField> Q;

    for (label i = 0; i < cols; i++)
    {
        Q.append(fluidPatchPointsV_[cols-1-i]);
    }

    for (label i = 0; i < cols; i++)
    {
        // Normalize column i
        R[i][i] = Foam::sqrt(sum(Q[i] & Q[i]));
        Q[i] /= R[i][i];

        // Orthogonalize columns to the right of column i
        for (label j = i+1; j < cols; j++)
        {
            R[i][j] = sum(Q[i] & Q[j]);
            Q[j] -= R[i][j]*Q[i];
        }

        // Project minus the residual vector on the Q
        C[i][0] = sum
        (
            Q[i]
          & (
                fluidZonePointsDispl()
              - solidZonePointsDispl()
            )
        );
    }

    // Solve the upper triangular system
    for (label j = 0; j < cols; j++)
    {
        Rcolsum[0][j] = 0.0;
        for (label i = 0; i < j+1; i++)
        {
            Rcolsum[0][j] += cmptMag(R[i][j]);
        }
    }
    scalar epsilon = 1.0E-10*max(Rcolsum);
    for (label i = 0; i < cols; i++)
    {
        if (cmptMag(R[i][i]) > epsilon)
        {
            for (label j = i+1; j < cols; j++)
            {
                R[i][j] /= R[i][i];
            }
            C[i][0] /= R[i][i];
            R[i][i] = 1.0;
        }
    }
    for (label j = cols-1; j >= 0; j--)
    {
        if (cmptMag(R[j][j]) > epsilon)
        {
            for (label i = 0; i < j; i++)
            {
                C[i][0] -= C[j][0]*R[i][j];
            }
        }
        else
        {
            C[j][0] = 0.0;
        }
    }

    fluidZonePointsDisplPrev() = fluidZonePointsDispl();

    fluidZonePointsDispl() = solidZonePointsDispl();

    for (label i = 0; i < cols; i++)
    {
        fluidZonePointsDispl() += fluidPatchPointsW_[i]*C[cols-1-i][0];
    }
}

void Foam::fluidStructureInterface::updateDisplacement2UsingIQNILS()
{
    // Consider fluidPatchPointsV as a matrix V
    // with as columns the items
    // in the DynamicList and calculate the QR-decomposition of V
    // with modified Gram-Schmidt
    label cols = fluid2PatchPointsV_.size();
    RectangularMatrix<scalar> R(cols, cols, 0.0);
    RectangularMatrix<scalar> C(cols, 1);
    RectangularMatrix<scalar> Rcolsum(1, cols);
    DynamicList<vectorField> Q;

    for (label i = 0; i < cols; i++)
    {
        Q.append(fluid2PatchPointsV_[cols-1-i]);
    }

    for (label i = 0; i < cols; i++)
    {
        // Normalize column i
        R[i][i] = Foam::sqrt(sum(Q[i] & Q[i]));
        Q[i] /= R[i][i];

        // Orthogonalize columns to the right of column i
        for (label j = i+1; j < cols; j++)
        {
            R[i][j] = sum(Q[i] & Q[j]);
            Q[j] -= R[i][j]*Q[i];
        }

        // Project minus the residual vector on the Q
        C[i][0] = sum
        (
            Q[i]
          & (
                fluid2ZonePointsDispl()
              - solidZone2PointsDispl()
            )
        );
    }

    // Solve the upper triangular system
    for (label j = 0; j < cols; j++)
    {
        Rcolsum[0][j] = 0.0;
        for (label i = 0; i < j+1; i++)
        {
            Rcolsum[0][j] += cmptMag(R[i][j]);
        }
    }
    scalar epsilon = 1.0E-10*max(Rcolsum);
    for (label i = 0; i < cols; i++)
    {
        if (cmptMag(R[i][i]) > epsilon)
        {
            for (label j = i+1; j < cols; j++)
            {
                R[i][j] /= R[i][i];
            }
            C[i][0] /= R[i][i];
            R[i][i] = 1.0;
        }
    }
    for (label j = cols-1; j >= 0; j--)
    {
        if (cmptMag(R[j][j]) > epsilon)
        {
            for (label i = 0; i < j; i++)
            {
                C[i][0] -= C[j][0]*R[i][j];
            }
        }
        else
        {
            C[j][0] = 0.0;
        }
    }

    fluid2ZonePointsDisplPrev() = fluid2ZonePointsDispl();

    fluid2ZonePointsDispl() = solidZone2PointsDispl();

    for (label i = 0; i < cols; i++)
    {
        fluid2ZonePointsDispl() += fluid2PatchPointsW_[i]*C[cols-1-i][0];
    }
}

void Foam::fluidStructureInterface::moveFluidMesh()
{
    // Get fluid patch displacement from fluid zone displacement

    vectorField fluidPatchPointsDispl
    (
        fluidMesh().boundaryMesh()[fluidPatchIndex()].nPoints(),
        vector::zero
    );

    vectorField fluidPatchPointsDisplPrev
    (
        fluidMesh().boundaryMesh()[fluidPatchIndex()].nPoints(),
        vector::zero
    );

    const labelList& fluidPatchMeshPoints =
        fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

    forAll(fluidPatchPointsDispl, pointI)
    {
        label curMeshPointID = fluidPatchMeshPoints[pointI];

        label curFluidZonePointID =
            fluidMesh().faceZones()[fluidZoneIndex()]()
           .whichPoint(curMeshPointID);

        fluidPatchPointsDispl[pointI] =
            fluidZonePointsDispl()[curFluidZonePointID];

        fluidPatchPointsDisplPrev[pointI] =
            fluidZonePointsDisplPrev()[curFluidZonePointID];
    }

    // Move fluid mesh
    const vectorField& n =
        fluidMesh().boundaryMesh()[fluidPatchIndex()].pointNormals();

    primitivePatchInterpolation patchInterpolator
    (
        fluidMesh().boundaryMesh()[fluidPatchIndex()]
    );

    scalarField pointDeltaCoeffs =
        patchInterpolator.faceToPointInterpolate
        (
            fluidMesh().boundary()[fluidPatchIndex()].deltaCoeffs()
        );

    scalar delta =
        gMax
        (
            mag
            (
                n
              & (
                    accumulatedFluidInterfaceDisplacement()
                  + fluidPatchPointsDispl
                  - fluidPatchPointsDisplPrev
                )
            )
           *pointDeltaCoeffs
        );

    Info << "Maximal accumulated displacement of interface points: "
        << delta << endl;

    if (delta < interfaceDeformationLimit())
    {
        // Move only interface points
        pointField newPoints = fluidMesh().allPoints();

        const labelList& meshPoints =
            fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

        forAll (fluidPatchPointsDispl, pointI)
        {
            newPoints[meshPoints[pointI]] +=
                fluidPatchPointsDispl[pointI]
              - fluidPatchPointsDisplPrev[pointI];
        }

        twoDPointCorrector twoDCorrector(fluidMesh());

        twoDCorrector.correctPoints(newPoints);

        fluidMesh_.movePoints(newPoints);

        // Accumulate interface points displacement
        accumulatedFluidInterfaceDisplacement() +=
            fluidPatchPointsDispl
          - fluidPatchPointsDisplPrev;
    }
    else
    {
        // Move whole fluid mesh
        pointField newPoints = fluidMesh().allPoints();

        const labelList& meshPoints =
            fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

        forAll (accumulatedFluidInterfaceDisplacement(), pointI)
        {
            newPoints[meshPoints[pointI]] -=
                accumulatedFluidInterfaceDisplacement()[pointI];
        }

        twoDPointCorrector twoDCorrector(fluidMesh());

        twoDCorrector.correctPoints(newPoints);

        fluidMesh_.movePoints(newPoints);

        accumulatedFluidInterfaceDisplacement() +=
            fluidPatchPointsDispl
          - fluidPatchPointsDisplPrev;

        // Check mesh motion solver type
        bool feMotionSolver =
            fluidMesh().objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );

        bool fvMotionSolver =
            fluidMesh().objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );

        bool rbfMotionSolver =
            fluidMesh().objectRegistry::foundObject<RBFMotionSolver>
            (
                "dynamicMeshDict"
            );

        if (rbfMotionSolver)
        {
            // Grab RBF motion solver
            RBFMotionSolver& ms =
                const_cast<RBFMotionSolver&>
                (
                    fluidMesh().objectRegistry::lookupObject<RBFMotionSolver>
                    (
                        "dynamicMeshDict"
                    )
                );

            Info << "RBF mesh motion" << endl;

            const labelList& movingMeshPoints = ms.movingIDs();

            vectorField motion(movingMeshPoints.size(), vector::zero);

            vectorField fluidPatchDisplacement =
                accumulatedFluidInterfaceDisplacement();
//                /flow().runTime().deltaT().value();

            const labelList& meshPoints =
                fluidMesh().boundaryMesh()[fluidPatchIndex()].meshPoints();

            forAll(meshPoints, pointI)
            {
                label curMovingPoint =
                    findIndex(movingMeshPoints, meshPoints[pointI]);

                if (curMovingPoint != -1)
                {
                    motion[curMovingPoint] = fluidPatchDisplacement[pointI];
                }
            }

            ms.setMotion(motion);

//             FatalErrorIn("fluidStructureInterface::moveFluidMesh()")
//                 << "Problem with fluid mesh motion solver selection"
//                     << abort(FatalError);
        }
        else if (feMotionSolver)
        {
            tetPointVectorField& motionU =
                const_cast<tetPointVectorField&>
                (
                    fluidMesh().objectRegistry::
                    lookupObject<tetPointVectorField>
                    (
                        "motionU"
                    )
                );

            fixedValueTetPolyPatchVectorField& motionUFluidPatch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[fluidPatchIndex()]
                );

            tetPolyPatchInterpolation tppi
            (
                refCast<const faceTetPolyPatch>(motionUFluidPatch.patch())
            );

            motionUFluidPatch ==
                tppi.pointToPointInterpolate
                (
                    accumulatedFluidInterfaceDisplacement()
                   /flow().runTime().deltaT().value()
                );
        }
        else if (fvMotionSolver)
        {
            pointVectorField& motionU =
                const_cast<pointVectorField&>
                (
                    fluidMesh().objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            fixedValuePointPatchVectorField& motionUFluidPatch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[fluidPatchIndex()]
                );

            motionUFluidPatch ==
                accumulatedFluidInterfaceDisplacement()
               /flow().runTime().deltaT().value();
        }
        else
        {
            FatalErrorIn("fluidStructureInterface::moveFluidMesh()")
                << "Problem with fluid mesh motion solver selection"
                    << abort(FatalError);
        }

        fluidMesh_.update();

        accumulatedFluidInterfaceDisplacement() =
            vectorField
            (
                accumulatedFluidInterfaceDisplacement().size(),
                vector::zero
            );
    }


    // Move unused fluid mesh points
    {
        vectorField newPoints = fluidMesh().allPoints();

        const labelList& fluidZoneMeshPoints =
            fluidMesh().faceZones()[fluidZoneIndex()]().meshPoints();

        forAll(fluidZonePointsDispl(), pointI)
        {
            if (fluidZoneMeshPoints[pointI] >= fluidMesh().nPoints())
            {
                newPoints[fluidZoneMeshPoints[pointI]] +=
                    fluidZonePointsDispl()[pointI]
                  - fluidZonePointsDisplPrev()[pointI];
            }
        }

        twoDPointCorrector twoDCorrector(fluidMesh());

        twoDCorrector.correctPoints(newPoints);

        fluidMesh_.movePoints(newPoints);
    }
}

void Foam::fluidStructureInterface::moveFluidMesh2()
{
    // Get fluid patch displacement from fluid zone displacement

    vectorField fluid2PatchPointsDispl
    (
        fluidMesh2().boundaryMesh()[fluid2PatchIndex()].nPoints(),
        vector::zero
    );

    vectorField fluid2PatchPointsDisplPrev
    (
        fluidMesh2().boundaryMesh()[fluid2PatchIndex()].nPoints(),
        vector::zero
    );

    const labelList& fluid2PatchMeshPoints =
        fluidMesh2().boundaryMesh()[fluid2PatchIndex()].meshPoints();

    forAll(fluid2PatchPointsDispl, pointI)
    {
        label curMeshPointID = fluid2PatchMeshPoints[pointI];

        label curFluid2ZonePointID =
            fluidMesh2().faceZones()[fluid2ZoneIndex()]()
           .whichPoint(curMeshPointID);

        fluid2PatchPointsDispl[pointI] =
            fluid2ZonePointsDispl()[curFluid2ZonePointID];

        fluid2PatchPointsDisplPrev[pointI] =
            fluid2ZonePointsDisplPrev()[curFluid2ZonePointID];
    }

    // Move fluid2 mesh
    const vectorField& n =
        fluidMesh2().boundaryMesh()[fluid2PatchIndex()].pointNormals();

    primitivePatchInterpolation patchInterpolator
    (
        fluidMesh2().boundaryMesh()[fluid2PatchIndex()]
    );

    scalarField pointDeltaCoeffs =
        patchInterpolator.faceToPointInterpolate
        (
            fluidMesh2().boundary()[fluid2PatchIndex()].deltaCoeffs()
        );

    scalar delta =
        gMax
        (
            mag
            (
                n
              & (
                    accumulatedFluid2InterfaceDisplacement()
                  + fluid2PatchPointsDispl
                  - fluid2PatchPointsDisplPrev
                )
            )
           *pointDeltaCoeffs
        );

    Info << "Maximal accumulated displacement of interface points: "
        << delta << endl;

    if (delta < interfaceDeformationLimit())
    {
        // Move only interface points
        pointField newPoints = fluidMesh2().allPoints();

        const labelList& meshPoints =
            fluidMesh2().boundaryMesh()[fluid2PatchIndex()].meshPoints();

        forAll (fluid2PatchPointsDispl, pointI)
        {
            newPoints[meshPoints[pointI]] +=
                fluid2PatchPointsDispl[pointI]
              - fluid2PatchPointsDisplPrev[pointI];
        }

        twoDPointCorrector twoDCorrector(fluidMesh2());

        twoDCorrector.correctPoints(newPoints);

        fluidMesh2_.movePoints(newPoints);

        // Accumulate interface points displacement
        accumulatedFluid2InterfaceDisplacement() +=
            fluid2PatchPointsDispl
          - fluid2PatchPointsDisplPrev;
    }
    else
    {
        // Move whole fluid mesh
        pointField newPoints = fluidMesh2().allPoints();

        const labelList& meshPoints =
            fluidMesh2().boundaryMesh()[fluid2PatchIndex()].meshPoints();

        forAll (accumulatedFluid2InterfaceDisplacement(), pointI)
        {
            newPoints[meshPoints[pointI]] -=
                accumulatedFluid2InterfaceDisplacement()[pointI];
        }

        twoDPointCorrector twoDCorrector(fluidMesh2());

        twoDCorrector.correctPoints(newPoints);

        fluidMesh2_.movePoints(newPoints);

        accumulatedFluid2InterfaceDisplacement() +=
            fluid2PatchPointsDispl
          - fluid2PatchPointsDisplPrev;

        // Check mesh motion solver type
        bool feMotionSolver =
            fluidMesh2().objectRegistry::foundObject<tetPointVectorField>
            (
                "motionU"
            );

        bool fvMotionSolver =
            fluidMesh2().objectRegistry::foundObject<pointVectorField>
            (
                "pointMotionU"
            );

        bool rbfMotionSolver =
            fluidMesh2().objectRegistry::foundObject<RBFMotionSolver>
            (
                "dynamicMeshDict"
            );

        if (rbfMotionSolver)
        {
            // Grab RBF motion solver
            RBFMotionSolver& ms =
                const_cast<RBFMotionSolver&>
                (
                    fluidMesh2().objectRegistry::lookupObject<RBFMotionSolver>
                    (
                        "dynamicMeshDict"
                    )
                );

            Info << "RBF mesh2 motion" << endl;

            const labelList& movingMesh2Points = ms.movingIDs();

            vectorField motion(movingMesh2Points.size(), vector::zero);

            vectorField fluid2PatchDisplacement =
                accumulatedFluid2InterfaceDisplacement();
//                /flow().runTime().deltaT().value();

            const labelList& meshPoints =
                fluidMesh2().boundaryMesh()[fluid2PatchIndex()].meshPoints();

            forAll(meshPoints, pointI)
            {
                label curMovingPoint =
                    findIndex(movingMesh2Points, meshPoints[pointI]);

                if (curMovingPoint != -1)
                {
                    motion[curMovingPoint] = fluid2PatchDisplacement[pointI];
                }
            }

            ms.setMotion(motion);

//             FatalErrorIn("fluidStructureInterface::moveFluidMesh()")
//                 << "Problem with fluid mesh motion solver selection"
//                     << abort(FatalError);
        }
        else if (feMotionSolver)
        {
            tetPointVectorField& motionU =
                const_cast<tetPointVectorField&>
                (
                    fluidMesh2().objectRegistry::
                    lookupObject<tetPointVectorField>
                    (
                        "motionU"
                    )
                );

            fixedValueTetPolyPatchVectorField& motionUFluid2Patch =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[fluid2PatchIndex()]
                );

            tetPolyPatchInterpolation tppi
            (
                refCast<const faceTetPolyPatch>(motionUFluid2Patch.patch())
            );

            motionUFluid2Patch ==
                tppi.pointToPointInterpolate
                (
                    accumulatedFluid2InterfaceDisplacement()
                   /flow2().runTime().deltaT().value()
                );
        }
        else if (fvMotionSolver)
        {
            pointVectorField& motionU =
                const_cast<pointVectorField&>
                (
                    fluidMesh2().objectRegistry::
                    lookupObject<pointVectorField>
                    (
                        "pointMotionU"
                    )
                );

            fixedValuePointPatchVectorField& motionUFluid2Patch =
                refCast<fixedValuePointPatchVectorField>
                (
                    motionU.boundaryField()[fluid2PatchIndex()]
                );

            motionUFluid2Patch ==
                accumulatedFluid2InterfaceDisplacement()
               /flow2().runTime().deltaT().value();
        }
        else
        {
            FatalErrorIn("fluidStructureInterface::moveFluidMesh2()")
                << "Problem with fluid2 mesh motion solver selection"
                    << abort(FatalError);
        }

        fluidMesh2_.update();

        accumulatedFluid2InterfaceDisplacement() =
            vectorField
            (
                accumulatedFluid2InterfaceDisplacement().size(),
                vector::zero
            );
    }


    // Move unused fluid mesh points
    {
        vectorField newPoints = fluidMesh2().allPoints();

        const labelList& fluid2ZoneMeshPoints =
            fluidMesh2().faceZones()[fluid2ZoneIndex()]().meshPoints();

        forAll(fluid2ZonePointsDispl(), pointI)
        {
            if (fluid2ZoneMeshPoints[pointI] >= fluidMesh2().nPoints())
            {
                newPoints[fluid2ZoneMeshPoints[pointI]] +=
                    fluid2ZonePointsDispl()[pointI]
                  - fluid2ZonePointsDisplPrev()[pointI];
            }
        }

        twoDPointCorrector twoDCorrector(fluidMesh2());

        twoDCorrector.correctPoints(newPoints);

        fluidMesh2_.movePoints(newPoints);
    }
}



void Foam::fluidStructureInterface::updateForce()
{
    Info << "Setting traction on solid patch" << endl;

    vectorField fluidZoneTraction =
        flow().faceZoneViscousForce
        (
            fluidZoneIndex(),
            fluidPatchIndex()
        );

    scalarField fluidZonePressure =
        flow().faceZonePressureForce(fluidZoneIndex(), fluidPatchIndex());


    // Fluid zone face normals
    const vectorField& p =
        fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();
    const faceList& f =
        fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

    vectorField n(f.size(), vector::zero);
    forAll(n, faceI)
    {
        n[faceI] = f[faceI].normal(p);
        n[faceI] /= mag(n[faceI]);
    }

    // Fluid zone total traction
    vectorField fluidZoneTotalTraction =
        fluidZoneTraction - fluidZonePressure*n;


    vectorField solidZoneTraction =
        ggiInterpolator().masterToSlave
        (
           -fluidZoneTraction
        );

    vectorField solidZoneTotalTraction =
        ggiInterpolator().masterToSlave
        (
            -fluidZoneTotalTraction
        );

    scalarField solidZoneMuEff =
        ggiInterpolator().masterToSlave
        (
            flow().faceZoneMuEff(fluidZoneIndex(), fluidPatchIndex())
        );

    tensorField solidZoneSurfaceGradientOfVelocity =
        stress().faceZoneSurfaceGradientOfVelocity
        (
            solidZoneIndex(),
            solidPatchIndex()
        );

    vectorField solidZoneNormal =
        stress().faceZoneNormal
        (
            solidZoneIndex(),
            solidPatchIndex()
        );

    // Add part of the viscous force present only
    // at the deforming and moving walls
    solidZoneTraction +=
        solidZoneMuEff
       *(
           -2*tr(solidZoneSurfaceGradientOfVelocity)*solidZoneNormal
          + (solidZoneSurfaceGradientOfVelocity&solidZoneNormal)
        );

    vectorField movingTraction =
        solidZoneMuEff
       *(
           -2*tr(solidZoneSurfaceGradientOfVelocity)*solidZoneNormal
          + (solidZoneSurfaceGradientOfVelocity&solidZoneNormal)
        );

//     Info << "Movig wall traction, max: "
//         << max(mag(movingTraction))
//             << ", avg: " << average(mag(movingTraction))
//             << ", min: " << min(mag(movingTraction)) << endl;

//     vector Fm =
//         sum
//         (
//             movingTraction
//            *solidMesh().magSf().boundaryField()[solidPatchIndex()]
//         );

//     Info << "Movig wall total force: " << Fm << endl;

    solidZonePressure_ =
//         fluidToSolid().faceInterpolate
//         ggiFluidToSolid().masterToSlave
        ggiInterpolator().masterToSlave
        (
            fluidZonePressure
        );


    if (coupled())
    {
//         stress().setPressure
//         (
//             solidPatchIndex(),
//             solidZoneIndex(),
//             solidZonePressure_
//         );

        stress().setTraction
        (
            solidPatchIndex(),
            solidZoneIndex(),
            solidZoneTotalTraction
//             solidZoneTraction
        );
    }

    // Total force at the fluid side of the interface
    {
        const vectorField& p =
            fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

        const faceList& f =
            fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
        }

        vector totalTractionForce = sum(fluidZoneTotalTraction*mag(S));

//         vector totalPressureForce = sum(fluidZonePressure*S);

        Info << "Total force (fluid) = "
            << totalTractionForce << endl;

//         Info << "Total pressure force (fluid) = "
//             << totalPressureForce << endl;
    }

    // Totla force at the solid side of the interface
    {
        const vectorField& p =
            solidMesh().faceZones()[solidZoneIndex_]().localPoints();

        const faceList& f =
            solidMesh().faceZones()[solidZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
        }

        vector totalTractionForce =
            sum(solidZoneTotalTraction*mag(S));

//         vector totalPressureForce =
//             sum(solidZonePressure_*S);

        Info << "Total force (solid) = "
            << totalTractionForce << endl;

//         Info << "Total pressure force (solid) = "
//             << totalPressureForce << endl;
    }
}


void Foam::fluidStructureInterface::updateForce2()
{
    Info << "Setting traction on solid patch" << endl;

    vectorField fluid2ZoneTraction =
        flow2().faceZoneViscousForce
        (
            fluid2ZoneIndex(),
            fluid2PatchIndex()
        );

    scalarField fluid2ZonePressure =
        flow2().faceZonePressureForce(fluid2ZoneIndex(), fluid2PatchIndex());


    // Fluid zone face normals
    const vectorField& p =
        fluidMesh2().faceZones()[fluid2ZoneIndex_]().localPoints();
    const faceList& f =
        fluidMesh2().faceZones()[fluid2ZoneIndex_]().localFaces();

    vectorField n(f.size(), vector::zero);
    forAll(n, faceI)
    {
        n[faceI] = f[faceI].normal(p);
        n[faceI] /= mag(n[faceI]);
    }

    // Fluid zone total traction
    vectorField fluid2ZoneTotalTraction =
        fluid2ZoneTraction - fluid2ZonePressure*n;


    vectorField solidZone2Traction =
        ggiInterpolator2().masterToSlave
        (
           -fluid2ZoneTraction
        );

    vectorField solidZone2TotalTraction =
        ggiInterpolator2().masterToSlave
        (
            -fluid2ZoneTotalTraction
        );

    scalarField solidZone2MuEff =
        ggiInterpolator2().masterToSlave
        (
            flow2().faceZoneMuEff(fluid2ZoneIndex(), fluid2PatchIndex())
        );

    tensorField solidZone2SurfaceGradientOfVelocity =
        stress().faceZoneSurfaceGradientOfVelocity
        (
            solidZone2Index(),
            solidPatch2Index()
        );

    vectorField solidZone2Normal =
        stress().faceZoneNormal
        (
            solidZone2Index(),
            solidPatch2Index()
        );

    // Add part of the viscous force present only
    // at the deforming and moving walls
    solidZone2Traction +=
        solidZone2MuEff
       *(
           -2*tr(solidZone2SurfaceGradientOfVelocity)*solidZone2Normal
          + (solidZone2SurfaceGradientOfVelocity&solidZone2Normal)
        );

    vectorField movingTraction =
        solidZone2MuEff
       *(
           -2*tr(solidZone2SurfaceGradientOfVelocity)*solidZone2Normal
          + (solidZone2SurfaceGradientOfVelocity&solidZone2Normal)
        );

//     Info << "Movig wall traction, max: "
//         << max(mag(movingTraction))
//             << ", avg: " << average(mag(movingTraction))
//             << ", min: " << min(mag(movingTraction)) << endl;

//     vector Fm =
//         sum
//         (
//             movingTraction
//            *solidMesh().magSf().boundaryField()[solidPatchIndex()]
//         );

//     Info << "Movig wall total force: " << Fm << endl;

    solidZone2Pressure_ =
//         fluidToSolid().faceInterpolate
//         ggiFluidToSolid().masterToSlave
        ggiInterpolator2().masterToSlave
        (
            fluid2ZonePressure
        );


    if (coupled())
    {
//         stress().setPressure
//         (
//             solidPatchIndex(),
//             solidZoneIndex(),
//             solidZonePressure_
//         );

        stress().setTraction
        (
            solidPatch2Index(),
            solidZone2Index(),
            solidZone2TotalTraction
//             solidZoneTraction
        );
    }

    // Total force at the fluid side of the interface
    {
        const vectorField& p =
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().localPoints();

        const faceList& f =
            fluidMesh2().faceZones()[fluid2ZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
        }

        vector totalTractionForce = sum(fluid2ZoneTotalTraction*mag(S));

//         vector totalPressureForce = sum(fluidZonePressure*S);

        Info << "Total force (fluid2) = "
            << totalTractionForce << endl;

//         Info << "Total pressure force (fluid) = "
//             << totalPressureForce << endl;
    }

    // Totla force at the solid side of the interface
    {
        const vectorField& p =
            solidMesh().faceZones()[solidZone2Index_]().localPoints();

        const faceList& f =
            solidMesh().faceZones()[solidZone2Index_]().localFaces();

        vectorField S(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
        }

        vector totalTractionForce =
            sum(solidZone2TotalTraction*mag(S));

//         vector totalPressureForce =
//             sum(solidZonePressure_*S);

        Info << "Total force (solid) = "
            << totalTractionForce << endl;

//         Info << "Total pressure force (solid) = "
//             << totalPressureForce << endl;
    }
}

void Foam::fluidStructureInterface::updateWeakForce()
{
    if (coupled())
    {
        Info << "Setting weak traction on solid patch" << endl;

        predictedSolidZoneTraction_ =
            2*solidZoneTraction_ - solidZoneTractionPrev_;

        stress().setTraction
        (
            solidPatchIndex(),
            solidZoneIndex(),
            predictedSolidZoneTraction_
        );
    }
}


void Foam::fluidStructureInterface::updateWeakTraction()
{
    Info << "Update weak traction on solid patch" << endl;

    solidZoneTractionPrev_ = solidZoneTraction_;

    // Calc fluid traction

    const vectorField& p =
        fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();
    const faceList& f =
        fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

    vectorField n(f.size(), vector::zero);
    forAll(n, faceI)
    {
        n[faceI] = f[faceI].normal(p);
        n[faceI] /= mag(n[faceI]);
    }

    vectorField fluidZoneTraction =
        flow().faceZoneViscousForce
        (
            fluidZoneIndex(),
            fluidPatchIndex()
        )
      - flow().faceZonePressureForce(fluidZoneIndex(), fluidPatchIndex())*n;

    vectorField fluidZoneTractionAtSolid =
        ggiInterpolator().masterToSlave
        (
            -fluidZoneTraction
        );

    scalar beta_ = relaxationFactor_;

    solidZoneTraction_ =
        beta_*fluidZoneTractionAtSolid
      + (1.0-beta_)*predictedSolidZoneTraction_;


    // Total force at the fluid side of the interface
    {
        const vectorField& p =
            fluidMesh().faceZones()[fluidZoneIndex_]().localPoints();

        const faceList& f =
            fluidMesh().faceZones()[fluidZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
        }

        vector totalTractionForce = sum(fluidZoneTraction*mag(S));

        Info << "Total force (fluid) = "
            << totalTractionForce << endl;
    }

    // Totla force at the solid side of the interface
    {
        const vectorField& p =
            solidMesh().faceZones()[solidZoneIndex_]().localPoints();

        const faceList& f =
            solidMesh().faceZones()[solidZoneIndex_]().localFaces();

        vectorField S(f.size(), vector::zero);

        forAll(S, faceI)
        {
            S[faceI] = f[faceI].normal(p);
        }

        vector totalTractionForce =
            sum(fluidZoneTractionAtSolid*mag(S));

        Info << "Total force (solid) = "
            << totalTractionForce << endl;
    }
}


void Foam::fluidStructureInterface::predictAndUpdateForce()
{
    if (coupled())
    {
        Info << "Setting traction on solid patch using prediction" << endl;

        stress().setPressure
        (
            solidPatchIndex(),
            solidZoneIndex(),
            stress().predictPressure(solidPatchIndex(), solidZoneIndex())
        );

        stress().setTraction
        (
            solidPatchIndex(),
            solidZoneIndex(),
            stress().predictTraction(solidPatchIndex(), solidZoneIndex())
        );
    }
}


void Foam::fluidStructureInterface::evolveStress()
{
//     if (closedFluidDomain())
//     {
//         DynamicList<scalar> p0;
//         DynamicList<scalar> dV;

//         scalar requiredVolumeIncrement = 0;
//         forAll(flow().U().boundaryField(), patchI)
//         {
//             if (patchI != fluidPatchIndex())
//             {
//                 requiredVolumeIncrement +=
//                     sum
//                     (
//                         flow().U().boundaryField()[patchI]
//                       & fluidMesh().Sf().boundaryField()[patchI]
//                     )
//                    *runTime().deltaT().value();
//             }
//         }

//         scalar volumeIncrement = 0;

//         label nIter = 0;

//         do
//         {
//             // Calc ref. pressure increment

//             if (dV.size() == 1)
//             {
//                 refPressureIncrement_ += compressibility_*dV[0];
//             }
//             else if (dV.size() >= 2)
//             {
//                 label i = p0.size() - 1;

//                 refPressureIncrement_ =
//                     p0[i-1]
//                   - dV[i-1]*(p0[i] - p0[i-1])
//                    /(dV[i] - dV[i-1] + SMALL);
//             }

//             p0.append(refPressureIncrement_);

//             stress().setPressure
//             (
//                 solidPatchIndex(),
//                 solidZoneIndex(),
//                 solidZonePressure_ + refPressure() + refPressureIncrement()
//             );


//             // Solve solid
//             stress().evolve();

//             // Calculate volume increment

//             const labelList& meshPoints =
//                 solidMesh().faceZones()[solidZoneIndex()]().meshPoints();

//             const faceList& localFaces =
//                 solidMesh().faceZones()[solidZoneIndex()]().localFaces();

//             const vectorField& localPoints =
//                 solidMesh().faceZones()[solidZoneIndex()]().localPoints();


//             vectorField oldLocalPoints =
//                 localPoints
//               + vectorField
//                 (
//                     stress().pointD().oldTime(),
//                     meshPoints
//                 );

//             vectorField newLocalPoints =
//                 localPoints
//               + vectorField
//                 (
//                     stress().pointD(),
//                     meshPoints
//                 );

//             volumeIncrement = 0;

//             forAll(localFaces, faceI)
//             {
//                 volumeIncrement +=
//                     localFaces[faceI].sweptVol
//                     (
//                         oldLocalPoints,
//                         newLocalPoints
//                     );
//             }

//             volumeIncrement -= requiredVolumeIncrement;

//             dV.append(volumeIncrement);
//         }
//         while(mag(volumeIncrement) > SMALL*10 && ++nIter < 10);


//         Info << "Solid volume increment: " << volumeIncrement << endl;
//         Info << "Ref pressure: " << refPressure_ << endl;
//         Info << "Ref pressure increment: " << refPressureIncrement_ << endl;
//         Info << "Calculated compressibility = "
//             << (refPressureIncrement() - p0[0])/(dV[0] + SMALL) << endl;
//     }
//     else
//     {
//         stress().evolve();
//     }

    stress().evolve();
}


Foam::scalar Foam::fluidStructureInterface::updateResidual()
{
    vectorField solidZonePointsDisplAtSolid =
        stress().faceZonePointDisplacementIncrement(solidZoneIndex());

    solidZonePointsDispl() =
//         solidToFluid().pointInterpolate
        ggiInterpolator().slaveToMasterPointInterpolate
        (
            solidZonePointsDisplAtSolid
        );

    residualPrev() = residual();

    residual() = solidZonePointsDispl() - fluidZonePointsDispl();

//     const scalarField& minLe = minEdgeLength();

//     scalar residualNorm = gMax(mag(residual())/(minLe + SMALL));

    scalar residualNorm = ::sqrt(sum(magSqr(residual())));

    vectorField solidZone2PointsDisplAtSolid =
        stress().faceZonePointDisplacementIncrement(solidZone2Index());

    solidZone2PointsDispl() =
//         solidToFluid().pointInterpolate
        ggiInterpolator2().slaveToMasterPointInterpolate
        (
            solidZone2PointsDisplAtSolid
        );

    residual2Prev() = residual2();

    residual2() = solidZone2PointsDispl() - fluid2ZonePointsDispl();

//     const scalarField& minLe = minEdgeLength();

//     scalar residualNorm = gMax(mag(residual())/(minLe + SMALL));

    scalar residual2Norm = ::sqrt(sum(magSqr(residual2())));

//     Info << "Current fsi residual norm: " << residualNorm << endl;
    if (residualNorm < residual2Norm)
    {
         residualNorm =residual2Norm;
    }

    if (residualNorm > maxResidualNorm_)
    {
        maxResidualNorm_ = residualNorm;
    }

//     Info << "Current fsi max residual norm: " << maxResidualNorm_ << endl;

    residualNorm /= maxResidualNorm_ + SMALL;

    Info << "Current fsi relative residual norm: " << residualNorm << endl;

    return residualNorm;
}


// bool Foam::fluidStructureInterface::read()
// {
//     if (regIOobject::read())
//     {
//         flowProperties_ = subDict(type() + "Coeffs");

//         return true;
//     }
//     else
//     {
//         return false;
//     }
// }


// ************************************************************************* //
