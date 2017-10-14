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

Class
    fluidStructureInterface

Description
    Fluid-structure interface class

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    fluidStructureInterface.C

\*---------------------------------------------------------------------------*/

#ifndef fluidStructureInterface_H
#define fluidStructureInterface_H

#include "flowModel.H"
#include "stressModel.H"
#include "IOdictionary.H"
#include "patchToPatchInterpolation.H"
#include "dynamicFvMesh.H"
#include "ggiInterpolation.H"
#include "extendedGgiInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidStructureInterface Declaration
\*---------------------------------------------------------------------------*/

class fluidStructureInterface
:
    public IOdictionary
{
    // Private data

        //- Fluid (flow) mesh
        dynamicFvMesh& fluidMesh_;

        //- Flow model
        autoPtr<flowModel> flow_;

        //- Solid (stress) mesh
        fvMesh& solidMesh_;

        //- Stress model
        autoPtr<stressModel> stress_;

        //- Solid patch index
        label solidPatchIndex_;

        //- Solid face zone index
        label solidZoneIndex_;

        //- Fluid patch index
        label fluidPatchIndex_;

        //- Fluid face zone index
        label fluidZoneIndex_;

        //- Solid face zone current points
        mutable vectorField* currentSolidZonePointsPtr_;

        //- Solid face zone primitive patch in current configuration
        mutable PrimitivePatch<face, List, const pointField&>* 
        currentSolidZonePatchPtr_;

        //- Fluid zone to solid zone interpolator
        mutable zoneToZoneInterpolation* fluidToSolidPtr_;

        //- Ggi zone-to-zone interpolation
        mutable ggiZoneInterpolation* ggiFluidToSolidPtr_;

        //- Ggi zone-to-zone interpolation
        mutable extendedGgiZoneInterpolation* ggiInterpolatorPtr_;

        //- Solid zone to fluid zone interpolator
        mutable zoneToZoneInterpolation* solidToFluidPtr_;

        //- Coupling scheme
        word couplingScheme_;

        //- Relaxation factor
        scalar relaxationFactor_;

        //- Relaxation factor
        scalar aitkenRelaxationFactor_;

        //- Outer correction loop stoping tolerance
        scalar outerCorrTolerance_;
    
        //- Maximal number of outer correctors
        label nOuterCorr_;

        //- Fsi coupling on/off
        Switch coupled_;

        //- Predictor on/off
        Switch predictor_;

        //- Coupling reuse
        label couplingReuse_;

        //- Interface deformation limit
        scalar interfaceDeformationLimit_;

        //- Fluid zone point displacement
        vectorField fluidZonePointsDispl_;

        //- Fluid zone ref. point displacement
        vectorField fluidZonePointsDisplRef_;

        //- Fluid zone previous point displacement
        vectorField fluidZonePointsDisplPrev_;

        //- Solid zone point displacement
        vectorField solidZonePointsDispl_;

        //- Solid zone ref. point displacement
        vectorField solidZonePointsDisplRef_;

        //- Solid zone pressure
        scalarField solidZonePressure_;

        //- Solid zone traction (pressure + vicous)
        vectorField solidZoneTraction_;

        //- Solid zone traction (pressure + vicous)
        vectorField solidZoneTractionPrev_;

        //- Solid zone traction (pressure + vicous)
        vectorField predictedSolidZoneTraction_;

        //- Current fsi residual
        vectorField residual_;

        //- Previous fsi residual
        vectorField residualPrev_;

        //- Maximal resudual norm
        scalar maxResidualNorm_;

        //- Outer corrector
        label outerCorr_;

//         //- Is it fluid domain pure Dirichlet (witout outlets)?
//         Switch closedFluidDomain_;

//         //- Reference pressure 
//         scalar refPressure_;

//         //- Reference pressure 
//         scalar refPressureIncrement_;

//         //- Current time index
//         label timeIndex_;

//         //- Copressibility
//         scalar compressibility_;

        //- Interpolator update frequency
        label interpolatorUpdateFrequency_;


        //- IQN-ILS coupling fields
        DynamicList<vectorField> fluidPatchPointsV_;
        DynamicList<vectorField> fluidPatchPointsW_;
        DynamicList<scalar> fluidPatchPointsT_;

        //- Accumulated fluid side interface displacement
        mutable vectorIOField* accumulatedFluidInterfaceDisplacementPtr_;

        //- Min edge length for interface points at fluid side
        mutable scalarField* minEdgeLengthPtr_;

    // Private Member Functions

        //- Calculate current solid zone points
        void calcCurrentSolidZonePoints() const;

        //- Calculate current solid zone primitive patch
        void calcCurrentSolidZonePatch() const;

        //- Calculate fluid to solid interpolator
        void calcFluidToSolidInterpolator() const;

        //- Calculate fluid to solid ggi interpolator
        void calcGgiFluidToSolidInterpolator() const;

        //- Calculate fluid to solid ggi interpolator
        void calcGgiInterpolator() const;

        //- Calculate fluid to solid interpolator
        void calcSolidToFluidInterpolator() const;

        //- Accumulated fluid interface displacement
        void calcAccumulatedFluidInterfaceDisplacement() const;

        //- Calculate minimal edge lengths
        void calcMinEdgeLength() const;

        //- Return accumulated interface displacement
        vectorIOField& accumulatedFluidInterfaceDisplacement();

        //- Return minimal edge length
        const scalarField& minEdgeLength() const;

        //- Disallow default bitwise copy construct
        fluidStructureInterface(const fluidStructureInterface&);

        //- Disallow default bitwise assignment
        void operator=(const fluidStructureInterface&);


protected:

    // Protected member functions

public:

    //- Runtime type information
    TypeName("fluidStructureInterface");


    // Declare run-time constructor selection table


    // Constructors

        //- Construct from components
        fluidStructureInterface
        (
            dynamicFvMesh& fluidMesh,
            fvMesh& solidMesh
        );


    // Selectors


    // Destructor

        virtual ~fluidStructureInterface();


    // Member Functions

        // Access

            //- Return fluid mesh
            const dynamicFvMesh& fluidMesh() const
            {
                return fluidMesh_;
            }

            //- Return solid mesh
            const fvMesh& solidMesh() const
            {
                return solidMesh_;
            }

            //- Return time
            const Time& runTime() const
            {
                return fluidMesh_.time();
            }

            //- Return flow model
            const flowModel& flow() const
            {
                return flow_();
            }

            //- Return flow model
            flowModel& flow()
            {
                return flow_();
            }

            //- Return stress model
            const stressModel& stress() const
            {
                return stress_();
            }

            //- Return stress model
            stressModel& stress()
            {
                return stress_();
            }

            //- Return solid patch index
            label solidPatchIndex() const
            {
                return solidPatchIndex_;
            }

            //- Return solid face zone index
            label solidZoneIndex() const
            {
                return solidZoneIndex_;
            }

            //- Return fluid patch index
            label fluidPatchIndex() const
            {
                return fluidPatchIndex_;
            }

            //- Return fluid face zone index
            label fluidZoneIndex() const
            {
                return fluidZoneIndex_;
            }

            //- Return current solid zone points
            const vectorField& currentSolidZonePoints() const;

            //- Return current solid zone patch
            const PrimitivePatch<face, List, const pointField&>& 
            currentSolidZonePatch() const;

            //- Return fluid to solid interpolator
            const zoneToZoneInterpolation& fluidToSolid() const;

            //- Return fluid to solid interpolator
            const ggiZoneInterpolation& ggiFluidToSolid() const;

            //- Return fluid to solid interpolator
            const extendedGgiZoneInterpolation& ggiInterpolator() const;

            //- Return fluid to solid interpolator
            const zoneToZoneInterpolation& solidToFluid() const;

            //- Return coupling scheme
            const word& couplingScheme() const
            {
                return couplingScheme_;
            }

            //- Return relaxation factor
            scalar relaxationFactor() const
            {
                return relaxationFactor_;
            }

            //- Return relaxation factor
            scalar& aitkenRelaxationFactor()
            {
                return aitkenRelaxationFactor_;
            }

            //- Return outer corrector loop tolerance
            scalar outerCorrTolerance() const
            {
                return outerCorrTolerance_;
            }

            //- Return max numter of outer correctors
            label nOuterCorr() const
            {
                return nOuterCorr_;
            }

            //- Is it fluid and structure coupled
            const Switch& coupled() const
            {
                return coupled_;
            }

            //- Is it fluid and structure coupled
            const Switch& predictor() const
            {
                return predictor_;
            }

            //- Is it fluid and structure coupled
            Switch& coupled()
            {
                return coupled_;
            }

            //- Is it fluid and structure coupled
            label couplingReuse() const
            {
                return couplingReuse_;
            }

            //- Return relaxation factor
            scalar interfaceDeformationLimit() const
            {
                return interfaceDeformationLimit_;
            }

            //- Return fluid zone point displacement
            vectorField& fluidZonePointsDispl()
            {
                return fluidZonePointsDispl_;
            }

            //- Return fluid zone ref. point displacement
            vectorField& fluidZonePointsDisplRef()
            {
                return fluidZonePointsDisplRef_;
            }

            //- Return fluid zone previous point displacement
            vectorField& fluidZonePointsDisplPrev()
            {
                return fluidZonePointsDisplPrev_;
            }

            //- Return solid zone point displacement
            vectorField& solidZonePointsDispl()
            {
                return solidZonePointsDispl_;
            }

            //- Return solid zone ref. point displacement
            vectorField& solidZonePointsDisplRef()
            {
                return solidZonePointsDisplRef_;
            }

            //- Return solid zone pressure
            scalarField& solidZonePressure()
            {
                return solidZonePressure_;
            }

            //- Return solid zone total traction
            vectorField& solidZoneTraction()
            {
                return solidZoneTraction_;
            }

            //- Return solid zone total traction
            vectorField& solidZoneTractionPrev()
            {
                return solidZoneTractionPrev_;
            }

            //- Return current fsi residual
            vectorField& residual()
            {
                return residual_;
            }

            //- Return previous fsi residual
            vectorField& residualPrev()
            {
                return residualPrev_;
            }

            //- Return current outer iteration
            label& outerCorr()
            {
                return outerCorr_;
            }

            //- Return current outer iteration
            const label& outerCorr() const
            {
                return outerCorr_;
            }

//             //- Is it fluid domain pure Dirichlet
//             Switch closedFluidDomain() const
//             {
//                 return closedFluidDomain_;
//             }

//             //- Ref. pressure used in case of closed fluid domain
//             scalar refPressure()
//             {
//                 return refPressure_;
//             }

//             //- Ref. pressure used in case of closed fluid domain
//             scalar refPressureIncrement()
//             {
//                 return refPressureIncrement_;
//             }

//             //- Compressibility used in case of closed fluid domain
//             scalar compressibility()
//             {
//                 return compressibility_;
//             }

        // Edit
    
            //- Initialize fields
            void initializeFields();

            //- Initialize fields
            void updateInterpolator();

            //- Calculate interface displacement
            void updateDisplacement();

            //- Calculate interface displacement
            void updateWeakDisplacement();

            //- Calculate interface displacement
            void updateDisplacementUsingIQNILS();

            //- Move fluid mesh
            void moveFluidMesh();

            //- Update interface force
            void updateForce();

            //- Update interface force
            void updateWeakForce();

            //- Update interface force
            void updateWeakTraction();

            //- Update interface force
            void predictAndUpdateForce();

            //- Solve solid
            void evolveStress();

            //- Update interface force
            scalar updateResidual();

//             //- Evolve the flow model
//             virtual void evolve() = 0;

//             //- Read dictionary
//             virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
