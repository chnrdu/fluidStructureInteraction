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
    stressModel

Description
    Virtual base class for stress analysis models

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    stressModel.C
    newStressModel.C

\*---------------------------------------------------------------------------*/

#ifndef stressModel_H
#define stressModel_H

#include "fvMesh.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "constitutiveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class stressModel Declaration
\*---------------------------------------------------------------------------*/

class stressModel
:
    public IOdictionary
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Flow properties dictionary
        dictionary stressProperties_;

        //- Global face zones
        mutable labelList* globalFaceZonesPtr_;

        //- Global to local face zone poin maping
        mutable labelListList* globalToLocalFaceZonePointMapPtr_;

    // Private Member Functions

        //- Find global face zones
        void calcGlobalFaceZones() const;

        //- Find global to local point addressing
        void calcGlobalToLocalFaceZonePointMap() const;
    
        //- Disallow default bitwise copy construct
        stressModel(const stressModel&);

        //- Disallow default bitwise assignment
        void operator=(const stressModel&);


protected:

    // Protected member functions

        //- Return flow properties dictionary
        dictionary& stressProperties()
        {
            return stressProperties_;
        }

public:

    //- Runtime type information
    TypeName("stressModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            stressModel,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from components
        stressModel
        (
            const word& type,
            const fvMesh& mesh
        );


    // Selectors

        //- Select constructed from mesh
        static autoPtr<stressModel> New(const fvMesh& mesh);


    // Destructor

        virtual ~stressModel();


    // Member Functions

        // Access

            //- Return mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return time
            const Time& runTime() const
            {
                return mesh_.time();
            }

            //- Return flow properties dictionary
            const dictionary& stressProperties() const
            {
                return stressProperties_;
            }

            //- Thermal stress
            bool thermalStress() const;

            //- Thermal stress
            virtual const constitutiveModel& rheology() const = 0;

            //- Return point displacement field
            virtual const pointVectorField& pointD() const = 0;

            //- Return point velocity field
            virtual vector pointU(label pointID) const = 0;

            //- Patch point displacement increment
            virtual tmp<vectorField> patchPointDisplacementIncrement
            (
                const label patchID
            ) const = 0;

            //- Face zone point displacement increment
            virtual tmp<vectorField> faceZonePointDisplacementIncrement
            (
                const label zoneID
            ) const = 0;

            //- Face zone velocity gradient
            virtual tmp<tensorField> faceZoneSurfaceGradientOfVelocity
            (
                const label zoneID,
                const label patchID
            ) const = 0;

            //- Face zone current normal
            virtual tmp<vectorField> faceZoneNormal
            (
                const label zoneID,
                const label patchID
            ) const = 0;

            //- Face zone primitive patch in current configuration
            virtual tmp<vectorField> 
            currentFaceZonePoints(const label zoneID) const = 0;

            //- Get global face zones list
            const labelList& globalFaceZones() const;

            //- Get global to local face zone point addressing
            const labelListList& globalToLocalFaceZonePointMap() const;

        // Edit

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label patchID, 
                const vectorField& traction
            ) = 0;

            //- Set traction at specified patch
            virtual void setPressure
            (
                const label patchID, 
                const scalarField& pressure
            ) = 0;

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label patchID,
                const label zoneID,
                const vectorField& faceZoneTraction
            ) = 0;

            //- Set traction at specified patch
            virtual void setPressure
            (
                const label patchID,
                const label zoneID,
                const scalarField& faceZonePressure
            ) = 0;

            //- Predict traction at specified zone
            virtual tmp<vectorField> predictTraction
            (
                const label patchID,
                const label zoneID
            ) = 0;

            //- Predict pressure at specified zone
            virtual tmp<scalarField> predictPressure
            (
                const label patchID,
                const label zoneID
            ) = 0;
    
            //- Evolve the stress model
            virtual bool evolve() = 0;

            //- Predict the stress model
            virtual void predict() = 0;

            //- Smooth solution
            virtual scalar smooth
            (
                const volVectorField& source,
                label nCorrectors = 0
            ) = 0;

            //- Smooth solution
            virtual scalar smooth
            (
                const volVectorField& source,
                const volVectorField& refSolution,
                label nCorrectors = 0
            ) = 0;

            //- Return reference to current solution
            virtual const volVectorField& solution() const = 0;

            //- Return reference to current solution field
            virtual volVectorField& solution() = 0;

            //- Initialise solution
            virtual void initialise
            (
                const volVectorField& sol,
                bool consistenBoundaryField = false
            ) = 0;

            //- Return current residual
            virtual tmp<volVectorField> residual
            (
                const volVectorField& source
            ) = 0;

            //- Update total fields
            virtual void updateTotalFields()
            {}

            //- Read flowProperties dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
