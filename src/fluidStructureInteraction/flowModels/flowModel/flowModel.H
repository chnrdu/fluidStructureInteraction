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
    flowModel

Description
    Virtual base class for flow models

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    flowModel.C
    newFlowModel.C

\*---------------------------------------------------------------------------*/

#ifndef flowModel_H
#define flowModel_H

#include "fvMesh.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class flowModel Declaration
\*---------------------------------------------------------------------------*/

class flowModel
:
    public IOdictionary
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Flow properties dictionary
        dictionary flowProperties_;

        //- Global face zones
        mutable labelList* globalFaceZonesPtr_;

        //- Global to local face zone poin maping
        mutable labelListList* globalToLocalFaceZonePointMapPtr_;

    // Private Member Functions

        //- Find global face zones
        void calcGlobalFaceZones() const;

        //- Calculate global to local point addressing
        void calcGlobalToLocalFaceZonePointMap() const;

        //- Disallow default bitwise copy construct
        flowModel(const flowModel&);

        //- Disallow default bitwise assignment
        void operator=(const flowModel&);


protected:

    // Protected member functions

        //- Return flow properties dictionary
        const dictionary& flowProperties() const
        {
            return flowProperties_;
        }


public:

    //- Runtime type information
    TypeName("flowModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            flowModel,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from components
        flowModel
        (
            const word& type,
            const fvMesh& mesh
        );


    // Selectors

        //- Select constructed from mesh
        static autoPtr<flowModel> New(const fvMesh& mesh);


    // Destructor

        virtual ~flowModel();


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

            //- Return velocity field
            virtual const volVectorField& U() const = 0;

            //- Return pressure field
            virtual const volScalarField& p() const = 0;

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const = 0;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const = 0;

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> faceZoneViscousForce
            (
                const label zoneID,
                const label patchID
            ) const = 0;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> faceZonePressureForce
            (
                const label zoneID,
                const label patchID
            ) const = 0;

            //- Face zone effective dynamic viscosity
            virtual tmp<scalarField> faceZoneMuEff
            (
                const label zoneID,
                const label patchID                
            ) const = 0;

            //- Get global face zones list
            const labelList& globalFaceZones() const;

            //- Get global to local face zone point addressing
            const labelListList& globalToLocalFaceZonePointMap() const;

        // Edit

            //- Evolve the flow model
            virtual void evolve() = 0;

            //- Read flowProperties dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
