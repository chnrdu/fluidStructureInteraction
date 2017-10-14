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
    icoFlow

Description
    icoFoam flow model

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    icoFlow.C

\*---------------------------------------------------------------------------*/

#ifndef icoFlow_H
#define icoFlow_H

#include "flowModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flowModels
{

/*---------------------------------------------------------------------------*\
                        Class icoFlow Declaration
\*---------------------------------------------------------------------------*/

class icoFlow
:
    public flowModel
{
    // Private data

        //- Velocity field
        volVectorField U_;

        //- Pressure field
        volScalarField p_;

        //- Pressure field
        volVectorField gradp_;

        //- Flux field
        surfaceScalarField phi_;

        // Transport properties dictionary
        IOdictionary transportProperties_;

        //- Kinematic viscosity
        dimensionedScalar nu_;

        //- Density
        dimensionedScalar rho_;

//         //-
//         scalarField nGradUn_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        icoFlow(const icoFlow&);

        //- Disallow default bitwise assignment
        void operator=(const icoFlow&);


public:

    //- Runtime type information
    TypeName("icoFlow");

    // Constructors

        //- Construct from components
        icoFlow(const fvMesh& mesh);

//         //- Construct from components
//         icoFlow
//         (
//             const word& type,
//             const fvMesh& mesh
//         );

    // Destructor

        virtual ~icoFlow()
        {}


    // Member Functions

        // Access

            //- Return velocity field
            virtual const volVectorField& U() const;

            //- Return velocity field
            volVectorField& U()
            {
                return U_;
            }

            //- Return pressure field
            virtual const volScalarField& p() const;

            //- Return pressure field
            volScalarField& p()
            {
                return p_;
            }

            //- Return pressure gradient
            volVectorField& gradp()
            {
                return gradp_;
            }

            //- Return flux field
            surfaceScalarField& phi()
            {
                return phi_;
            }

            //- Return kinematic viscosity
            const dimensionedScalar& nu() const
            {
                return nu_;
            }

            //- Density
            const dimensionedScalar& rho()
            {
                return rho_;
            }

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const;

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> faceZoneViscousForce
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> faceZonePressureForce
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone effective dynamic viscosity
            virtual tmp<scalarField> faceZoneMuEff
            (
                const label zoneID,
                const label patchID                
            ) const;

        // Edit

            //- Evolve the flow model
            virtual void evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
