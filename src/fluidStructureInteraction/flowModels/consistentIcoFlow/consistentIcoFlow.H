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
    consistentIcoFlow

Description
    icoFoam flow model

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    consistentIcoFlow.C

\*---------------------------------------------------------------------------*/

#ifndef consistentIcoFlow_H
#define consistentIcoFlow_H

#include "icoFlow.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flowModels
{

/*---------------------------------------------------------------------------*\
                        Class consistentIcoFlow Declaration
\*---------------------------------------------------------------------------*/

class consistentIcoFlow
:
    public icoFlow
{
    // Private data

        //- Face area vectors
        mutable surfaceVectorField* SfPtr_;
    
    // Private Member Functions

        //- Make face area vectors
        void makeSf() const;

        //- Make face area vectors
        void updateSf();

        //- Disallow default bitwise copy construct
        consistentIcoFlow(const consistentIcoFlow&);

        //- Disallow default bitwise assignment
        void operator=(const consistentIcoFlow&);


public:

    //- Runtime type information
    TypeName("consistentIcoFlow");

    // Constructors

        //- Construct from components
        consistentIcoFlow(const fvMesh& mesh);


    // Destructor

        virtual ~consistentIcoFlow()
        {}


    // Member Functions

        // Access

            //- Return face area vectors
            const surfaceVectorField& Sf() const;

            //- Return face area vectors
            surfaceVectorField& Sf();

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
