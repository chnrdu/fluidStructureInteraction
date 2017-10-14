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

stressModels/pRveUnsIncrTotalLagrangianStress/pRveUnsIncrTotalLagrangianStress.C    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    pRveUnsTotalLagrangianStress

Description
    

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    pRveUnsTotalLagrangianStress.C

\*---------------------------------------------------------------------------*/

#ifndef pRveUnsTotalLagrangianStress_H
#define pRveUnsTotalLagrangianStress_H

#include "unsTotalLagrangianStress.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stressModels
{

/*---------------------------------------------------------------------------*\
                        Class pRveUnsTotalLagrangianStress Declaration
\*---------------------------------------------------------------------------*/

class pRveUnsTotalLagrangianStress
:
    public unsTotalLagrangianStress
{
    // Private data

        //- Average deformation gradient
        tensor avgDeformationGradient_;

        //- Total point displacement field
        pointVectorField totPointD_;

        //- Total stress field
        volSymmTensorField totSigma_;

        //- Total strain field
        volSymmTensorField totEpsilon_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pRveUnsTotalLagrangianStress(const pRveUnsTotalLagrangianStress&);

        //- Disallow default bitwise assignment
        void operator=(const pRveUnsTotalLagrangianStress&);


public:

    //- Runtime type information
    TypeName("pRveUnsTotalLagrangianStress");

    // Constructors

        //- Construct from components
        pRveUnsTotalLagrangianStress(const fvMesh& mesh);


    // Destructor

        virtual ~pRveUnsTotalLagrangianStress()
        {}


    // Member Functions

        // Access

            //- Return average deformation gradient
            const tensor& avgDeformationGradient() const
            {
                return avgDeformationGradient_;
            }

            //- Return average deformation gradient
            tensor& avgDeformationGradient()
            {
                return avgDeformationGradient_;
            }

        // Edit

            //- Update total fields
            virtual void updateTotalFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
