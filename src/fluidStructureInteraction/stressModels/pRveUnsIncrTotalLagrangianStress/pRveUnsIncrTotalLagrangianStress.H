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
    pRveUnsIncrTotalLagrangianStress

Description
    

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    pRveUnsIncrTotalLagrangianStress.C

\*---------------------------------------------------------------------------*/

#ifndef pRveUnsIncrTotalLagrangianStress_H
#define pRveUnsIncrTotalLagrangianStress_H

#include "unsIncrTotalLagrangianStress.H"

// #include "volFields.H"
// #include "surfaceFields.H"
// #include "pointFields.H"
// #include "constitutiveModel.H"
// #include "leastSquaresVolPointInterpolation.H"
// #include "ITLMaterialInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace stressModels
{

/*---------------------------------------------------------------------------*\
                        Class pRveUnsIncrTotalLagrangianStress Declaration
\*---------------------------------------------------------------------------*/

class pRveUnsIncrTotalLagrangianStress
:
    public unsIncrTotalLagrangianStress
{
    // Private data

        //- Average strain increment
        symmTensor avgDEpsilon_;

        //- Average strain
        symmTensor avgEpsilon_;

//         //- Average deformation gradient increment
//         tensor avgDeformationGradientIncrement_;

//         //- Average deformation gradient increment
//         tensor avgDeformationGradient_;

        //- Total displacement
        volVectorField totD_;

        //- Total displacement
        pointVectorField totPointD_;

        //- Total stress field
        volSymmTensorField totSigma_;

        //- Total stress field
        surfaceSymmTensorField totSigmaf_;

    // Private Member Functions
    
        //- Disallow default bitwise copy construct
        pRveUnsIncrTotalLagrangianStress
        (
            const pRveUnsIncrTotalLagrangianStress&
        );

        //- Disallow default bitwise assignment
        void operator=(const pRveUnsIncrTotalLagrangianStress&);

public:

    //- Runtime type information
    TypeName("pRveUnsIncrTotalLagrangianStress");

    // Constructors

        //- Construct from components
        pRveUnsIncrTotalLagrangianStress(const fvMesh& mesh);


    // Destructor

        virtual ~pRveUnsIncrTotalLagrangianStress()
        {}


    // Member Functions

        // Access

//             //- Average deformation gradient increment
//             const tensor& avgDeformationGradientIncrement() const
//             {
//                 return avgDeformationGradientIncrement_;
//             }

//             //- Average deformation gradient increment
//             const tensor& avgDeformationGradient() const
//             {
//                 return avgDeformationGradient_;
//             }

            //- Average strain increment
            symmTensor avgDEpsilon() const
            {
                return avgDEpsilon_;
            }

            //- Average strain increment
            symmTensor avgEpsilon() const
            {
                return avgEpsilon_;
            }

        // Edit

//             //- Evolve the flow model
//             virtual bool evolve();

            //- Update total fields
            virtual void updateTotalFields();

        // Write

            //- Write using given format, version and compression
            virtual bool writeObject
            (
                IOstream::streamFormat,
                IOstream::versionNumber,
                IOstream::compressionType
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
