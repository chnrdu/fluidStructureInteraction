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
    Foam::quadraticReconstruction

Description
    Foam::quadraticReconstruction

SourceFiles
    quadraticReconstruction.C
    quadraticReconstructionTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticReconstruction_H
#define quadraticReconstruction_H

#include "regIOobject.H"
#include "PtrList.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "scalarList.H"
#include "scalarMatrices.H"
#include "pointFieldsFwd.H"

#include "FieldField.H"
#include "IOPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class quadraticReconstruction Declaration
\*---------------------------------------------------------------------------*/

class quadraticReconstruction
:
    public regIOobject
{
    // Private data

        //- Reference to finite volume mesh
        const fvMesh& mesh_;

        //- Cell cell addressig
        mutable labelListList* cellCellsPtr_;

        //- Cell boundary face addressig
        mutable labelListList* cellFacesPtr_;

        //- Cell boundary face addressig
        mutable labelListList* cellConstrainedFacesPtr_;

        //- Interpolation matrices
        mutable PtrList<scalarRectangularMatrix> invLsMatrices_;

        //- Cell based reference length
        mutable scalarField* refLenghtsPtr_;

        //- Interpolation vol-to-point interpolation 
        //  weighting factor array
        mutable scalarListList* pointWeightsPtr_;


    // Private member functions

        //- Make cell boundary face addressing
        void makeCellCells() const;

        //- Make cell boundary face addressing
        void makeCellFaces() const;

        //- Make cell boundary face addressing
        void makeCellConstrainedFaces() const;

        //- Make interpolation matrix
        void makeInvLsMatrices() const;

        //- Make reference lengths
        void makeRefLenghts() const;

        //- Make weights
        void makePointWeights() const;

        //- Disallow default bitwise copy construct
        quadraticReconstruction
        (
            const quadraticReconstruction&
        );

        //- Disallow default bitwise assignment
        void operator=(const quadraticReconstruction&);


public:

    // Declare name of the class and its debug switch
    ClassName("quadraticReconstruction");


    // Constructors

        //- Constructor given fvMesh
        explicit quadraticReconstruction(const fvMesh&);


    // Destructor

        ~quadraticReconstruction();


    // Member functions

        // Access

            //- Get point boundary face addressing
            const labelListList& cellCells() const;

            //- Get point boundary face addressing
            const labelListList& cellFaces() const;

            //- Get point boundary face addressing
            const labelListList& cellConstrainedFaces() const;

            //- Get interpolation matrices
            const PtrList<scalarRectangularMatrix>& invLsMatrices() const;

            //- Get reference lengths
            const scalarField& refLenghts() const;

            //- Get weights
            const scalarListList& pointWeights() const;

        // Edit

            //- Correct reconstruction
            bool movePoints();

        // Writing

            //- Virtual writaData function.
            virtual bool writeData(Ostream&) const
            {
                return true;
            }

    // Interpolation functions

        //- Calc reconstruction
        template<class Type>
        Type reconstruct
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const label cell,
            const vector& point
        ) const;

        //- Calc reconstruction
        template<class Type>
        Type reconstruct
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const label cell,
            const vector& point,
            const PtrList<Field<Type> >& coeffs
        ) const;

        //- Calc reconstruction
        template<class Type>
        tmp<Field<Type> > reconstruct
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const labelList&,
            const Field<vector>&
        ) const;

        //- Calc difference
        template<class Type>
        tmp<Field<Type> > difference
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const labelList&,
            const Field<vector>&
        ) const;

        //- Calc difference
        template<class Type>
        tmp<Field<Type> > difference
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const labelList&,
            const Field<vector>&,
            const IOPtrList<Field<Type> >&
        ) const;

        //- Calc derivative
        template<class Type>
        tmp<Field<Type> > derivative
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const labelList&,
            const Field<vector>&,
            const Field<vector>&
        ) const;

        //- Calc coefficients
        template<class Type>
        tmp<FieldField<Field, Type> > coeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Calc coefficients
        template<class Type>
        void calcCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            IOPtrList<Field<Type> >&
        ) const;

        //- Cell to face interpolation
        template<class Type>
        tmp
        < 
            GeometricField<Type, fvsPatchField, surfaceMesh>
        > 
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<volVectorField> grad(const volScalarField& vf) const;
        tmp<volTensorField> grad(const volVectorField& vf) const;
        tmp<volTensorField> grad
        (
            const volVectorField& vf,
            const IOPtrList<Field<vector> >& Coeffs
        ) const;

        tmp<surfaceTensorField> sGrad(const volVectorField& vf) const;
        tmp<surfaceTensorField> sGrad
        (
            const volVectorField& vf,
            const IOPtrList<Field<vector> >& Coeffs
        ) const;
        tmp<surfaceVectorField> sGrad
        (
            const volScalarField& vf,
            const PtrList<Field<scalar> >& Coeffs
        ) const;

        void calcPointGrad
        (
            const volVectorField& vf,
            const IOPtrList<Field<vector> >& Coeffs,
            pointTensorField& pGradVf
        ) const;

        //- Interpolate from volField to pointField
        //  using inverse distance weighting
        template<class Type>
        void volToPointInterpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField<Type, pointPatchField, pointMesh>&,
            const IOPtrList<Field<Type> >& coeffs
        ) const;

        //- Interpolate from volField to pointField
        //  using inverse distance weighting
        template<class Type>
        void volToPointInterpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "quadraticReconstructionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
