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

#include "quadraticReconstruction.H"
#include "volFields.H"
#include "volSurfaceMapping.H"
#include "zeroGradientFvPatchFields.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Type quadraticReconstruction::reconstruct
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const label cell,
    const vector& point
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::reconstruct("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "const label, "
            << "const vector&) : "
            << "reconstruct field in provided points"
            << endl;
    }

    Type reconValue = pTraits<Type>::zero;

    const Field<Type>& vfI = vf.internalField();

    label nCoeffs = 5;
    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();
    const labelListList& cCells = cellCells();
    const labelListList& cFaces = cellFaces();
    const labelListList& ccFaces = cellConstrainedFaces();

    const vectorField& C = mesh_.cellCentres();

    const scalarField& refL = refLenghts();

//     forAll(cells, cellI)
    {
        label curCell = cell;

        const scalarRectangularMatrix& curInvMatrix = invMatrices[curCell];

        const labelList& interpCells = cCells[curCell];
        const labelList& interpFaces = cFaces[curCell];
        const labelList& constrFaces = ccFaces[curCell];

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source
        (
            interpCells.size() 
          + interpFaces.size()
          + constrFaces.size(),
            pTraits<Type>::zero
        );

        label pointID = 0;

        for (label i=0; i<interpCells.size(); i++)
        {
            source[pointID++] = vfI[interpCells[i]] - vfI[curCell];
        }

        for (label i=0; i<interpFaces.size(); i++)
        {
            label faceID = interpFaces[i];
            label patchID = 
                mesh_.boundaryMesh().whichPatch(faceID);

            label start = mesh_.boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID++] = 
                vf.boundaryField()[patchID][localFaceID] - vfI[curCell];
        }

        for (label i=0; i<constrFaces.size(); i++)
        {
            label faceID = constrFaces[i];
            label patchID = 
                mesh_.boundaryMesh().whichPatch(faceID);

            label start = mesh_.boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID++] = 
                vf.boundaryField()[patchID][localFaceID] - vfI[curCell];
        }

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        vector r = (point - C[curCell])/refL[curCell];

        reconValue = 
            vfI[curCell]
          + coeffs[0]*r.x()
          + coeffs[1]*r.y()
          + coeffs[2]*r.x()*r.y()
          + coeffs[3]*sqr(r.x())
          + coeffs[4]*sqr(r.y());

        if (mesh_.nGeometricD() == 3)
        {
            reconValue += 
                coeffs[5]*r.z()
              + coeffs[6]*r.x()*r.z()
              + coeffs[7]*r.y()*r.z()
              + coeffs[8]*sqr(r.z());
        }
    }

    return reconValue;
}


template<class Type>
Type quadraticReconstruction::reconstruct
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const label cell,
    const vector& point,
    const PtrList<Field<Type> >& coeffs
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::reconstruct("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "const label, "
            << "const vector&, "
            << "const FieldField<Field, Type>&) : "
            << "reconstruct field in provided point"
            << endl;
    }

    Type reconValue = pTraits<Type>::zero;

    const Field<Type>& vfI = vf.internalField();

    label nCoeffs = 5;
    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    const vectorField& C = mesh_.cellCentres();

    const scalarField& refL = refLenghts();

    vector r = (point - C[cell])/refL[cell];

    reconValue = 
        vfI[cell]
      + coeffs[cell][0]*r.x()
      + coeffs[cell][1]*r.y()
      + coeffs[cell][2]*r.x()*r.y()
      + coeffs[cell][3]*sqr(r.x())
      + coeffs[cell][4]*sqr(r.y());

    if (mesh_.nGeometricD() == 3)
    {
        reconValue += 
            coeffs[cell][5]*r.z()
          + coeffs[cell][6]*r.x()*r.z()
          + coeffs[cell][7]*r.y()*r.z()
          + coeffs[cell][8]*sqr(r.z());
    }

    return reconValue;
}


template<class Type>
tmp<Field<Type> > quadraticReconstruction::reconstruct
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells,
    const vectorField& points
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::reconstruct("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "const labelList&, "
            << "const vectorField&) : "
            << "reconstruct field in provided points"
            << endl;
    }

    tmp<Field<Type> > tReconValues
    (
        new Field<Type>(cells.size(), pTraits<Type>::zero)
    );
    Field<Type>& reconValues = tReconValues();

    const Field<Type>& vfI = vf.internalField();

    label nCoeffs = 5;
    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();
    const labelListList& cCells = cellCells();
    const labelListList& cFaces = cellFaces();
    const labelListList& ccFaces = cellConstrainedFaces();

    const vectorField& C = mesh_.cellCentres();

    const scalarField& refL = refLenghts();

    forAll(cells, cellI)
    {
        label curCell = cells[cellI];

        const scalarRectangularMatrix& curInvMatrix = invMatrices[curCell];

        const labelList& interpCells = cCells[curCell];
        const labelList& interpFaces = cFaces[curCell];
        const labelList& constrFaces = ccFaces[curCell];

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source
        (
            interpCells.size() 
          + interpFaces.size()
          + constrFaces.size(),
            pTraits<Type>::zero
        );

        label pointID = 0;

        for (label i=0; i<interpCells.size(); i++)
        {
            source[pointID++] = vfI[interpCells[i]] - vfI[curCell];
        }

        for (label i=0; i<interpFaces.size(); i++)
        {
            label faceID = interpFaces[i];
            label patchID = 
                mesh_.boundaryMesh().whichPatch(faceID);

            label start = mesh_.boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID++] = 
                vf.boundaryField()[patchID][localFaceID] - vfI[curCell];
        }

        for (label i=0; i<constrFaces.size(); i++)
        {
            label faceID = constrFaces[i];
            label patchID = 
                mesh_.boundaryMesh().whichPatch(faceID);

            label start = mesh_.boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID++] = 
                vf.boundaryField()[patchID][localFaceID] - vfI[curCell];
        }

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        vector r = (points[cellI] - C[curCell])/refL[curCell];

        reconValues[cellI] =
            vfI[curCell]
          + coeffs[0]*r.x()
          + coeffs[1]*r.y()
          + coeffs[2]*r.x()*r.y()
          + coeffs[3]*sqr(r.x())
          + coeffs[4]*sqr(r.y());

        if (mesh_.nGeometricD() == 3)
        {
            reconValues[cellI] +=
                coeffs[5]*r.z()
              + coeffs[6]*r.x()*r.z()
              + coeffs[7]*r.y()*r.z()
              + coeffs[8]*sqr(r.z());
        }
    }

    return tReconValues;
}

template<class Type>
tmp<Field<Type> > quadraticReconstruction::difference
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells,
    const vectorField& delta
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::reconstruct("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "const labelList&, "
            << "const vectorField&) : "
            << "reconstruct field in provided points"
            << endl;
    }

    tmp<Field<Type> > tDiff
    (
        new Field<Type>(cells.size(), pTraits<Type>::zero)
    );
    Field<Type>& diff = tDiff();

    const Field<Type>& vfI = vf.internalField();

    label nCoeffs = 5;
    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();
    const labelListList& cCells = cellCells();
    const labelListList& cFaces = cellFaces();
    const labelListList& ccFaces = cellConstrainedFaces();

    const scalarField& refL = refLenghts();

    forAll(cells, cellI)
    {
        label curCell = cells[cellI];

        const scalarRectangularMatrix& curInvMatrix = invMatrices[curCell];

        const labelList& interpCells = cCells[curCell];
        const labelList& interpFaces = cFaces[curCell];
        const labelList& constrFaces = ccFaces[curCell];

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source
        (
            interpCells.size() 
          + interpFaces.size()
          + constrFaces.size(),
            pTraits<Type>::zero
        );

        label pointID = 0;

        for (label i=0; i<interpCells.size(); i++)
        {
            source[pointID++] = vfI[interpCells[i]] - vfI[curCell];
        }

        for (label i=0; i<interpFaces.size(); i++)
        {
            label faceID = interpFaces[i];
            label patchID = 
                mesh_.boundaryMesh().whichPatch(faceID);

            label start = mesh_.boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID++] = 
                vf.boundaryField()[patchID][localFaceID] - vfI[curCell];
        }

        for (label i=0; i<constrFaces.size(); i++)
        {
            label faceID = constrFaces[i];
            label patchID = 
                mesh_.boundaryMesh().whichPatch(faceID);

            label start = mesh_.boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID++] = 
                vf.boundaryField()[patchID][localFaceID] - vfI[curCell];
        }

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        vector r = delta[cellI]/refL[curCell];

        diff[cellI] = 
            coeffs[0]*r.x()
          + coeffs[1]*r.y()
          + coeffs[2]*r.x()*r.y()
          + coeffs[3]*sqr(r.x())
          + coeffs[4]*sqr(r.y());

        if (mesh_.nGeometricD() == 3)
        {
            diff[cellI] += 
                coeffs[5]*r.z()
              + coeffs[6]*r.x()*r.z()
              + coeffs[7]*r.y()*r.z()
              + coeffs[8]*sqr(r.z());
        }
    }

    return tDiff;
}

template<class Type>
tmp<Field<Type> > quadraticReconstruction::difference
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells,
    const vectorField& delta,
    const IOPtrList<Field<Type> >& coeffs
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::reconstruct("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "const labelList&, "
            << "const vectorField&) : "
            << "reconstruct field in provided points"
            << endl;
    }

    tmp<Field<Type> > tDiff
    (
        new Field<Type>(cells.size(), pTraits<Type>::zero)
    );
    Field<Type>& diff = tDiff();

//     const Field<Type>& vfI = vf.internalField();

    label nCoeffs = 5;
    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

//     const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();
//     const labelListList& cCells = cellCells();
//     const labelListList& cFaces = cellFaces();

    const scalarField& refL = refLenghts();

    forAll(cells, cellI)
    {
        label curCell = cells[cellI];

        vector r = delta[cellI]/refL[curCell];

        diff[cellI] = 
            coeffs[curCell][0]*r.x()
          + coeffs[curCell][1]*r.y()
          + coeffs[curCell][2]*r.x()*r.y()
          + coeffs[curCell][3]*sqr(r.x())
          + coeffs[curCell][4]*sqr(r.y());

        if (mesh_.nGeometricD() == 3)
        {
            diff[cellI] += 
                coeffs[curCell][5]*r.z()
              + coeffs[curCell][6]*r.x()*r.z()
              + coeffs[curCell][7]*r.y()*r.z()
              + coeffs[curCell][8]*sqr(r.z());
        }
    }

    return tDiff;
}


template<class Type>
tmp<Field<Type> > quadraticReconstruction::derivative
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const labelList& cells,
    const vectorField& delta,
    const vectorField& direction
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::reconstruct("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "const labelList&, "
            << "const vectorField&) : "
            << "reconstruct field in provided points"
            << endl;
    }

    tmp<Field<Type> > tDrv
    (
        new Field<Type>(cells.size(), pTraits<Type>::zero)
    );
    Field<Type>& drv = tDrv();

    const Field<Type>& vfI = vf.internalField();

    label nCoeffs = 5;
    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();
    const labelListList& cCells = cellCells();
    const labelListList& cFaces = cellFaces();
    const labelListList& ccFaces = cellConstrainedFaces();

    const scalarField& refL = refLenghts();

    forAll(cells, cellI)
    {
        label curCell = cells[cellI];

        const scalarRectangularMatrix& curInvMatrix = invMatrices[curCell];

        const labelList& interpCells = cCells[curCell];
        const labelList& interpFaces = cFaces[curCell];
        const labelList& constrFaces = ccFaces[curCell];

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source
        (
            interpCells.size() 
          + interpFaces.size()
          + constrFaces.size(),
            pTraits<Type>::zero
        );

        label pointID = 0;

        for (label i=0; i<interpCells.size(); i++)
        {
            source[pointID++] = vfI[interpCells[i]] - vfI[curCell];
        }

        for (label i=0; i<interpFaces.size(); i++)
        {
            label faceID = interpFaces[i];
            label patchID = 
                mesh_.boundaryMesh().whichPatch(faceID);

            label start = mesh_.boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID++] = 
                vf.boundaryField()[patchID][localFaceID] - vfI[curCell];
        }

        for (label i=0; i<constrFaces.size(); i++)
        {
            label faceID = constrFaces[i];
            label patchID = 
                mesh_.boundaryMesh().whichPatch(faceID);

            label start = mesh_.boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID++] = 
                vf.boundaryField()[patchID][localFaceID] - vfI[curCell];
        }

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        vector r = delta[cellI]/refL[curCell];

        drv[cellI] =
            direction[cellI].x()
           *(
                coeffs[0] + coeffs[2]*r.y() + 2*coeffs[3]*r.x()
            )	  
          + direction[cellI].y()
           *(
                coeffs[1] + coeffs[2]*r.x() + 2*coeffs[4]*r.y()
            );

        if (mesh_.nGeometricD() == 3)
        {
            drv[cellI] +=
                direction[cellI].z()
               *(
                    coeffs[5] 
                  + coeffs[6]*r.x() 
                  + coeffs[7]*r.y() 
                  + 2*coeffs[8]*r.z()
                )
              + direction[cellI].x()*coeffs[6]*r.z()
              + direction[cellI].y()*coeffs[7]*r.z();
        }

        drv[cellI] /= refL[curCell];
    }

    return tDrv;
}

template<class Type>
tmp<FieldField<Field, Type> > quadraticReconstruction::coeffs
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::coeffs("
            << "const GeometricField<Type, fvPatchField, volMesh>& ) : "
            << "calc coefficients of reconstraction"
            << endl;
    }

    const vectorField& C = mesh_.cellCentres();

    tmp<FieldField<Field, Type> > tCoeffs
    (
        new FieldField<Field, Type>(C.size())
    );
    FieldField<Field, Type>& Coeffs = tCoeffs();

//     tmp<IOPtrList<Field<Type> > > tCoeffs
//     (
//         new IOPtrList<Field<Type> >
//         (
//             IOobject
//             (
//                 "qCoeffs("+vf.name()+')',
//                 vf.instance(),
//                 mesh_,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             C.size()
//         )
//     );
//     IOPtrList<Field<Type> >& Coeffs = tCoeffs();

    const Field<Type>& vfI = vf.internalField();

    label nCoeffs = 5;
    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();
    const labelListList& cCells = cellCells();
    const labelListList& cFaces = cellFaces();
    const labelListList& ccFaces = cellConstrainedFaces();

    forAll(C, cellI)
    {
        label curCell = cellI;

        const scalarRectangularMatrix& curInvMatrix = invMatrices[curCell];

        const labelList& interpCells = cCells[curCell];
        const labelList& interpFaces = cFaces[curCell];
        const labelList& constrFaces = ccFaces[curCell];

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source
        (
            interpCells.size() 
          + interpFaces.size()
          + constrFaces.size(),
            pTraits<Type>::zero
        );

        label pointID = 0;

        for (label i=0; i<interpCells.size(); i++)
        {
            source[pointID++] = vfI[interpCells[i]] - vfI[curCell];
        }

        for (label i=0; i<interpFaces.size(); i++)
        {
            label faceID = interpFaces[i];
            label patchID = 
                mesh_.boundaryMesh().whichPatch(faceID);

            label start = mesh_.boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID++] = 
                vf.boundaryField()[patchID][localFaceID] - vfI[curCell];
        }

        for (label i=0; i<constrFaces.size(); i++)
        {
            label faceID = constrFaces[i];
            label patchID = 
                mesh_.boundaryMesh().whichPatch(faceID);

            label start = mesh_.boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID++] = 
                vf.boundaryField()[patchID][localFaceID] - vfI[curCell];
        }

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        Coeffs.set(cellI, new Field<Type>(coeffs));
    }

    return tCoeffs;
}


template<class Type>
void quadraticReconstruction::calcCoeffs
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    IOPtrList<Field<Type> >& qCoeffs
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::coeffs("
            << "const GeometricField<Type, fvPatchField, volMesh>& ) : "
            << "calc coefficients of reconstraction"
            << endl;
    }

    const vectorField& C = mesh_.cellCentres();

    qCoeffs.clear();
    qCoeffs.setSize(C.size());
//     qCoeffs.rename("qCoeffs("+vf.name()+')');

    const Field<Type>& vfI = vf.internalField();

    label nCoeffs = 5;
    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 4;
    }

    const PtrList<scalarRectangularMatrix>& invMatrices = invLsMatrices();
    const labelListList& cCells = cellCells();
    const labelListList& cFaces = cellFaces();
    const labelListList& ccFaces = cellConstrainedFaces();

    forAll(C, cellI)
    {
        label curCell = cellI;

        const scalarRectangularMatrix& curInvMatrix = invMatrices[curCell];

        const labelList& interpCells = cCells[curCell];
        const labelList& interpFaces = cFaces[curCell];
        const labelList& constrFaces = ccFaces[curCell];

        Field<Type> coeffs(nCoeffs, pTraits<Type>::zero);
        Field<Type> source
        (
            interpCells.size() 
          + interpFaces.size()
          + constrFaces.size(),
            pTraits<Type>::zero
        );

        label pointID = 0;

        for (label i=0; i<interpCells.size(); i++)
        {
            source[pointID++] = vfI[interpCells[i]] - vfI[curCell];
        }

        for (label i=0; i<interpFaces.size(); i++)
        {
            label faceID = interpFaces[i];
            label patchID = 
                mesh_.boundaryMesh().whichPatch(faceID);

            label start = mesh_.boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID++] = 
                vf.boundaryField()[patchID][localFaceID] - vfI[curCell];
        }

        for (label i=0; i<constrFaces.size(); i++)
        {
            label faceID = constrFaces[i];
            label patchID = 
                mesh_.boundaryMesh().whichPatch(faceID);

            label start = mesh_.boundaryMesh()[patchID].start();
            label localFaceID = faceID - start;

            source[pointID++] = 
                vf.boundaryField()[patchID][localFaceID] - vfI[curCell];
        }

        for (label i=0; i<nCoeffs; i++)
        {
            for (label j=0; j<source.size(); j++)
            {
                coeffs[i] += curInvMatrix[i][j]*source[j];
            }
        }

        qCoeffs.set(cellI, new Field<Type>(coeffs));
    }
}


template<class Type>
tmp
< 
    GeometricField<Type, fvsPatchField, surfaceMesh>
> 
quadraticReconstruction::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::interpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>& ) : "
            << "calc interpolation"
            << endl;
    }

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<Type>
            (
                "0",
                vf.dimensions(),
                pTraits<Type>::zero
            )
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();
   
    Field<Type>& sfI = sf.internalField();
 
    const Field<Type>& vfI = vf.internalField();

    FieldField<Field, Type> Coeffs = coeffs(vf);

    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();

    const unallocLabelList& own = mesh_.owner(); 
    const unallocLabelList& ngb = mesh_.neighbour();
    
    const scalarField& w = mesh_.weights().internalField();

    const scalarField& refL = refLenghts();

    forAll(sfI, faceI)
    {
        label ownCell = own[faceI];

        vector r = (Cf[faceI] - C[ownCell])/refL[ownCell];

        Type ownValue =
            vfI[ownCell]
          + Coeffs[ownCell][0]*r.x()
          + Coeffs[ownCell][1]*r.y()
          + Coeffs[ownCell][2]*r.x()*r.y()
          + Coeffs[ownCell][3]*sqr(r.x())
          + Coeffs[ownCell][4]*sqr(r.y());

        if (mesh_.nGeometricD() == 3)
        {
            ownValue += 
                Coeffs[ownCell][5]*r.z()
              + Coeffs[ownCell][6]*r.x()*r.z()
              + Coeffs[ownCell][7]*r.y()*r.z()
              + Coeffs[ownCell][8]*sqr(r.z());
        }

        label ngbCell = ngb[faceI];

        r = (Cf[faceI] - C[ngbCell])/refL[ngbCell];

        Type ngbValue =
            vfI[ngbCell]
          + Coeffs[ngbCell][0]*r.x()
          + Coeffs[ngbCell][1]*r.y()
          + Coeffs[ngbCell][2]*r.x()*r.y()
          + Coeffs[ngbCell][3]*sqr(r.x())
          + Coeffs[ngbCell][4]*sqr(r.y());

        if (mesh_.nGeometricD() == 3)
        {
            ngbValue += 
                Coeffs[ngbCell][5]*r.z()
              + Coeffs[ngbCell][6]*r.x()*r.z()
              + Coeffs[ngbCell][7]*r.y()*r.z()
              + Coeffs[ngbCell][8]*sqr(r.z());
        }
        
        sfI[faceI] = w[faceI]*(ownValue - ngbValue) + ngbValue;
    }

    forAll(sf.boundaryField(), patchI)
    {
        sf.boundaryField()[patchI] = vf.boundaryField()[patchI];
    }

    return tsf;
}


template<class Type>
void quadraticReconstruction::volToPointInterpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf,
    const IOPtrList<Field<Type> >& coeffs
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::volToPointInterpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    const labelListList& pointCells = mesh_.pointCells();
    const pointField& points = mesh_.points();
    const vectorField& C = mesh_.cellCentres();
    const scalarListList& pw = pointWeights();

    const scalarField& refL = refLenghts();

    forAll(pointCells, pointI)
    {
        const scalarList& pWeights = pw[pointI];
        const labelList& pCells = pointCells[pointI];

        pf[pointI] = pTraits<Type>::zero;

        forAll(pCells, cellI)
        {
            label curCell = pCells[cellI];

            vector curR = (points[pointI] - C[curCell])/refL[curCell];

            Type reconPointValue =
                vf[curCell]
              + coeffs[curCell][0]*curR.x()
              + coeffs[curCell][1]*curR.y()
              + coeffs[curCell][2]*curR.x()*curR.y()
              + coeffs[curCell][3]*sqr(curR.x())
              + coeffs[curCell][4]*sqr(curR.y());

            if (mesh_.nGeometricD() == 3)
            {
                reconPointValue += 
                    coeffs[curCell][5]*curR.z()
                  + coeffs[curCell][6]*curR.x()*curR.z()
                  + coeffs[curCell][7]*curR.y()*curR.z()
                  + coeffs[curCell][8]*sqr(curR.z());
            }

            pf[pointI] += pWeights[cellI]*reconPointValue;
        }
    }

    forAll(pf.boundaryField(), patchi)
    {
        if (pf.boundaryField()[patchi].coupled())
        {
            pf.boundaryField()[patchi].initAddField();
        }
    }

    forAll(pf.boundaryField(), patchi)
    {
        if (pf.boundaryField()[patchi].coupled())
        {
            pf.boundaryField()[patchi].addField
            (
                pf.internalField()
            );
        }
    }

    pf.correctBoundaryConditions();
}


template<class Type>
void quadraticReconstruction::volToPointInterpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::volToPointInterpolate("
            << "const GeometricField<Type, fvPatchField, volMesh>&, "
            << "GeometricField<Type, pointPatchField, pointMesh>&) : "
            << "interpolating field from cells to points"
            << endl;
    }

    const labelListList& pointCells = mesh_.pointCells();
    const pointField& points = mesh_.points();
    const vectorField& C = mesh_.cellCentres();
    const scalarListList& pw = pointWeights();

    FieldField<Field, Type> Coeffs = coeffs(vf);
    
    const scalarField& refL = refLenghts();

    forAll(pointCells, pointI)
    {
        const scalarList& pWeights = pw[pointI];
        const labelList& pCells = pointCells[pointI];

        pf[pointI] = pTraits<Type>::zero;

        forAll(pCells, cellI)
        {
            label curCell = pCells[cellI];

            vector curR = (points[pointI] - C[curCell])/refL[curCell];

            Type reconPointValue =
                vf[curCell]
              + Coeffs[curCell][0]*curR.x()
              + Coeffs[curCell][1]*curR.y()
              + Coeffs[curCell][2]*curR.x()*curR.y()
              + Coeffs[curCell][3]*sqr(curR.x())
              + Coeffs[curCell][4]*sqr(curR.y());

            if (mesh_.nGeometricD() == 3)
            {
                reconPointValue += 
                    Coeffs[curCell][5]*curR.z()
                  + Coeffs[curCell][6]*curR.x()*curR.z()
                  + Coeffs[curCell][7]*curR.y()*curR.z()
                  + Coeffs[curCell][8]*sqr(curR.z());
            }

            pf[pointI] += pWeights[cellI]*reconPointValue;
        }
    }

    forAll(pf.boundaryField(), patchi)
    {
        if (pf.boundaryField()[patchi].coupled())
        {
            pf.boundaryField()[patchi].initAddField();
        }
    }

    forAll(pf.boundaryField(), patchi)
    {
        if (pf.boundaryField()[patchi].coupled())
        {
            pf.boundaryField()[patchi].addField
            (
                pf.internalField()
            );
        }
    }

    pf.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
