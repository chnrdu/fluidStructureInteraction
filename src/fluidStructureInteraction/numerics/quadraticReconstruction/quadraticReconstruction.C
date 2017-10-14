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

#include "fvMesh.H"
#include "volFields.H"
#include "demandDrivenData.H"
#include "emptyPolyPatch.H"
#include "fvc.H"
#include "pointFields.H"
#include "symmetryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(quadraticReconstruction, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void quadraticReconstruction::makeCellCells() const
{
    if (debug)
    {
        InfoIn("void quadraticReconstruction::makeCellCells() const")
            << "create extended stencil"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellCellsPtr_)
    {
        FatalErrorIn
        (
            "quadraticReconstruction::makeCellCells() const"
        )   << "cell-cell addressing already exists"
            << abort(FatalError);
    }

    const vectorField& C = mesh_.cellCentres();
    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointCells = mesh_.pointCells();

    cellCellsPtr_ = new labelListList(C.size());
    labelListList& cellCells = *cellCellsPtr_;

    forAll (cellCells, cellI)
    {
        labelHashSet cellSet;

        const labelList& curCellPoints = cellPoints[cellI];

        forAll (curCellPoints, pointI)
        {
            label curPoint = curCellPoints[pointI];
            const labelList& curPointCells = pointCells[curPoint];

            forAll (curPointCells, cI)
            {
                if (!cellSet.found(curPointCells[cI]))
                {
                    cellSet.insert(curPointCells[cI]);
                }
            }
        }

        cellSet.erase(cellI);

        cellCells[cellI] = cellSet.toc();
    }
}


void quadraticReconstruction::makeCellFaces() const
{
    if (debug)
    {
        InfoIn("void quadraticReconstruction::makeCellFaces() const")
            << "create extended cell-cell stencil"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellFacesPtr_)
    {
        FatalErrorIn
        (
            "quadraticReconstruction::makeCellFaces() const"
        )   << "cell-faces addressing already exists"
            << abort(FatalError);
    }

    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointFaces = mesh_.pointFaces();

    const unallocLabelList& own = mesh_.owner(); 

    cellFacesPtr_ = new labelListList(cellPoints.size());
    labelListList& cellFaces = *cellFacesPtr_;

    forAll (cellFaces, cellI)
    {
        labelHashSet faceSet;

        const labelList& curCellPoints = cellPoints[cellI];

        forAll (curCellPoints, pointI)
        {
            label curPoint = curCellPoints[pointI];
            const labelList& curPointFaces = pointFaces[curPoint];

            forAll (curPointFaces, fI)
            {
                if (curPointFaces[fI] >= mesh_.nInternalFaces())
                {
                    label patchID = 
                        mesh_.boundaryMesh().whichPatch(curPointFaces[fI]);

                    if
                    (
                        mesh_.boundaryMesh()[patchID].type()
                     != emptyPolyPatch::typeName
                    )
                    {
                        if (!faceSet.found(curPointFaces[fI]))
                        {
                            if (own[curPointFaces[fI]] != cellI)
                            {
                                faceSet.insert(curPointFaces[fI]);
                            }
                        }
                    }
                }
            }
        }

        cellFaces[cellI] = faceSet.toc();
    }
}


void quadraticReconstruction::makeCellConstrainedFaces() const
{
    if (debug)
    {
        InfoIn("void quadraticReconstruction::makeCellConstrainedFaces() const")
            << "create extended cell-cell stencil"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellConstrainedFacesPtr_)
    {
        FatalErrorIn
        (
            "quadraticReconstruction::makeCellConstrainedFaces() const"
        )   << "cell-constrainedFaces addressing already exists"
            << abort(FatalError);
    }

    const labelListList& cellPoints = mesh_.cellPoints();
//     const labelListList& pointFaces = mesh_.pointFaces();

    cellConstrainedFacesPtr_ = new labelListList(cellPoints.size());
    labelListList& cellConstrainedFaces = *cellConstrainedFacesPtr_;

    forAll(cellConstrainedFaces, cellI)
    {
        labelHashSet faceSet;

        const labelList& curCellFaces = mesh_.cells()[cellI];

        forAll(curCellFaces, faceI)
        {
            if (curCellFaces[faceI] >= mesh_.nInternalFaces())
            {
                label patchID = 
                    mesh_.boundaryMesh().whichPatch(curCellFaces[faceI]);

                if
                (
                    mesh_.boundaryMesh()[patchID].type()
                 != emptyPolyPatch::typeName
                )
                {
                    if (!faceSet.found(curCellFaces[faceI]))
                    {
                        faceSet.insert(curCellFaces[faceI]);
                    }
                }
            }
        }

        cellConstrainedFaces[cellI] = faceSet.toc();
    }
}


void quadraticReconstruction::makeInvLsMatrices() const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::makeInvLsMatrices() : "
            << "making least squares linear interpolation matrices"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (invLsMatrices_.size() != 0)
    {
        FatalErrorIn("quadraticReconstruction::makeInvLsMatrices()")
            << "least square quadratic inerpolation matrices already exist"
            << abort(FatalError);
    }

    invLsMatrices_.setSize(mesh_.cells().size());

    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();

    const labelListList& cCells = cellCells();
    const labelListList& cFaces = cellFaces();
    const labelListList& ccFaces = cellConstrainedFaces();
        
    label nCoeffs = 5;
    if (mesh_.nGeometricD() == 3)
    {
        nCoeffs += 3;
    }

    scalarField conditionNumber(C.size(), 0.0);

    const scalarField& refL = refLenghts();

    forAll(invLsMatrices_, cellI)
    {
        const labelList& interpCells = cCells[cellI];
        const labelList& interpFaces = cFaces[cellI];
        const labelList& constrFaces = ccFaces[cellI];

        vectorField allPoints
        (
            interpCells.size()
          + interpFaces.size()
          + constrFaces.size(),
            vector::zero
        );

        if (allPoints.size() < nCoeffs)
        {
            FatalErrorIn
            (
                "quadraticReconstruction::makeInvLsMatrices()"
            )
                << "allPoints.size() < " << nCoeffs << " : "
                    << allPoints.size() << abort(FatalError);
        }

        // Scale points

        label pointID = 0;

        // Cells
        for (label i=0; i<interpCells.size(); i++)
        {
            allPoints[pointID++] = C[interpCells[i]];
        }

        // Faces
        for (label i=0; i<interpFaces.size(); i++)
        {
            allPoints[pointID++] = Cf[interpFaces[i]];
        }

        // Constrained faces
        for (label i=0; i<constrFaces.size(); i++)
        {
            allPoints[pointID++] = Cf[constrFaces[i]];
        }

        // Local origin
        vector origin = C[cellI];

        // Weights
        scalarField W(allPoints.size(), 1.0);

        for (label i=0; i<allPoints.size(); i++)
        {
            scalar curR =  mag(allPoints[i] - origin)/refL[cellI];
            W[i] = 1.0/sqr(curR);
        }

//         // Correct weights for constrained faces
//         for (label i=0; i<constrFaces.size(); i++)
//         {
//             W[constrFaces.size() - 1 - i] = 10/sqr(0.5);
//         }
        

        invLsMatrices_.set
        (
            cellI, 
            new scalarRectangularMatrix
            (
                nCoeffs,
                allPoints.size(),
                0.0
            ) 
        );
        scalarRectangularMatrix& curMatrix = invLsMatrices_[cellI];

        scalarRectangularMatrix M
        (
            allPoints.size(),
            nCoeffs,
            0.0
        );

        for(label i=0; i<allPoints.size(); i++)
        {
            scalar X = allPoints[i].x() - origin.x();
            scalar Y = allPoints[i].y() - origin.y();

            X /= refL[cellI];
            Y /= refL[cellI];

            label coeff = 0;
            M[i][coeff++] = X;
            M[i][coeff++] = Y;
            M[i][coeff++] = X*Y;
            M[i][coeff++] = sqr(X);
            M[i][coeff++] = sqr(Y);

            if (mesh_.nGeometricD() == 3)
            {
                scalar Z = allPoints[i].z() - origin.z();
                Z /= refL[cellI];

                M[i][coeff++] = Z;
                M[i][coeff++] = X*Z;
                M[i][coeff++] = Y*Z;
                M[i][coeff++] = sqr(Z);
            }
        }

        for (label i=0; i<M.n(); i++)
        {
            for (label j=0; j<M.m(); j++)
            {
                M[i][j] *= W[i];
            }
        }

        scalarSquareMatrix lsM(nCoeffs, 0.0);

        for (label i=0; i<lsM.n(); i++)
        {
            for (label j=0; j<lsM.m(); j++)
            {
                for (label k=0; k<M.n(); k++)
                {
                    lsM[i][j] += M[k][i]*M[k][j];
                }
            }
        }

        // Calculate matrix norm
        scalar maxRowSum = 0.0;
        for (label i=0; i<lsM.n(); i++)
        {
            scalar curRowSum = 0.0;

            for (label j=0; j<lsM.m(); j++)
            {
                curRowSum += lsM[i][j];
            }
            if(curRowSum > maxRowSum)
            {
                maxRowSum = curRowSum;
            }
        }
        conditionNumber[cellI] = maxRowSum;

        // Calculate inverse
        scalarSquareMatrix invLsM = lsM.LUinvert();

        for (label i=0; i<lsM.n(); i++)
        {
            for (label j=0; j<M.n(); j++)
            {
                for (label k=0; k<lsM.n(); k++)
                {
                    curMatrix[i][j] += invLsM[i][k]*M[j][k]*W[j];
                }
            }
        }

        // Calculate condition number
        maxRowSum = 0.0;
        for (label i = 0; i < lsM.n(); i++)
        {
            scalar curRowSum = 0.0;

            for (label j = 0; j < lsM.m(); j++)
            {
                curRowSum += invLsM[i][j];
            }
            if (curRowSum > maxRowSum)
            {
                maxRowSum = curRowSum;
            }
        }
        conditionNumber[cellI] *= maxRowSum;
    }

    Info << "Max matrix condition number: "
        << gMax(conditionNumber) << endl;
}



void quadraticReconstruction::makeRefLenghts() const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::makeRefLenghts() : "
            << "making cell based reference lenghts"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (refLenghtsPtr_)
    {
        FatalErrorIn("quadraticReconstruction::makeRefLenghts()")
            << "reference lenghts already exist"
            << abort(FatalError);
    }

    refLenghtsPtr_ = new scalarField(mesh_.nCells(), 0.0);
    scalarField& delta = *refLenghtsPtr_;

    if(mesh_.nGeometricD() == 3)
    {
        scalarField V(mesh_.V().field());

        delta = Foam::pow(V, 1.0/3.0);
    }
    else
    {
        scalar thickness = 0.0;
        const Vector<label>& directions = mesh_.geometricD();
        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh_.bounds().span()[dir];
                break;
            }
        }

        delta = sqrt(scalarField(mesh_.V().field())/thickness);
    }

    delta = average(delta);

//     delta = 1;

//     boundBox box(mesh_.points());
//     scalar refL = mag(box.max() - box.min())/2;
   
//     refL = 1;
 
//     Info << refL << endl;

//     delta = refL;
}



void quadraticReconstruction::makePointWeights() const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::makePointWeights() : "
            << "making vol-to-point interpolation weights"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (pointWeightsPtr_)
    {
        FatalErrorIn("quadraticReconstruction::makePointWeights()")
            << "weights already exist"
            << abort(FatalError);
    }

    const pointField& points = mesh_.points();
    const labelListList& pointCells = mesh_.pointCells();
    const vectorField& cellCentres = mesh_.cellCentres();

    pointWeightsPtr_ = new scalarListList(points.size());
    scalarListList& pointWeights = *pointWeightsPtr_;


    forAll(pointWeights, pointi)
    {
        pointWeights[pointi].setSize(pointCells[pointi].size());
    }

    pointScalarField sumWeights
    (
        IOobject
        (
            "volPointSumWeights",
            mesh_.polyMesh::instance(),
            mesh_
        ),
        pointMesh::New(mesh_),
        dimensionedScalar("zero", dimless, 0)
    );

    // Calculate inverse distances between cell centres and points
    // and store in weighting factor array
    forAll(points, pointi)
    {
        scalarList& pw = pointWeights[pointi];
        const labelList& pcp = pointCells[pointi];

        forAll(pcp, pointCelli)
        {
            pw[pointCelli] =
                1.0/magSqr(points[pointi] - cellCentres[pcp[pointCelli]]);

            sumWeights[pointi] += pw[pointCelli];
        }
    }

    forAll(sumWeights.boundaryField(), patchi)
    {
        if (sumWeights.boundaryField()[patchi].coupled())
        {
            sumWeights.boundaryField()[patchi].initAddField();
        }
    }

    forAll(sumWeights.boundaryField(), patchi)
    {
        if (sumWeights.boundaryField()[patchi].coupled())
        {
            sumWeights.boundaryField()[patchi].addField
            (
                sumWeights.internalField()
            );
        }
    }

    forAll(points, pointi)
    {
        scalarList& pw = pointWeights[pointi];

        forAll(pw, pointCelli)
        {
            pw[pointCelli] /= sumWeights[pointi];
        }
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

quadraticReconstruction::quadraticReconstruction
(
    const fvMesh& vm
)
:
    regIOobject
    (
        IOobject
        (
            "quadraticReconstruction",
            vm.time().constant(),
            vm,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(vm),
    cellCellsPtr_(NULL),
    cellFacesPtr_(NULL),
    cellConstrainedFacesPtr_(NULL),
    invLsMatrices_(),
    refLenghtsPtr_(NULL),
    pointWeightsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

quadraticReconstruction::~quadraticReconstruction()
{
    deleteDemandDrivenData(cellCellsPtr_);
    deleteDemandDrivenData(cellFacesPtr_);
    deleteDemandDrivenData(cellConstrainedFacesPtr_);
    invLsMatrices_.clear();
    deleteDemandDrivenData(refLenghtsPtr_);
    deleteDemandDrivenData(pointWeightsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool quadraticReconstruction::movePoints()
{
    invLsMatrices_.clear();

    return true;
}

const Foam::labelListList& quadraticReconstruction::cellCells() const
{
    if (!cellCellsPtr_)
    {
        makeCellCells();
    }

    return *cellCellsPtr_;
}

const Foam::labelListList& quadraticReconstruction::cellFaces() const
{
    if (!cellFacesPtr_)
    {
        makeCellFaces();
    }

    return *cellFacesPtr_;
}

const Foam::labelListList& quadraticReconstruction::
cellConstrainedFaces() const
{
    if (!cellConstrainedFacesPtr_)
    {
        makeCellConstrainedFaces();
    }

    return *cellConstrainedFacesPtr_;
}

const PtrList<scalarRectangularMatrix>&
quadraticReconstruction::invLsMatrices() const
{
    label size = invLsMatrices_.size();

    reduce(size, maxOp<label>());

    if (size == 0)
    {
        makeInvLsMatrices();
    }

    return invLsMatrices_;
}


const scalarField& quadraticReconstruction::refLenghts() const
{
    if (!refLenghtsPtr_)
    {
        makeRefLenghts();
    }

    return *refLenghtsPtr_;
}


const scalarListList& quadraticReconstruction::pointWeights() const
{
    if (!pointWeightsPtr_)
    {
        makePointWeights();
    }

    return *pointWeightsPtr_;
}


tmp<volVectorField> quadraticReconstruction::grad
(
    const volScalarField& vf
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::grad("
            << "const volScalarField& ) : "
            << "calc cell centre gradient"
            << endl;
    }

    tmp<volVectorField> tGrad
    (
        new volVectorField
        (
            IOobject
            (
                "grad("+vf.name()+')',
                vf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<vector>
            (
                "0",
                vf.dimensions()/dimLength,
                vector::zero
            ),
            zeroGradientFvPatchField<vector>::typeName
        )
    );
    volVectorField& Grad = tGrad();

    Field<vector>& GradI = Grad.internalField();

    FieldField<Field, scalar> Coeffs = coeffs(vf);
//     IOPtrList<Field<scalar> > Coeffs = coeffs(vf);

    const scalarField& refL = refLenghts();

    forAll(GradI, cellI)
    {
        GradI[cellI].x() = Coeffs[cellI][0];
        GradI[cellI].y() = Coeffs[cellI][1];

        if (mesh_.nGeometricD() == 3)
        {
            GradI[cellI].z() = Coeffs[cellI][5];
        }
    }

    GradI /= refL;

    Grad.correctBoundaryConditions();

    fv::gaussGrad<scalar>(mesh_).correctBoundaryConditions(vf, Grad);

    return tGrad;
}


tmp<volTensorField> quadraticReconstruction::grad
(
    const volVectorField& vf
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::grad("
            << "const volScalarField& ) : "
            << "calc cell centre gradient"
            << endl;
    }

    tmp<volTensorField> tGrad
    (
        new volTensorField
        (
            IOobject
            (
                "grad("+vf.name()+')',
                vf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<tensor>
            (
                "0",
                vf.dimensions()/dimLength,
                tensor::zero
            ),
            zeroGradientFvPatchField<tensor>::typeName
        )
    );
    volTensorField& Grad = tGrad();

    Field<tensor>& GradI = Grad.internalField();

    FieldField<Field, vector> Coeffs = coeffs(vf);
//     IOPtrList<Field<vector> > Coeffs = coeffs(vf);

    const scalarField& refL = refLenghts();

    forAll(GradI, cellI)
    {
        GradI[cellI].xx() = Coeffs[cellI][0].x();
        GradI[cellI].xy() = Coeffs[cellI][0].y();
        GradI[cellI].xz() = 0;

        GradI[cellI].yx() = Coeffs[cellI][1].x();
        GradI[cellI].yy() = Coeffs[cellI][1].y();
        GradI[cellI].yz() = 0;

        if (mesh_.nGeometricD() == 3)
        {
            GradI[cellI].xz() = Coeffs[cellI][0].z();
            GradI[cellI].yz() = Coeffs[cellI][1].z();
            
            GradI[cellI].zx() = Coeffs[cellI][5].x();
            GradI[cellI].zy() = Coeffs[cellI][5].y();
            GradI[cellI].zz() = Coeffs[cellI][5].z();
        }
    }

    GradI /= refL;

    Grad.correctBoundaryConditions();

    fv::gaussGrad<vector>(mesh_).correctBoundaryConditions(vf, Grad);

    return tGrad;
}


tmp<volTensorField> quadraticReconstruction::grad
(
    const volVectorField& vf,
    const IOPtrList<vectorField>& Coeffs
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::grad("
            << "const volScalarField& ) : "
            << "calc cell centre gradient"
            << endl;
    }

    tmp<volTensorField> tGrad
    (
        new volTensorField
        (
            IOobject
            (
                "grad("+vf.name()+')',
                vf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<tensor>
            (
                "0",
                vf.dimensions()/dimLength,
                tensor::zero
            ),
            zeroGradientFvPatchField<tensor>::typeName
        )
    );
    volTensorField& Grad = tGrad();

    Field<tensor>& GradI = Grad.internalField();

//     FieldField<Field, vector> Coeffs = coeffs(vf);
//     IOPtrList<Field<vector> > Coeffs = coeffs(vf);

    const scalarField& refL = refLenghts();

    forAll(GradI, cellI)
    {
        GradI[cellI].xx() = Coeffs[cellI][0].x();
        GradI[cellI].xy() = Coeffs[cellI][0].y();
        GradI[cellI].xz() = 0;

        GradI[cellI].yx() = Coeffs[cellI][1].x();
        GradI[cellI].yy() = Coeffs[cellI][1].y();
        GradI[cellI].yz() = 0;

        if (mesh_.nGeometricD() == 3)
        {
            GradI[cellI].xz() = Coeffs[cellI][0].z();
            GradI[cellI].yz() = Coeffs[cellI][1].z();
            
            GradI[cellI].zx() = Coeffs[cellI][5].x();
            GradI[cellI].zy() = Coeffs[cellI][5].y();
            GradI[cellI].zz() = Coeffs[cellI][5].z();
        }
    }

    GradI /= refL;

    Grad.correctBoundaryConditions();

    fv::gaussGrad<vector>(mesh_).correctBoundaryConditions(vf, Grad);

    return tGrad;
}


tmp<surfaceTensorField> quadraticReconstruction::sGrad
(
    const volVectorField& vf
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::sGrad("
            << "const volVectorField& ) : "
            << "calc surface gradient"
            << endl;
    }

    tmp<surfaceTensorField> tSGrad
    (
        new surfaceTensorField
        (
            IOobject
            (
                "grad"+vf.name()+'f',
                vf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<tensor>
            (
                "0",
                vf.dimensions()/dimLength,
                pTraits<tensor>::zero
            )
        )
    );
    surfaceTensorField& sGrad = tSGrad();
   
    tensorField& sGradI = sGrad.internalField();
 
    FieldField<Field, vector> Coeffs = coeffs(vf);
//     IOPtrList<Field<vector> > Coeffs = coeffs(vf);

    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();

    const surfaceVectorField n = mesh_.Sf()/mesh_.magSf();
//     const vectorField& nI = n.internalField();

    const unallocLabelList& own = mesh_.owner(); 
    const unallocLabelList& ngb = mesh_.neighbour();
    
    const scalarField& w = mesh_.weights().internalField();

    const scalarField& refL = refLenghts();

    forAll(sGradI, faceI)
    {
        // Face grad using owner cell reconstruction
        label ownCell = own[faceI];

        vector r = (Cf[faceI] - C[ownCell])/refL[ownCell];

        tensor ownGrad = tensor::zero;

        ownGrad.xx() = 
            Coeffs[ownCell][0].x() 
          + Coeffs[ownCell][2].x()*r.y()
          + 2*Coeffs[ownCell][3].x()*r.x();
        ownGrad.xy() =
            Coeffs[ownCell][0].y() 
          + Coeffs[ownCell][2].y()*r.y()
          + 2*Coeffs[ownCell][3].y()*r.x();
        ownGrad.xz() = 0;

        ownGrad.yx() = 
            Coeffs[ownCell][1].x() 
          + Coeffs[ownCell][2].x()*r.x()
          + 2*Coeffs[ownCell][4].x()*r.y();
        ownGrad.yy() =
            Coeffs[ownCell][1].y() 
          + Coeffs[ownCell][2].y()*r.x()
          + 2*Coeffs[ownCell][4].y()*r.y();
        ownGrad.yz() = 0;

        if (mesh_.nGeometricD() == 3)
        {
            ownGrad.xx() +=
                Coeffs[ownCell][6].x()*r.z();
            ownGrad.xy() +=
                Coeffs[ownCell][6].y()*r.z();
            ownGrad.xz() =
                Coeffs[ownCell][0].z() 
              + Coeffs[ownCell][2].z()*r.y()
              + 2*Coeffs[ownCell][3].z()*r.x()
              + Coeffs[ownCell][6].z()*r.z();

            ownGrad.yx() +=
                Coeffs[ownCell][7].x()*r.z();
            ownGrad.yy() +=
                Coeffs[ownCell][7].y()*r.z();
            ownGrad.yz() =
                Coeffs[ownCell][1].z() 
              + Coeffs[ownCell][2].z()*r.x()
              + 2*Coeffs[ownCell][4].z()*r.y()
              + Coeffs[ownCell][7].z()*r.z();

            ownGrad.zx() = 
                Coeffs[ownCell][5].x()
              + Coeffs[ownCell][6].x()*r.x()
              + Coeffs[ownCell][7].x()*r.y()
              + 2*Coeffs[ownCell][8].x()*r.z();
            ownGrad.zy() =
                Coeffs[ownCell][5].y()
              + Coeffs[ownCell][6].y()*r.x()
              + Coeffs[ownCell][7].y()*r.y()
              + 2*Coeffs[ownCell][8].y()*r.z();
            ownGrad.zz() =
                Coeffs[ownCell][5].z()
              + Coeffs[ownCell][6].z()*r.x()
              + Coeffs[ownCell][7].z()*r.y()
              + 2*Coeffs[ownCell][8].z()*r.z();
        }

        ownGrad /= refL[ownCell];


        // Face grad using neighbour cell reconstruction
        label ngbCell = ngb[faceI];

        r = (Cf[faceI] - C[ngbCell])/refL[ngbCell];

        tensor ngbGrad = tensor::zero;

        ngbGrad.xx() = 
            Coeffs[ngbCell][0].x() 
          + Coeffs[ngbCell][2].x()*r.y()
          + 2*Coeffs[ngbCell][3].x()*r.x();
        ngbGrad.xy() =
            Coeffs[ngbCell][0].y() 
          + Coeffs[ngbCell][2].y()*r.y()
          + 2*Coeffs[ngbCell][3].y()*r.x();
        ngbGrad.xz() = 0;

        ngbGrad.yx() = 
            Coeffs[ngbCell][1].x() 
          + Coeffs[ngbCell][2].x()*r.x()
          + 2*Coeffs[ngbCell][4].x()*r.y();
        ngbGrad.yy() =
            Coeffs[ngbCell][1].y() 
          + Coeffs[ngbCell][2].y()*r.x()
          + 2*Coeffs[ngbCell][4].y()*r.y();
        ngbGrad.yz() = 0;

        if (mesh_.nGeometricD() == 3)
        {
            ngbGrad.xx() +=
                Coeffs[ngbCell][6].x()*r.z();
            ngbGrad.xy() +=
                Coeffs[ngbCell][6].y()*r.z();
            ngbGrad.xz() =
                Coeffs[ngbCell][0].z() 
              + Coeffs[ngbCell][2].z()*r.y()
              + 2*Coeffs[ngbCell][3].z()*r.x()
              + Coeffs[ngbCell][6].z()*r.z();

            ngbGrad.yx() +=
                Coeffs[ngbCell][7].x()*r.z();
            ngbGrad.yy() +=
                Coeffs[ngbCell][7].y()*r.z();
            ngbGrad.yz() =
                Coeffs[ngbCell][1].z() 
              + Coeffs[ngbCell][2].z()*r.x()
              + 2*Coeffs[ngbCell][4].z()*r.y()
              + Coeffs[ngbCell][7].z()*r.z();

            ngbGrad.zx() = 
                Coeffs[ngbCell][5].x()
              + Coeffs[ngbCell][6].x()*r.x()
              + Coeffs[ngbCell][7].x()*r.y()
              + 2*Coeffs[ngbCell][8].x()*r.z();
            ngbGrad.zy() =
                Coeffs[ngbCell][5].y()
              + Coeffs[ngbCell][6].y()*r.x()
              + Coeffs[ngbCell][7].y()*r.y()
              + 2*Coeffs[ngbCell][8].y()*r.z();
            ngbGrad.zz() =
                Coeffs[ngbCell][5].z()
              + Coeffs[ngbCell][6].z()*r.x()
              + Coeffs[ngbCell][7].z()*r.y()
              + 2*Coeffs[ngbCell][8].z()*r.z();
        }

        ngbGrad /= refL[ngbCell];


        // Averaging
        sGradI[faceI] = w[faceI]*(ownGrad - ngbGrad) + ngbGrad;
    }

    forAll(sGrad.boundaryField(), patchI)
    {
        const unallocLabelList& faceCells = 
            mesh_.boundary()[patchI].faceCells();

        const vectorField& pCf = mesh_.boundary()[patchI].Cf();

        Field<tensor>& pSGrad = sGrad.boundaryField()[patchI];

        forAll(pSGrad, faceI)
        {
            label ownCell = faceCells[faceI];

            vector r = (pCf[faceI] - C[ownCell])/refL[ownCell];

            tensor ownGrad = tensor::zero;

            ownGrad.xx() = 
                Coeffs[ownCell][0].x() 
              + Coeffs[ownCell][2].x()*r.y()
              + 2*Coeffs[ownCell][3].x()*r.x();
            ownGrad.xy() =
                Coeffs[ownCell][0].y() 
              + Coeffs[ownCell][2].y()*r.y()
              + 2*Coeffs[ownCell][3].y()*r.x();
            ownGrad.xz() = 0;

            ownGrad.yx() = 
                Coeffs[ownCell][1].x() 
              + Coeffs[ownCell][2].x()*r.x()
              + 2*Coeffs[ownCell][4].x()*r.y();
            ownGrad.yy() =
                Coeffs[ownCell][1].y() 
              + Coeffs[ownCell][2].y()*r.x()
              + 2*Coeffs[ownCell][4].y()*r.y();
            ownGrad.yz() = 0;
            
            if (mesh_.nGeometricD() == 3)
            {
                ownGrad.xx() +=
                    Coeffs[ownCell][6].x()*r.z();
                ownGrad.xy() +=
                    Coeffs[ownCell][6].y()*r.z();
                ownGrad.xz() =
                    Coeffs[ownCell][0].z() 
                  + Coeffs[ownCell][2].z()*r.y()
                  + 2*Coeffs[ownCell][3].z()*r.x()
                  + Coeffs[ownCell][6].z()*r.z();

                ownGrad.yx() +=
                    Coeffs[ownCell][7].x()*r.z();
                ownGrad.yy() +=
                    Coeffs[ownCell][7].y()*r.z();
                ownGrad.yz() =
                    Coeffs[ownCell][1].z() 
                  + Coeffs[ownCell][2].z()*r.x()
                  + 2*Coeffs[ownCell][4].z()*r.y()
                  + Coeffs[ownCell][7].z()*r.z();
                
                ownGrad.zx() = 
                    Coeffs[ownCell][5].x()
                  + Coeffs[ownCell][6].x()*r.x()
                  + Coeffs[ownCell][7].x()*r.y()
                  + 2*Coeffs[ownCell][8].x()*r.z();
                ownGrad.zy() =
                    Coeffs[ownCell][5].y()
                  + Coeffs[ownCell][6].y()*r.x()
                  + Coeffs[ownCell][7].y()*r.y()
                  + 2*Coeffs[ownCell][8].y()*r.z();
                ownGrad.zz() =
                    Coeffs[ownCell][5].z()
                  + Coeffs[ownCell][6].z()*r.x()
                  + Coeffs[ownCell][7].z()*r.y()
                  + 2*Coeffs[ownCell][8].z()*r.z();
            }

            ownGrad /= refL[ownCell];

            pSGrad[faceI] = ownGrad;
        }
    }

    sGrad = ((I-n*n)&sGrad);
    sGrad += n*fvc::snGrad(vf);

    return tSGrad;
}


tmp<surfaceTensorField> quadraticReconstruction::sGrad
(
    const volVectorField& vf,
    const IOPtrList<vectorField>& Coeffs
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::sGrad("
            << "const volVectorField& ) : "
            << "calc surface gradient"
            << endl;
    }

    tmp<surfaceTensorField> tSGrad
    (
        new surfaceTensorField
        (
            IOobject
            (
                "grad"+vf.name()+'f',
                vf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<tensor>
            (
                "0",
                vf.dimensions()/dimLength,
                pTraits<tensor>::zero
            )
        )
    );
    surfaceTensorField& sGrad = tSGrad();
   
    tensorField& sGradI = sGrad.internalField();
 
//     FieldField<Field, vector> Coeffs = coeffs(vf);
//     IOPtrList<Field<vector> >& Coeffs = coeffs(vf);

    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();

    const surfaceVectorField n = mesh_.Sf()/mesh_.magSf();

    const unallocLabelList& own = mesh_.owner(); 
    const unallocLabelList& ngb = mesh_.neighbour();
    
    const scalarField& w = mesh_.weights().internalField();

    const scalarField& refL = refLenghts();

    forAll(sGradI, faceI)
    {
        // Face grad using owner cell reconstruction
        label ownCell = own[faceI];

        vector r = (Cf[faceI] - C[ownCell])/refL[ownCell];

        tensor ownGrad = tensor::zero;

        ownGrad.xx() = 
            Coeffs[ownCell][0].x() 
          + Coeffs[ownCell][2].x()*r.y()
          + 2*Coeffs[ownCell][3].x()*r.x();
        ownGrad.xy() =
            Coeffs[ownCell][0].y() 
          + Coeffs[ownCell][2].y()*r.y()
          + 2*Coeffs[ownCell][3].y()*r.x();
        ownGrad.xz() = 0;

        ownGrad.yx() = 
            Coeffs[ownCell][1].x() 
          + Coeffs[ownCell][2].x()*r.x()
          + 2*Coeffs[ownCell][4].x()*r.y();
        ownGrad.yy() =
            Coeffs[ownCell][1].y() 
          + Coeffs[ownCell][2].y()*r.x()
          + 2*Coeffs[ownCell][4].y()*r.y();
        ownGrad.yz() = 0;

        if (mesh_.nGeometricD() == 3)
        {
            ownGrad.xx() +=
                Coeffs[ownCell][6].x()*r.z();
            ownGrad.xy() +=
                Coeffs[ownCell][6].y()*r.z();
            ownGrad.xz() =
                Coeffs[ownCell][0].z() 
              + Coeffs[ownCell][2].z()*r.y()
              + 2*Coeffs[ownCell][3].z()*r.x()
              + Coeffs[ownCell][6].z()*r.z();

            ownGrad.yx() +=
                Coeffs[ownCell][7].x()*r.z();
            ownGrad.yy() +=
                Coeffs[ownCell][7].y()*r.z();
            ownGrad.yz() =
                Coeffs[ownCell][1].z() 
              + Coeffs[ownCell][2].z()*r.x()
              + 2*Coeffs[ownCell][4].z()*r.y()
              + Coeffs[ownCell][7].z()*r.z();

            ownGrad.zx() = 
                Coeffs[ownCell][5].x()
              + Coeffs[ownCell][6].x()*r.x()
              + Coeffs[ownCell][7].x()*r.y()
              + 2*Coeffs[ownCell][8].x()*r.z();
            ownGrad.zy() =
                Coeffs[ownCell][5].y()
              + Coeffs[ownCell][6].y()*r.x()
              + Coeffs[ownCell][7].y()*r.y()
              + 2*Coeffs[ownCell][8].y()*r.z();
            ownGrad.zz() =
                Coeffs[ownCell][5].z()
              + Coeffs[ownCell][6].z()*r.x()
              + Coeffs[ownCell][7].z()*r.y()
              + 2*Coeffs[ownCell][8].z()*r.z();
        }

        ownGrad /= refL[ownCell];


        // Face grad using neighbour cell reconstruction
        label ngbCell = ngb[faceI];

        r = (Cf[faceI] - C[ngbCell])/refL[ngbCell];

        tensor ngbGrad = tensor::zero;

        ngbGrad.xx() = 
            Coeffs[ngbCell][0].x() 
          + Coeffs[ngbCell][2].x()*r.y()
          + 2*Coeffs[ngbCell][3].x()*r.x();
        ngbGrad.xy() =
            Coeffs[ngbCell][0].y() 
          + Coeffs[ngbCell][2].y()*r.y()
          + 2*Coeffs[ngbCell][3].y()*r.x();
        ngbGrad.xz() = 0;

        ngbGrad.yx() = 
            Coeffs[ngbCell][1].x() 
          + Coeffs[ngbCell][2].x()*r.x()
          + 2*Coeffs[ngbCell][4].x()*r.y();
        ngbGrad.yy() =
            Coeffs[ngbCell][1].y() 
          + Coeffs[ngbCell][2].y()*r.x()
          + 2*Coeffs[ngbCell][4].y()*r.y();
        ngbGrad.yz() = 0;

        if (mesh_.nGeometricD() == 3)
        {
            ngbGrad.xx() +=
                Coeffs[ngbCell][6].x()*r.z();
            ngbGrad.xy() +=
                Coeffs[ngbCell][6].y()*r.z();
            ngbGrad.xz() =
                Coeffs[ngbCell][0].z() 
              + Coeffs[ngbCell][2].z()*r.y()
              + 2*Coeffs[ngbCell][3].z()*r.x()
              + Coeffs[ngbCell][6].z()*r.z();

            ngbGrad.yx() +=
                Coeffs[ngbCell][7].x()*r.z();
            ngbGrad.yy() +=
                Coeffs[ngbCell][7].y()*r.z();
            ngbGrad.yz() =
                Coeffs[ngbCell][1].z() 
              + Coeffs[ngbCell][2].z()*r.x()
              + 2*Coeffs[ngbCell][4].z()*r.y()
              + Coeffs[ngbCell][7].z()*r.z();

            ngbGrad.zx() = 
                Coeffs[ngbCell][5].x()
              + Coeffs[ngbCell][6].x()*r.x()
              + Coeffs[ngbCell][7].x()*r.y()
              + 2*Coeffs[ngbCell][8].x()*r.z();
            ngbGrad.zy() =
                Coeffs[ngbCell][5].y()
              + Coeffs[ngbCell][6].y()*r.x()
              + Coeffs[ngbCell][7].y()*r.y()
              + 2*Coeffs[ngbCell][8].y()*r.z();
            ngbGrad.zz() =
                Coeffs[ngbCell][5].z()
              + Coeffs[ngbCell][6].z()*r.x()
              + Coeffs[ngbCell][7].z()*r.y()
              + 2*Coeffs[ngbCell][8].z()*r.z();
        }

        ngbGrad /= refL[ngbCell];

        // Averaging
        sGradI[faceI] = w[faceI]*(ownGrad - ngbGrad) + ngbGrad;
    }

    forAll(sGrad.boundaryField(), patchI)
    {
        const unallocLabelList& faceCells = 
            mesh_.boundary()[patchI].faceCells();

        const vectorField& pCf = mesh_.boundary()[patchI].Cf();

        Field<tensor>& pSGrad = sGrad.boundaryField()[patchI];

        forAll(pSGrad, faceI)
        {
            label ownCell = faceCells[faceI];

            vector r = (pCf[faceI] - C[ownCell])/refL[ownCell];

            tensor ownGrad = tensor::zero;

            ownGrad.xx() = 
                Coeffs[ownCell][0].x() 
              + Coeffs[ownCell][2].x()*r.y()
              + 2*Coeffs[ownCell][3].x()*r.x();
            ownGrad.xy() =
                Coeffs[ownCell][0].y() 
              + Coeffs[ownCell][2].y()*r.y()
              + 2*Coeffs[ownCell][3].y()*r.x();
            ownGrad.xz() = 0;

            ownGrad.yx() = 
                Coeffs[ownCell][1].x() 
              + Coeffs[ownCell][2].x()*r.x()
              + 2*Coeffs[ownCell][4].x()*r.y();
            ownGrad.yy() =
                Coeffs[ownCell][1].y() 
              + Coeffs[ownCell][2].y()*r.x()
              + 2*Coeffs[ownCell][4].y()*r.y();
            ownGrad.yz() = 0;

            if (mesh_.nGeometricD() == 3)
            {
                ownGrad.xx() +=
                    Coeffs[ownCell][6].x()*r.z();
                ownGrad.xy() +=
                    Coeffs[ownCell][6].y()*r.z();
                ownGrad.xz() =
                    Coeffs[ownCell][0].z() 
                  + Coeffs[ownCell][2].z()*r.y()
                  + 2*Coeffs[ownCell][3].z()*r.x()
                  + Coeffs[ownCell][6].z()*r.z();

                ownGrad.yx() +=
                    Coeffs[ownCell][7].x()*r.z();
                ownGrad.yy() +=
                    Coeffs[ownCell][7].y()*r.z();
                ownGrad.yz() =
                    Coeffs[ownCell][1].z() 
                  + Coeffs[ownCell][2].z()*r.x()
                  + 2*Coeffs[ownCell][4].z()*r.y()
                  + Coeffs[ownCell][7].z()*r.z();

                ownGrad.zx() = 
                    Coeffs[ownCell][5].x()
                  + Coeffs[ownCell][6].x()*r.x()
                  + Coeffs[ownCell][7].x()*r.y()
                  + 2*Coeffs[ownCell][8].x()*r.z();
                ownGrad.zy() =
                    Coeffs[ownCell][5].y()
                  + Coeffs[ownCell][6].y()*r.x()
                  + Coeffs[ownCell][7].y()*r.y()
                  + 2*Coeffs[ownCell][8].y()*r.z();
                ownGrad.zz() =
                    Coeffs[ownCell][5].z()
                  + Coeffs[ownCell][6].z()*r.x()
                  + Coeffs[ownCell][7].z()*r.y()
                  + 2*Coeffs[ownCell][8].z()*r.z();
            }

            ownGrad /= refL[ownCell];

            pSGrad[faceI] = ownGrad;
        }
    }

    sGrad = ((I-n*n)&sGrad);
    sGrad += n*fvc::snGrad(vf);

    return tSGrad;
}


tmp<surfaceVectorField> quadraticReconstruction::sGrad
(
    const volScalarField& vf,
    const PtrList<scalarField>& Coeffs
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::sGrad("
            << "const volScalarField& ) : "
            << "calc surface gradient"
            << endl;
    }

    tmp<surfaceVectorField> tSGrad
    (
        new surfaceVectorField
        (
            IOobject
            (
                "grad"+vf.name()+'f',
                vf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<vector>
            (
                "0",
                vf.dimensions()/dimLength,
                pTraits<vector>::zero
            )
        )
    );
    surfaceVectorField& sGrad = tSGrad();
    vectorField& sGradI = sGrad.internalField();

//     FieldField<Field, vector> Coeffs = coeffs(vf);
//     IOPtrList<Field<vector> >& Coeffs = coeffs(vf);

    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();

    const surfaceVectorField n = mesh_.Sf()/mesh_.magSf();

    const unallocLabelList& own = mesh_.owner(); 
    const unallocLabelList& ngb = mesh_.neighbour();
    
    const scalarField& w = mesh_.weights().internalField();

    const scalarField& refL = refLenghts();

    forAll(sGradI, faceI)
    {
        // Face grad using owner cell reconstruction
        label ownCell = own[faceI];

        vector r = (Cf[faceI] - C[ownCell])/refL[ownCell];

        vector ownGrad = vector::zero;

        ownGrad.x() =
            Coeffs[ownCell][0]
          + Coeffs[ownCell][2]*r.y()
          + 2*Coeffs[ownCell][3]*r.x();

        ownGrad.y() =
            Coeffs[ownCell][1]
          + Coeffs[ownCell][2]*r.x()
          + 2*Coeffs[ownCell][4]*r.y();

        ownGrad.z() = 0;

        if (mesh_.nGeometricD() == 3)
        {
            ownGrad.x() +=
                Coeffs[ownCell][6]*r.z();

            ownGrad.y() +=
                Coeffs[ownCell][7]*r.z();

            ownGrad.z() =
                Coeffs[ownCell][5]
              + Coeffs[ownCell][6]*r.x()
              + Coeffs[ownCell][7]*r.y()
              + 2*Coeffs[ownCell][8]*r.z();
        }

        ownGrad /= refL[ownCell];


        // Face grad using neighbour cell reconstruction
        label ngbCell = ngb[faceI];

        r = (Cf[faceI] - C[ngbCell])/refL[ngbCell];

        vector ngbGrad = vector::zero;

        ngbGrad.x() =
            Coeffs[ngbCell][0]
          + Coeffs[ngbCell][2]*r.y()
          + 2*Coeffs[ngbCell][3]*r.x();

        ngbGrad.y() =
            Coeffs[ngbCell][1]
          + Coeffs[ngbCell][2]*r.x()
          + 2*Coeffs[ngbCell][4]*r.y();

        ngbGrad.z() = 0;

        if (mesh_.nGeometricD() == 3)
        {
            ngbGrad.x() +=
                Coeffs[ngbCell][6]*r.z();

            ngbGrad.y() +=
                Coeffs[ngbCell][7]*r.z();

            ngbGrad.z() =
                Coeffs[ngbCell][5]
              + Coeffs[ngbCell][6]*r.x()
              + Coeffs[ngbCell][7]*r.y()
              + 2*Coeffs[ngbCell][8]*r.z();
        }

        ngbGrad /= refL[ngbCell];

        // Averaging
        sGradI[faceI] = w[faceI]*(ownGrad - ngbGrad) + ngbGrad;
    }

    forAll(sGrad.boundaryField(), patchI)
    {
        const unallocLabelList& faceCells = 
            mesh_.boundary()[patchI].faceCells();

        const vectorField& pCf = mesh_.boundary()[patchI].Cf();

        Field<vector>& pSGrad = sGrad.boundaryField()[patchI];

        forAll(pSGrad, faceI)
        {
            label ownCell = faceCells[faceI];

            vector r = (pCf[faceI] - C[ownCell])/refL[ownCell];

            vector ownGrad = vector::zero;

            ownGrad.x() =
                Coeffs[ownCell][0]
              + Coeffs[ownCell][2]*r.y()
              + 2*Coeffs[ownCell][3]*r.x();

            ownGrad.y() =
                Coeffs[ownCell][1]
              + Coeffs[ownCell][2]*r.x()
              + 2*Coeffs[ownCell][4]*r.y();

            ownGrad.z() = 0;

            if (mesh_.nGeometricD() == 3)
            {
                ownGrad.x() +=
                    Coeffs[ownCell][6]*r.z();

                ownGrad.y() +=
                    Coeffs[ownCell][7]*r.z();

                ownGrad.z() =
                    Coeffs[ownCell][5]
                  + Coeffs[ownCell][6]*r.x()
                  + Coeffs[ownCell][7]*r.y()
                  + 2*Coeffs[ownCell][8]*r.z();
            }

            ownGrad /= refL[ownCell];

            pSGrad[faceI] = ownGrad;
        }
    }

    sGrad = ((I-n*n)&sGrad);
    sGrad += n*fvc::snGrad(vf);

    return tSGrad;
}


void quadraticReconstruction::calcPointGrad
(
    const volVectorField& vf,
    const IOPtrList<vectorField>& Coeffs,
    pointTensorField& pGradVf
) const
{
    if (debug)
    {
        Info<< "quadraticReconstruction::calcPointGrad("
            << "const volVectorField& ) : "
            << "calc point gradient"
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

        pGradVf[pointI] = tensor::zero;

        forAll(pCells, cellI)
        {
            label curCell = pCells[cellI];

            vector r = (points[pointI] - C[curCell])/refL[curCell];

            tensor curGrad = tensor::zero;
            
            curGrad.xx() = 
                Coeffs[curCell][0].x() 
              + Coeffs[curCell][2].x()*r.y()
              + 2*Coeffs[curCell][3].x()*r.x();
            curGrad.xy() =
                Coeffs[curCell][0].y() 
              + Coeffs[curCell][2].y()*r.y()
              + 2*Coeffs[curCell][3].y()*r.x();
            curGrad.xz() = 0;

            curGrad.yx() = 
                Coeffs[curCell][1].x() 
              + Coeffs[curCell][2].x()*r.x()
              + 2*Coeffs[curCell][4].x()*r.y();
            curGrad.yy() =
                Coeffs[curCell][1].y() 
              + Coeffs[curCell][2].y()*r.x()
              + 2*Coeffs[curCell][4].y()*r.y();
            curGrad.yz() = 0;

            if (mesh_.nGeometricD() == 3)
            {
                curGrad.xx() +=
                    Coeffs[curCell][6].x()*r.z();
                curGrad.xy() +=
                    Coeffs[curCell][6].y()*r.z();
                curGrad.xz() =
                    Coeffs[curCell][0].z() 
                  + Coeffs[curCell][2].z()*r.y()
                  + 2*Coeffs[curCell][3].z()*r.x()
                  + Coeffs[curCell][6].z()*r.z();

                curGrad.yx() +=
                    Coeffs[curCell][7].x()*r.z();
                curGrad.yy() +=
                    Coeffs[curCell][7].y()*r.z();
                curGrad.yz() =
                    Coeffs[curCell][1].z() 
                  + Coeffs[curCell][2].z()*r.x()
                  + 2*Coeffs[curCell][4].z()*r.y()
                  + Coeffs[curCell][7].z()*r.z();

                curGrad.zx() = 
                    Coeffs[curCell][5].x()
                  + Coeffs[curCell][6].x()*r.x()
                  + Coeffs[curCell][7].x()*r.y()
                  + 2*Coeffs[curCell][8].x()*r.z();
                curGrad.zy() =
                    Coeffs[curCell][5].y()
                  + Coeffs[curCell][6].y()*r.x()
                  + Coeffs[curCell][7].y()*r.y()
                  + 2*Coeffs[curCell][8].y()*r.z();
                curGrad.zz() =
                    Coeffs[curCell][5].z()
                  + Coeffs[curCell][6].z()*r.x()
                  + Coeffs[curCell][7].z()*r.y()
                  + 2*Coeffs[curCell][8].z()*r.z();
            }

            curGrad /= refL[curCell];

            pGradVf[pointI] += pWeights[cellI]*curGrad;
        }
    }

//     // Correct gradiet for ponts at symmetryPlane
//     forAll(mesh_.boundary(), patchI)
//     {
//         if ( isA<symmetryFvPatch>(mesh_.boundary()[patchI]) )
// 	{
//             pointPatchField<tensor>& ppGradVf = 
//                 pGradVf.boundaryField()[patchI];

// 	    const vectorField& n = ppGradVf.patch().pointNormals();

// 	    tensorField nnGradVft = ppGradVf.patchInternalField();
// 	    nnGradVft = (nnGradVft&(I-sqr(n)));
// 	    nnGradVft = -(sqr(n)&nnGradVft);
	    
//             ppGradVf.addToInternalField
//             (
//                 pGradVf.internalField(),
// 		nnGradVft
//             );
// 	}
//     }

//     forAll(pGradVf.boundaryField(), patchi)
//     {
//         if (pGradVf.boundaryField()[patchi].coupled())
//         {
//             pGradVf.boundaryField()[patchi].initAddField();
//         }
//     }

//     forAll(pGradVf.boundaryField(), patchi)
//     {
//         if (pGradVf.boundaryField()[patchi].coupled())
//         {
//             pGradVf.boundaryField()[patchi].addField
//             (
//                 pGradVf.internalField()
//             );
//         }
//     }

//     pGradVf.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
