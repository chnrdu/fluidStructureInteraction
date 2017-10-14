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

#include "strain.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "fvcGradf.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace calcTypes
    {
        defineTypeNameAndDebug(strain, 0);
        addToRunTimeSelectionTable(calcType, strain, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcTypes::strain::strain()
:
    calcType()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcTypes::strain::~strain()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcTypes::strain::init()
{
    argList::validArgs.append("strain");
    argList::validArgs.append("fieldName");
}


void Foam::calcTypes::strain::preCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{}


void Foam::calcTypes::strain::calc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    const word& fieldName = args.additionalArgs()[1];

    IOobject fieldHeader
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject pointFieldHeader
    (
        "point"+fieldName,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check field exists
    if (fieldHeader.headerOk() && pointFieldHeader.headerOk())
    {
        bool processed = false;

        calcStrain(fieldHeader, pointFieldHeader, mesh, processed);

        if (!processed)
        {
            FatalError
                << "Unable to process " << fieldName << nl
                << "No call to mag for fields of type "
                << fieldHeader.headerClassName() << nl << nl
                << exit(FatalError);
        }
    }
    else
    {
        Info<< "    No " << fieldName << " or point" << fieldName << endl; 
    }
}

void Foam::calcTypes::strain::calcStrain
(
    const IOobject& fieldHeader,
    const IOobject& pointFieldHeader,
    const fvMesh& mesh,
    bool& processed
)
{
    if (fieldHeader.headerClassName() == volVectorField::typeName)
    {
        Info<< "    Reading " << fieldHeader.name() << endl;
        volVectorField D(fieldHeader, mesh);

        Info<< "    Reading " << pointFieldHeader.name() << endl;
        pointMesh pMesh(mesh);
        pointVectorField pointD(pointFieldHeader, pMesh);

        volTensorField gradD = fvc::grad(D, pointD);

        volSymmTensorField epsilon = 
            symm(gradD) + 0.5*symm(gradD & gradD.T());

        Info<< "    Calculating equivalent stress" << endl;
        volScalarField epsilonEq
        (
            IOobject
            (
                "epsilonEq",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            sqrt((3.0/2.0)*magSqr(dev(epsilon)))
        );
        epsilonEq.write();

        Info<< "Max epsilonEq = " << max(epsilonEq).value()
            << endl;
        
        processed = true;
    }
}

// ************************************************************************* //

