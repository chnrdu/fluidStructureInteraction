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

#include "equivalent.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "fvcGradf.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace calcTypes
    {
        defineTypeNameAndDebug(equivalent, 0);
        addToRunTimeSelectionTable(calcType, equivalent, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcTypes::equivalent::equivalent()
:
    calcType()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcTypes::equivalent::~equivalent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcTypes::equivalent::init()
{
    argList::validArgs.append("equivalent");
    argList::validArgs.append("fieldName");
}


void Foam::calcTypes::equivalent::preCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{}


void Foam::calcTypes::equivalent::calc
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

    // Check field exists
    if (fieldHeader.headerOk())
    {
        bool processed = false;

        calcEquivalent(fieldHeader, mesh, processed);

        if (!processed)
        {
            FatalError
                << "Unable to process " << fieldName << nl
                << "No equivalent for fields of type "
                << fieldHeader.headerClassName() << nl << nl
                << exit(FatalError);
        }
    }
    else
    {
        Info<< "    No " << fieldName << endl; 
    }
}

void Foam::calcTypes::equivalent::calcEquivalent
(
    const IOobject& fieldHeader,
    const fvMesh& mesh,
    bool& processed
)
{
    if (fieldHeader.headerClassName() == volSymmTensorField::typeName)
    {
        Info<< "    Reading " << fieldHeader.name() << endl;
        volSymmTensorField field(fieldHeader, mesh);

        Info<< "    Calculating equivalent " << fieldHeader.name() << endl;
        volScalarField fieldEq
        (
            IOobject
            (
                fieldHeader.name() + "Eq",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            sqrt((3.0/2.0)*magSqr(dev(field)))
        );
        fieldEq.write();

        Info<< "Max " << fieldHeader.name() << "Eq = " 
            << max(fieldEq).value() << endl;

        processed = true;
    }
}

// ************************************************************************* //

