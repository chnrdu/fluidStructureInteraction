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

#include "avgConstraint.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace calcTypes
    {
        defineTypeNameAndDebug(avgConstraint, 0);
        addToRunTimeSelectionTable(calcType, avgConstraint, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcTypes::avgConstraint::avgConstraint()
:
    calcType()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcTypes::avgConstraint::~avgConstraint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcTypes::avgConstraint::init()
{
    argList::validArgs.append("avgConstraint");
    argList::validArgs.append("fieldName");
}


void Foam::calcTypes::avgConstraint::preCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{}


void Foam::calcTypes::avgConstraint::calc
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

        calcAvgConstraint(fieldHeader, mesh, processed);

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
        Info<< "    No " << fieldName << endl;
    }
}

void Foam::calcTypes::avgConstraint::calcAvgConstraint
(
    const IOobject& header,
    const fvMesh& mesh,
    bool& processed
)
{
    if (header.headerClassName() == volSymmTensorField::typeName)
    {
        Info<< "    Reading " << header.name() << endl;
        volSymmTensorField sigma(header, mesh);

        volScalarField sigmaHyd = tr(sigma)/3;

        volScalarField sigmaEq = sqrt((3.0/2.0)*magSqr(dev(sigma)));
        
        scalar avgSigmaHyd = 
            gSum(sigmaHyd.internalField()*mesh.V().field())
           /gSum(mesh.V().field());

        scalar avgSigmaEq = 
            gSum(sigmaEq.internalField()*mesh.V().field())
           /gSum(mesh.V().field());

        Info<< "    Calculating average constrait for field "
            << header.name() << nl
            << "    " << avgSigmaHyd/(avgSigmaEq+SMALL)
            << endl;

        processed = true;
    }
}

// ************************************************************************* //

