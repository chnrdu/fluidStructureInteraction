/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "extendedPointToCell.H"
#include "polyMesh.H"
#include "pointSet.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(extendedPointToCell, 0);

addToRunTimeSelectionTable(topoSetSource, extendedPointToCell, word);

addToRunTimeSelectionTable(topoSetSource, extendedPointToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::extendedPointToCell::usage_
(
    extendedPointToCell::typeName,
    "\n    Usage: extendedPointToCell <pointSet> any\n\n"
    "    Select all cells with any point in the pointSet\n\n"
);

template<>
const char* Foam::NamedEnum<Foam::extendedPointToCell::pointAction, 1>::names[] =
{
    "any"
};


const Foam::NamedEnum<Foam::extendedPointToCell::pointAction, 1>
    Foam::extendedPointToCell::pointActionNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::extendedPointToCell::combine(topoSet& set, const bool add) const
{
    // Load the set
    pointSet loadedSet(mesh_, setName_);

    const labelListList& cellCells = mesh_.cellCells();

    // Handle any selection
    if (option_ == ANY)
    {
        for
        (
            pointSet::const_iterator iter = loadedSet.begin();
            iter != loadedSet.end();
            ++iter
        )
        {
            label pointI = iter.key();

            const labelList& pCells = mesh_.pointCells()[pointI];

            labelHashSet cellSet;
            
            forAll(pCells, cI)
            {
                if (!cellSet.found(pCells[cI]))
                {
                    cellSet.insert(pCells[cI]);
                }
            }

            labelList tmpCells(cellSet.toc());

            forAll(tmpCells, cellI)
            {
                const labelList& curCells = cellCells[tmpCells[cellI]];

                forAll(curCells, cI)
                {
                    if (!cellSet.found(curCells[cI]))
                    {
                        cellSet.insert(curCells[cI]);
                    }
                }
            }

            labelList cellLabels(cellSet.toc());

            forAll(cellLabels, cellI)
            {
                addOrDelete(set, cellLabels[cellI], add);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::extendedPointToCell::extendedPointToCell
(
    const polyMesh& mesh,
    const word& setName,
    const pointAction option
)
:
    topoSetSource(mesh),
    setName_(setName),
    option_(option)
{}


// Construct from dictionary
Foam::extendedPointToCell::extendedPointToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("set")),
    option_(pointActionNames_.read(dict.lookup("option")))
{}


// Construct from Istream
Foam::extendedPointToCell::extendedPointToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is)),
    option_(pointActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extendedPointToCell::~extendedPointToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::extendedPointToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells according to pointSet " << setName_
            << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells according to pointSet " << setName_
            << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
