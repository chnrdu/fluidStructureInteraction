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

#include "patchPointToCell.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(patchPointToCell, 0);

addToRunTimeSelectionTable(topoSetSource, patchPointToCell, word);

addToRunTimeSelectionTable(topoSetSource, patchPointToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::patchPointToCell::usage_
(
    patchPointToCell::typeName,
    "\n    Usage: patchPointToCell zone master|slave\n\n"
    "    Select master or slave side of the facePoint."
    " Note:accepts wildcards for zone.\n\n"
);

// template<>
// const char* Foam::NamedEnum<Foam::patchPointToCell::faceAction, 2>::names[] =
// {
//     "master",
//     "slave"
// };

// const Foam::NamedEnum<Foam::patchPointToCell::faceAction, 2>
//     Foam::patchPointToCell::faceActionNames_;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchPointToCell::combine(topoSet& set, const bool add) const
{
    bool hasMatched = false;

    forAll(mesh_.boundaryMesh(), i)
    {
//         const faceZone& zone = mesh_.faceZones()[i];
        const polyPatch& patch = mesh_.boundaryMesh()[i];

        if (patchName_.match(patch.name()))
        {
            labelHashSet cellSet;

            const labelList& patchPoints =
                mesh_.boundaryMesh()[i].meshPoints();

            const labelListList& pointCells = mesh_.pointCells();

            forAll(patchPoints, pointI)
            {
                const labelList& curCells = 
                    pointCells[patchPoints[pointI]];

                forAll(curCells, cellI)
                {
                    if (!cellSet.found(curCells[cellI]))
                    {
                        cellSet.insert(curCells[cellI]);
                    }
                }
            }

            labelList cellLabels(cellSet.toc());

            

//             const labelList& cellLabels =
//             (
//                 option_ == MASTER
//               ? zone.masterCells()
//               : zone.slaveCells()
//             );

//             Info<< "    Found matching zone " << zone.name()
//                 << " with " << cellLabels.size() << " cells on selected side."
//                 << endl;

//             hasMatched = true;

            forAll(cellLabels, i)
            {
                // Only do active cells
                if (cellLabels[i] < mesh_.nCells())
                {
                    addOrDelete(set, cellLabels[i], add);
                }
            }
        }
    }

    if (!hasMatched)
    {
        WarningIn("patchPointToCell::combine(topoSet&, const bool)")
            << "Cannot find any patch named " << patchName_ << endl;
//             << "Valid names are " << mesh_.faceZones().names() << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::patchPointToCell::patchPointToCell
(
    const polyMesh& mesh,
    const word& patchName
//     const faceAction option
)
:
    topoSetSource(mesh),
    patchName_(patchName)
//     option_(option)
{}


// Construct from dictionary
Foam::patchPointToCell::patchPointToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    patchName_(dict.lookup("name"))
//     option_(faceActionNames_.read(dict.lookup("option")))
{}


// Construct from Istream
Foam::patchPointToCell::patchPointToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    patchName_(checkIs(is))
//     option_(faceActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchPointToCell::~patchPointToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchPointToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
//         Info<< "    Adding all " << faceActionNames_[option_]
        Info<< "    Adding all cells of patch " 
            << patchName_ << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
//         Info<< "    Removing all " << faceActionNames_[option_]
        Info<< "    Removing all cells of patch " 
            << patchName_ << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
