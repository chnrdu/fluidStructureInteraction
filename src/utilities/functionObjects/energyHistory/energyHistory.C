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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "energyHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "boundBox.H"
#include "polyPatchID.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(energyHistory, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        energyHistory,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::energyHistory::writeData()
{
    Info << "Writing energy data" << endl;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const volScalarField& rho = 
        mesh.lookupObject<volScalarField>("rho");

    const volVectorField& U = 
        mesh.lookupObject<volVectorField>(velocityName_);

    const volVectorField& D = 
        mesh.lookupObject<volVectorField>(displacementName_);

    const volSymmTensorField& sigma = 
        mesh.lookupObject<volSymmTensorField>("sigma");


    if (mesh.foundObject<volVectorField>(displacementIncrementName_))
    {
        const volVectorField& DD = 
            mesh.lookupObject<volVectorField>(displacementIncrementName_);

        const volSymmTensorField& DSigma = 
            mesh.lookupObject<volSymmTensorField>("DSigma");

        const volSymmTensorField& DEpsilon = 
            mesh.lookupObject<volSymmTensorField>("DEpsilon");


        // External force work done
        scalarField DW(DD.boundaryField().size(), 0);
        forAll(DD.boundaryField(), patchI)
        {
            if (!D.boundaryField()[patchI].coupled())
            {
                const vectorField& patchDD =
                    DD.boundaryField()[patchI];

                const symmTensorField& patchSigma = 
                    sigma.boundaryField()[patchI];

                const symmTensorField& patchDSigma = 
                    DSigma.boundaryField()[patchI];

                const vectorField& patchS = 
                    mesh.Sf().boundaryField()[patchI];
        
                symmTensorField patchSigmaM = 
                    patchSigma - 0.5*patchDSigma;

                DW[patchI] = gSum((patchS & patchSigmaM) & patchDD);
                
                patchWork_[patchI] += DW[patchI];
            }
        }


        // Kinetic energy

        scalar kineticEnergyOld = kineticEnergy_;

        kineticEnergy_ =
            gSum
            (
                0.5*rho.internalField()
               *(U.internalField() & U.internalField())
               *mesh.V().field()
            );

        scalar DEk = kineticEnergy_ - kineticEnergyOld;

        // Elastic energy
        
        scalar DEs =
            gSum
            (
                (
                    (
                        sigma.internalField() 
                      - 0.5*DSigma.internalField()
                    )
                 && DEpsilon.internalField()
                )
               *mesh.V().field()
            );

        elasticEnergy_ += DEs;


        // Plastic energy

        scalar DEp = 0;

        if (mesh.foundObject<volSymmTensorField>("DEpsilonP"))
        {
            const volSymmTensorField& DEpsilonP = 
                mesh.lookupObject<volSymmTensorField>("DEpsilonP");

            DEp =
                gSum
                (
                    (
                        (
                            sigma.internalField() 
                          - 0.5*DSigma.internalField()
                        )
                     && DEpsilonP.internalField()
                    )
                   *mesh.V().field()
                );
        }

        plasticEnergy_ += DEp;

        Info << "DW: " << sum(DW) << ", DEk: " << DEk << ", DEs: " << DEs 
            << ", DEp: " << DEp << ", sum(DEk,DEs,DEp): " 
            << DEk + DEs + DEp << endl;
    }
    else
    {
        const volSymmTensorField& epsilon = 
            mesh.lookupObject<volSymmTensorField>("epsilon");

        // External force work done
        scalarField DW(D.boundaryField().size(), 0);

        forAll(D.boundaryField(), patchI)
        {
            if (!D.boundaryField()[patchI].coupled())
            {
                const vectorField patchDD =
                    D.boundaryField()[patchI]
                  - D.oldTime().boundaryField()[patchI];

                const symmTensorField& patchSigma = 
                    sigma.boundaryField()[patchI];

                const symmTensorField& patchSigmaOld = 
                    sigma.oldTime().boundaryField()[patchI];

                const vectorField& patchS = 
                    mesh.Sf().boundaryField()[patchI];

                symmTensorField patchSigmaM = 
                    0.5*(patchSigma + patchSigmaOld);

                DW[patchI] = gSum((patchS&patchSigmaM) & patchDD);

                patchWork_[patchI] += DW[patchI];
            }
        }


        // Kinetic energy

        scalar kineticEnergyOld = kineticEnergy_;

        kineticEnergy_ =
            gSum
            (
                0.5*rho.internalField()
               *(U.internalField() & U.internalField())
               *mesh.V()
            );

        scalar DEk = kineticEnergy_ - kineticEnergyOld;


        // Elastic energy
        
        scalar DEs =
            gSum
            (
                0.5
               *(
                    (
                        sigma.internalField() 
                      + sigma.oldTime().internalField()
                    )
                 && (
                        epsilon.internalField()
                      - epsilon.oldTime().internalField()
                    )
                )
               *mesh.V()
            );

        elasticEnergy_ += DEs;

        Info << "DW: " << sum(DW) << ", DEk: " << DEk << ", DEs: " << DEs  
            << ", sum(DEk,DEs): " << DEk + DEs << endl;
    }

    Info << "W: " << sum(patchWork_) 
        << ", Ek: " << kineticEnergy_
        << ", Es: " << elasticEnergy_
        << ", DEp: " << plasticEnergy_
        << ", sum(Ek,Es,Ep): " 
        << kineticEnergy_ + elasticEnergy_ + plasticEnergy_ 
        << endl << endl;

    if (Pstream::master())
    {
        historyFilePtr_()
            << mesh.time().value() << tab
                << tab << kineticEnergy_
                << tab << elasticEnergy_
                << tab << plasticEnergy_
                << tab << sum(patchWork_);

        forAll(mesh.boundary(), patchI)
        {
            if (!mesh.boundary()[patchI].coupled())
            {
                historyFilePtr_() 
                    << tab << patchWork_[patchI];
            }
        }

        historyFilePtr_() << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energyHistory::energyHistory
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    historyFilePtr_(NULL),
    displacementName_("D"),
    displacementIncrementName_("DD"),
    velocityName_("U"),
    kineticEnergy_(0),
    elasticEnergy_(0),
    plasticEnergy_(0),
    patchWork_()
{
    Info << "Creating functio object " << name_ << "\n" << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    if (dict.found("displacementName"))
    {
        dict.lookup("displacementName") >> displacementName_;
    }

    if (dict.found("displacementIncrementName"))
    {
        dict.lookup("displacementIncrementName") >> displacementIncrementName_;
    }

    if (dict.found("velocityName"))
    {
        dict.lookup("velocityName") >> velocityName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    patchWork_.setSize(mesh.boundaryMesh().size(), 0);

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                mesh.time().timeName(mesh.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"history"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"history"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            // Open new file at start up
            historyFilePtr_.reset(new OFstream(historyDir/"energy.dat"));

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time"
                        << tab << "kineticEnergy" << tab
                        << tab << "elasticEnergy" << tab
                        << tab << "plasticEnergy" << tab
                        << tab << "externalWork";

                forAll(mesh.boundary(), patchI)
                {
                    if (!mesh.boundary()[patchI].coupled())
                    {
                        historyFilePtr_()
                            << tab 
                            << word("work_") 
                          + mesh.boundaryMesh()[patchI].name();
                    }
                }

                historyFilePtr_() << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::energyHistory::start()
{
    return writeData();
}


bool Foam::energyHistory::execute()
{
    return writeData();
}


bool Foam::energyHistory::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
