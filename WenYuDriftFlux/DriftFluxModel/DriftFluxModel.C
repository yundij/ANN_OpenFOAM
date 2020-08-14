/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include "IOdictionary.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "phasePair.H"
#include "DriftFluxModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

   namespace dragModels
   {

      namespace driftFluxModels
      {

	  defineTypeNameAndDebug(DriftFluxModel, 0);
	  defineRunTimeSelectionTable(DriftFluxModel, dictionary);
      }

   }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::driftFluxModels::DriftFluxModel::DriftFluxModel
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    driftFlux_
    (
	IOobject
	(
            "driftFlux",
            pair.dispersed().U().time().timeName(),
            pair.dispersed().U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
	),
	pair.dispersed().U().mesh()
    ),
    pair_( pair ) 
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::driftFluxModels::DriftFluxModel::~DriftFluxModel()
{}


Foam::autoPtr< Foam::dragModels::driftFluxModels::DriftFluxModel >
Foam::dragModels::driftFluxModels::DriftFluxModel::New
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting driftFluxModel: " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown drift flux model type " << modelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid damping model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return
        autoPtr<Foam::dragModels::driftFluxModels::DriftFluxModel>
        (
            cstrIter()(dict,pair,registerObject)
        );
}

Foam::volVectorField& Foam::dragModels::driftFluxModels::DriftFluxModel::driftFlux()
{
   return driftFlux_;
}

void Foam::dragModels::driftFluxModels::DriftFluxModel::evolve()
{
   Info<< "Evolving drift flux ... " << endl;
}

bool Foam::dragModels::driftFluxModels::DriftFluxModel::writeData(Ostream& os) const
{
    return os.good();
}


