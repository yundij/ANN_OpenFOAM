/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "volFields.H"
#include "DriftFluxModel.H"
#include "ConstantDriftFlux.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

   namespace dragModels
   {

      namespace driftFluxModels
      {
	  defineTypeNameAndDebug(ConstantDriftFlux, 0);
	  addToRunTimeSelectionTable(DriftFluxModel, ConstantDriftFlux, dictionary);
      }

   }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::driftFluxModels::ConstantDriftFlux::ConstantDriftFlux
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject 
)
:
  Foam::dragModels::driftFluxModels::DriftFluxModel( dict, pair, registerObject ),
   uSlip
   (
      pair_.continuous().U() - pair_.dispersed().U()
    ),
  solidVolumeFraction 
   (
       max( scalar(1) - pair_.continuous(), pair_.dispersed().residualAlpha() )
   )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::dragModels::driftFluxModels::ConstantDriftFlux::~ConstantDriftFlux()
{}



void Foam::dragModels::driftFluxModels::ConstantDriftFlux::evolve()
{
   dragModels::driftFluxModels::DriftFluxModel::evolve();
   uSlip = pair_.continuous().U() - pair_.dispersed().U();
   solidVolumeFraction = max(max( scalar(1) - pair_.continuous(), pair_.dispersed().residualAlpha() ),1e-3);
   Info<< "Setting drift flux to a Multiple of Uslip ... " << endl;
   forAll(driftFlux_, iCell )
     {
       driftFlux_[iCell].z() = -0.7* uSlip[iCell].z();
       driftFlux_[iCell].x() = 0;
       driftFlux_[iCell].y() = 0;

     }
   driftFlux_.correctBoundaryConditions();
}
