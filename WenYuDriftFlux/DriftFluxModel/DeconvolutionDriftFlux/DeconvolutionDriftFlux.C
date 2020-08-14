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
#include "DeconvolutionDriftFlux.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "LESfilter.H"
#include "fvc.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

   namespace dragModels
   {

      namespace driftFluxModels
      {
	  defineTypeNameAndDebug(DeconvolutionDriftFlux, 0);
	  addToRunTimeSelectionTable(DriftFluxModel, DeconvolutionDriftFlux, dictionary);
      }

   }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::driftFluxModels::DeconvolutionDriftFlux::DeconvolutionDriftFlux
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject 
)
:
  Foam::dragModels::driftFluxModels::DriftFluxModel( dict, pair, registerObject ),
  filterPtr_(NULL),
  phiSolidResolved_
  (
      pair_.dispersed()	 
  ),
  phiGasResolved_
  (
      pair_.continuous()	 
  ),  
  uGasResolved_
  (
      pair_.continuous().U()
  ),
  numberOfFilters_
  (
      readScalar(dict.subDict( type() + "Props" ).lookup("NumberOfIterations"))
  )   
{
    Info << "Using ADM "<<numberOfFilters_<<"th order model for drag correction" << endl;

    //- LES filter
    filterPtr_ = LESfilter::New
    			    ( 
    				pair_.continuous().mesh(), 
				dict.subDict( type() + "Props" ).subDict("filterProps")
			    );

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::dragModels::driftFluxModels::DeconvolutionDriftFlux::~DeconvolutionDriftFlux()
{}

//performs one iteration of deconvolution 
void Foam::dragModels::driftFluxModels::DeconvolutionDriftFlux::filterFields
( 
    const volScalarField& phiSolid,
    const volScalarField& phiGas,  
    const volVectorField& ugas  
)
{

   tmp<volScalarField> phiGasResFiltered = phiGasResolved_;   
   tmp<volScalarField> phiSolidResFiltered = phiSolidResolved_;
   
   // G * ( alpha1 * ug ) at current iteration
   tmp<volVectorField> ugResFiltered = phiGasResolved_ * uGasResolved_;
   
   //apply filter to the fields
   applyFilter( phiSolidResFiltered );
   applyFilter( phiGasResFiltered );
   applyFilter( ugResFiltered );
   
   //update resolved solid volume fraction to next iteration
   phiSolidResolved_ = phiSolidResolved_ + ( phiSolid - phiSolidResFiltered() );
   phiGasResolved_   = phiGasResolved_   + ( phiGas   - phiGasResFiltered()   );
   
   //update gas velocity to next iteration
   uGasResolved_ = uGasResolved_ + ( ugas - ( ugResFiltered()/ phiGasResolved_ ) );
   
}

void Foam::dragModels::driftFluxModels::DeconvolutionDriftFlux::applyFilter( tmp<volScalarField>& field )
{
    volScalarField& fieldRef = field();  
    fieldRef = max( filterPtr_()( fieldRef ),1e-12 );
}

void Foam::dragModels::driftFluxModels::DeconvolutionDriftFlux::applyFilter( tmp<volVectorField>& field )
{
    volVectorField& fieldRef = field(); 
    fieldRef = filterPtr_()( fieldRef );
}

void Foam::dragModels::driftFluxModels::DeconvolutionDriftFlux::evolve()
{
   
   dragModels::driftFluxModels::DriftFluxModel::evolve();
   
   const volScalarField& phiSolid = pair_.dispersed();
   const volScalarField& phiGas = pair_.continuous();
   const volVectorField& ugas = pair_.continuous().U();
   
   //set resolved fields to course fields
   phiGasResolved_ = phiGas;
   phiSolidResolved_ = phiSolid;
   uGasResolved_ = ugas;
   
   //reconstruct the fine grid velocity values by iterative deconvolution
   for( int i = 0; i < numberOfFilters_; ++i )
      filterFields( phiSolid, phiGas, ugas );
   
   //compute drift flux
   tmp<volVectorField> driftTmp = phiSolidResolved_ * uGasResolved_;   
   applyFilter( driftTmp );
   
//   driftFlux_ = phiSolid * ( driftTmp()/max(phiSolid,1e-6) - ugas );
   driftFlux_ = driftTmp() - phiSolid * ugas;
  forAll( driftFlux_, iCell )
  {    
    
    driftFlux_[iCell].x()=0.0;
    driftFlux_[iCell].y()=0.0;

  }
   
   driftFlux_.correctBoundaryConditions();
   
}
