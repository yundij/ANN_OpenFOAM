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
#include "DriftFluxStatic3Par.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
//#include "fvOption.H"
#include "phasePair.H"
#include "extendedCentredCellToCellStencil.H"
#include "MeshObject.H"
#include "keras_model.H"
#include <vector>



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace dragModels
  {
    namespace driftFluxModels
    {
      defineTypeNameAndDebug(DriftFluxStatic3Par,0);
      addToRunTimeSelectionTable(DriftFluxModel, DriftFluxStatic3Par, dictionary);
    }
  }
}

using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::driftFluxModels::DriftFluxStatic3Par::DriftFluxStatic3Par
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject 
)
:
   Foam::dragModels::driftFluxModels::DriftFluxModel( dict, pair, registerObject ),
   mesh_( pair.dispersed().mesh() ),
   gradP
   (
	IOobject
	(
            "gradP",
            pair.dispersed().U().time().timeName(),
            pair.dispersed().U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	pair.dispersed().U().mesh(),
	dimensionedVector("gradP",dimensionSet(1,-2,-2,0,0),Foam::vector(0,0,0))
   ),  
   gradAlphaP
   (
	IOobject
	(
            "gradAlphaP",
            pair.dispersed().U().time().timeName(),
            pair.dispersed().U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	pair.dispersed().U().mesh(),
	dimensionedVector("gradAlphaP",dimensionSet(0,-1,0,0,0),Foam::vector(0,0,0))
   ),  
   DFnnModel_
   (
       "DFkerasParameters.nnet",
       false
   ),
   uTerminal
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("terminalVelocity"))
   ),
   alphaMax
   (
       readScalar(dict.subDict("driftFluxModelProps").lookup("MaximumSolidVolumeFraction"))
   ),
   /*
   filterSize_
   (
       dict.subDict("driftFluxModelProps").lookup("filterSize")   
   ),
   */
   g_
   (
       pair.g()
   ), 
   uSlip
   (
      pair_.continuous().U() - pair_.dispersed().U()
   ),
   uSlipVdrift
   ( 
      pair_.continuous().U() - pair_.dispersed().U() + driftFlux_
   ),
   solidVolumeFraction 
   (
       max( scalar(1) - pair_.continuous(), pair_.dispersed().residualAlpha() )
   ),
   rhoParticle_
   (
	pair_.dispersed().rho()[0]
    ),
   rhoFluid_
   (
	pair_.continuous().rho()[0]
    ),
   DFNFeatures
   (
       int( readScalar(dict.subDict("driftFluxModelProps").lookup("DF_NumberOfFeatures")) )  
   ),
   DFmeans_(NULL),
   DFvars_(NULL)
{
   
  //Info<<"Number of features: "<< DFNFeatures <<endl;
   
  DFmeans_ = new double[DFNFeatures];
  DFvars_ = new double[DFNFeatures];

  
     
    std::ifstream reader;

    reader.open( "DF_mean.csv" );
    
    if( !reader.is_open() )
    {
        Info<<"Failed to open file for DF_mean!"<<endl;
	
	for( int i = 0; i < DFNFeatures; ++i )
	    DFmeans_[i] = 0.0;
	
    }else
    {
	for( int i = 0; i < DFNFeatures; ++i )
            reader>>DFmeans_[i];
    }
    
    reader.close();

    reader.open( "DF_std.csv" );

    if( !reader.is_open() )
    {
        Info<<"Failed to open file for DF_std!"<<endl;
	
	for( int i = 0; i < DFNFeatures; ++i )
	    DFvars_[i] = 1.0;
	
    }else{
    
	for( int i = 0; i < DFNFeatures; ++i )
            reader>>DFvars_[i];
	    
    }
    
    reader.close();

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
Foam::dragModels::driftFluxModels::DriftFluxStatic3Par::~DriftFluxStatic3Par()
{
  if( DFmeans_ )
  {
     delete[] DFmeans_;
  }
  if( DFvars_ )
  {
     delete[] DFvars_;
  }
}



void Foam::dragModels::driftFluxModels::DriftFluxStatic3Par::evolve()
{

  //fluxes
  
  dragModels::driftFluxModels::DriftFluxModel::evolve();
  Info<< "Solving for drift flux transport ... " << endl;
  solidVolumeFraction = max( scalar(1) - pair_.continuous(), pair_.dispersed().residualAlpha() );
  uSlip = pair_.continuous().U() - pair_.dispersed().U();
  //  gradAlphaP = fvc::grad( solidVolumeFraction );
  gradP = fvc::grad( pair_.continuous().thermo().p());
  //driftFlux_min = -solidVolumeFraction*uSlip;
  //  tmp<volVectorField> gradGradPz = fvc::grad(gradP.component(2));
  //tmp<volScalarField> laplaceGradPz = fvc::laplacian(gradP.component(2));
  //tmp<volScalarField> laplaceAlphaP = fvc::laplacian(solidVolumeFraction);

  /******* Neural Network for DriftFlux********/


  //const extendedCentredCellToCellStencil& stencil = this->stencil();
  
  #include "DF_z.H"



  
  forAll( driftFlux_, iCell )
  {    
    //  if( solidVolumeFraction[iCell] <= 1e-3 )
    if(( solidVolumeFraction[iCell] <= 1e-2 ) || ( solidVolumeFraction[iCell] >= 0.55 ))
     {
         driftFlux_[iCell] = Foam::vector( 0.0, 0.0, 0.0 );
     }
  }
  driftFlux_.correctBoundaryConditions();

}

