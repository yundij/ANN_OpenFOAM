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

Class
    Foam::dragModels::driftFluxModels::DriftFluxModel::DeconvolutionDriftFlux

Description

SourceFiles
    DeconvolutionDriftFlux.C

\*---------------------------------------------------------------------------*/

#ifndef DeconvolutionDriftFlux_H
#define DeconvolutionDriftFlux_H

#include "volFields.H"
#include "DriftFluxModel.H"
#include "LESfilter.H"
#include "fvc.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace dragModels
{

namespace driftFluxModels
{

/*---------------------------------------------------------------------------*\
                         Class DeconvolutionDriftFlux Declaration
\*---------------------------------------------------------------------------*/

class DeconvolutionDriftFlux
:
   public DriftFluxModel
{

protected:
   
   int numberOfFilters_;
   
   mutable volScalarField phiGasResolved_;
   mutable volScalarField phiSolidResolved_;
   mutable volVectorField uGasResolved_;
   
   void applyFilter( tmp<volScalarField>& field );
   void applyFilter( tmp<volVectorField>& field );
   
   void filterFields( const volScalarField&, const volScalarField&, const volVectorField& );

   //- Pointer for filter
   autoPtr<LESfilter> filterPtr_;   

public:

  //- Runtime type information
  TypeName("DeconvolutionDriftFlux");

  //Constructor
  DeconvolutionDriftFlux
  ( 
     const dictionary& dict,
     const phasePair& pair,
     const bool registerObject
  );

  //- Destructor
  virtual ~DeconvolutionDriftFlux();
  
  virtual void evolve();
  
};

}
}
}

#endif
