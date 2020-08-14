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

#include "WenYuDriftFlux.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
   namespace dragModels
   {
       defineTypeNameAndDebug(WenYuDriftFlux, 0);
       addToRunTimeSelectionTable(dragModel, WenYuDriftFlux, dictionary);
   }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::WenYuDriftFlux::WenYuDriftFlux
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict),
    Hscale_("correction_factor", dimless, dict),

    driftFluxModel_
    (
       Foam::dragModels::driftFluxModels::DriftFluxModel::New
       (
           dict,
	   pair,
	   registerObject
       ).ptr()
    ),
    DragCorr_
    (
	IOobject
	(
            "dragCorrection",
            pair.dispersed().U().time().timeName(),
            pair.dispersed().U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	pair.dispersed().U().mesh(),
	dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 1.0)    
    ),
    resolvedDragCoefficient_
    (
	IOobject
	(
            "resolvedDragCoefficient",
            pair.dispersed().U().time().timeName(),
            pair.dispersed().U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
	),
	pair.dispersed().U().mesh(),
	dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 1.0)    
    ),
    debug_
    (
        Switch( dict.lookup("debug") )
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::WenYuDriftFlux::~WenYuDriftFlux()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::WenYuDriftFlux::CdRe() const
{
    volScalarField alpha1
    (
        max(pair_.dispersed(), pair_.dispersed().residualAlpha())
    );
    volScalarField alpha2
    (
        max(scalar(1) - pair_.dispersed(), pair_.continuous().residualAlpha())
    );

    volScalarField Res(alpha2*pair_.Re());
    //volScalarField Hscale(Hscale_);
    volScalarField CdsRes
    (
        neg(Res - 1000)*24.0*(1.0 + 0.15*pow(Res, 0.687))
      + pos(Res - 1000)*0.44*max(Res, residualRe_)
    );
    
    resolvedDragCoefficient_ = CdsRes * pow(alpha2, -3.65) * max(pair_.continuous(), pair_.continuous().residualAlpha());
   
    // ***** drag correction and drift flux evolution ***** //
    
    driftFluxModel_().evolve();
    
    volVectorField& driftFlux = driftFluxModel_().driftFlux();
    
    volVectorField slipVelocity
    (
    	pair_.continuous().U() - pair_.dispersed().U()
    );

    


    forAll( DragCorr_, iCell )
    {
    
      DragCorr_[iCell] = min(
			    max(
				1.0 + ( driftFlux[iCell] & slipVelocity[iCell] ) /alpha1[iCell] / max( 1e-6, magSqr( slipVelocity[iCell] ) ),
				//1.0 + ( driftFlux[iCell].component(2)  ) /alpha1[iCell] / max( 1e-6,  slipVelocity[iCell].component(2) ),
				0.001
				),
			        1.5
			    );
    }


    // ***** drag correction and drift flux evolution ***** //
    Info<< "DragCorr is " << DragCorr_[0]<< endl;
    
    if( debug_ )
    {
      return resolvedDragCoefficient_;
    }else
        DragCorr_.correctBoundaryConditions();
        return Hscale_* DragCorr_ * resolvedDragCoefficient_;
    
    DragCorr_.correctBoundaryConditions();
    return Hscale_ * DragCorr_ * resolvedDragCoefficient_;
        
}


// ************************************************************************* //
