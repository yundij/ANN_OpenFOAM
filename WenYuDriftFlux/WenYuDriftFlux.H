/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::dragModels::WenYu

Description
    Wen and Yu drag model

    Reference:
    \verbatim
        "Eulerian Two-Phase Flow Theory Applied to Fluidization"
        Enwald, H., Peirano, E., Almstedt, A-E.,
        Int. J. Multiphase Flow, Vol. 22, Suppl, 1996, pp. 21-66
        Eq. 86-87, p. 40
    \endverbatim

SourceFiles
    WenYu.C

\*---------------------------------------------------------------------------*/

#ifndef WenYuDriftFlux_H
#define WenYuDriftFlux_H

#include "DriftFluxModel.H"
#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class WenYu Declaration
\*---------------------------------------------------------------------------*/

class WenYuDriftFlux
:
    public dragModel
{
    // Private data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;
        const dimensionedScalar Hscale_;

	
	mutable autoPtr< driftFluxModels::DriftFluxModel > driftFluxModel_;
	
	mutable volScalarField DragCorr_;
	mutable volScalarField resolvedDragCoefficient_;
	
	const bool debug_;

public:

    //- Runtime type information
    TypeName("WenYuDriftFlux");


    // Constructors

        //- Construct from a dictionary and a phase pair
        WenYuDriftFlux
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~WenYuDriftFlux();


    // Member Functions
	
	const volScalarField& DragCorr(){ return DragCorr_; } 
	const volScalarField& resolvedDragCoefficient(){ return resolvedDragCoefficient_; }
		
        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
