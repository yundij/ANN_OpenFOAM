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
    DriftFluxModel.C

\*---------------------------------------------------------------------------*/

#ifndef DriftFluxModel_H
#define DriftFluxModel_H

#include "IOdictionary.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "phasePair.H"

namespace Foam
{

namespace dragModels
{

namespace driftFluxModels
{

class DriftFluxModel
:
    public regIOobject
{
    
    protected:
    
       volVectorField driftFlux_;
       const phasePair& pair_;

    public:

    //- Runtime type information
    TypeName("driftFluxModel");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            DriftFluxModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair,
                const bool registerObject
            ),
            (dict, pair, registerObject)
        );

    	
     DriftFluxModel
     ( 
	const dictionary& dict,
	const phasePair& pair,
	const bool registerObject
     );	

    //- Destructor
    virtual ~DriftFluxModel();	

    //- Selector
    static autoPtr<DriftFluxModel> New
    (
        const dictionary& dict,
	const phasePair& pair,
        const bool registerObject
    );
    
    volVectorField& driftFlux();
    virtual void evolve();
    
    //- Dummy write for regIOobject
    bool writeData(Ostream& os) const;
    
};

}

}

}

#endif



