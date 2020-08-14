/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "WenYuANNDrag.H"
#include "addToRunTimeSelectionTable.H"
#include "averagingModel.H"
#include <fstream>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(WenYuDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    WenYuDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
WenYuDrag::WenYuDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    phi_(propsDict_.lookupOrDefault<scalar>("phi",1.)),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    switchingVoidfraction_(0.1),
    //---------Added for NN------------//
    Resys_(propsDict_.lookupOrDefault<scalar>("Resys",5.0)),
    fs_(propsDict_.lookupOrDefault<scalar>("filter_size",0.001)),
    pFieldName_(propsDict_.lookup("pFieldName")),
    p_(sm.mesh().lookupObject<volScalarField> (pFieldName_)),
    useReNNCorr_(false),
    DFnnModel_
    (
     "DFkerasParameters.nnet",
     false
     ),
    uTerminal(propsDict_.lookupOrDefault<scalar>("terminal_velocity",0.2)),
    alphaMax(propsDict_.lookupOrDefault<scalar>("MaximumSolidVolumeFraction",0.64)),
    rho_s(propsDict_.lookupOrDefault<scalar>("solidDensity",1500)),
    DFNFeatures(propsDict_.lookupOrDefault<scalar>("DF_numberOfFeatures", 5)),
    DFmeans_(NULL),
    DFvars_(NULL),
    driftFlux_
    (
        IOobject
        (
            "driftFlux",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sm.mesh()
    )
{
    //Append the field names to be probed
    // suppress particle probe
    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, typeName+".logDat");
        particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must  be the force
        particleCloud_.probeM().vectorFields_.append("Urel");
        particleCloud_.probeM().scalarFields_.append("Rep");
        particleCloud_.probeM().scalarFields_.append("betaP");
        particleCloud_.probeM().scalarFields_.append("voidfraction");
        particleCloud_.probeM().writeHeader();
    }

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(2,true); // activate implDEM switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(true);

    if (propsDict_.found("switchingVoidfraction"))
        switchingVoidfraction_ = readScalar(propsDict_.lookup("switchingVoidfraction"));
    //-------------------Added for NN--------------//
    if (propsDict_.found("useReNNCorr")) useReNNCorr_=true;

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

WenYuDrag::~WenYuDrag()
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void WenYuDrag::setForce() const
{

    const volScalarField& rhoField = forceSubM(0).rhoField();
    #if defined(version24Dev)
       // there seems to have been a change in the return value of 
       // particleCloud_.turbulence().nu() used by forceSubM(0).nuField();
       const volScalarField& nufField = particleCloud_.turbulence().nu();
    #else
        const volScalarField& nufField = forceSubM(0).nuField();
    #endif


    //update force submodels to prepare for loop
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).preParticleLoop(forceSubM(iFSub).verbose());

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI=0;

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar dParcel(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Vs(0);
    scalar localPhiP(0);

    scalar CdMagUrLag(0);       //Cd of the very particle
    scalar betaP(0);             //momentum exchange of the very particle

    vector dragExplicit(0,0,0);
    scalar dragCoefficient(0);

    //-------------Added for NN-------//
    volVectorField gradP_ = fvc::grad(p_);
    vector gradP;
    scalar dragCorrection(1.);
    scalar driftFlux(0.0);
    
    #include "resetVoidfractionInterpolator.H"
    #include "resetUInterpolator.H"
    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {
        //if(mask[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector(0,0,0);
            dragExplicit = vector(0,0,0);
            betaP = 0;
            Vs = 0;
            Ufluid =vector(0,0,0);
            voidfraction=0;
            dragCoefficient = 0;

            if (cellI > -1) // particle Found
            {

                if( forceSubM(0).interpolation() )
                {
	            position     = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_().interpolate(position,cellI);
                    Ufluid       = UInterpolator_().interpolate(position,cellI);
		    gradP = gradPInterpolator_().interpolate(position,cellI); //Added for NN
		    //Ensure interpolated void fraction to be meaningful
                    // Info << " --> voidfraction: " << voidfraction << endl;
                    if(voidfraction>1.00) voidfraction = 1.0;
                    if(voidfraction<0.10) voidfraction = 0.10;
                }
                else
                {
		    voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
		    gradP = gradP_[cellI]; // Added for NN
                }

                Us = particleCloud_.velocity(index);
                ds = 2*particleCloud_.radius(index);
                dParcel = ds;
                forceSubM(0).scaleDia(ds,index); //caution: this fct will scale ds!
                rho = rhoField[cellI];
                nuf = nufField[cellI];

                //Update any scalar or vector quantity
                for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                      forceSubM(iFSub).update(  index, 
                                                cellI,
                                                ds,
                                                Ufluid, 
                                                Us, 
                                                nuf,
                                                rho,
                                                forceSubM(0).verbose()
                                             );

                Vs = ds*ds*ds*M_PI/6;
                Ur = Ufluid-Us;
                magUr = mag(Ur);
                Rep=0.0;
                localPhiP = 1.0f-voidfraction+SMALL;

                // calc particle's drag coefficient (i.e., Force per unit slip velocity and per m³ PARTICLE)
                
                
                Rep=ds*voidfraction*magUr/(nuf+SMALL);
                CdMagUrLag = (24.0*nuf/(ds*voidfraction)) //1/magUr missing here, but compensated in expression for betaP!
                                 *(scalar(1.0)+0.15*Foam::pow(Rep, 0.687));

                betaP = 0.75*(                                  //this is betaP = beta / localPhiP!
                               rho*voidfraction*CdMagUrLag
                                /
                               (ds*Foam::pow(voidfraction,2.65))
                             );
                
		// ------ Added for NN----------//
		// -------calc drag correction with eulerian variables ----- //
		if(useReNNCorr_)
		{
		  if ( ((1 - voidfraction) > 1e-2) && ((1-voidfraction) < 0.55) ) {
		    #include "DF_z.H"
		  }
		  driftFlux_[cellI] = driftFlux;
		  dragCorrection = min(
				       max(
					   1.0 + driftFlux * Ur.z() / (1-voidfraction) / max(1e-6, magUr*magUr),
					   0.001
					   ),
				       1.5
				       );
		}else
		{
		  dragCorrection = 1.0;
		}



                // calc particle's drag
                dragCoefficient = dragCorrection * Vs*betaP;     // added for NN
                if (modelType_=="B")
                    dragCoefficient /= voidfraction;

                forceSubM(0).scaleCoeff(dragCoefficient,dParcel,index);
                drag = dragCoefficient * Ur;

                // explicitCorr
                for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                    forceSubM(iFSub).explicitCorr( drag, 
                                                   dragExplicit,
                                                   dragCoefficient,
                                                   Ufluid, U_[cellI], Us, UsField_[cellI],
                                                   forceSubM(iFSub).verbose()
                                                 );
                if(forceSubM(0).verbose() && index >=0 && index <2)
                {
                    Pout << "cellI = " << cellI << endl;
                    Pout << "index = " << index << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "dprim = " << ds << endl;
                    Pout << "phi = " << phi_ << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "betaP = " << betaP << endl;
                    Pout << "drag = " << drag << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    // Note: for other than ext one could use vValues.append(x)
                    // instead of setSize
                    vValues.setSize(vValues.size()+1, drag);           //first entry must the be the force
                    vValues.setSize(vValues.size()+1, Ur);
                    sValues.setSize(sValues.size()+1, Rep); 
                    sValues.setSize(sValues.size()+1, betaP);
                    sValues.setSize(sValues.size()+1, voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,drag,dragExplicit,Ufluid,dragCoefficient);

        //}// end if mask
    }// end loop particles
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
