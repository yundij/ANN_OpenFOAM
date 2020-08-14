keras::DataChunkFlat DFdataChunk( DFNFeatures, 0.0 );
std::vector<float>& DFfeatures = DFdataChunk.get_1d_rw();

forAll( driftFlux_, iCell )
{

  float DFtarget = 1.0e-6;
  unsigned int index = 0;
  DFfeatures.clear();
   //solid volume fraction at the cell 
   DFfeatures.push_back( ( solidVolumeFraction[iCell]/alphaMax - DFmeans_[index])/DFvars_[index] );
   ++index;

   DFfeatures.push_back((gradP[iCell].z()/9.8/rhoParticle_-DFmeans_[index])/DFvars_[index]);
   ++index;

   //uSlip.z at the center 
   DFfeatures.push_back((uSlip[iCell].z()/uTerminal-DFmeans_[index])/DFvars_[index]);
   ++index;
   //evaluate the neural network model for the current point
   std::vector<float> DFresponse = DFnnModel_.compute_output( &DFdataChunk );
   DFtarget = DFresponse.at(0);
   if(iCell==215)
     {
	Info<<"icell:"<<iCell<<endl;
	Info << "DFfeature values:" << endl;
	for(int i =0; i<DFNFeatures; ++i) Info << DFfeatures[i] << " ";
	Info << "DFresponse values:" << DFresponse.at(0) << endl;
     }
   driftFlux_[iCell].x()=0.0;
   driftFlux_[iCell].y()=0.0;
   driftFlux_[iCell].z()=DFtarget*alphaMax*uTerminal;

}
