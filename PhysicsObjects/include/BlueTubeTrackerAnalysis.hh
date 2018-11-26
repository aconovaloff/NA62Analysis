#ifndef BLUETUBETRACKERANALYSIS_h
#define BLUETUBETRACKERANALYSIS_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "BlueTubeTracker.hh"
#include "UserMethods.hh"

class BlueTubeTrackerAnalysis
{
  public :
        BlueTubeTrackerAnalysis(NA62Analysis::Core::BaseAnalysis *);
        void Set(BlueTubeTracker*, TVector3, TVector3, TVector3, TVector3);
//	void SaveAllPlots();
//	void BookHistos_Set(); 
        ~BlueTubeTrackerAnalysis();
	UserMethods* fUserMethods;
 	TLorentzVector pionMomentumCorr;
	TLorentzVector GetMomentumCorrection() {return pionMomentumCorr;};
	TVector3 vertexCorr;
	TVector3 GetVertexCorrection() {return vertexCorr;};
  private:
        BlueTubeTracker *Tracker;
	TVector3 K_pos;
	TVector3 DecayVertex;
        TVector3 PionPosition;
	TVector3 PionMomentum; 
};

#endif

