#ifndef CalorimeterAnalysis_h
#define CalorimeterAnalysis_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "CHODAnalysis.hh"
#include "StrawAnalysis.hh"
#include "CHANTIAnalysis.hh"
#include "LKrAnalysis.hh"
#include "RICHAnalysis_G.hh"
#include "CedarAnalysis.hh"
#include "GigaTrackerAnalysis.hh"
#include "BlueTubeTrackerAnalysis.hh"
#include "UserMethods.hh"
#include "MUV1Analysis.hh"
#include "MUV2Analysis.hh"
#include "MUV3Analysis.hh"

class CalorimeterAnalysis
{
  public :
        CalorimeterAnalysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoMUV1Event*, TRecoMUV2Event*, TRecoMUV3Event*, TRecoLKrCandidate*, Double_t, Double_t, TVector3, TVector3, Bool_t);
	Int_t Rejection();
	Int_t MuonSelection();
	Int_t MUV3Efficiency();  
	void SaveAllPlots();
	void BookHistos();
	void Phil(); 
        ~CalorimeterAnalysis();
	UserMethods* fUserMethods;
  private:
	TRecoMUV1Event* MUV1Event;
	TRecoMUV2Event* MUV2Event;
	TRecoMUV3Event* MUV3Event;
	MUV1Analysis* MUV1a;
	MUV2Analysis* MUV2a;
	MUV3Analysis* MUV3a;
	TRecoLKrCandidate* LKrCand; 	
	Double_t Reft;
	Double_t mmsq;
	TVector3 PionPosition;
	TVector3 PionMomentum;
	Bool_t MCflag;
//	TLorentzVector ComputeMomentum(TRecoSpectrometerCandidate*, Double_t);
//	TVector3 ComputePosition(TRecoSpectrometerCandidate*, Double_t); 
};

#endif

