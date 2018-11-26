#ifndef PhotonAnalysis_h
#define PhotonAnalysis_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "CHODAnalysis.hh"
#include "StrawAnalysis.hh"
#include "CHANTIAnalysis.hh"
#include "LKrAnalysis.hh"
#include "LKrAnalysis_G.hh"
#include "SAVAnalysis.hh"
#include "RICHAnalysis_G.hh"
#include "CedarAnalysis.hh"
#include "GigaTrackerAnalysis.hh"
#include "BlueTubeTrackerAnalysis.hh"
#include "UserMethods.hh"
#include "MUV1Analysis.hh"
#include "MUV2Analysis.hh"
#include "MUV3Analysis.hh"
#include "SAVMatching.hh"
#include "LAVMatching.hh"
#include "LAVAnalysis.hh"
#include "IRCAnalysis.hh"
#include "SACAnalysis.hh"

class PhotonAnalysis
{
  public :
        PhotonAnalysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoLKrEvent*, TRecoLAVEvent*, TRecoIRCEvent*, TRecoSACEvent*,TRecoSAVEvent*, TRecoLKrCandidate*, Double_t, Double_t, Int_t, TVector3, Bool_t);
	Int_t Rejection();
	void SaveAllPlots();
	void BookHistos();
        ~PhotonAnalysis();
	UserMethods* fUserMethods;
	LKrAnalysis_G* LKra_G;
	SAVAnalysis* SAVa; 
	LAVAnalysis* LAVa;
	SACAnalysis* SACa;
	IRCAnalysis* IRCa; 
	Int_t OtherVetoes(TRecoLAVEvent*, TRecoIRCEvent*, TRecoSACEvent*, TRecoSAVEvent*, Double_t, Bool_t);  
  private:
//	LKrAnalysis_G* LKra_G;
//	SAVAnalysis* SAVa; 
	TRecoLAVEvent* LAVEvent;
	TRecoIRCEvent* IRCEvent;
	TRecoSACEvent* SACEvent;
	TRecoSAVEvent* SAVEvent;
	TRecoLKrEvent* LKrEvent; 
	TRecoLKrCandidate* LKrCand; 	
	TVector3 LKr_pos;  
	Double_t Reft;
	Double_t mmsq;
	Int_t UnmatchedClusters;
	Bool_t MCflag;
//	LAVMatching* fLAVMatching;
//	SAVMatching* fSAVMatching;
//	TLorentzVector ComputeMomentum(TRecoSpectrometerCandidate*, Double_t);
//	TVector3 ComputePosition(TRecoSpectrometerCandidate*, Double_t); 
};

#endif

