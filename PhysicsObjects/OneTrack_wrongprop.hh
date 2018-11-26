#ifndef OneTrack_h
#define OneTrack_h
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

class OneTrack
{
  public :
        OneTrack(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoSpectrometerEvent*, TRecoCHODEvent*, TRecoCedarEvent*, TRecoGigaTrackerEvent*, TRecoCHANTIEvent*, TRecoLKrEvent*, TRecoRICHEvent*, L0TPData*, Bool_t, Int_t);
	Int_t OneTrackSelection();
	Int_t BackgroundSelection(); 
	void SaveAllPlots();
	void BookHistos();
        ~OneTrack();
	UserMethods* fUserMethods;
	Double_t finalmmsq;
	Double_t Getmmsq(){return finalmmsq;};
	Double_t Reft;
	Double_t GetReft(){return Reft;};
	Double_t RICHmass;
	Double_t GetRICHmass(){return RICHmass;};
   	Double_t RingRadius;
   	Double_t GetRingRadius() {return RingRadius;};
	TRecoLKrCandidate* LKrCand;
	TRecoLKrCandidate* GetLKrCand(){return LKrCand;};
	TVector3 PionPosition;
	TVector3 GetPionPosition(){return PionPosition;};
	TVector3 PionMomentum; 
	TVector3 GetPionMomentum(){return PionMomentum;};
	TVector3 LKr_pos;
        TVector3 GetLKrPosition(){return LKr_pos;};
	LKrAnalysis* LKra; 
	GigaTrackerAnalysis* Giga;
	StrawAnalysis* Strawa;
	Int_t MatchedRICHCandidate(TRecoRICHEvent*, TRecoSpectrometerCandidate*, TVector3, TVector3, Double_t);
  private:
//	StrawAnalysis* Strawa;
	CHODAnalysis* CHODa;
	CedarAnalysis* Cedara;
//	GigaTrackerAnalysis* Giga;
	CHANTIAnalysis* CHANTIa;
//	LKrAnalysis* LKra;
	BlueTubeTrackerAnalysis* Bluea;
	RICHAnalysis_G* RICHa_G;
	TRecoSpectrometerEvent* SpectrometerEvent;
	TRecoCHODEvent* CHODEvent;
	TRecoCedarEvent* CedarEvent;
	TRecoGigaTrackerEvent* GigaTrackerEvent;
	TRecoCHANTIEvent* CHANTIEvent;
	TRecoLKrEvent* LKrEvent;
	TRecoRICHEvent* RICHEvent;
	L0TPData* L0;
	Bool_t MCflag;
	Int_t RunNumber;
        TString acc[6];
	TString pinn_mmsq[4]; 
//	TH2I *acc[6];  
//	TLorentzVector ComputeMomentum(TRecoSpectrometerCandidate*, Double_t);
//	TVector3 ComputePosition(TRecoSpectrometerCandidate*, Double_t); 
};

#endif

