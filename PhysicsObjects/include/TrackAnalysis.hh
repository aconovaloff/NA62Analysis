#ifndef TrackAnalysis_h
#define TrackAnalysis_h
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
#include "PhotonAnalysis.hh"
#include "CalorimeterAnalysis.hh"
#include "OneTrack.hh"
#include "CalorimeterPlotGenerator.hh"

class TrackAnalysis
{
  public :
        TrackAnalysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoCedarEvent*, TRecoGigaTrackerEvent*, TRecoCHANTIEvent*, TRecoSpectrometerEvent*, TRecoCHODEvent*, TRecoLKrEvent*, TRecoRICHEvent*, TRecoMUV1Event*, TRecoMUV2Event*, TRecoMUV3Event*,TRecoLAVEvent*,  TRecoIRCEvent*, TRecoSACEvent*, TRecoSAVEvent*, L0TPData*, Bool_t, Int_t);
	Int_t PhotonRejectionFactor();
	void MUV3Efficiency();
	void Kmu2Normalization();  
	void CalorimeterRejectionFactor();
	void RICHRejectionFactor();
	void PionKineFactor();
	void MuonKineFactor();
	void RICHAnalysisPion();
	void RICHAnalysisMuon();
	void MuonAccidentalPhotonFactor();
	void BeamAccidentalPhotonFactor();  
	void SaveAllPlots();
	void BookHistos();
	void CalorimeterPlots();
	void ep(); 
        ~TrackAnalysis();
	UserMethods* fUserMethods;
/*
	Double_t finalmmsq;
	Double_t Getmmsq(){return finalmmsq;};
	Double_t Reft;
	Double_t GetReft(){return Reft;};
	Double_t RICHmass;
	Double_t GetRICHmass(){return RICHmass;};
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
	Int_t MatchedRICHCandidate(TRecoRICHEvent*, TRecoSpectrometerCandidate*, TVector3, TVector3, Double_t);
*/	
	OneTrack* One;
	CalorimeterAnalysis* Calo;
	PhotonAnalysis* Photon;
  private:
	TString mp_bin[6];
	TString pp_bin[6];
        TString pp_bin_post[6];
	TString mp_bin_post[6];
	LKrAnalysis* LKra;  
	MUV3Analysis* MUV3a;
//	StrawAnalysis* Strawa;
//	CHODAnalysis* CHODa;
//	CedarAnalysis* Cedara;
//	GigaTrackerAnalysis* Giga;
//	CHANTIAnalysis* CHANTIa;
//	LKrAnalysis* LKra;
//	BlueTubeTrackerAnalysis* Bluea;
//	RICHAnalysis_G* RICHa_G;
	TRecoSpectrometerEvent* SpectrometerEvent;
	TRecoCHODEvent* CHODEvent;
	TRecoCedarEvent* CedarEvent;
	TRecoGigaTrackerEvent* GigaTrackerEvent;
	TRecoCHANTIEvent* CHANTIEvent;
	TRecoLKrEvent* LKrEvent;
	TRecoRICHEvent* RICHEvent;
	TRecoMUV1Event* MUV1Event;
	TRecoMUV2Event* MUV2Event;
	TRecoMUV3Event* MUV3Event;
	TRecoLAVEvent* LAVEvent;
	TRecoIRCEvent* IRCEvent;
	TRecoSACEvent* SACEvent;
	TRecoSAVEvent* SAVEvent;	
	L0TPData* L0;
	Bool_t MCflag;
	Int_t RunNumber; 
	CalorimeterPlotGenerator* CaloPlot;
	TString epCal[6];
	TString epRICH[6];
	TString epCal_post[6];
	TString epRICH_post[6];  
//	TLorentzVector ComputeMomentum(TRecoSpectrometerCandidate*, Double_t);
//	TVector3 ComputePosition(TRecoSpectrometerCandidate*, Double_t); 

};

#endif

