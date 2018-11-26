#ifndef StrawAnalysis_h
#define StrawAnalysis_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TRecoSpectrometerEvent.hh"
#include "UserMethods.hh"
class StrawAnalysis
{
  public :
        StrawAnalysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoSpectrometerEvent*, Double_t, Double_t, Double_t, bool);
	void SaveAllPlots();
	void BookHistos_MatchedTime();
	void BookHistos_SingleTrack();
	Int_t MatchedTime(Double_t);  
        Int_t SingleTrack();
	TVector3 GetPionMomentum();
	TVector3 GetPionPosition();
	TVector3 GetPionMomentumProp();
        TVector3 GetPionPositionProp();
	TLorentzVector GetPionFourMomentum();  
        ~StrawAnalysis();
	UserMethods* fUserMethods;
	Double_t GetTime(){return ftime;}; 
 	TLorentzVector ComputeMomentum(TRecoSpectrometerCandidate*, Double_t);
 	TLorentzVector ComputeMomentumAfterMagnet(TRecoSpectrometerCandidate*, Double_t);
      	TVector3 ComputePosition(TRecoSpectrometerCandidate*, Double_t); 

  private:
	Double_t ftime;
        TRecoSpectrometerEvent *event;
        Double_t reftime;
        Double_t tLow;
        Double_t tHi;
	TLorentzVector Pi4;
	TVector3 Pi3;
	TVector3 PiPos;
        TLorentzVector Pi4Prop;
        TVector3 Pi3Prop;
        TVector3 PiPosProp;
	bool MCflag; 
//	TLorentzVector ComputeMomentum(TRecoSpectrometerCandidate*, Double_t);
//	TVector3 ComputePosition(TRecoSpectrometerCandidate*, Double_t); 
};

#endif

