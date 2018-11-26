#ifndef GigaTrackerAnalysis_h
#define GigaTrackerAnalysis_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TRecoGigaTrackerEvent.hh"
#include "UserMethods.hh"

class GigaTrackerAnalysis
{
  public :
	GigaTrackerAnalysis(NA62Analysis::Core::BaseAnalysis *); 
        void Set(TRecoGigaTrackerEvent*, Double_t, Double_t, Double_t, bool);
	void SaveAllPlots();
	void BookHistos_MatchedTrack();
	void BookHistos_Set();  
        Int_t MatchedTrack(TVector3, TVector3, Int_t);
	TVector3 GetKaonMomentum();
	TVector3 GetKaonPosition();
	TLorentzVector GetKaonFourMomentum();  
        ~GigaTrackerAnalysis();
	UserMethods* fUserMethods;
  private:
        TRecoGigaTrackerEvent *event;
        Double_t reftime;
        Double_t tLow;
        Double_t tHi;
	TLorentzVector K4;
	TVector3 K3;
	TVector3 KPos;
	bool MCflag;
};

#endif

