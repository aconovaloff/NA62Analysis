#ifndef RICHAnalysis_h
#define RICHAnalysis_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TRecoRICHEvent.hh"
#include "UserMethods.hh"

class RICHAnalysis
{
  public :
        RICHAnalysis(NA62Analysis::Core::BaseAnalysis *);
        void Set(TRecoRICHEvent*, Double_t, Double_t, Double_t, bool);
	void BookHistos_PionID();
	void BookHistos_MatchedTrack();
	void BookHistos_Set();  
	void SaveAllPlots();
        Int_t PionID(TVector3, TVector3);
	Int_t MuonID(TVector3, TVector3);
	Int_t MatchedTrack(TVector3, TLorentzVector, TLorentzVector);
	Double_t mass_RICH;
	Double_t Getmass_RICH(){return mass_RICH;};
	Double_t pick_RICH_mmsq;
	Double_t Get_RICH_mmsq(){return pick_RICH_mmsq;}; 
        ~RICHAnalysis();
	UserMethods* fUserMethods; 
  private:
        TRecoRICHEvent *event;
        Double_t reftime;
        Double_t tLow;
        Double_t tHi;
	bool MCflag; 
};

#endif
