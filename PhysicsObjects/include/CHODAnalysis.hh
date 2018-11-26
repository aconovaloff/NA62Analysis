#ifndef CHODAnalysis_h
#define CHODAnalysis_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TRecoCHODEvent.hh"
#include "UserMethods.hh"

class CHODAnalysis
{
  public :
        CHODAnalysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoCHODEvent*, Double_t, Double_t, Double_t, bool);
	void BookHistos_MatchedTrack();
	void SaveAllPlots();
	void BookHistos_Set();
	void BookHistos_MinDisc();
	Int_t MinDisc(TVector3, Double_t); 
        Int_t MatchedTrack(TVector3);
	Int_t MatchedTrack_t(TVector3);
	Double_t GetTime();
        ~CHODAnalysis();
	UserMethods* fUserMethods; 
	Int_t MinDiscID;
	Int_t GetMinDiscID(){return MinDiscID;}; 
  private:
        TRecoCHODEvent *event;
        Double_t reftime;
        Double_t tLow;
        Double_t tHi;
	Int_t minDCandCHOD;
	bool MCflag;
};

#endif
