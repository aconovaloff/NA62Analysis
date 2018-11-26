#ifndef CedarAnalysis_h
#define CedarAnalysis_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TRecoCedarEvent.hh"
#include "UserMethods.hh"
class CedarAnalysis
{
  public :
        CedarAnalysis(NA62Analysis::Core::BaseAnalysis *);
        void Set(TRecoCedarEvent*, Double_t, Double_t, Double_t, bool);
	void SaveAllPlots();
	void BookHistos_MatchedKaon();
	void BookHistos_Set();  
        Int_t MatchedKaon();
	Int_t MatchedKaon_Straw(Double_t);
	Int_t GetnSectors();
        ~CedarAnalysis();
	UserMethods* fUserMethods;
	Double_t GetTime(){return fTime;}; 
	Double_t fTime;
  private:
        TRecoCedarEvent *event;
        Double_t reftime;
        Double_t tLow;
        Double_t tHi;
	Int_t nSectors;
	bool MCflag; 
};

#endif

