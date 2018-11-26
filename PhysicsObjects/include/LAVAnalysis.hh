#ifndef LAVAnalysis_h
#define LAVAnalysis_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TRecoLAVEvent.hh"
#include "UserMethods.hh"

class LAVAnalysis
{
  public :
        LAVAnalysis(NA62Analysis::Core::BaseAnalysis *);
        void Set(TRecoLAVEvent*, Double_t, Double_t, Double_t, bool);
	void SaveAllPlots();
        Int_t Veto();
	void BookHistos_Set(); 
        ~LAVAnalysis();
	UserMethods* fUserMethods;
  private:
        TRecoLAVEvent *event;
        Double_t reftime;
        Double_t tLow;
        Double_t tHi;
	bool MCflag;
};

#endif

