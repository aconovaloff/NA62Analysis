#ifndef CHANTIAnalysis_h
#define CHANTIAnalysis_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TRecoCHANTIEvent.hh"
#include "UserMethods.hh"

class CHANTIAnalysis
{
  public :
        CHANTIAnalysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoCHANTIEvent*, Double_t, Double_t, Double_t, bool);
	void SaveAllPlots();
	void BookHistos_Set();
        Int_t Veto();
        ~CHANTIAnalysis();
	UserMethods* fUserMethods;
  private:
        TRecoCHANTIEvent *event;
        Double_t reftime;
        Double_t tLow;
        Double_t tHi;
	bool MCflag;
};

#endif
