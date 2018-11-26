#ifndef IRCANALYSIS_h
#define IRCANALYSIS_h

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "UserMethods.hh"
#include "TRecoIRCEvent.hh"

class IRCAnalysis 
{
  public :
 	IRCAnalysis(NA62Analysis::Core::BaseAnalysis *);
	void Set(TRecoIRCEvent*, Double_t, Double_t, Double_t, bool); 
     	Int_t Veto();
        void BookHistos_Set();
	void SaveAllPlots();  
	UserMethods* fUserMethods;
  	~IRCAnalysis();
  private:
    	TRecoIRCEvent *event;
    	Double_t reftime;
    	Double_t tLow;
    	Double_t tHi;  
	bool MCflag;
};

#endif
