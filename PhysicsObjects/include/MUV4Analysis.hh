#ifndef MUV4Analysis_h
#define MUV4Analysis_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TRecoMUV2Event.hh"
#include "UserMethods.hh"
class MUV4Analysis
{
  public :
        MUV4Analysis(NA62Analysis::Core::BaseAnalysis *);
	void BookHistos(); 
        void FillHistos();
	void SaveHistos(); 
    UserMethods *fUserMethods;
        ~MUV4Analysis();
  private:
};

#endif

