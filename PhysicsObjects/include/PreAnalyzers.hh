#ifndef PREANALYZERS_h
#define PREANALYZERS_h
#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TRecoLKrEvent.hh"
#include "UserMethods.hh"
class PreAnalyzers
{
  public :
        PreAnalyzers();
        void LKrCorrections(TRecoLKrEvent*&);
	void LKrMToG(TRecoLKrEvent*&);
	void StrawCorrections(TRecoSpectrometerEvent*&); 
        ~PreAnalyzers();
  private:
	TRecoLKrEvent *LKrEvent;
	TRecoSpectrometerEvent *SpectrometerEvent; 
};

#endif

