#ifndef DUMMY_HH
#define DUMMY_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include "MCSimple.hh"
#include "DetectorAcceptance.hh"
#include <TCanvas.h>

class TH1I;
class TH2F;
class TGraph;
class TTree;


class dummy : public NA62Analysis::Analyzer
{
	public:
		dummy(NA62Analysis::Core::BaseAnalysis *ba);
		~dummy();
		void InitHist();
		void InitOutput();
		void DefineMCSimple();
		void ProcessSpecialTriggerUser(int iEvent, unsigned int triggerType);
		void Process(int iEvent);
		void PostProcess();
		void StartOfBurstUser();
		void EndOfBurstUser();
		void StartOfRunUser();
		void EndOfRunUser();
        void EndOfJobUser();
		void DrawPlot();
	protected:


};
#endif
