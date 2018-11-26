#ifndef IRCAnalysis_h
#define IRCAnalysis_h

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>

#include "UserMethods.hh"
#include "TRecoIRCEvent.hh"
#include "PhotonVetoCandidate.hh"

//class Parameters;
//class AnalysisTools;

class IRCAnalysis2 
{
  public :
    IRCAnalysis2(NA62Analysis::Core::BaseAnalysis *);
     ~IRCAnalysis2();
    UserMethods *GetUserMethods() { return fUserMethods;};
    void StartBurst(Int_t);
    void Clear();
    PhotonVetoCandidate *GetPhotonCandidate(Int_t j) {return fPhotonCandidate[j];};

  public:
     Int_t MakeCandidate(Double_t,TRecoIRCEvent*);

  private:
//    Parameters *fPar;
//    AnalysisTools *fTools;
    UserMethods *fUserMethods;
    PhotonVetoCandidate *fPhotonCandidate[20];

  private:
    Int_t fIRCPriorityMask[16];

  private:
};

#endif
