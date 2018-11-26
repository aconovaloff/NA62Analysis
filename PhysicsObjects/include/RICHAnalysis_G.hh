#ifndef RICHAnalysis_G_h
#define RICHAnalysis_G_h

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TVector3.h"
#include "AnalysisTools_G.hh"
#include "UserMethods.hh"
#include "TRecoRICHEvent.hh"
#include "TRecoSpectrometerCandidate.hh"
#include "RICHCandidate.hh"

//class AnalysisTools;
class RICHImprovedRing;

class RICHAnalysis_G 
{
  public :
    RICHAnalysis_G(NA62Analysis::Core::BaseAnalysis *);
    virtual ~RICHAnalysis_G();
    UserMethods *GetUserMethods() { return fUserMethods;};
    void Clear();
    Int_t TrackMultiRingMatching(Double_t,TRecoRICHEvent*,TRecoSpectrometerCandidate*, TVector3, TVector3);
    Int_t TrackSingleRingMatching(Double_t,TRecoRICHEvent*,TRecoSpectrometerCandidate*, TVector3, TVector3);
    RICHCandidate *GetRICHMultiCandidate() {return fRICHMultiCandidate;};
    RICHCandidate *GetRICHSingleCandidate() {return fRICHSingleCandidate;};
///new 
   Bool_t RICHSingleCandidate(RICHCandidate, TRecoSpectrometerCandidate*, TVector3, TVector3);
   Bool_t RICHMultiCandidate(RICHCandidate, TRecoSpectrometerCandidate*, TVector3, TVector3);
   Double_t DistCheck;
   Double_t GetDistCheck() {return DistCheck;}; 
   Double_t singlerichmass;
   Double_t multirichmass;
   Double_t Getsinglerichmass() {return singlerichmass;};
   Double_t Getmultirichmass() {return multirichmass;};
   void SaveAllPlots(); 
  public:

  private:
   AnalysisTools_G *fTools;
    UserMethods *fUserMethods;
    RICHCandidate *fRICHMultiCandidate;
    RICHCandidate *fRICHSingleCandidate;
    RICHImprovedRing *fImprovedRing;

  private:
    Double_t **fDeltaMisa; 
    Double_t ***fMirrorPos;
 
  private:
    TRecoRICHEvent *fRICHEvent;
    void MirrorVector(TVector3*,TVector3*);
    Int_t MirrorSurface(Double_t,Double_t,Double_t,Bool_t);
};

#endif
