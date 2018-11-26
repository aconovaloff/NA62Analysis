#ifndef RICHImprovedRing_H
#define RICHImprovedRing_H 1

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMinuit.h"
#include "AnalysisTools_G.hh"
using namespace std;

class AnalysisTools_G;
class TRecoRICHEvent;
class TRecoRICHCandidate;
class TRecoRICHHit;

class RICHImprovedRing
{
  public:
    RICHImprovedRing();
    ~RICHImprovedRing();
    Bool_t Chi2Fit(TRecoRICHEvent*,TRecoRICHCandidate*,TVector3,TVector3,Double_t **);
    static void RingChi2FCN(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);
    Double_t RingChi2(Double_t *);
    Double_t GetDeltaMisa(Int_t i, Int_t j) {return fDeltaMisa[i][j];};

  private:
    AnalysisTools_G *tools;
    Bool_t UpdateRingFit(TRecoRICHCandidate*,TVector3,TVector3);
    void InitMirrorPositions(); 
    Bool_t CorrectHitPosition(TRecoRICHHit*,TVector3,TVector3);
    TVector2 GetChPosAngCorr(Int_t);
    TVector3 ComputeCorrection(Int_t,Int_t,Double_t,Double_t);

  private:
    Int_t fNPars;
    TMinuit *fFitter;
//    Double_t fDeltaMisa[100][2];
    Double_t **fDeltaMisa;
    Double_t p[25][6][2];
    Double_t v[25][25][25][2];
    Double_t ac[25][28];
    Double_t bcost[25][28];
    Int_t fakemirror[25][25];
};

#endif
