#ifndef AnalysisTools_G_G_H
#define AnalysisTools_G_G_H 1

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
#include "TVector3.h"
#include "TLorentzVector.h"
using namespace std;

//class BTubeField;
//class BlueTubeTracker;
class TRecoSpectrometerCandidate;
class TRecoLKrCandidate;
//class Particle;

class AnalysisTools_G 
{
  public:
  AnalysisTools_G();
  static AnalysisTools_G* GetInstance();
  
  private:
  static AnalysisTools_G* fInstance;

  public:
  TVector3 SingleTrackVertex(TVector3*,TVector3*,TVector3*,TVector3*,Double_t*);
  TVector3 MultiTrackVertexSimple(Int_t,TLorentzVector*,TVector3*,Double_t*); 
  TLorentzVector Get4Momentum(Double_t*);
  TVector3 GetPositionAtZ(TRecoSpectrometerCandidate*, Double_t);
  void ClusterCorrections(Int_t,TRecoLKrCandidate *);
  TVector3 Propagate(Int_t,TVector3 *,TVector3 *,Double_t *,Double_t);
  void BTubeCorrection(TVector3*,Double_t*,Double_t*,Double_t);
 // TVector3 BlueFieldCorrection(Particle*,Double_t,Double_t); 
  TVector3 BlueFieldCorrection(TVector3*,TVector3,Int_t,Double_t); 
  Bool_t LKrAcceptance(Double_t,Double_t,Double_t,Double_t);
  void MirrorVector(TVector3*,TVector3*);
  Int_t MirrorSurface(Double_t,Double_t,Double_t,Bool_t);
  Int_t GetLAVStation(Double_t);
  Int_t TriggerMask(Int_t);
  
  public:
//  BTubeField *fBTubeField;
//  BlueTubeTracker *fTracker;
  Double_t fMirrorPos[25][6][2];

  private:
  void InitializeMirrorCoordinate();  
  Bool_t fIsClusterCorrected;
};  
    
#endif
