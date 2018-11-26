#ifndef GIGATRACKEREVTRECO2_HH
#define GIGATRACKEREVTRECO2_HH

// This analyzer can only be run as a pre-analyzer
#pragma pre-analyzer

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
class TAxis;


class GigaTrackerEvtReco2 : public NA62Analysis::Analyzer
{
public:
  GigaTrackerEvtReco2(NA62Analysis::Core::BaseAnalysis *ba);
//Int_t fEventNumber;
  void InitHist();
  //  void InitHist2();
  void InitOutput();
  void DefineMCSimple();
  void Process(int iEvent);
  void StartOfBurstUser();
  void EndOfBurstUser();
  void StartOfRunUser();
  void EndOfJobUser();
  void EndOfRunUser();
  void PostProcess();
  void DrawPlot();

  template<typename Order>
  void Clusterize(vector< int >, vector< vector<int> >&, double minDist, Order order, int minCont);
  void BuildCandidate(TRecoGigaTrackerCandidate* cand);
  void LinearLeastSquareFit(Double_t *x, Double_t *y, Int_t Nsample, Double_t *sigma, Double_t &a, Double_t &b, Double_t &rho, Double_t &chi2);
  void FillHisto(TRecoGigaTrackerCandidate* cand);

  //--------------
  struct TimeOrder{
    GigaTrackerEvtReco2* fGTKReco;
    TimeOrder(GigaTrackerEvtReco2* p) : fGTKReco(p) {};
    bool operator() ( int i, int j ){
      TRecoGigaTrackerHit* h1 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(i);
      if(h1 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;}
      TRecoGigaTrackerHit* h2 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(j);
      if(h2 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;};
      return h1->GetTime() < h2->GetTime();
    }
    double dist(int i, int j){
      TRecoGigaTrackerHit* h1 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(i);
      if(h1 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;}
      TRecoGigaTrackerHit* h2 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(j);
      if(h2 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;};
      return abs(h1->GetTime() - h2->GetTime());
    }
  };

  //-------------
  struct XOrder{
    GigaTrackerEvtReco2* fGTKReco;
    XOrder(GigaTrackerEvtReco2* p) : fGTKReco(p) {};
    bool operator() ( int i, int j ){
      TRecoGigaTrackerHit* h1 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(i);
      if(h1 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;}
      TRecoGigaTrackerHit* h2 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(j);
      if(h2 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;};
      return h1->GetPosition().X() < h2->GetPosition().X();
    }
    double dist(int i, int j){
      TRecoGigaTrackerHit* h1 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(i);
      if(h1 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;}
      TRecoGigaTrackerHit* h2 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(j);
      if(h2 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;};
      return abs(h1->GetPosition().X() - h2->GetPosition().X());
    }

  };

  //-------------
  struct YOrder{
    GigaTrackerEvtReco2* fGTKReco;
    YOrder(GigaTrackerEvtReco2* p) : fGTKReco(p) {};
    bool operator() ( int i, int j ){
      TRecoGigaTrackerHit* h1 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(i);
      if(h1 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;}
      TRecoGigaTrackerHit* h2 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(j);
      if(h2 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;};
      return h1->GetPosition().Y() < h2->GetPosition().Y();
    }
    double dist(int i, int j){
      TRecoGigaTrackerHit* h1 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(i);
      if(h1 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;}
      TRecoGigaTrackerHit* h2 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(j);
      if(h2 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;};
      return abs(h1->GetPosition().Y() - h2->GetPosition().Y());
    }

  };

  //-------------
  struct StationOrder{
    GigaTrackerEvtReco2* fGTKReco;
    StationOrder(GigaTrackerEvtReco2* p) : fGTKReco(p) {};
    bool operator() ( int i, int j ){
      TRecoGigaTrackerHit* h1 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(i);
      if(h1 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;}
      TRecoGigaTrackerHit* h2 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(j);
      if(h2 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;};
      return h1->GetStationNo()<h2->GetStationNo();

    }
    double dist(int i, int j){
      TRecoGigaTrackerHit* h1 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(i);
      if(h1 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;}
      TRecoGigaTrackerHit* h2 = (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(j);
      if(h2 == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return 0;};
      return abs(h1->GetStationNo()-h2->GetStationNo());

    }

  };

  //-------------
  struct Cluster{
    double X,Y,Z,T;
    int N,S;
    vector<int> hits;
    GigaTrackerEvtReco2* fGTKReco;
    Cluster(GigaTrackerEvtReco2* p):X(0),Y(0),Z(0),T(0),N(0), fGTKReco(p) { };
    
    void add(int i){
      TRecoGigaTrackerHit* h= (TRecoGigaTrackerHit*) fGTKReco->fGigaTrackerEvent->GetHit(i);
      if(h == NULL)  {std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl; return;}
      if(N>0 && h->GetStationNo()!=S){ std::cout<<"Cluster Error: trying to add hit from different stations"<<endl; return;}
      if(N==0) S = h->GetStationNo();
      X = (N*X+h->GetPosition().X())/(N+1);
      Y = (N*Y+h->GetPosition().Y())/(N+1);
      Z = (N*Z+h->GetPosition().Z())/(N+1);
      T = (N*T+h->GetTime()        )/(N+1);
      N = N+1;
      hits.push_back(i);
  }

    ~Cluster(){
      hits.clear();
    }

  };

protected:
//  TH1F* EventType; 
  bool fNoisy[3][18000];
  double fTOffsetChip[3][10];

  TH2F* hit_map;
  TH2F* t_block_map;
  TH2F* cluster_map;

  TH2F* hit_map2;
  TH2F* t_block_map2;
  TH2F* cluster_map2;

  TH1D* fHNHitPerPixel[3];
  double fBins[100];

  TH1D* fHDt[3];
  TH2I* fHMap[3];
  TH2I* fHMapNoise[3];
  TH2D* fHDtToT[3];

  TH1D* fHMomentum[4];
  TH1D* fHThX[4];
  TH1D* fHThY[4];
  TH1D* fHChi2[4];
  TH1D* fHChi2X[4];
  TH1D* fHChi2Y[4];
  TH1D* fHChi2T[4];

  TH1D* fHDX_12[4];
  TH1D* fHDX_13[4];
  TH1D* fHDX_23[4];

  TH1D* fHDY_12[4];
  TH1D* fHDY_13[4];
  TH1D* fHDY_23[4];

  TH1D* fHDt_12[4];
  TH1D* fHDt_13[4];
  TH1D* fHDt_23[4];

  TH2D* fHDX_X_12[4];
  TH2D* fHDY_Y_12[4];

  TH2D* fHDX_X_13[4];
  TH2D* fHDY_Y_13[4];

  TH2D* fHDX_X_23[4];
  TH2D* fHDY_Y_23[4];


  TRecoGigaTrackerEvent* fGigaTrackerEvent;
  double fTimeWindow;
  double fTimeWindowTrigger;
  double fXWindow;
  double fYWindow;
  double fChi2X;
  double fChi2Y;
  double fChi2T;
  bool   fRemoveHit;
  bool   fRefineTime;

  TString fTzeroFolder;
  TString fTWalkFolder;
  TString fXYFolder;


  Double_t fBLbend    ; //- Integrated magnetic mield mor bending magnets [T * mm]
  Double_t fDeltaBend ; //- Distance between bending magnets doublets [mm]
  Double_t fBeta      ; //- [mm * MeV/c]
  Double_t fBLtrim    ; //- Integrated magnetic mield mor trim magnet [T * mm]
  Double_t fDeltaTrim ; //- Distance between trim magnet and GTK3 [mm]
  Double_t fKaonMass  ; //- Kaon mass in MeV/c
  Double_t fDelta12 ;   //- Distance between GTK1 and GTK2 22800 [mm]
  Double_t fDelta13 ;   //- Distance between GTK1 and GTK3 13200 [mm]
  Double_t fDelta23 ;   //- Distance between GTK2 and GTK3 [mm]
  Double_t fAlpha   ;   //- [Adimensional] 
  Double_t fDeltaZ  ;   //- shift in Z between station 1,3 and 2
  Double_t fShiftTrim;
  double   fClight;

  int fMissingChips[30];
  int fNMissingChips;
  float fT0[3][18000];
  float fTW[3][10][410];
  TAxis fTWbins;
  TVector3 fPosOff[3];
  bool fRedoTimeCorr;
  bool fRedoXYCorr;

  Long64_t fStartDraw;
  int fFileNumber;
  Long64_t fNDraw;
  bool fNewBurst;

  static int fNPrints;
  TCanvas* fCv;
private:
  Bool_t fWarnOnce;
};
#endif
