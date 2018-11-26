#ifndef LKrAnalysis_G_h
#define LKrAnalysis_G_h

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>

#include "UserMethods.hh"
#include "TRecoLKrEvent.hh"
#include "TRecoSpectrometerCandidate.hh"
#include "LKrCandidate.hh"
#include "L0TPData.hh"

class AnalysisTools_G;

class LKrAnalysis_G 
{
  public :
    LKrAnalysis_G(Int_t,NA62Analysis::Core::BaseAnalysis *);
    virtual ~LKrAnalysis_G();
    void StartBurst(Int_t,Int_t);
    UserMethods *GetUserMethods() { return fUserMethods;};
    void Clear(Int_t);

    // 2-photon analysis
    Int_t PhotonAnalysis(Bool_t,Int_t,TRecoLKrEvent*,L0TPData*);
    Int_t GetNEMDoublets() {return fNEMDoublets;};
    Int_t *GetTwoEMClustersID() {return fTwoEMClustersID;};

    // Track cluster matching for pinunu analysis
//    Int_t TrackClusterMatching(Bool_t,Int_t,Double_t,TRecoLKrEvent*,TRecoSpectrometerCandidate*);
//    Int_t TrackCellMatching(Bool_t,Int_t,Double_t,TRecoLKrEvent*,TRecoSpectrometerCandidate*);
    Int_t TrackClusterMatching(Bool_t,Int_t,Double_t,Double_t,TRecoLKrEvent*,TVector3*);
    Int_t TrackCellMatching(Bool_t,Int_t,Double_t,Double_t,TRecoLKrEvent*,TVector3*);
    LKrCandidate *GetLKrClusterCandidate() {return fLKrClusterCandidate;};
    LKrCandidate *GetLKrCellCandidate() {return fLKrCellCandidate;};

    // Photon search for pinunu analysis
    Int_t ExtraClusters(Double_t,TRecoLKrEvent*,TVector2*,TVector3*,Int_t);
    LKrCandidate *GetPhotonCandidate(Int_t j) {return fPhotonCandidate[j];};
    Int_t FindNewClusters(Double_t,TRecoLKrEvent*,TVector2*,TVector3*);
    LKrCandidate *GetNewClusterCandidate(Int_t j) {return fNewClusterCandidate[j];};

  public:

  private:
    AnalysisTools_G *fTools;
    UserMethods *fUserMethods;
    Double_t **fLKrCellT0;
    TRecoLKrEvent *fLKrEvent;
    Int_t fFlag;
    L0TPData *fL0Data;
    Int_t fYear;

  private:  
    Int_t fNEMClusters;
    Int_t fNMIPClusters;
    Int_t fNOtherClusters;
    Int_t fEMClusterID[10];          
    Int_t fMIPClusterID[10];
    Int_t fOtherClusterID[10];
    Int_t fNEMDoublets;
    Int_t *fTwoEMClustersID; 

    LKrCandidate *fLKrClusterCandidate;
    LKrCandidate *fLKrCellCandidate;

    LKrCandidate *fPhotonCandidate[10];

    LKrCandidate *fNewClusterCandidate[5];

  private:
    void ClusterAnalysis(Bool_t,Int_t,TRecoLKrEvent*);
    Bool_t SelectEMCluster(Int_t,TRecoLKrCandidate*);
    Bool_t SelectMIPCluster(Int_t,TRecoLKrCandidate*);
    Bool_t SelectOtherCluster(Int_t,TRecoLKrCandidate*,Bool_t,Bool_t);
    inline void AddEMClusters() {fNEMClusters++;};
    inline void AddMIPClusters() {fNMIPClusters++;};
    inline void AddOtherClusters() {fNOtherClusters++;};
    inline void SetEMClusterID(Int_t val) {fEMClusterID[fNEMClusters]=val;};
    inline void SetMIPClusterID(Int_t val) {fMIPClusterID[fNMIPClusters]=val;};
    inline void SetOtherClusterID(Int_t val) {fOtherClusterID[fNOtherClusters]=val;};
    Bool_t TwoPhotonCandidate(TRecoLKrEvent*);
    Bool_t CloseCluster(Int_t,Double_t,Double_t,TRecoLKrEvent*);
     
    Bool_t SelectTiming(Double_t,Double_t); 
  
   // Alternative reconstruction
   private: 
   Int_t Nearby(Int_t,Int_t);
   void ReadCell(Double_t,TRecoLKrEvent*,TVector2*,TVector3*);
   void OrderRawColumn();
   Int_t Ncell;
   Double_t Tcells[14000],Ecells[14000],Xcells[14000],Ycells[14000];
   Int_t xcellcoor[14000],ycellcoor[14000];
 
};

#endif
