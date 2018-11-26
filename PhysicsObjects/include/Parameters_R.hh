#ifndef PARAMETERS_R_HH
#define PARAMETERS_R_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
using namespace std;

class Parameters_R
{

public:
  Parameters_R();
  static Parameters_R* GetInstance();

public:
  int LoadParameters_R(TString);
  double GetValue(TString); 

public: // Detector related parameters
  void StoreCHODParameters_R();
  Double_t *GetCHODPosV() {return fCHODPosV;};
  Double_t *GetCHODPosH() {return fCHODPosH;};
  Double_t *GetCHODLengthV() {return fCHODLengthV;};
  Double_t *GetCHODLengthH() {return fCHODLengthV;};
  Double_t *GetCHODT0() {return fCHODT0;};
  Double_t *GetCHODVelocity() {return fCHODVelocity;};
  Double_t **GetCHODAllT0() {return fCHODAllT0;};
  Double_t **GetCHODAllSlewSlope() {return fCHODAllSlewSlope;};
  Double_t **GetCHODAllSlewConst() {return fCHODAllSlewConst;};

  void StoreRICHParameters_R();
  Double_t **GetRICHDeltaMisa() {return fRICHDeltaMisa;};
  Double_t ***GetRICHPos() {return fRICHPos;};

  void StoreLKrParameters_R();
  Double_t **GetLKrCellT0() {return fLKrCellT0;};

  void StoreGigaTrackerParameters_R();
  Bool_t **GetGigaTrackerBadPixel() {return fGigaTrackerBadPixel;};

  void StoreLAVParameters_R();
  Bool_t *GetLAVBadChannel() {return fLAVBadChannel;};
  Double_t *GetLAVChannelT0() {return fLAVChannelT0;};

  void StoreSACParameters_R();
  Double_t *GetSACChannelT0() {return fSACChannelT0;};

private:
  static Parameters_R* fInstance;

private:
  NA62Analysis::NA62Map<TString,double>::type fParameters_R;

private: // Detector related parameters
  Double_t *fCHODPosV;
  Double_t *fCHODPosH;
  Double_t *fCHODLengthV;
  Double_t *fCHODLengthH;
  Double_t *fCHODT0;
  Double_t *fCHODVelocity;
  Double_t **fCHODAllT0;
  Double_t **fCHODAllSlewSlope;
  Double_t **fCHODAllSlewConst;
  Bool_t fIsCHODStored;

  Double_t **fRICHDeltaMisa;
  Double_t ***fRICHPos;
  Bool_t fIsRICHStored;

  Double_t **fLKrCellT0;
  Bool_t fIsLKrStored;

  Bool_t **fGigaTrackerBadPixel;
  Bool_t fIsGTKStored;

  Bool_t *fLAVBadChannel;
  Double_t *fLAVChannelT0;
  Bool_t fIsLAVStored;

  Double_t *fSACChannelT0;
  Bool_t fIsSACStored;
};
#endif
