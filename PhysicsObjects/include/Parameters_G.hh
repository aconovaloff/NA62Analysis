#ifndef PARAMETERS_G_HH
#define PARAMETERS_G_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include <TCanvas.h>
using namespace std;

class Parameters_G
{

public:
  Parameters_G();
  static Parameters_G* GetInstance();

public:
  int LoadParameters_G(TString);
  double GetValue(TString); 

public: // Detector related parameters
  void StoreCHODParameters_G();
  Double_t *GetCHODPosV() {return fCHODPosV;};
  Double_t *GetCHODPosH() {return fCHODPosH;};
  Double_t *GetCHODLengthV() {return fCHODLengthV;};
  Double_t *GetCHODLengthH() {return fCHODLengthV;};
  Double_t *GetCHODT0() {return fCHODT0;};
  Double_t *GetCHODVelocity() {return fCHODVelocity;};
  Double_t **GetCHODAllT0() {return fCHODAllT0;};
  Double_t **GetCHODAllFineT0() {return fCHODAllFineT0;};
  Double_t **GetCHODAllSlewSlope() {return fCHODAllSlewSlope;};
  Double_t **GetCHODAllSlewConst() {return fCHODAllSlewConst;};

  void StoreRICHParameters_G();
  Double_t **GetRICHDeltaMisa() {return fRICHDeltaMisa;};
  Double_t ***GetRICHPos() {return fRICHPos;};

  void StoreLKrParameters_G();
  Double_t **GetLKrCellT0() {return fLKrCellT0;};

  void StoreGigaTrackerParameters_G();
  Bool_t **GetGigaTrackerBadPixel() {return fGigaTrackerBadPixel;};

  void StoreLAVParameters_G();
  Bool_t *GetLAVBadChannel() {return fLAVBadChannel;};
  Double_t *GetLAVChannelT0() {return fLAVChannelT0;};

  void StoreSACParameters_G();
  Double_t *GetSACChannelT0() {return fSACChannelT0;};

  void StoreSAVParameters_G();
  Double_t *GetSAVChannelT0() {return fSAVChannelT0;};

private:
  static Parameters_G* fInstance;

private:
  NA62Analysis::NA62Map<TString,double>::type fParameters_G;

private: // Detector related parameters
  Double_t *fCHODPosV;
  Double_t *fCHODPosH;
  Double_t *fCHODLengthV;
  Double_t *fCHODLengthH;
  Double_t *fCHODT0;
  Double_t *fCHODVelocity;
  Double_t **fCHODAllT0;
  Double_t **fCHODAllFineT0;
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

  Double_t *fSAVChannelT0;
  Bool_t fIsSAVStored;
};
#endif
