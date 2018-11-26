#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "Parameters_R.hh"
using namespace std;
using namespace NA62Constants;

Parameters_R* Parameters_R::fInstance = 0;
Parameters_R::Parameters_R()
{ 
  // CHOD parameters
  fCHODPosV = new Double_t[64];
  fCHODPosH = new Double_t[64];
  fCHODLengthV = new Double_t[64];
  fCHODLengthH = new Double_t[64];
  fCHODT0 = new Double_t[128];
  fCHODVelocity = new Double_t[128];
  fCHODAllT0 = new Double_t*[128];
  fCHODAllSlewSlope = new Double_t*[128];
  fCHODAllSlewConst = new Double_t*[128];
  for (Int_t k=0; k<128; k++) {
    fCHODAllT0[k] = new Double_t[16];
    fCHODAllSlewSlope[k] = new Double_t[16];
    fCHODAllSlewConst[k] = new Double_t[16];
  }
  fIsCHODStored = 0;

  // RICH parameters
  fRICHDeltaMisa = new Double_t*[25];  
  fRICHPos = new Double_t**[25];  
  for (Int_t k=0; k<25; k++) {
    fRICHDeltaMisa[k] = new Double_t[2];
    fRICHPos[k] = new Double_t*[6];
    for (Int_t j=0; j<6; j++) {
      fRICHPos[k][j] = new Double_t[2];
    }
  }
  fIsRICHStored = 0;

  // LKr parameters
  fLKrCellT0 = new Double_t*[128];
  for (Int_t k=0; k<128; k++) {
    fLKrCellT0[k] = new Double_t[128];
  }
  fIsLKrStored = 0;

  // GTK parameters
  fGigaTrackerBadPixel = new Bool_t*[3];
  for (Int_t k=0; k<3; k++) {
    fGigaTrackerBadPixel[k] = new Bool_t[18000];
  }  
  fIsGTKStored = 0;

  // LAV parameters
  fLAVBadChannel = new Bool_t[125000]; 
  fLAVChannelT0 = new Double_t[125000]; 
  fIsLAVStored = 0;

  // SAC parameters
  fSACChannelT0 = new Double_t[4];
  fIsSACStored = 0;
}

int Parameters_R::LoadParameters_R(TString fileparname) {
  ifstream parFile(fileparname.Data());
  if ( !parFile.is_open()) return 0; 
  TString Line;
  while(Line.ReadLine(parFile)) {
    if(Line.BeginsWith("#")) continue;
    TObjArray * l = Line.Tokenize(" ");
    TString name = ((TObjString*)(l->At(0)))->GetString();
    double value = ((TObjString*)(l->At(1)))->GetString().Atof();
    fParameters_R.insert(std::pair<TString,double>(name,value));    
  }
  parFile.close();
  return 1;
}


Parameters_R *Parameters_R::GetInstance() {
  if (fInstance==0) {
    fInstance = new Parameters_R();
  }

  return fInstance;
}

double Parameters_R::GetValue(TString name) {
  NA62Analysis::NA62Map<TString,double>::type::const_iterator ptr;
  if ((ptr=fParameters_R.find(name))!=fParameters_R.end()) return ptr->second;
  return -9999999.;
}

// Detector related parameters repository //

void Parameters_R::StoreCHODParameters_R() {
  if (fIsCHODStored) return;

  Double_t pos[16];
  Double_t slabLength[16];
  for (Int_t jslab=0; jslab<16; jslab++) {
   TString namevar1;  namevar1.Form("chod_slab%d_pos",jslab);
   TString namevar2;  namevar2.Form("chod_slab%d_width",jslab);
   pos[jslab] = GetValue(namevar1.Data()); 
   slabLength[jslab] =GetValue(namevar2.Data())*2.;
  }
  for (Int_t j=0; j<16; j++)  fCHODPosV[j] = -pos[j];
  for (Int_t j=16; j<32; j++) fCHODPosV[j] = -pos[31-j];
  for (Int_t j=32; j<48; j++) fCHODPosV[j] = pos[j-32];
  for (Int_t j=48; j<64; j++) fCHODPosV[j] = pos[63-j];
  for (Int_t j=0; j<16; j++)  fCHODPosH[j] = pos[15-j];
  for (Int_t j=16; j<32; j++) fCHODPosH[j] = -pos[j-16];
  for (Int_t j=32; j<48; j++) fCHODPosH[j] = -pos[47-j];
  for (Int_t j=48; j<64; j++) fCHODPosH[j] = pos[j-48];
  for (Int_t j=0; j<16; j++)  fCHODLengthV[j] = slabLength[j];
  for (Int_t j=16; j<32; j++) fCHODLengthV[j] = slabLength[31-j];
  for (Int_t j=32; j<48; j++) fCHODLengthV[j] = slabLength[j-32];
  for (Int_t j=48; j<64; j++) fCHODLengthV[j] = slabLength[63-j];
  for (Int_t j=0; j<16; j++)  fCHODLengthH[j] = slabLength[15-j];
  for (Int_t j=16; j<32; j++) fCHODLengthH[j] = slabLength[j-16];
  for (Int_t j=32; j<48; j++) fCHODLengthH[j] = slabLength[47-j];
  for (Int_t j=48; j<64; j++) fCHODLengthH[j] = slabLength[j-48];
  for (Int_t jchannel=0; jchannel<128; jchannel++) {
   TString namevar1;  namevar1.Form("chod_channel%d_t0",jchannel);
   TString namevar2;  namevar2.Form("chod_channel%d_velocity",jchannel);
   fCHODT0[jchannel] = GetValue(namevar1.Data()); 
   fCHODVelocity[jchannel] = GetValue(namevar2.Data());
  }

  // T0 taking into account the impact point correction
  for (Int_t jpair=0; jpair<2048; jpair++) {
    TString namevar1;  namevar1.Form("pairt0%d",jpair);
    TString namevar2;  namevar2.Form("slewslope%d",jpair);
    TString namevar3;  namevar3.Form("slewconst%d",jpair);
    Int_t islab = jpair/16;
    Int_t ip = jpair%16;
    Double_t t0 = GetValue(namevar1.Data());
    Double_t sslope = GetValue(namevar2.Data());
    Double_t sconst = GetValue(namevar3.Data());
    if (fabs(t0)>999) t0 = 0.0;
    fCHODAllT0[islab][ip] = t0;
    fCHODAllSlewSlope[islab][ip] = sslope;
    fCHODAllSlewConst[islab][ip] = sconst;
  }    
  fIsCHODStored = 1;
}

void Parameters_R::StoreRICHParameters_R() {
  if (fIsRICHStored) return;
  for (Int_t jmirror=0; jmirror<25; jmirror++) {
    for (Int_t jxy=0; jxy<2; jxy++) {
      TString namevar1;  namevar1.Form("mirroroff%d_%d",jmirror,jxy);
      Double_t mis = GetValue(namevar1.Data()); 
      fRICHDeltaMisa[jmirror][jxy] = mis;
      for (Int_t jcorn=0; jcorn<6; jcorn++) {
        TString namevar2;  namevar2.Form("mirrorpos%d_%d_%d",jmirror,jcorn,jxy);
        Double_t pos = GetValue(namevar2.Data()); 
        fRICHPos[jmirror][jcorn][jxy] = pos;
      }
    }
  } 
  fIsRICHStored = 1;
}

void Parameters_R::StoreLKrParameters_R() {
  if (fIsLKrStored) return;
  for (Int_t jxcell=0; jxcell<128; jxcell++) {
    for (Int_t jycell=0; jycell<128; jycell++) {
      TString namevar1; namevar1.Form("cell%d_%d",jxcell,jycell);
      Double_t t0 = GetValue(namevar1.Data());
      fLKrCellT0[jxcell][jycell] = t0;
    }
  }
  fIsLKrStored = 1;
}

void Parameters_R::StoreGigaTrackerParameters_R() {
  if (fIsGTKStored) return;
  for (Int_t jst=0; jst<3; jst++) {
    for (Int_t jpix=0; jpix<18000; jpix++) {
      TString namevar; namevar.Form("pixel_%d_%d",jst,jpix);
      fGigaTrackerBadPixel[jst][jpix] = 0;
      if (GetValue(namevar.Data())==1) fGigaTrackerBadPixel[jst][jpix] = 1; 
    }
  }
  fIsGTKStored = 1;
}

void Parameters_R::StoreLAVParameters_R() {
  if (fIsLAVStored) return;
  for (Int_t jch=0; jch<125000; jch++) {
    TString namevar1; namevar1.Form("channel_%d",jch);
    fLAVBadChannel[jch] = 0;
    if (GetValue(namevar1.Data())==1) fLAVBadChannel[jch] = 1; 
    TString namevar2; namevar2.Form("channelt0_%d",jch);
    Double_t t0 = GetValue(namevar2.Data());
    if (t0==-9999999.) t0 = 0;
    fLAVChannelT0[jch] = t0; 
  }
  fIsLAVStored = 1;
}

void Parameters_R::StoreSACParameters_R() {
  if (fIsSACStored) return;
  for (Int_t jch=0; jch<4; jch++) {
    TString namevar; namevar.Form("sac_channelt0_%d",jch);
    Double_t t0 = GetValue(namevar.Data());
    if (t0==-9999999.) t0 = 0;
    fSACChannelT0[jch] = t0; 
  }
  fIsSACStored = 1;
}
