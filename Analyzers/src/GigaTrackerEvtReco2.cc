// ---------------------------------------------------------------
//
// History:
//
// Created by Mathieu Perrin-Terrin (mathieu.perrin-terrin@cern.ch) 01/03/2016
//
// ---------------------------------------------------------------

/// \class GigaTrackerEvtReco2
/// \Brief
/// Rebuild the TRecoGigaTrackerEvent
/// \EndBrief
/// \Detailed
/// Rebuild the TRecoGigaTrackerEvent. Noisy pixels are identified at the
/// beginning of each burst and masked. The possibility to apply time and
/// space corrections is maintained but disabled by default. For all data
/// reprocessed with a revision below or equal to 1299 these corrections 
/// must be enabled by setting RedoCorr to 1.
/// Hits out of time (>5ns) from the trigger are removed if RemoveHit flag is set to 1.
/// The new TRecoGigaTrackerEvent is filled with TRecoGigaTrackerCandidates
/// that can be of type 123, 12, 13, 23, depending if they are made from hit
/// in GTK1, GTK2 and GTK3. (NOW ONLY 123)
/// Usage: TRecoGigaTrackerEvent* gigaTrackerEvent = (TRecoGigaTrackerEvent*)GetOutput("GigaTrackerEvtReco2.GigaTrackerEvent", isEvtValid);
/// \author  Mathieu Perrin-Terrin (mathieu.perrin-terrin@cern.ch)
/// \EndDetailed

#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include <TF1.h>
#include <TKey.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "GigaTrackerEvtReco2.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "BaseAnalysis.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

int GigaTrackerEvtReco2::fNPrints = 0;

//-------------------------------
GigaTrackerEvtReco2::GigaTrackerEvtReco2(Core::BaseAnalysis *ba) : Analyzer(ba, "GigaTrackerEvtReco2"), fWarnOnce(false){

  fNewBurst = 1;
  fStartDraw = 0;
  fNDraw = 30000;
  fFileNumber = -1;

  fGigaTrackerEvent = new TRecoGigaTrackerEvent;
  RequestTree("GigaTracker", fGigaTrackerEvent,"Reco");

  for(int iS(0); iS<3; iS++){
    for(int iP(0); iP<18000; iP++){
      fNoisy[iS][iP]=0;
    }
    for(int iC(0); iC<10; iC++){
      fTOffsetChip[iS][iC]=0;
    }
  }


  // Basic Geometrical Parameters for the Reconstruction
  double	magnetFieldStrength0 = -1.6678;//Tm
  double	magnetFieldStrength4 = -0.7505;//Tm
  double	magnetZLength1       = 2.5 * 1.e3;
  double	magnetZLength2       = 0.4 * 1.e3;
  double	detectorZPosition    = 0.5 * (79.580 + 102.420) * 1.e3;
  double	magnet1PositionZ     = 82.980 * 1.e3 - 0.5 * magnetZLength1 - detectorZPosition;
  double	magnet2PositionZ     = 86.580 * 1.e3 - 0.5 * magnetZLength1 - detectorZPosition;
  double	magnet5PositionZ     = 102.000 * 1.e3 - 0.5 * magnetZLength2 - detectorZPosition;
  double	stationZLength	     = 200e-6*1e3;
  double	station1PositionZ    = 79.600 * 1.e3 - 0.5 * stationZLength - detectorZPosition;
  double	station2PositionZ    = 92.800 * 1.e3 - 0.5 * stationZLength - detectorZPosition;
  double	station3PositionZ    = 102.400 * 1.e3 - 0.5 * stationZLength - detectorZPosition;
  double	station2PositionY    = -60.e-3 * 1.e3;


  // Derived Parameters
  fClight                 = TMath::C();
  fBLbend	       	  = -1.0*magnetFieldStrength0*magnetZLength1;
  fDeltaBend		  = magnet2PositionZ - magnet1PositionZ;
  fBeta			  = 1e-3*fClight * fBLbend * fDeltaBend ;
  fBLtrim		  = magnetFieldStrength4*magnetZLength2;
  fDeltaTrim		  = station3PositionZ - magnet5PositionZ;
  fKaonMass		  = 493.667;
  fDelta12		  = station2PositionZ - station1PositionZ;
  fDelta13		  = station3PositionZ - station1PositionZ;
  fDelta23		  = station3PositionZ - station2PositionZ;
  fAlpha		  = fDelta12 / fDelta13;
  fDeltaZ		  = -station2PositionY;
  fShiftTrim              = -((1e-3*fClight * fBLtrim) / 75e9) * fDeltaTrim;

  // Algorithm Paramters
  AddParam("TimeWindow"	                , &fTimeWindow                      ,  0.7	); //ns
  AddParam("TimeWindowTrigger"          , &fTimeWindowTrigger               ,  3	); //ns
  AddParam("XWindow"	                , &fXWindow                         ,  10	); //mm
  AddParam("YWindow"	                , &fYWindow                         ,  10	); //mm
  AddParam("Chi2X" 	                , &fChi2X                           ,  20	); //chi2 cut
  AddParam("Chi2Y" 	                , &fChi2Y                           ,  20	); //chi2 cut
  AddParam("Chi2T" 	                , &fChi2T                           ,  20	); //chi2 cut
  AddParam("RemoveHit"                  , &fRemoveHit                       ,  0        );
  AddParam("RefineTime"                 , &fRefineTime                      ,  1        );

  AddParam("RedoTimeCorr"               , &fRedoTimeCorr                    ,0  );
  AddParam("RedoXYCorr"                 , &fRedoXYCorr                      ,0  );
  AddParam("TZeroFolder"                , &fTzeroFolder                     ,Form("%s/../NA62Analysis/config/gtk_time-corrections/",getenv ("NA62MCSOURCE")));
  AddParam("TWalkFolder"                , &fTWalkFolder                     ,Form("%s/../NA62Analysis/config/gtk_time-corrections/",getenv ("NA62MCSOURCE")));
  AddParam("XYFolder"                   , &fXYFolder                        ,Form("%s/../NA62Analysis/config/gtk_xy-corrections/",getenv ("NA62MCSOURCE")));

}

//-------------------------------
void GigaTrackerEvtReco2::InitOutput(){
  RegisterOutput("GigaTrackerEvent", fGigaTrackerEvent);
//  RegisterOutput("EventNumber", fEventNumber);
}

//-------------------------------
void GigaTrackerEvtReco2::InitHist(){

  for (int i(0); i< 98;i++) fBins[i] = i*3;
  fBins[98] = 1e3;
  fBins[99] = 10e3;
//  EventType = new TH1F("EventType","EventType",0,0,200);
//  BookHisto(EventType);
  for(int iS(0); iS<3; iS++){
    fHDt[iS]       = new TH1D(Form("hDt_Station%d",iS),Form("Station %d; t_{hit}-t_{trigger};count",iS),160*4,-40,40);
    fHDtToT[iS]    = new TH2D(Form("mDtToT_Station%d",iS),Form("Station %d;t_{1} - t_{2};ToT [ns]",iS),160*4,-40,+40,256,0,2*ClockPeriod);
    fHMap[iS]      = new TH2I(Form("hMap_Station%d",iS),Form("Station %d; X Channel; Y Channel",iS),200,0,200,90,0,90);
    fHMapNoise[iS] = new TH2I(Form("hMapNoise_Station%d",iS),Form("Noise Station %d; X Channel; Y Channel",iS),200,0,200,90,0,90);    

    fHNHitPerPixel[iS] = new TH1D(Form("mHNHitPerPixelStation%d",iS),Form("Station %d;NHit;NPixel",iS),250,0,500);

    BookHisto(fHDt[iS]);
    BookHisto(fHMap[iS]);
    BookHisto(fHMapNoise[iS]);
    BookHisto(fHDtToT[iS]);
    BookHisto(fHNHitPerPixel[iS]);

    
  }

  int tTypes[4]={123,12,13,23};
  for(int iT(0); iT<4; iT++){
    fHMomentum[iT]   = new TH1D(Form("mHMomentum_%d",tTypes[iT]),Form("Type %d ;p [GeV];count",tTypes[iT]),200,60,90);
    fHThX[iT]        = new TH1D(Form("mHThX_%d",tTypes[iT]),Form("Type %d;#theta_{X} - 1.2 [mrad];count",tTypes[iT]),100,-1,1);
    fHThY[iT]        = new TH1D(Form("mHThY_%d",tTypes[iT]),Form("Type %d;#theta_{Y} [mrad];count",tTypes[iT]),100,-1,1);
    fHChi2[iT]       = new TH1D(Form("mHChi2_%d",tTypes[iT]),Form("Type %d;#chi^{2};count",tTypes[iT]),200,0,500);
    fHChi2X[iT]      = new TH1D(Form("mHChi2X_%d",tTypes[iT]),Form("Type %d;#chi^{2} X;count",tTypes[iT]),200,0,500);
    fHChi2Y[iT]      = new TH1D(Form("mHChi2Y_%d",tTypes[iT]),Form("Type %d;#chi^{2} Y;count",tTypes[iT]),200,0,500);
    fHChi2T[iT]      = new TH1D(Form("mHChi2T_%d",tTypes[iT]),Form("Type %d;#chi^{2} T;count",tTypes[iT]),200,0,500);

    fHDX_12[iT]   = new TH1D(Form("mDX_12_%d",tTypes[iT]),Form("Type %d;X_{1} - X_{2};count",tTypes[iT]),200,-25,25);
    fHDX_13[iT]   = new TH1D(Form("mDX_13_%d",tTypes[iT]),Form("Type %d;X_{1} - X_{3};count",tTypes[iT]),200,-25,25);
    fHDX_23[iT]   = new TH1D(Form("mDX_23_%d",tTypes[iT]),Form("Type %d;X_{2} - X_{3};count",tTypes[iT]),200,-25,25);

    fHDY_12[iT]   = new TH1D(Form("mDY_12_%d",tTypes[iT]),Form("Type %d;Y_{1} - Y_{2};count",tTypes[iT]),200,-25,25);
    fHDY_13[iT]   = new TH1D(Form("mDY_13_%d",tTypes[iT]),Form("Type %d;Y_{1} - Y_{3};count",tTypes[iT]),200,-25,25);
    fHDY_23[iT]   = new TH1D(Form("mDY_23_%d",tTypes[iT]),Form("Type %d;Y_{2} - Y_{3};count",tTypes[iT]),200,-25,25);

    fHDt_12[iT]   = new TH1D(Form("mDt_12_%d",tTypes[iT]),Form("Type %d;t_{1} - t_{2};count",tTypes[iT]),200,-5,+5);
    fHDt_13[iT]   = new TH1D(Form("mDt_13_%d",tTypes[iT]),Form("Type %d;t_{1} - t_{3};count",tTypes[iT]),200,-5,+5);
    fHDt_23[iT]   = new TH1D(Form("mDt_23_%d",tTypes[iT]),Form("Type %d;t_{2} - t_{3};count",tTypes[iT]),200,-5,+5);


    fHDX_X_12[iT]   = new TH2D(Form("mDX_X_12_%d",tTypes[iT]),Form("Type %d;X_{1};X_{1} - X_{2}",tTypes[iT]),200,-40,40,200,-25,25);
    fHDX_X_13[iT]   = new TH2D(Form("mDX_X_13_%d",tTypes[iT]),Form("Type %d;X_{1};X_{1} - X_{3}",tTypes[iT]),200,-40,40,200,-25,25);
    fHDX_X_23[iT]   = new TH2D(Form("mDX_X_23_%d",tTypes[iT]),Form("Type %d;X_{2};X_{2} - X_{3}",tTypes[iT]),200,-40,40,200,-25,25);

    fHDY_Y_12[iT]   = new TH2D(Form("mDY_Y_12_%d",tTypes[iT]),Form("Type %d;Y_{1};Y_{1} - Y_{2}",tTypes[iT]),200,-25,25,200,-25,25);
    fHDY_Y_13[iT]   = new TH2D(Form("mDY_Y_13_%d",tTypes[iT]),Form("Type %d;Y_{1};Y_{1} - Y_{3}",tTypes[iT]),200,-25,25,200,-25,25);
    fHDY_Y_23[iT]   = new TH2D(Form("mDY_Y_23_%d",tTypes[iT]),Form("Type %d;Y_{2};Y_{2} - Y_{3}",tTypes[iT]),200,-25,25,200,-25,25);



    BookHisto(fHMomentum[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHThX[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHThY[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHChi2[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHChi2X[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHChi2Y[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHChi2T[iT],0,Form("Type%d",tTypes[iT]));
    
    BookHisto(fHDX_12[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHDX_13[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHDX_23[iT],0,Form("Type%d",tTypes[iT]));
    
    BookHisto(fHDY_12[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHDY_13[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHDY_23[iT],0,Form("Type%d",tTypes[iT]));
    
    BookHisto(fHDt_12[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHDt_13[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHDt_23[iT],0,Form("Type%d",tTypes[iT]));

    BookHisto(fHDX_X_12[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHDX_X_13[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHDX_X_23[iT],0,Form("Type%d",tTypes[iT]));

    BookHisto(fHDY_Y_12[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHDY_Y_13[iT],0,Form("Type%d",tTypes[iT]));
    BookHisto(fHDY_Y_23[iT],0,Form("Type%d",tTypes[iT]));

  }
}

//-------------------------------
void GigaTrackerEvtReco2::DefineMCSimple(){
}

void GigaTrackerEvtReco2::StartOfRunUser(){

  // Protection against multiple application of the corrections
  if (IsAnalyzerInHistory(GetAnalyzerName())) {
    if (!fWarnOnce) std::cout << user() << "Corrections already applied" << std::endl;
    fWarnOnce = true;
    return;
  }

  if (!GetIsTree()) return;
  if (GetWithMC()) return; // no corrections for MC

  EventHeader* rawHeader = (EventHeader*) GetEventHeader();
  int runId = rawHeader->GetRunID();
  FILE *fp;

  if(fRedoTimeCorr == 1 ){

    cout<<"Load Time Corrections Lookup Table"<<endl;

    //LOAD LOOKUP-TABLE T0
    TString tmp = Form("%s/switchboard",fTzeroFolder.Data());
    fp = fopen(tmp.Data(),"r");
    if(fp == NULL ) {
      cout<<"Cannot load the lookup table file "<<tmp<<endl;
      exit(1);  
    }
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    bool foundCorr = 0;
    string corrTzeroVersion;
    while ((read = getline(&line, &len, fp)) != -1) {
      char *line2 = line;
      char *p = strtok(line2, " ");
      string v = string(p);
      p = strtok(NULL, " ");
      int run1 = atoi(p);
      p = strtok(NULL, " ");
      int run2 = atoi(p);

      if(runId>run1 && runId<run2){
	foundCorr = 1;
	corrTzeroVersion = v;
      }
    }
    fclose(fp);

    //LOAD LOOKUP-TABLE TWalk
    tmp = Form("%s/switchboard",fTWalkFolder.Data());
    fp = fopen(tmp.Data(),"r");
    if(fp == NULL ) {
      cout<<"Cannot load the twalk lookup table file "<<tmp<<endl;
      exit(1);  
    }
    line = NULL;
    len = 0;
    foundCorr = 0;
    string corrTWalkVersion;
    while ((read = getline(&line, &len, fp)) != -1) {
      char *line2 = line;
      char *p = strtok(line2, " ");
      string v = string(p);
      p = strtok(NULL, " ");
      int run1 = atoi(p);
      p = strtok(NULL, " ");
      int run2 = atoi(p);

      if(runId>run1 && runId<run2){
	foundCorr = 1;
	corrTWalkVersion = v;
      }
    }
     

 
    if (foundCorr == 0) exit(1);
    cout<<"Run "<<runId<<" -> Taking Corrections from folders :"<<endl;
    cout<<"  T0: "<<fTzeroFolder.Data()<<corrTzeroVersion.c_str()<<endl;
    cout<<"  TW: "<<fTWalkFolder.Data()<<corrTWalkVersion.c_str()<<endl;


    for(int gtk(0); gtk<3; gtk ++){
      for(int chip(0); chip<10; chip++){
      
	//TZERO
	tmp = Form("%s/%s/t-zero_gtk%d-chip%d",fTzeroFolder.Data(),corrTzeroVersion.c_str(),gtk,chip); 
	fp = fopen(tmp.Data(),"r");
	if(fp == NULL ) {
	  cout<<"Cannot load the lookup table file "<<tmp<<endl;
	  exit(1);  
	}
      
	line = NULL;
	len = 0;
	while ((read = getline(&line, &len, fp)) != -1) {
	  char *line2 = line;
	  char *p = strtok(line2, " ");
	  int uid = atoi(p);
	  p = strtok(NULL, " ");
	  float t0 = atof(p);
	  fT0[gtk][uid] = t0;
	
	}
	fclose(fp);
      

	//TWALK
	tmp = Form("%s/%s/time-walk_gtk%d-chip%d",fTWalkFolder.Data(),corrTWalkVersion.c_str(),gtk,chip); 
	fp = fopen(tmp.Data(),"r");
	if(fp == NULL )  {
	  cout<<"Cannot load the lookup table file "<<tmp<<endl;
	  exit(1);  
	}

	line = NULL;
	len = 0;
	int ibin(0);
	while ((read = getline(&line, &len, fp)) != -1) {
	  float tw = atof(line);
	  fTW[gtk][chip][ibin] = tw;
	  ibin++;
	}
	fclose(fp);
	if (line)  free(line);
      }
    
    }

  
    // TWALK BINS
    tmp = Form("%s/%s/time-walk_binning",fTWalkFolder.Data(),corrTWalkVersion.c_str());
    fp = fopen(tmp.Data(),"r");
    if(fp == NULL ) {
      cout<<"Cannot load the lookup table file "<<tmp<<endl;
      exit(1);  
    }
  
    line = NULL;
    len = 0;
    int ibin(0);
    double twBins[411];
    while ((read = getline(&line, &len, fp)) != -1) {
      float bin = atof(line);
      twBins[ibin] = bin;
      ibin++;
    }
    fTWbins = TAxis(410,twBins);
    fclose(fp);
    if (line)  free(line);
  }

  if(fRedoXYCorr == 1 ){

    cout<<"Load XY Corrections Lookup Table"<<endl;

    //LOAD LOOKUP-TABLE XY
    TString tmp = Form("%s/switchboard",fXYFolder.Data());
    fp = fopen(tmp.Data(),"r");
    if(fp == NULL ) {
      cout<<"Cannot load the xy lookup table file "<<tmp<<endl;
      exit(1);  
    }

    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    bool foundCorr = 0;
    string corrXYVersion;
    while ((read = getline(&line, &len, fp)) != -1) {
      char *line2 = line;
      char *p = strtok(line2, " ");
      string v = string(p);
      p = strtok(NULL, " ");
      int run1 = atoi(p);
      p = strtok(NULL, " ");
      int run2 = atoi(p);

      if(runId>run1 && runId<run2){
	foundCorr = 1;
	corrXYVersion = v;
      }
    }

    fclose(fp);

    // XY Correction
    tmp = Form("%s/%s/xy_offset.txt",fXYFolder.Data(),corrXYVersion.c_str());
    fp = fopen(tmp.Data(),"r");
    if(fp == NULL )  exit(1);  
  
    line = NULL;
    len = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
      char *p = strtok(line, " ");
      int gtk = atoi(p);
      p = strtok(NULL, " ");
      float xOffset = atof(p);
      p = strtok(NULL, " ");
      float yOffset = atof(p);
      fPosOff[gtk].SetXYZ(xOffset,yOffset,0);

    }
    fclose(fp);

    if (foundCorr == 0) exit(1);
    cout<<"Run "<<runId<<" -> Taking Corrections from folders :"<<endl;
    cout<<"  XY: "<<fXYFolder.Data()<<corrXYVersion.c_str()<<endl;
    if (line)  free(line);
  }

  

}

//-------------------------------
void GigaTrackerEvtReco2::StartOfBurstUser(){
cout<<"right program"<<endl;
  // Protection against multiple application of the corrections
  if (IsAnalyzerInHistory(GetAnalyzerName())) {
    if (!fWarnOnce) std::cout << user() << "Corrections already applied" << std::endl;
    fWarnOnce = true;
    return;
  }

  if (!GetIsTree()) return;

  // no corrections for MC
  if (GetWithMC()) return;

  if(  fParent->GetIOHandler()->GetCurrentFileNumber() != fFileNumber){
    //cout<<"New File! "<<endl;
    fStartDraw = 0;
    fFileNumber = fParent->GetIOHandler()->GetCurrentFileNumber();
  }

  // Get Burst Id
  EventHeader* rawHeader = (EventHeader*) GetEventHeader();
  int burstID = rawHeader->GetBurstID();
  printf("---------- Burst %3d ----------\n",burstID);

  //Preparing for drawing
  TChain* ttree = GetTree("Reco");
  TString variables, condi;

  // Check all chips are here 
  TH1I *hchip = new TH1I("hchip","",30,0,30);
  variables = Form("GigaTracker.fHits.fChipID + 10*GigaTracker.fHits.fStationNo>>hchip");
  ttree->GetTree()->SetEstimate(1e+07);
  ttree->GetTree()->Draw(variables,Form("EventHeader.fBurstID == %d",burstID),"goff",fNDraw,fStartDraw);
  fNMissingChips = 0; //reset
  for( int c(0);  c<30;c++){
    fMissingChips[c] = 0; //reset
    if(hchip->GetBinContent(c+1) == 0){
      fMissingChips[c] = 1;
      fNMissingChips = fNMissingChips+1;
    }
  }
  for(int iS(0); iS<3;iS++){
    if(iS==0) printf(" ________________Missing Chips: %2d ______________\n",fNMissingChips);
    else      printf(" _________________________________________________\n");
    printf("  %6d ⁵| %6d ⁶| %6d ⁷| %6d ⁸| %6d ⁹\n",fMissingChips[iS*10+5],fMissingChips[iS*10+6],fMissingChips[iS*10+7],fMissingChips[iS*10+8],fMissingChips[iS*10+9]);
    printf("%d-------------------------------------------------\n",iS+1);
    printf("  %6d ⁰| %6d ¹| %6d ²| %6d ³| %6d ⁴\n",fMissingChips[iS*10+0],fMissingChips[iS*10+1],fMissingChips[iS*10+2],fMissingChips[iS*10+3],fMissingChips[iS*10+4]);
  }
  hchip->Delete();

  // Check global trigger time offset
  //double triggeroffset = 0;
  TF1  *f1 = new TF1("f1","gaus",-100,+100);
  double tdcCalib = 24.951059536/256.0;
  TH1F *hktag = new TH1F("hktag","",500,-10,10);
  if(  ttree->GetListOfBranches()->Contains("Cedar") ){
    variables = Form("Cedar.fCandidates.fTime-EventHeader.fFineTime*%f",tdcCalib);
    cout<<"[GigaTrackerEvtReco2] Correcting for trigger jumps using Cedar"<<endl;
  }
  else{
    variables = Form("RICH.fCandidates.fTime-EventHeader.fFineTime*%f",tdcCalib);
    cout<<"[GigaTrackerEvtReco2] Correcting for trigger jumps using RICH"<<endl;
  }

  condi = Form("L0TP.fDataType<16 && EventHeader.fBurstID == %d",burstID);
  ttree->GetTree()->SetEstimate(1e+07);
  ttree->GetTree()->Draw(variables,condi,"goff",fNDraw,fStartDraw);


  int nevt = ttree->GetTree()->GetSelectedRows();
  double *dtimektag = ttree->GetTree()->GetV1(); 
  for (int iEv(0); iEv<nevt; iEv++) hktag->Fill(dtimektag[iEv]);
  double tav = hktag->GetXaxis()->GetBinCenter(hktag->GetMaximumBin());
  hktag->Fit(f1,"q0","",tav-0.5,tav+0.5);
  double mean = f1->GetParameter(1);
  double sigm = f1->GetParameter(2);
  cout <<"Trigger offset: "<< mean << " +/- " << sigm << endl;
  //triggeroffset = mean;

  //Check time offset and noisy pixels
  std::vector<double> v_med;
  double med = 0;


  TH1F *h[3][10];
  TH1F* hm[3];
  for(int iS(0); iS<3;iS++){
    for (int iC(0); iC<10; iC++) {
      h[iS][iC] = new TH1F(Form("histgtk%d-%d",iS,iC),"",200,-10,10);
    }
    hm[iS] = new TH1F(Form("hm%d",iS),"",18000,0,18000);
  }



  if(fRedoTimeCorr == 1){
    variables = Form("(GigaTracker.fHits.fRawTime-EventHeader.fFineTime*%f):GigaTracker.fHits.fChannelID:GigaTracker.fHits.fToT",tdcCalib);
    condi = Form("GigaTracker.fHits.fRawTime&& EventHeader.fBurstID == %d",burstID);
  }
  else{
    variables = Form("(GigaTracker.fHits.fTime-EventHeader.fFineTime*%f):GigaTracker.fHits.fChannelID:GigaTracker.fHits.fToT",tdcCalib);
    condi = Form("GigaTracker.fHits.fTime && EventHeader.fBurstID == %d",burstID);
    
  }


  TChain* reco = GetTree("Reco");
  reco->GetTree()->Draw(variables.Data(),condi.Data(),"goff",fNDraw,fStartDraw);
  int nentries = reco->GetTree()->GetSelectedRows();
  cout << "Total number of hits: " << nentries << endl;
  double *dtime = reco->GetTree()->GetV1(); 
  double *pixel = reco->GetTree()->GetV2(); 
  double *tovth = reco->GetTree()->GetV3(); 
  for (int iEv(0); iEv<nentries; iEv++) {
    int gtk =   int(pixel[iEv])/100000;
    int chip = (int(pixel[iEv])%100000)/10000;
    int pix =   int(pixel[iEv])%10000;
    int x = pix%40 + (chip%5)*40;
    int y = pix/40 + (chip/5)*45;
    int uid = x + y*200;
    if(fRedoTimeCorr == 1){
      int totBin    = fTWbins.FindBin(tovth[iEv]) -1;
      double dtcorr = dtime[iEv]- fT0[gtk][uid] - fTW[gtk][chip][totBin];
      h[gtk][chip]->Fill(dtcorr);
      hm[gtk]->Fill(uid);
    }
    else{
      h[gtk][chip]->Fill(dtime[iEv]);
      hm[gtk]->Fill(uid);
    }
  }

  //time correction
  for(int iS(0); iS<3;iS++){
    for (int iC(0); iC<10; iC++) {
      fTOffsetChip[iS][iC] = 0.;
      double xc = h[iS][iC]->GetXaxis()->GetBinCenter(h[iS][iC]->GetMaximumBin()); 
      int nen = h[iS][iC]->Integral();
      if (nen<250){
	cout <<"GTK "<< iS << " Chip " << iC << ": " << nen << " entries, skip " << endl;
	continue;
      }
      h[iS][iC]->Fit(f1,"q0","",xc-0.6,xc+0.6);
      double tmean = f1->GetParameter(1);
      double tsigm = f1->GetParameter(2);
      //      cout <<"GTK "<< iS << " Chip " << iC << ": " << nen << " entries, mean time " << tmean << " reso " << tsigm << endl;
      //      if (fabs(tsigm)<0.7) fTOffsetChip[iS][iC] = tmean-triggeroffset;//-triggeroffset (0.1)?;
      if(fabs(tsigm)<0.6 && fabs(tmean)>2) fTOffsetChip[iS][iC] = (tmean/fabs(tmean))*24.951059536/8.;
    }
  }
  for(int iS(0); iS<3;iS++){
    if(iS==0) printf(" ________________Chip Time Offsets________________\n");
    else      printf(" _________________________________________________\n");
    printf("  %6.3f ⁵| %6.3f ⁶| %6.3f ⁷| %6.3f ⁸| %6.3f ⁹\n",fTOffsetChip[iS][5],fTOffsetChip[iS][6],fTOffsetChip[iS][7],fTOffsetChip[iS][8],fTOffsetChip[iS][9]);
    printf("%d-------------------------------------------------\n",iS+1);
    printf("  %6.3f ⁰| %6.3f ¹| %6.3f ²| %6.3f ³| %6.3f ⁴\n",fTOffsetChip[iS][0],fTOffsetChip[iS][1],fTOffsetChip[iS][2],fTOffsetChip[iS][3],fTOffsetChip[iS][4]);
  }
  

  //noise
  for(int iS(0); iS<3;iS++){
    if(hm[iS]->Integral()<100) continue;
    v_med.clear();
    for(int iX(1); iX<=18000; iX++){
      double iC = hm[iS]->GetBinContent(iX);
      if(iC > 0 ){
	v_med.push_back(iC);
	fHNHitPerPixel[iS]->Fill(iC);
      }
    }
    std::sort(v_med.begin(),v_med.end());


    med = v_med[8*v_med.size()/10];
    int nbMasked = 0;
    for(int iX(1); iX<=18000; iX++){
	double iC = hm[iS]->GetBinContent(iX);
	if (iC>5*med){
	  fNoisy[iS][iX-1] = 1;
	  nbMasked++;
      }

    }
    printf("GTK Reco: %d pixels with more than %f hits in this burst have been declared noisy for station %d\n",nbMasked,20*med,iS);
  }

  //CLEANING
  for(int iS(0); iS<3;iS++){
    delete hm[iS];
    for (int iC(0); iC<10; iC++) {
      delete h[iS][iC];
    }
  }
  delete f1;
  delete hktag;
  fNewBurst = 1;

  return;
}

//-------------------------------
void GigaTrackerEvtReco2::Process(int iEE) {

  // Protection against multiple application of the corrections
//  if (IsAnalyzerInHistory(GetAnalyzerName())) {
 //   if (!fWarnOnce) std::cout << user() << "Corrections already applied" << std::endl;
//    fWarnOnce = true;
//    return;
//  }

  if (!GetIsTree()) return;

  fStartDraw++;

  if (fNewBurst) {
    //cout<<"New Burst!"<<endl;
    GetTree("Reco")->GetEntry(iEE);  
    fNewBurst = 0;
  }

  if(fGigaTrackerEvent==NULL) {
    cout<<"Empty GTK Event!"<<endl;
    return;
  }
//fEventType = 0;
 SetOutputState("GigaTrackerEvent", kOValid);
 //SetOutputState("EventNumber", kOInvalid);
  double triggerTime(0);
/*
  if (!GetWithMC()) {
    EventHeader* rawHeader = (EventHeader*) GetEventHeader();
    if(rawHeader == NULL) std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl;
    triggerTime =  rawHeader->GetFineTime()*ClockPeriod/256.;
  }
*/
  double eventTime = 0;
TRecoGigaTrackerHit       *gigatrackerHit;
/*
  //Set Error Mask if Missing Chips
  if(fNMissingChips>0) {
    ULong64_t errorMask = fGigaTrackerEvent->GetErrorMask();
    errorMask = errorMask | 0x8000000000000000;
    fGigaTrackerEvent->SetErrorMask(errorMask);
  }

  // REMOVE HIT FROM NOISY PIXELS OUT OF TIME HIT AND CLEAN CANDIDATES
  TRecoGigaTrackerHit       *gigatrackerHit;
  int pixelID,stationNo,chip;
  double deltaT, ToT;
  int iH(0), totBin(0);
  while(iH<fGigaTrackerEvent->GetNHits()){
    gigatrackerHit = (TRecoGigaTrackerHit*)fGigaTrackerEvent->GetHit(iH);
    if(gigatrackerHit == NULL) {
      std::cout<<__FILE__<<" at line: "<<__LINE__<<"null pointer"<<endl;
      cout<<"iH "<<iH<<"/"<<fGigaTrackerEvent->GetNHits()<<endl;
      SetOutputState("GigaTrackerEvent", kOInvalid);
      return;
    }
    
    pixelID   = gigatrackerHit->GetPixelID();
    stationNo = gigatrackerHit->GetStationNo();
    chip      = gigatrackerHit->GetChipID();
    ToT       = gigatrackerHit->GetToT();
    totBin    = fTWbins.FindBin(ToT)-1;

    if(stationNo>=3 || pixelID >=18000 ) std::cout<<__FILE__<<" at line: "<<__LINE__<<"invalid stationNo or pixelID ("<<stationNo <<", "<<pixelID<<")" <<endl;

    //correct for time and position
    if(fRedoTimeCorr == 1){
      gigatrackerHit->SetTime(gigatrackerHit->GetRawTime() - fT0[stationNo][pixelID] - fTW[stationNo][chip][totBin] - fTOffsetChip[stationNo][chip] );
    }
    else  gigatrackerHit->SetTime(gigatrackerHit->GetTime() - fTOffsetChip[stationNo][chip]);

    if(fRedoXYCorr == 1){
      gigatrackerHit->SetPosition(gigatrackerHit->GetRawPosition() - fPosOff[stationNo]);
    }
    else  gigatrackerHit->SetPosition(gigatrackerHit->GetPosition() );
    
    deltaT    = gigatrackerHit->GetTime() - triggerTime;
       
    if(fNoisy[stationNo][pixelID] == 1)  fHMapNoise[stationNo]->Fill(pixelID%200, pixelID/200);
    else                                 fHMap[stationNo]->Fill(pixelID%200, pixelID/200);

    // clean out of time hits and hits from noisy pixels
    if( fNoisy[stationNo][pixelID] == 1 ||  (abs(deltaT) > fTimeWindowTrigger && fRemoveHit) ) fGigaTrackerEvent->RemoveHit(iH);
    else {
      iH++;
      fHDt[stationNo]->Fill(deltaT);
      fHDtToT[stationNo]->Fill(deltaT,ToT);
    }
    
  }
*/
//good for seg
  //REMOVE ALL EXISTING CANDIDATE
  int nCand = fGigaTrackerEvent->GetNCandidates();
  for(int iC(0);iC<nCand; iC++) fGigaTrackerEvent->RemoveCandidate(0); //remove candidate


  // BUILD CANDIDATE
  vector<int> vT;
  for(int i(0); i<fGigaTrackerEvent->GetNHits(); i++){
    gigatrackerHit = (TRecoGigaTrackerHit*)fGigaTrackerEvent->GetHit(i);
    if(abs(gigatrackerHit->GetTime() - triggerTime)<=fTimeWindowTrigger) vT.push_back(i);
  }

  TimeOrder    to(this);
  XOrder       xo(this);
  YOrder       yo(this);
  StationOrder so(this);


  // Regroup hit per Block
  // i.e. hits that could come from the same particle
  vector< vector<int> > vBt, vBtx, vBtxy;
  vector< vector<int> >::iterator iBt, iBtx, iBtxy;
  Clusterize(vT, vBt, fTimeWindow, to, 1); //make [time] block

  //for(iBt  = vBt.begin() ; iBt !=vBt.end() ; iBt++ )  Clusterize(*iBt , vBtx , fXWindow   , xo, 1);  //make [time, x] block
  //for(iBtx = vBtx.begin(); iBtx!=vBtx.end(); iBtx++)  Clusterize(*iBtx, vBtxy, fYWindow   , yo, 1);  //make [time, x, y] block

/*  
  //Debug: Plot the hits for the event and the clusters
//  hit_map->Reset();
//  t_block_map->Reset();
//  cluster_map->Reset();

//  hit_map2->Reset();
//  t_block_map2->Reset();
//  cluster_map2->Reset();

  vector<int>::iterator iT;
  for(iT = vT.begin(); iT!=vT.end(); iT++) { 
     gigatrackerHit = (TRecoGigaTrackerHit*)fGigaTrackerEvent->GetHit(*iT);
     int x = gigatrackerHit->GetPixelID()%200;
     int y = gigatrackerHit->GetPixelID()/200;
     double t = gigatrackerHit->GetTime();
     hit_map->Fill(x,y,t+10);
     hit_map2->Fill(x,y,gigatrackerHit->GetStationNo()+1);
  }

  double iClusterTemp(0);
  int iC(1);
  //  for(iBt = vBt.begin(); iBt!=vBt.end(); iBt++) { 
  for(iBt = vBt.begin(); iBt!=vBt.end(); iBt++) { 
    iClusterTemp = iClusterTemp+10;
    for(iT = iBt->begin(); iT!=iBt->end(); iT++) { 
      gigatrackerHit = (TRecoGigaTrackerHit*)fGigaTrackerEvent->GetHit(*iT);
      int x = gigatrackerHit->GetPixelID()%200;
      int y = gigatrackerHit->GetPixelID()/200;
      t_block_map->Fill(x,y,iClusterTemp);
      t_block_map2->Fill(x,y,iC);
    }
    iC++;
  }
  fCv->cd(1);
  gStyle->SetOptStat("e");
  hit_map2->SetTitle("Hit Map | Time + 10 ns | Station No; x [pixel]; y [pixel]");
  hit_map2->GetZaxis()->SetRangeUser(hit_map2->GetBinContent(hit_map2->GetMinimumBin()),hit_map2->GetBinContent(hit_map2->GetMaximumBin()));
  hit_map2->Draw("TEXT");
  hit_map->Draw("COLZ SAME");

  fCv->cd(2);
  t_block_map2->SetTitle("Time Block Map | Time + 10 ns | Block Id; x [pixel]; y [pixel]");
  t_block_map2->GetZaxis()->SetRangeUser(hit_map2->GetBinContent(hit_map2->GetMinimumBin()),hit_map2->GetBinContent(hit_map2->GetMaximumBin()));
  t_block_map2->Draw("TEXT");
  t_block_map->Draw("COLZ SAME");
  //Might need to move this up to above BUILD CANDIDATE (who knows?)

*/
  // Building Candidate from Block
  multimap<int, Cluster > mC;
  vector< vector<int> > aB, vCx, vCxy;
  vector< vector<int> >::iterator iS, iCx, iCxy;
  vector<int> ::iterator iHit;
  for(iBtxy = vBt.begin(); iBtxy!=vBt.end(); iBtxy++) { // { h1,h0,h0,h2... }
      //for(iBtxy = vBtxy.begin(); iBtxy!=vBtxy.end(); iBtxy++) { // { h1,h0,h0,h2... }
    mC.clear();

    // Clustering Hits
    aB.clear();
    // Split block in station{ (h0,h0..), (h1,h1..), (h2,h2..) }
    Clusterize(*iBtxy, aB, 0, so, 0);
    for(iS = aB.begin(); iS!=aB.end(); iS++) { //iS contains hit from one station (h0,h0...)

      // Split iS in space (x then y) clusters
      vCx.clear();      vCxy.clear();
      Clusterize(*iS , vCx , 0.410   , xo, 0);
      for(iCx = vCx.begin(); iCx!=vCx.end(); iCx++)  Clusterize(*iCx, vCxy,0.410   , yo, 0);

      // Merge hits in Cluster
      for(iCxy = vCxy.begin(); iCxy!=vCxy.end(); iCxy++) {
	Cluster cluster(this);
	for(iHit= iCxy->begin(); iHit!=iCxy->end(); iHit++){
	  cluster.add(*iHit);
	}
	mC.insert( pair<const int,Cluster>(cluster.S, cluster) );
      }
    }


    //Building Candidate
    // Four Cases: GTK123, GTK13, GTK12, GTK23
//    int NC0 = mC.count(0);
    int NC1 = mC.count(1);
    int NC2 = mC.count(2);
    int NC0 = 0; 
//int NC1 = 0; 
//    int NC2 = 0; 
    int cType = -1;
    if(NC0>0  && NC1>0  && NC2>0)  cType = 123;
    if(NC0>0  && NC1>0  && NC2==0) cType =  12;
    if(NC0>0  && NC1==0 && NC2>0)  cType =  13;
    if(NC0==0 && NC1>0  && NC2>0)  cType =  23;
    TRecoGigaTrackerCandidate* newCand;
cout<<cType<<endl;
if(cType!=23) {
SetOutputState("GigaTrackerEvent", kOInvalid);
return;
}
 
//	fEventNumber = cType;
//SetOutputState("EventNumber", kOInvalid); 
//	EventType->Fill(cType); 
    pair <multimap<int, Cluster>::iterator, multimap<int, Cluster>::iterator> ii0,ii1,ii2;
    multimap<int,Cluster>::iterator i0,i1,i2;
    switch (cType) {
      
    case 123: //############## GTK 123 ##############

      ii0=mC.equal_range(0);	ii1=mC.equal_range(1);	ii2=mC.equal_range(2);
      for(i0=ii0.first; i0!=ii0.second; ++i0){
	for(i1=ii1.first; i1!=ii1.second; ++i1){
	  for(i2=ii2.first; i2!=ii2.second; ++i2){
	    Cluster cs[3] = {i0->second, i1->second, i2->second};
	    //cout<<"Ncand: "<<fGigaTrackerEvent->GetNCandidates();
	    newCand = (TRecoGigaTrackerCandidate*)fGigaTrackerEvent->AddCandidate();
	    //cout<<" -> "<<fGigaTrackerEvent->GetNCandidates()<<endl;

	    for(int s(0); s<3; s++){
	      TVector3 pos(cs[s].X,cs[s].Y,cs[s].Z);
	      newCand->SetType(cType);
	      newCand->SetPosition(s,pos);
	      newCand->SetTimeStation(s,cs[s].T - eventTime);
	      for(iHit = cs[s].hits.begin();iHit != cs[s].hits.end(); iHit++)  newCand->AddHit(*iHit);
	    }
	    BuildCandidate(newCand);

	    //tests a la giuseppe
	    double pkaon = newCand->GetMomentum().Mag();
	    double dxdz = newCand->GetMomentum().X()/newCand->GetMomentum().Z();
	    double dydz = newCand->GetMomentum().Y()/newCand->GetMomentum().Z();
	    //double chi2X = newCand->GetChi2X();
	    //double chi2Y = newCand->GetChi2Y();
	    //double chi2Time = newCand->GetChi2Time();

	    if (pkaon<72000 || pkaon>78000 ||
		dxdz>0.0016 || dxdz<0.0009||
		dydz>0.0004 || dydz<-0.0003){
	      //		chi2X>20||chi2Y>20|| chi2Time>30){
	      fGigaTrackerEvent->RemoveCandidate(fGigaTrackerEvent->GetNCandidates()-1);
	    }
	    else{
	      FillHisto(newCand);
	    }
	  }
	}
      }
      break;


       
    case 12: //############## GTK 12  ##############
      ii0=mC.equal_range(0);	ii1=mC.equal_range(1);
      for(i0=ii0.first; i0!=ii0.second; ++i0){
	for(i1=ii1.first; i1!=ii1.second; ++i1){
	  //derived c3 position
	  Cluster mC(this);
	  mC.S = 2;
	  mC.X = i0->second.X + fDelta13 * (i1->second.X-i0->second.X) / fDelta12 + fShiftTrim;
	  mC.Y = i0->second.Y + fDelta13 * (i1->second.Y-i0->second.Y) / fDelta12;
	  mC.T = 0.5*(i0->second.T + i1->second.T);
	  Cluster cs[3] ={ i0->second, i1->second, mC};

	  newCand = (TRecoGigaTrackerCandidate*)fGigaTrackerEvent->AddCandidate();
	  for(int s(0); s<3; s++){
	    TVector3 pos(cs[s].X,cs[s].Y,cs[s].Z);
	    newCand->SetType(cType);
	    newCand->SetPosition(s,pos);
	    newCand->SetTimeStation(s,cs[s].T - eventTime);
	    for(iHit = cs[0].hits.begin();iHit != cs[0].hits.end(); iHit++)  newCand->AddHit(*iHit);
	  }
	  BuildCandidate(newCand);
	  FillHisto(newCand);
	}
      }
      break;


    case 13:  //############## GTK 13  ##############
      ii0=mC.equal_range(0);	ii2=mC.equal_range(2);
      for(i0=ii0.first; i0!=ii0.second; ++i0){
	for(i2=ii2.first; i2!=ii2.second; ++i2){
	  //derived c2 position
	  Cluster mC(this);
	  mC.S = 1;
	  mC.X = i0->second.X + fDelta12 * (i2->second.X - fShiftTrim - i0->second.X) / fDelta13 ;
	  mC.Y = i0->second.Y + fDelta12 * (i2->second.Y-i0->second.Y) / fDelta13;
	  mC.T = 0.5*(i0->second.T + i2->second.T);
	  Cluster cs[3] ={ i0->second, mC, i2->second};

	  newCand = (TRecoGigaTrackerCandidate*)fGigaTrackerEvent->AddCandidate();
	  for(int s(0); s<3; s++){
	    TVector3 pos(cs[s].X,cs[s].Y,cs[s].Z);
	    newCand->SetType(cType);
	    newCand->SetPosition(s,pos);
	    newCand->SetTimeStation(s,cs[s].T  - eventTime);
	    for(iHit = cs[0].hits.begin();iHit != cs[0].hits.end(); iHit++)  newCand->AddHit(*iHit);
	  }
	  BuildCandidate(newCand);
	  FillHisto(newCand);
	}
      }
      break;

    
    case 23:  //############## GTK 23  ##############
      ii1=mC.equal_range(1);	ii2=mC.equal_range(2);
      for(i1=ii1.first; i1!=ii1.second; ++i1){
	for(i2=ii2.first; i2!=ii2.second; ++i2){
	  //derived c1 position
	  Cluster mC(this);
	  mC.S = 0;
	  mC.X = i1->second.X - fDelta12 * (i2->second.X - fShiftTrim - i1->second.X) / fDelta23 ;
	  mC.Y = i1->second.Y - fDelta12 * (i2->second.Y-i1->second.Y) / fDelta23;
	  mC.T = 0.5*(i1->second.T + i2->second.T);
	  Cluster cs[3] ={mC, i1->second, i2->second};

	  newCand = (TRecoGigaTrackerCandidate*)fGigaTrackerEvent->AddCandidate();
	  for(int s(0); s<3; s++){
	    TVector3 pos(cs[s].X,cs[s].Y,cs[s].Z);
	    newCand->SetType(cType);
	    newCand->SetPosition(s,pos);
	    newCand->SetTimeStation(s,cs[s].T - eventTime);
	    for(iHit = cs[0].hits.begin();iHit != cs[0].hits.end(); iHit++)  newCand->AddHit(*iHit);
	  }
	  BuildCandidate(newCand);
	  FillHisto(newCand);
	}
      }
      break;
   

    
    }
    mC.clear();
  }

    /*
    //debug print clusters...
    multimap<int, Cluster >::iterator iC;
    for(iC = mC.begin(); iC!=mC.end(); iC ++){
      int x(0),y(0);
      for(vector<int>::iterator iH = iC->second.hits.begin(); iH != iC->second.hits.end(); iH++){
	gigatrackerHit = (TRecoGigaTrackerHit*)fGigaTrackerEvent->GetHit(*iH);
	x = gigatrackerHit->GetPixelID()%200;
	y = gigatrackerHit->GetPixelID()/200;
	cluster_map->Fill(x, y,std::distance(mC.begin(), iC)+1);
      }
      cluster_map2->Fill(x, y, std::distance(vBt.begin(), iBtxy)+1);
    }

  }
  fCv->cd(4);
  cluster_map2->SetTitle("Cluster Map | Time + 10 ns | Cluster Id; x [pixel]; y [pixel]");
  cluster_map2->Draw("TEXT");
  cluster_map->Draw("COLZ SAME");


  if(vT.size() >3 ){
    fNPrints++;
    if( fNPrints<10000 && fNPrints%100 == 1 ) fCv->Print("DebugCluster.pdf");
  }
    */

}

void GigaTrackerEvtReco2::PostProcess(){}

void GigaTrackerEvtReco2::EndOfBurstUser(){}

void GigaTrackerEvtReco2::EndOfRunUser(){}

void GigaTrackerEvtReco2::EndOfJobUser(){
  // fCv->Print("DebugCluster.pdf]");
  // SaveAllPlots();
}

void GigaTrackerEvtReco2::DrawPlot(){}

template<typename Order>
void GigaTrackerEvtReco2::Clusterize(vector< int > v, vector< vector<int> >& clusters, double minDist, Order order, int minCont){
  if(v.size()==0){
    return;
  }

  std::sort(v.begin(),v.end(),order); //sort first mpt 13/12/2016

  vector<int> aCluster;
  vector< int >::iterator iV = v.begin(); //sort first mpt 13/12/2016

  while(iV!=v.end()){
    //seed
    aCluster.clear();
    int seed = (*iV);
    //clusterize
    while(iV!=v.end() && order.dist(seed,(*iV))<=minDist ) {
      aCluster.push_back(*iV);
      seed = (*iV); // mpt 13/12/2016
      iV++;
    }
    //store clusters
    if(aCluster.size()>(UInt_t)minCont) clusters.push_back(aCluster);
    aCluster.clear();
  }
  return;
}

//-------------------------------
void GigaTrackerEvtReco2::BuildCandidate(TRecoGigaTrackerCandidate* cand){
  Double_t X[3];
  Double_t Xshift[3];
  Double_t Y[3];
  Double_t T[3];

  for(Int_t iStation=0; iStation < 3 ; iStation++){
    X[iStation] =  cand->GetPosition(iStation).X();
    Xshift[iStation] = X[iStation];
    Y[iStation] =  cand->GetPosition(iStation).Y();
    T[iStation] =  cand->GetTimeStation(iStation); //ns
  }
  Xshift[2] -= fShiftTrim ;   // Correct the TRIM5 effect

  // Compute kinematics and time
  Double_t p = fBeta / (Y[0] * (1.0 - fAlpha) - Y[1] + fDeltaZ + (fAlpha * Y[2]));
  Double_t dydz = (Y[2] - Y[0]) / fDelta13;
  Double_t dxdz = (X[2] - X[0]) / fDelta13 - (((1e-3*fClight * fBLtrim) / p) * (1. - (fDeltaTrim / fDelta13)));
  Double_t pz = p / TMath::Sqrt(1. + dxdz * dxdz + dydz * dydz);
  TVector3 Momentum;
  Momentum.SetXYZ(1e-6*pz * dxdz, 1e-6*pz * dydz, 1e-6*pz);
  Double_t Time = (T[0] + T[1] + T[2]) / 3.0; //ns

  // Track fitting to straight lines for horizontal view (X), with specific offset for GTK3 (trim effect corrected)
  Double_t a,b,rho,chi2X;
  Double_t sigma[3] = {0.0866,0.0866,0.220};
  Double_t z[3]     = {0,fDelta12,fDelta13};
  LinearLeastSquareFit(z,Xshift,3,sigma,a,b,rho,chi2X);

  // Constraints on vertical view (Y)
  Double_t sigmaY12 = 1.42 ;
  Double_t sigmaY23 = 1.20 ;
  Double_t chi2Y = TMath::Power((Y[1] - Y[0])/sigmaY12 , 2.) + TMath::Power((Y[2] - Y[1])/sigmaY23 , 2.) ;

  // Constraints on relative cluster times
  Double_t sigmaT = 0.250;
  Double_t chi2Time = TMath::Power((T[1] - T[0]  )/sigmaT , 2.) + TMath::Power((T[2] - T[1] )/sigmaT , 2.) ;

  // Global Chi2
  Double_t chi2 = chi2X + chi2Y + chi2Time ;

  // Candidate
  cand->SetMomentum(Momentum);
  cand->SetTime(Time);
  cand->SetChi2X(chi2X);
  cand->SetChi2Y(chi2Y);
  cand->SetChi2Time(chi2Time);
  cand->SetChi2(chi2);
}

//-------------------------------
void GigaTrackerEvtReco2::LinearLeastSquareFit(Double_t *x, Double_t *y, Int_t Nsample, Double_t *sigma, Double_t &a, Double_t &b, Double_t &rho, Double_t &chi2){
    // least square method applied to straight line (Y = a + b*X)
    Double_t xmean = TMath::Mean(Nsample,x);
    Double_t ymean = TMath::Mean(Nsample,y);
    Double_t varx = 0.;
    Double_t vary = 0.;
    Double_t covxy = 0.;

    for(Int_t i=0; i < Nsample ; i++){
      varx += (x[i] - xmean) * (x[i] - xmean) ;
      vary += (y[i] - ymean) * (y[i] - ymean) ;
      covxy += (x[i] - xmean) * (y[i] - ymean) ;
    }
    varx = varx / (Double_t)Nsample;
    vary = vary / (Double_t)Nsample;
    covxy = covxy / (Double_t)Nsample;
    b = covxy / varx;
    a = ymean - b * xmean;
    rho = covxy / (TMath::Sqrt(varx) * TMath::Sqrt(vary));

    chi2 = 0;
    for(Int_t i=0; i < Nsample ; i++){
      chi2 += ((y[i] - a - b*x[i])/sigma[i]) * ((y[i] - a - b*x[i])/sigma[i]) ;
    }

    return;
}

//-------------------------------
void GigaTrackerEvtReco2::FillHisto(TRecoGigaTrackerCandidate* cand){

  int type = cand->GetType();
  int iH = 0;
  if (type == 12) iH = 1;
  if (type == 13) iH = 2;
  if (type == 23) iH = 3;

  if(type == 123){
    fHChi2[iH]->Fill(cand->GetChi2());
    fHChi2X[iH]->Fill(cand->GetChi2X());
    fHChi2Y[iH]->Fill(cand->GetChi2Y());
    fHChi2T[iH]->Fill(cand->GetChi2Time());
    //if(cand->GetChi2X()>fChi2X) return;
    //if(cand->GetChi2Time()>fChi2T) return;
  }

  fHMomentum[iH]->Fill(1e-3*cand->GetMomentum().Mag());
  fHThX[iH]->Fill(1e3*cand->GetMomentum().X()/cand->GetMomentum().Z() - 1.2);
  fHThY[iH]->Fill(1e3*cand->GetMomentum().Y()/cand->GetMomentum().Z());


  fHDX_12[iH]->Fill(cand->GetPosition(0).X()-cand->GetPosition(1).X());
  fHDX_13[iH]->Fill(cand->GetPosition(0).X()-cand->GetPosition(2).X());
  fHDX_23[iH]->Fill(cand->GetPosition(1).X()-cand->GetPosition(2).X());

  fHDY_12[iH]->Fill(cand->GetPosition(0).Y()-cand->GetPosition(1).Y());
  fHDY_13[iH]->Fill(cand->GetPosition(0).Y()-cand->GetPosition(2).Y());
  fHDY_23[iH]->Fill(cand->GetPosition(1).Y()-cand->GetPosition(2).Y());

  fHDt_12[iH]->Fill(cand->GetTimeStation(0)-cand->GetTimeStation(1));
  fHDt_13[iH]->Fill(cand->GetTimeStation(0)-cand->GetTimeStation(2));
  fHDt_23[iH]->Fill(cand->GetTimeStation(1)-cand->GetTimeStation(2));

  fHDX_X_12[iH]->Fill(cand->GetPosition(0).X(), cand->GetPosition(0).X()-cand->GetPosition(1).X());
  fHDX_X_13[iH]->Fill(cand->GetPosition(0).X(), cand->GetPosition(0).X()-cand->GetPosition(2).X());
  fHDX_X_23[iH]->Fill(cand->GetPosition(1).X(), cand->GetPosition(1).X()-cand->GetPosition(2).X());

  fHDY_Y_12[iH]->Fill(cand->GetPosition(0).Y(), cand->GetPosition(0).Y()-cand->GetPosition(1).Y());
  fHDY_Y_13[iH]->Fill(cand->GetPosition(0).Y(), cand->GetPosition(0).Y()-cand->GetPosition(2).Y());
  fHDY_Y_23[iH]->Fill(cand->GetPosition(1).Y(), cand->GetPosition(1).Y()-cand->GetPosition(2).Y());

  return;

}
