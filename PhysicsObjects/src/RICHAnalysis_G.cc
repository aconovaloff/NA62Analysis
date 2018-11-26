#include "RICHAnalysis_G.hh"
#include "AnalysisTools.hh"
#include "AnalysisTools_G.hh"
#include "Parameters_R.hh"
#include "RICHImprovedRing.hh"
using namespace TMath; 
//changes marked by "change" search to find where
RICHAnalysis_G::RICHAnalysis_G(NA62Analysis::Core::BaseAnalysis *ba) { 
  Parameters_R *par = Parameters_R::GetInstance(); 
//  par->LoadParameters_R("/afs/cern.ch/user/r/ruggierg/workspace/public/run2015/database/richmirroroff.dat");
//  par->LoadParameters_R("/afs/cern.ch/user/r/ruggierg/workspace/public/run2015/database/richmirrorpos.dat");
par->LoadParameters_R("/afs/cern.ch/user/r/ruggierg/workspace/public/na62git/database/richmirroroff_2015.dat");
    par->LoadParameters_R("/afs/cern.ch/user/r/ruggierg/workspace/public/na62git/database/richmirrorpos_2015.dat");
  par->StoreRICHParameters_R();
  fDeltaMisa = (Double_t **)par->GetRICHDeltaMisa(); 
  fMirrorPos = (Double_t ***)par->GetRICHPos();

  fTools = AnalysisTools_G::GetInstance();
//change
  fUserMethods = new UserMethods(ba);
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_multiring_xy","",400,-100,100,400,-100,100));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_multiring_modvst","",200,-100,100,400,0,400));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_multiring_modvstchod","",400,-20,20,400,0,200));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_multiring_dtimevschi2rich","",200,0,50,200,-100,100));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_multiring_distvschi2rich","",200,0,50,200,0,400));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_multiring_minxy","",400,-100,100,400,-100,100));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_multiring_minmodvst","",200,-100,100,400,0,400));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_multiring_minmodvstchod","",400,-20,20,400,0,200));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_multiring_tchodvsxpos","",128,-1280,1280,200,-20,20));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_multiring_tchodvsypos","",128,-1280,1280,200,-20,20));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_multiring_radvsp","",200,0.,100.,300,0,300));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_singlering_xy","",400,-100,100,400,-100,100));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_singlering_modvst","",200,-100,100,400,0,400));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_singlering_dtimevschi2rich","",200,0,50,200,-100,100));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_singlering_distvschi2rich","",200,0,50,200,0,400));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_singlering_minxy","",400,-100,100,400,-100,100));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_singlering_minmodvst","",200,-100,100,400,0,400));
  fUserMethods->BookHisto(new TH2F("RICHAnalysis_G_singlering_radvsp","",200,0.,100.,300,0,300));
fUserMethods->BookHisto(new TH1I("RICHAnalysis_G_singlering_dist","",200,0,100));
fUserMethods->BookHisto(new TH1I("RICHAnalysis_G_singlering_time","",160,-20,20));
fUserMethods->BookHisto(new TH1I("RICHAnalysis_G_multiring_dist","",200,0,100));
fUserMethods->BookHisto(new TH1I("RICHAnalysis_G_multiring_time","",160,-20,20));
  fRICHMultiCandidate = new RICHCandidate();
  fRICHSingleCandidate = new RICHCandidate();
  fImprovedRing = new RICHImprovedRing();
}

RICHAnalysis_G::~RICHAnalysis_G()
{ }

void RICHAnalysis_G::Clear() {
  fRICHMultiCandidate->Clear();
  fRICHSingleCandidate->Clear();
}
//change, added pion momentum and used own funtion for position at mirror
Int_t RICHAnalysis_G::TrackMultiRingMatching(Double_t chodTime, TRecoRICHEvent* event, TRecoSpectrometerCandidate *thisTrack, TVector3 PionMomentum, TVector3 DecayVertex) {
  fRICHEvent = event;

  // Project the track onto the focal plane
 TVector3 posTrackAtMirror = fTools->GetPositionAtZ(thisTrack,236873.);
//  TVector3 posTrackAtMirror = prop(DecayVertex, PionMomentum*1000, 236873);//change
  Int_t mirrorid = MirrorSurface(posTrackAtMirror.X(),posTrackAtMirror.Y(),0.,0); // edge effect on mirrors (default: no) ?
  TVector3 deltaVers = (posTrackAtMirror-TVector3(0.,0.,202873.)).Unit(); 
  Double_t partrack[4] = {thisTrack->GetMomentum(),thisTrack->GetSlopeXAfterMagnet(),thisTrack->GetSlopeYAfterMagnet(),0.13957018};
  TVector3 versorTrack = fTools->Get4Momentum(partrack).Vect().Unit();
//  TVector3 versorTrack = PionMomentum.Unit(); //change
  MirrorVector(&versorTrack,&deltaVers);
  TVector3 pmtPos = posTrackAtMirror+(-17000./versorTrack.Z())*versorTrack;
  if (mirrorid==99) return -1; // No ring candidate if track is not hitting a mirror (no protection against edge effect) 
  if (mirrorid==2) return -1; // No ring candidate if track is hitting semi-hexagonal Jura

  // Associate tracks with a good ring
  Double_t ttime = thisTrack->GetTime();
  Double_t mindist = 9999999.;
  Double_t mintime = 9999999.;
  Double_t mintimechod = 9999999.;
  Int_t minidring = -1;
  Double_t minchi2 = 9999999.;
  TVector2 mindeltapos(-9999999.,-9999999.); 
  Int_t minNhits =-1;
  // Loop on multi rings (apply correction if they are not applied at reconstruction level) 
//   Double_t corrx = posTrackAtMirror.X()<0 ? +6.0 : -7.67;
//   Double_t corry = posTrackAtMirror.X()<0 ? +10.32 : -11.1;

  Double_t corrx = 0;
  Double_t corry = 0;
// corrx = posTrackAtMirror.X()<0 ? +8.09 : -7.00;
//    corry = posTrackAtMirror.X()<0 ? +10.71 : -11.92;
//    cout <<" pmtPos.X() "<<  pmtPos.X()  <<"  fDeltaMisa[mirrorid][0] "<<  fDeltaMisa[mirrorid][0]    <<" corrx   "<<  corrx<<endl;
//   cout <<" pmtPos.Y() "<<  pmtPos.Y()  <<"  fDeltaMisa[mirrorid][1] "<<  fDeltaMisa[mirrorid][1]    <<" corry   "<<  corry<<endl;
  pmtPos.SetX(pmtPos.X()+fDeltaMisa[mirrorid][0]+corrx);
  pmtPos.SetY(pmtPos.Y()+fDeltaMisa[mirrorid][1]+corry);
  for (Int_t jring=0; jring<fRICHEvent->GetNRingCandidates(); jring++) {
    TRecoRICHCandidate *ring = (TRecoRICHCandidate *)fRICHEvent->GetRingCandidate(jring);
    ring->SetEvent(fRICHEvent);
    if (ring->GetNHits()<4) continue;
    if (ring->GetRingChi2()>=10) continue; 
    TVector2 ringcenter = ring->GetRingCenter();
//      cout <<" ringcenter.X() "<<  ringcenter.X() <<"  pmtPos.X() "<<pmtPos.X() << endl;
//     cout <<" ringcenter.Y() "<<  ringcenter.Y() <<"  pmtPos.Y() "<<pmtPos.Y() << endl;
    TVector2 deltapos(ringcenter.X()-pmtPos.X(),ringcenter.Y()-pmtPos.Y());
//cout<<"center x: "<<ringcenter.X()<<endl;
//cout<<"center y: "<<ringcenter.Y()<<endl;
//cout<<"pmt x: "<<pmtPos.X()<<endl;
//cout<<"pmt y: "<<pmtPos.Y()<<endl;
//cout<<"total x: "<<ringcenter.X()-pmtPos.X()<<endl;
//cout<<"total y: "<<ringcenter.Y()-pmtPos.Y()<<endl;

    Double_t dist = deltapos.Mod();
    Double_t dtime = ring->GetTime()-0.6-ttime;
    Double_t dtimechod = ring->GetTime()-0.6-chodTime;

    //  cout <<" dist "<<   dist <<"   dtime "<< dtime   <<"  dtimechod "<< dtimechod  << endl;
    fUserMethods->FillHisto("RICHAnalysis_G_multiring_xy",deltapos.X(),deltapos.Y());
    fUserMethods->FillHisto("RICHAnalysis_G_multiring_modvst",dtime,dist);
    fUserMethods->FillHisto("RICHAnalysis_G_multiring_modvstchod",dtimechod,dist);
    Double_t chi2rich = dtime*dtime/(2*6*6)+dist*dist/(4*1.6*1.6); //how it was from roberta
//Double_t chi2rich = dtime*dtime/(2*6*6)+dist*dist/(4*16*16);//this makes it mm perhaps?
//cout<<"time: "<<dtime<<endl;
//cout<<"dist: "<<dist<<endl; 
    fUserMethods->FillHisto("RICHAnalysis_G_multiring_dtimevschi2rich",chi2rich,dtime);
    fUserMethods->FillHisto("RICHAnalysis_G_multiring_distvschi2rich",chi2rich,dist);
    fUserMethods->FillHisto("RICHAnalysis_G_multiring_dist",dist);
    fUserMethods->FillHisto("RICHAnalysis_G_multiring_time",dtime);
    if (chi2rich<minchi2) {
      minchi2 = chi2rich;
      mindist = dist;
      mintime = dtime;
      mintimechod = dtimechod;
      minidring = jring;
      mindeltapos = deltapos;
      minNhits = ring->GetNHits();
    } 
  }
  DistCheck = mindist; 
  if (minidring>-1) {
    fUserMethods->FillHisto("RICHAnalysis_G_multiring_minxy",mindeltapos.X(),mindeltapos.Y());
    fUserMethods->FillHisto("RICHAnalysis_G_multiring_minmodvst",mintime,mindist);
    fUserMethods->FillHisto("RICHAnalysis_G_multiring_minmodvstchod",mintimechod,mindist);
    fUserMethods->FillHisto("RICHAnalysis_G_multiring_tchodvsxpos",posTrackAtMirror.X(),mintimechod);
    fUserMethods->FillHisto("RICHAnalysis_G_multiring_tchodvsypos",posTrackAtMirror.Y(),mintimechod);
    TRecoRICHCandidate *minring = (TRecoRICHCandidate *)fRICHEvent->GetRingCandidate(minidring);
    if (mindist<20 && fabs(mintime)<30) fUserMethods->FillHisto("RICHAnalysis_G_multiring_radvsp",thisTrack->GetMomentum()/1000,minring->GetRingRadius()); 
    fRICHMultiCandidate->SetDiscriminant(minchi2);
//cout<<minchi2<<endl;
    fRICHMultiCandidate->SetDX(mindeltapos.X());
    fRICHMultiCandidate->SetDY(mindeltapos.Y());
    TVector2 rcenter = minring->GetRingCenter();
    fRICHMultiCandidate->SetXcenter(rcenter.X());
    fRICHMultiCandidate->SetYcenter(rcenter.Y());
    fRICHMultiCandidate->SetRadius(minring->GetRingRadius());
    fRICHMultiCandidate->SetDeltaTime(mintime);
    fRICHMultiCandidate->SetChi2(minring->GetRingChi2());
    fRICHMultiCandidate->SetMirrorID(mirrorid);
    // cout <<"inside multiring:  radius "<<minring->GetRingRadius() 
    //	 <<"  minNhits  "<<  minNhits 
    //	 <<" rcenter.X()  "<<  rcenter.X()  <<endl;
    fRICHMultiCandidate->SetNHits(minNhits);
    //  fRICHMultiCandidate->SetRingCandidateID(minidring);
  }
  return minidring;
}
//change, added pionmomentum and decayvertz as arguments
Int_t RICHAnalysis_G::TrackSingleRingMatching(Double_t chodTime, TRecoRICHEvent* event, TRecoSpectrometerCandidate *thisTrack, TVector3 PionMomentum, TVector3 DecayVertex) {
  fRICHEvent = event;

  // Project the track onto the focal plane
  TVector3 posTrackAtMirror = fTools->GetPositionAtZ(thisTrack,236873.);
//TVector3 posTrackAtMirror = prop(DecayVertex, PionMomentum*1000, 236837);//change
  // Int_t mirrorid = MirrorSurface(posTrackAtMirror.X(),posTrackAtMirror.Y(),0.,0); // edge effect on mirrors (default: no) ?
  Int_t mirrorid = MirrorSurface(posTrackAtMirror.X(),posTrackAtMirror.Y(),0.,1); // edge effect on mirrors (default: no) ?


  TVector3 deltaVers = (posTrackAtMirror-TVector3(0.,0.,202873.)).Unit(); 
  Double_t partrack[4] = {thisTrack->GetMomentum(),thisTrack->GetSlopeXAfterMagnet(),thisTrack->GetSlopeYAfterMagnet(),0.13957018};
  TVector3 versorTrack = fTools->Get4Momentum(partrack).Vect().Unit();
//  TVector3 versorTrack = PionMomentum.Unit(); //change
  MirrorVector(&versorTrack,&deltaVers);
  TVector3 pmtPos = posTrackAtMirror+(-17000./versorTrack.Z())*versorTrack;
  if (mirrorid==99) return -1; // No ring candidate if track is not hitting a mirror (no protection against edge effect) 
  if (mirrorid==2) return -1; // No ring candidate if track is hitting semi-hexagonal Jura

  // Associate tracks with a good ring
  Double_t ttime = thisTrack->GetTime();
  Double_t mindist = 9999999.;
  Double_t mintime = 9999999.;
  Int_t minidring = -1;
  Double_t minchi2 = 9999999.;
  TVector2 mindeltapos(-9999999.,-9999999.); 
  TVector2 minrcenter(-9999999.,-9999999.); 
  Double_t minradius = 0.; 
  Double_t minringchi2 = 99999999.;  
  Int_t minNhits =-1;
  // Loop on single rings 
  for (Int_t jring=0; jring<fRICHEvent->GetNTimeCandidates(); jring++) {
    TRecoRICHCandidate *ring = (TRecoRICHCandidate *)fRICHEvent->GetTimeCandidate(jring);
    Double_t goodFit = fImprovedRing->Chi2Fit(fRICHEvent,ring,pmtPos,posTrackAtMirror,fDeltaMisa);
    if (!goodFit) continue;
    TVector2 ringcenter = ring->GetRingCenter();
    TVector2 deltapos(ringcenter.X()-pmtPos.X(),ringcenter.Y()-pmtPos.Y());
    Double_t dist = deltapos.Mod();
    Double_t dtime = ring->GetTime()-0.6-ttime;
    fUserMethods->FillHisto("RICHAnalysis_G_singlering_xy",deltapos.X(),deltapos.Y());
    fUserMethods->FillHisto("RICHAnalysis_G_singlering_modvst",dtime,dist);
    Double_t chi2rich = dtime*dtime/(2*6*6)+dist*dist/(4*1.6*1.6);//how it was from roberta
// Double_t chi2rich = dtime*dtime/(2*6*6)+dist*dist/(4*25*25); //this makes it mm perhaps?? //set to 25 (instead 16) to deal with extra distance in multi
    fUserMethods->FillHisto("RICHAnalysis_G_singlering_dtimevschi2rich",chi2rich,dtime);
    fUserMethods->FillHisto("RICHAnalysis_G_singlering_distvschi2rich",chi2rich,dist);
    fUserMethods->FillHisto("RICHAnalysis_G_singlering_dist",dist);
    fUserMethods->FillHisto("RICHAnalysis_G_singlering_time",dtime);
    if (chi2rich<minchi2) {
      minchi2 = chi2rich;
      mindist = dist;
      mintime = dtime;
      mindeltapos = deltapos;
      minrcenter = ring->GetRingCenter();
      minradius = ring->GetRingRadius();
      minringchi2 = ring->GetRingChi2();
      minNhits = ring->GetNHits();
      minidring = jring;
    } 
  }
  if (minidring>-1) {
    fUserMethods->FillHisto("RICHAnalysis_G_singlering_minxy",mindeltapos.X(),mindeltapos.Y());
    fUserMethods->FillHisto("RICHAnalysis_G_singlering_minmodvst",mintime,mindist);
    if (mindist<20 && fabs(mintime)<30) fUserMethods->FillHisto("RICHAnalysis_G_singlering_radvsp",thisTrack->GetMomentum()/1000,minradius); 
    fRICHSingleCandidate->SetDiscriminant(minchi2);
    fRICHSingleCandidate->SetDX(mindeltapos.X());
    fRICHSingleCandidate->SetDY(mindeltapos.Y());
    fRICHSingleCandidate->SetXcenter(minrcenter.X());
    fRICHSingleCandidate->SetYcenter(minrcenter.Y());
    fRICHSingleCandidate->SetRadius(minradius);
    fRICHSingleCandidate->SetDeltaTime(mintime);
    fRICHSingleCandidate->SetChi2(minringchi2);
    fRICHSingleCandidate->SetMirrorID(mirrorid);
    fRICHSingleCandidate->SetNHits(minNhits);
    fRICHSingleCandidate->SetTimeCandInd(minidring);
  }
  return minidring;
}

///////////////////////
// Mirror reflection //
///////////////////////
void RICHAnalysis_G::MirrorVector(TVector3 *vec, TVector3 *axis) {
  TVector3 S = ((vec->Dot(*axis))/(axis->Mag2()))*(*axis);
  TVector3 d = S-*vec;
  TVector3 ret = S+d;
  *vec = ret.Unit();
}

////////////////////////
// Find the mirror ID //
////////////////////////
Int_t RICHAnalysis_G::MirrorSurface(Double_t xin, Double_t yin, Double_t SafeFactor, Bool_t flag) {
  Int_t imirror=0;
  for (imirror=0; imirror<25; imirror++) {
    if (imirror==0 || imirror==7 || imirror==18 || imirror==19 ) continue;  
    if (flag==1 && imirror==2) continue;

    Double_t ac1=(fMirrorPos[imirror][0][1]-fMirrorPos[imirror][1][1])/(fMirrorPos[imirror][0][0]-fMirrorPos[imirror][1][0]);
    Double_t bcost1=(fMirrorPos[imirror][1][1]*fMirrorPos[imirror][0][0]-fMirrorPos[imirror][1][0]*fMirrorPos[imirror][0][1])/(fMirrorPos[imirror][0][0]-fMirrorPos[imirror][1][0]);
    Double_t ac3=(fMirrorPos[imirror][2][1]-fMirrorPos[imirror][3][1])/(fMirrorPos[imirror][2][0]-fMirrorPos[imirror][3][0]);
    Double_t bcost3=(fMirrorPos[imirror][3][1]*fMirrorPos[imirror][2][0]-fMirrorPos[imirror][3][0]*fMirrorPos[imirror][2][1])/(fMirrorPos[imirror][2][0]-fMirrorPos[imirror][3][0]);
    Double_t ac4=(fMirrorPos[imirror][3][1]-fMirrorPos[imirror][4][1])/(fMirrorPos[imirror][3][0]-fMirrorPos[imirror][4][0]);
    Double_t bcost4=(fMirrorPos[imirror][4][1]*fMirrorPos[imirror][3][0]-fMirrorPos[imirror][4][0]*fMirrorPos[imirror][3][1])/(fMirrorPos[imirror][3][0]-fMirrorPos[imirror][4][0]);
    Double_t ac6=(fMirrorPos[imirror][5][1]-fMirrorPos[imirror][0][1])/(fMirrorPos[imirror][5][0]-fMirrorPos[imirror][0][0]);
    Double_t bcost6=(fMirrorPos[imirror][0][1]*fMirrorPos[imirror][5][0]-fMirrorPos[imirror][0][0]*fMirrorPos[imirror][5][1])/(fMirrorPos[imirror][5][0]-fMirrorPos[imirror][0][0]);  
    Double_t ConvSafeFactor1=SafeFactor*sqrt(1.+pow(ac1,2));
    Double_t ConvSafeFactor3=SafeFactor*sqrt(1.+pow(ac3,2));
    Double_t ConvSafeFactor4=SafeFactor*sqrt(1.+pow(ac4,2));
    Double_t ConvSafeFactor6=SafeFactor*sqrt(1.+pow(ac6,2));

    if (imirror!=2) {
      if ((yin<=ac1*xin+bcost1-ConvSafeFactor1) && 
          (xin<=(fMirrorPos[imirror][1][0]+fMirrorPos[imirror][2][0])/2.-SafeFactor) && 
          (yin>=ac3*xin+bcost3+ConvSafeFactor3) && 
          (yin>=ac4*xin+bcost4+ConvSafeFactor4) && 
          (xin>=(fMirrorPos[imirror][4][0]+fMirrorPos[imirror][5][0])/2.+SafeFactor) &&  
          (yin<=ac6*xin+bcost6-ConvSafeFactor6) && (imirror<23)) return imirror;
    } else {
      if ((yin<=ac1*xin+bcost1-ConvSafeFactor1) && 
          (xin<=(fMirrorPos[imirror][1][0]+fMirrorPos[imirror][2][0])/2.-SafeFactor) && 
          (yin>=ac3*xin+bcost3+ConvSafeFactor3) && 
          (yin>=ac4*xin+bcost4+ConvSafeFactor4) && 
          (xin>=0.+SafeFactor) &&  
          (yin<=ac6*xin+bcost6-ConvSafeFactor6) && (imirror<23)) return imirror;
    }

    if (flag==1) {
      if (imirror==23 || imirror==24){
        Double_t acs1=(fMirrorPos[imirror][1][1]-fMirrorPos[imirror][0][1])/(fMirrorPos[imirror][1][0]-fMirrorPos[imirror][0][0]);
        Double_t bcosts1=(fMirrorPos[imirror][0][1]*fMirrorPos[imirror][1][0]-fMirrorPos[imirror][0][0]*fMirrorPos[imirror][1][1])/(fMirrorPos[imirror][1][0]-fMirrorPos[imirror][0][0]);
        Double_t acs2=(fMirrorPos[imirror][3][1]-fMirrorPos[imirror][2][1])/(fMirrorPos[imirror][3][0]-fMirrorPos[imirror][2][0]);
        Double_t bcosts2=(fMirrorPos[imirror][2][1]*fMirrorPos[imirror][3][0]-fMirrorPos[imirror][2][0]*fMirrorPos[imirror][3][1])/(fMirrorPos[imirror][3][0]-fMirrorPos[imirror][2][0]);
        SafeFactor=-4.;
        Double_t ConvSafeFactors1=SafeFactor*sqrt(1.+pow(acs1,2));
        Double_t ConvSafeFactors2=SafeFactor*sqrt(1.+pow(acs2,2));
        if ((yin<=acs1*xin+bcosts1+ConvSafeFactors1) && (xin>=(fMirrorPos[imirror][0][0]+fMirrorPos[imirror][3][0])/2.-SafeFactor) && (yin>=acs2*xin+bcosts2-ConvSafeFactors2) && (xin<=(fMirrorPos[imirror][1][0]+fMirrorPos[imirror][2][0])/2.+SafeFactor) && (imirror==23)) return imirror;
        if ((yin<=acs1*xin+bcosts1+ConvSafeFactors1) && (xin<=(fMirrorPos[imirror][0][0]+fMirrorPos[imirror][3][0])/2.+SafeFactor) && (yin>=acs2*xin+bcosts2-ConvSafeFactors2) && (xin>=(fMirrorPos[imirror][1][0]+fMirrorPos[imirror][2][0])/2.-SafeFactor) && (imirror==24)) return imirror;
      }
    }

  }
  imirror=99;
  return imirror;  
}

//////////try copy/paste OneTrackEventAnalysis RICHMultiCand/RICHSingleRing Here

//Bool_t OneTrackEventAnalysis::RICHMultiCandidate(Int_t jTrack) {
Bool_t RICHAnalysis_G::RICHMultiCandidate(RICHCandidate RICHCand, TRecoSpectrometerCandidate *track, TVector3 PionMomentum, TVector3 DecayVertex) {
//  if (!fRICHMultiCandidate[jTrack].GetIsRICHCandidate()) return 0; // A least a ring associated to the track 
//  FillHisto("Track_multirich_chi2",fRICHMultiCandidate[jTrack].GetDiscriminant());
//  if (fRICHMultiCandidate[jTrack].GetDiscriminant()>50) return 0; // Minimum quality association required 
//cout<<"multi disc: "<<RICHCand.GetDiscriminant()<<endl;
if (RICHCand.GetDiscriminant()>35) return 0; // Minimum quality association required  
//  FillHisto("Track_multirich_ringchi2",fRICHMultiCandidate[jTrack].GetChi2());
//  fDownstreamTrack.SetRICHMultiIsCandidate(1);
//  TRecoSpectrometerCandidate *track = (TRecoSpectrometerCandidate *)fSpectrometerEvent->GetCandidate(jTrack);
 // Double_t richtime = fRICHMultiCandidate[jTrack].GetDeltaTime()+track->GetTime();   
//  Double_t chodtime = fCHODCandidate[jTrack].GetDeltaTime()+track->GetTime(); 
//  Double_t dtimechod = richtime-chodtime;
  Double_t dx = RICHCand.GetDX();
  Double_t dy = RICHCand.GetDY();
  Double_t mindist = sqrt(dx*dx+dy*dy);
//  FillHisto("Track_multirich_distvstime",dtimechod,mindist);
//  TVector3 posTrackAtMirror =ftools->GetPositionAtZ(track,236873.);
TVector3 posTrackAtMirror = prop(DecayVertex,PionMomentum*1000, 236873.);
  Double_t corrx = posTrackAtMirror.X()<0 ? +6.0 : -7.67;
  Double_t corry = posTrackAtMirror.X()<0 ? +10.32 : -11.1;
  Int_t mirrorid = RICHCand.GetMirrorID();
  Double_t slopex = (RICHCand.GetXcenter()-fDeltaMisa[mirrorid][0]-corrx-1.8)/17000.;
  Double_t slopey = (RICHCand.GetYcenter()-fDeltaMisa[mirrorid][1]-corry)/17000.;
  Double_t dslopex = slopex-track->GetSlopeXAfterMagnet();
  Double_t dslopey = slopey-track->GetSlopeYAfterMagnet();
//  FillHisto("Track_multirich_xyslope",dslopex,dslopey);
  Double_t radius = RICHCand.GetRadius();
//cout<<"rad: "<<radius<<endl; 
//  Double_t pmag = fDownstreamTrack.GetMomentum().P();
   Double_t pmag = PionMomentum.Mag();//////not sure if this is correct !!!!!!!!!!!!!!!
//  FillHisto("Track_multirich_radvsp",pmag,radius);
//  fDownstreamTrack.SetRICHMultiTime(richtime);
//  fDownstreamTrack.SetRICHMultiChi2(fRICHMultiCandidate[jTrack].GetChi2());
  Double_t changle = atan(radius/17000.);
  if(1.000062*1.000062*cos(changle)*cos(changle)<1)return 0; 
  multirichmass = pmag*sqrt(1.000062*1.000062*cos(changle)*cos(changle)-1);
//  FillHisto("Track_multirich_massvsp",pmag,richmass);
//  fDownstreamTrack.SetRICHMultiMass(richmass);
  Double_t d2ring = sqrt(fabs(190.2*190.2-radius*radius));   
  Double_t prich = 0.13957018*17000./d2ring;
  Double_t partrack[4] = {prich*1000,slopex+0.270/pmag,slopey,0.13957018};
//  FillHisto("Track_multirich_pvsp",pmag,prich);
//  fDownstreamTrack.SetRICHMultiMomentum(tools->Get4Momentum(partrack));
  return 1;
}

Bool_t RICHAnalysis_G::RICHSingleCandidate(RICHCandidate RICHCand, TRecoSpectrometerCandidate *track, TVector3 PionMomentum, TVector3 DecayVertex) {
//  if (!fRICHSingleCandidate[jTrack].GetIsRICHCandidate()) return 0; // A least a ring associated to the track 
//  FillHisto("Track_singlerich_chi2",fRICHSingleCandidate[jTrack].GetDiscriminant());
//cout<<"single disc: "<<RICHCand.GetDiscriminant()<<endl;
  if (RICHCand.GetDiscriminant()>50) return 0; // Minimum quality association required 
//  FillHisto("Track_singlerich_ringchi2",fRICHSingleCandidate[jTrack].GetChi2());
  if (RICHCand.GetChi2()>3) return 0; // Minimum ring quality required 
//  fDownstreamTrack.SetRICHSingleIsCandidate(1);
//  TRecoSpectrometerCandidate *track = (TRecoSpectrometerCandidate *)fSpectrometerEvent->GetCandidate(jTrack);
//  Double_t richtime = fRICHSingleCandidate[jTrack].GetDeltaTime()+track->GetTime();   
//  Double_t chodtime = fCHODCandidate[jTrack].GetDeltaTime()+track->GetTime(); 
//  Double_t dtimechod = richtime-chodtime;
  Double_t dx = RICHCand.GetDX();
  Double_t dy = RICHCand.GetDY();
  Double_t mindist = sqrt(dx*dx+dy*dy);
//  FillHisto("Track_singlerich_distvstime",dtimechod,mindist);
//  TVector3 posTrackAtMirror = ftools->GetPositionAtZ(track,236873.);
  TVector3 posTrackAtMirror = prop(DecayVertex,PionMomentum*1000, 236873.);
  Int_t mirrorid = RICHCand.GetMirrorID();
  Double_t slopex = (RICHCand.GetXcenter())/17000.;
  Double_t slopey = (RICHCand.GetYcenter())/17000.;
  Double_t dslopex = slopex-track->GetSlopeXAfterMagnet();
  Double_t dslopey = slopey-track->GetSlopeYAfterMagnet();
//  FillHisto("Track_singlerich_xyslope",dslopex,dslopey);
  Double_t radius = RICHCand.GetRadius();
//  Double_t pmag = fDownstreamTrack.GetMomentum().P();
  Double_t pmag = PionMomentum.Mag(); ///////////Don't know what this is!!!!!!!!!!!!!!
//  FillHisto("Track_singlerich_radvsp",pmag,radius);
//  fDownstreamTrack.SetRICHSingleTime(richtime);
//  fDownstreamTrack.SetRICHSingleChi2(fRICHSingleCandidate[jTrack].GetChi2());
  Double_t changle = atan(radius/17000.);
if(1.000062*1.000062*cos(changle)*cos(changle)<1)return 0;
  singlerichmass = pmag*sqrt(1.000062*1.000062*cos(changle)*cos(changle)-1);
//  FillHisto("Track_singlerich_massvsp",pmag,richmass);
//  fDownstreamTrack.SetRICHSingleMass(richmass);
  Double_t d2ring = sqrt(fabs(190.2*190.2-radius*radius));   
  Double_t prich = 0.13957018*17000./d2ring;
  Double_t trackptkick = track->GetSlopeXAfterMagnet()-track->GetSlopeXBeforeMagnet();
//  Double_t partrack[4] = {prich*1000,slopex+0.270/pmag,slopey,0.13957018};
//  Double_t partrack[4] = {prich*1000,slopex-trackptkick,slopey,0.13957018};
//  Double_t partrack[4] = {prich*1000,track->GetSlopeXBeforeMagnet(),track->GetSlopeYBeforeMagnet(),0.13957018};
//  FillHisto("Track_singlerich_pvsp",pmag,prich);
//  fDownstreamTrack.SetRICHSingleMomentum(tools->Get4Momentum(partrack));
  return 1;
          }


void RICHAnalysis_G::SaveAllPlots(){
fUserMethods->SaveAllPlots();
}
