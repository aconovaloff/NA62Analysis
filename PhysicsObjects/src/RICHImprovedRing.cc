#include "RICHImprovedRing.hh"
#include "AnalysisTools_G.hh"
#include "TRecoRICHEvent.hh"
#include "TRecoRICHCandidate.hh"
#include "TRecoRICHHit.hh"
#include "TMath.h"

using namespace NA62Constants;

std::vector<TVector2> hitPosForChiFit;


RICHImprovedRing::RICHImprovedRing() {
  InitMirrorPositions(); 
  fNPars = 3; // Parameters for fit
  fFitter = new TMinuit(fNPars);
  fFitter->SetFCN(RingChi2FCN);
  tools = AnalysisTools_G::GetInstance();
}

RICHImprovedRing::~RICHImprovedRing() {
}


/********************************************//**
 * Minuit minimization function
 ************************************************/
void RICHImprovedRing::RingChi2FCN(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag)
{
  if(1 || iflag == 4) {
    Int_t nHits = hitPosForChiFit.size();
    f = 0.;
    for (Int_t iHit = 0; iHit < nHits; iHit++) {
      TVector2 PMPosition = hitPosForChiFit[iHit];
      Double_t u = PMPosition.X() - par[0];
      Double_t v = PMPosition.Y() - par[1];
      Double_t d = TMath::Sqrt(u*u+v*v);
      Double_t dr = d - par[2];
      f += dr*dr/81.;
    }
  } // iflag=4 for Minuit
}

/********************************************//**
 * Ring Chi2
 ************************************************/
Double_t RICHImprovedRing::RingChi2(Double_t *par)
{
  Int_t nHits = hitPosForChiFit.size();
  Double_t f = 0.;
  for (Int_t iHit = 0; iHit < nHits; iHit++) {
    TVector2 PMPosition = hitPosForChiFit[iHit];
    Double_t u = PMPosition.X() - par[0];
    Double_t v = PMPosition.Y() - par[1];
    Double_t d = TMath::Sqrt(u*u+v*v);
    Double_t dr = par[2] - d;
    f += dr*dr/81.;
  }
  return f/(nHits-3);
}

/********************************************//**
 * Chi2 Fit
 ************************************************/
Bool_t RICHImprovedRing::Chi2Fit(TRecoRICHEvent* event, TRecoRICHCandidate *thisRing, TVector3 pmtpos, TVector3 posTrackAtMirror, Double_t **deltamirror) {
  Float_t Xcog = 0;
  Float_t Ycog = 0;
  Float_t Rmean = 0;
  fDeltaMisa = deltamirror;

  // Extract the hits belonging to the ring candidate and compute the CoG
  hitPosForChiFit.clear();
  Int_t nHits = thisRing->GetNHits();
  TClonesArray& Hits = (*(event->GetHits()));
  Int_t *hitIndex = (Int_t *)thisRing->GetHitsIndexes();     
  Int_t goodNHits = 0;
  for (Int_t jHit=0; jHit<nHits; jHit++) {
    TRecoRICHHit *richhit = (TRecoRICHHit*)Hits[hitIndex[jHit]];
    Bool_t corrected = 0;
    corrected = CorrectHitPosition(richhit,pmtpos,posTrackAtMirror); 
    Double_t xpos = richhit->GetFitPosition().X();
    Double_t ypos = richhit->GetFitPosition().Y();
    if (corrected) {
      xpos = richhit->GetPosition().X()-GetChPosAngCorr(richhit->GetChannelSeqID()).X();
      ypos = richhit->GetPosition().Y()-GetChPosAngCorr(richhit->GetChannelSeqID()).Y();
    }
    Double_t hdist2 = (xpos-pmtpos.X())*(xpos-pmtpos.X())+(ypos-pmtpos.Y())*(ypos-pmtpos.Y());
    if (hdist2>250.*250.) continue;
    if (hdist2<80.*80.) continue;
    if (richhit->GetOrSuperCellID()==1) continue;
    hitPosForChiFit.push_back(TVector2(xpos,ypos)); 
    Xcog += hitPosForChiFit[goodNHits].X();
    Ycog += hitPosForChiFit[goodNHits].Y();
    goodNHits++;
  }
  Xcog /= goodNHits;
  Ycog /= goodNHits;
  for (Int_t iHit=0; iHit<goodNHits; iHit++) {
    Rmean += sqrt(pow(Xcog-hitPosForChiFit[iHit].X(),2)+pow(Ycog-hitPosForChiFit[iHit].Y(),2));
  }
  Rmean /= goodNHits;
  if (goodNHits<4) return 0;

  // Perform the fit
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat,ierflag;
  Double_t pars[10],epars[10];
  Double_t arglist[1];
  arglist[0] = -1;
  fFitter->mnexcm("SET PRI", arglist, 1, ierflag); // Set MINUIT print level
  fFitter->mnexcm("SET NOW", arglist, 0, ierflag); // Set MINUIT warnings
  fFitter->mnparm(0, "x0", Xcog, 0.01, 0., 0., ierflag);
  fFitter->mnparm(1, "y0", Ycog, 0.01, 0., 0., ierflag);
  fFitter->mnparm(2, "R", Rmean, 0.01, 0., 0., ierflag);
  fFitter->mnparm(0, "x0", Xcog, 0.01, 0., 0., ierflag);
  fFitter->mnparm(1, "y0", Ycog, 0.01, 0., 0., ierflag);
  fFitter->mnparm(2, "R", Rmean, 0.01, 0., 0., ierflag);
  arglist[0] = 0;
  fFitter->mnexcm("MIGRAD", arglist, 1, ierflag); // Calls the minimization
  fFitter->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  for(Int_t iPar = 0; iPar < fNPars; iPar++) fFitter->GetParameter(iPar,pars[iPar],epars[iPar]);
  TVector2 CurrentRingCenter;
  Double_t CurrentRingRadius, CurrentRingChi2;

  CurrentRingCenter.Set(pars[0],pars[1]);
  CurrentRingRadius = pars[2];
  CurrentRingChi2 = RingChi2(pars);

  Double_t thetamax=-TMath::Pi();
  Double_t thetamin=TMath::Pi();
  for (Int_t iHit=0; iHit<hitPosForChiFit.size(); iHit++) {
    Double_t theta=atan2(hitPosForChiFit[iHit].X()-pars[0],hitPosForChiFit[iHit].Y()-pars[1]);
    if (thetamax<theta) thetamax=theta;
    if (thetamin>theta) thetamin=theta;
  }
  
  //cout <<" thetamax "<< thetamax   <<   "   thetamin "<<   thetamin   <<"  thetamax-thetamin "<<  (thetamax-thetamin)*180/TMath::Pi()  <<endl;
  
  
  if ((thetamax-thetamin)<185.*TMath::Pi()/180.) CurrentRingChi2=9.;
  

  // thisRing->SetHitsMaxAngle(thetamax-thetamin);
  thisRing->SetRingCenter(CurrentRingCenter);
  thisRing->SetRingRadius(CurrentRingRadius);
  thisRing->SetRingChi2(CurrentRingChi2);

  // Printout
//  cout << "New " << CurrentRingCenter.X() << " " << CurrentRingCenter.Y() << " " << CurrentRingRadius << " " << CurrentRingChi2 << endl;

  return 1;
}

Bool_t RICHImprovedRing::CorrectHitPosition(TRecoRICHHit *hit, TVector3 PMTPos, TVector3 posAtRICHMirrors) {

  Double_t p[25][6][2];
  Double_t v[25][25][25][2];
  Double_t ac[25][28];
  Double_t bcost[25][28];
  Int_t fakemirror[25][25];



  TVector3 correction(0,0,0);
  Double_t mircorr = hit->GetROChannelID()<1000 ? 127. : 177.;
  if (fabs(hit->GetPosition().X())>2000. || fabs(hit->GetPosition().Y())>2000.) return 0;
  Double_t alpha = atan2(hit->GetPosition().X()-PMTPos.X()-mircorr,hit->GetPosition().Y()-PMTPos.Y());
  if (alpha<0) alpha += TMath::Pi();
  if (alpha<0.000000000000001 || fabs(alpha-TMath::Pi())<0.000000000000001) alpha=0.;
  Double_t gammangle = atan2(18.,(2*sqrt(pow(hit->GetPosition().X()-PMTPos.X()-mircorr,2)+pow(hit->GetPosition().Y()-PMTPos.Y(),2))));
  Double_t alpha1 = alpha+gammangle;
  Double_t alpha2 = alpha-gammangle;
  Double_t Point1_X = posAtRICHMirrors.X();
  Double_t Point1_Y = posAtRICHMirrors.Y();
  Double_t Point2_X = posAtRICHMirrors.X()+sqrt(pow(hit->GetPosition().X()-PMTPos.X()-mircorr,2)+pow(hit->GetPosition().Y()-PMTPos.Y(),2))*sin(alpha1);
  Double_t Point2_Y = posAtRICHMirrors.Y()+sqrt(pow(hit->GetPosition().X()-PMTPos.X()-mircorr,2)+pow(hit->GetPosition().Y()-PMTPos.Y(),2))*cos(alpha1);
  Double_t Point3_X = posAtRICHMirrors.X()+sqrt(pow(hit->GetPosition().X()-PMTPos.X()-mircorr,2)+pow(hit->GetPosition().Y()-PMTPos.Y(),2))*sin(alpha2);
  Double_t Point3_Y = posAtRICHMirrors.Y()+sqrt(pow(hit->GetPosition().X()-PMTPos.X()-mircorr,2)+pow(hit->GetPosition().Y()-PMTPos.Y(),2))*cos(alpha2);

  Double_t vertex_X;
  Double_t vertex_Y;
  Int_t M1 = tools->MirrorSurface(Point1_X,Point1_Y,-1.,1);
  Int_t M2 = tools->MirrorSurface(Point2_X,Point2_Y,0.,1);
  Int_t M3 = tools->MirrorSurface(Point3_X,Point3_Y,0.,1);
  TVector3 new_hit_position; 
  if (M1==23 || M1 == 24 ) return 0;
  if (M1==99 && M2==99 && M3==99) return 0;
  if (M1<=0 || (M1>24 && M1!=99)) return 0;
  if (M2<=0 || (M2>24 && M2!=99)) return 0;
  if (M3<=0 || (M3>24 && M3!=99)) return 0;
  if (M1==M2 && M2==M3 && M1!=99 && M1!=23 && M1!=24) {
    new_hit_position.SetX(hit->GetPosition().X()-fDeltaMisa[M1][0]);
    new_hit_position.SetY(hit->GetPosition().Y()-fDeltaMisa[M1][1]); 
    new_hit_position.SetZ(hit->GetPosition().Z());
    hit->SetPosition(new_hit_position);
    return 1;
  } else if (M1==M2 && M1!=M3 &&  M1!=99 && M1!=23 && M1!=24 && (M3==99 || M3==23 || M3==24)){
    new_hit_position.SetX(hit->GetPosition().X()-fDeltaMisa[M1][0]);
    new_hit_position.SetY(hit->GetPosition().Y()-fDeltaMisa[M1][1]); 
    new_hit_position.SetZ(hit->GetPosition().Z());
    hit->SetPosition(new_hit_position);
    return 1;
  } else if (M1==M3 && M1!=M2 &&  M1!=99 && M1!=23 && M1!=24 && (M2==99 || M2==23 || M2==24)){
    new_hit_position.SetX(hit->GetPosition().X()-fDeltaMisa[M1][0]);
    new_hit_position.SetY(hit->GetPosition().Y()-fDeltaMisa[M1][1]); 
    new_hit_position.SetZ(hit->GetPosition().Z());
    hit->SetPosition(new_hit_position);
  } else if (M2==M3 && M1!=M2 &&  M2!=99 && M2!=23 && M2!=24 && (M1==99 || M1==23 || M1==24)){
    new_hit_position.SetX(hit->GetPosition().X()-fDeltaMisa[M2][0]);
    new_hit_position.SetY(hit->GetPosition().Y()-fDeltaMisa[M2][1]); 
    new_hit_position.SetZ(hit->GetPosition().Z());
    hit->SetPosition(new_hit_position);
    return 1;
  } else if ( M1!=99 && M1!=23 && M1!=24 && (M2==99 || M2==23 || M2==24) && (M3==99 || M3==23 || M3==24)){
    new_hit_position.SetX(hit->GetPosition().X()-fDeltaMisa[M1][0]);
    new_hit_position.SetY(hit->GetPosition().Y()-fDeltaMisa[M1][1]); 
    new_hit_position.SetZ(hit->GetPosition().Z());
    hit->SetPosition(new_hit_position);
    return 1;
  } else if ( M2!=99 && M2!=23 && M2!=24 && (M1==99 || M1==23 || M1==24) && (M3==99 || M3==23 || M3==24)){
    new_hit_position.SetX(hit->GetPosition().X()-fDeltaMisa[M2][0]);
    new_hit_position.SetY(hit->GetPosition().Y()-fDeltaMisa[M2][1]); 
    new_hit_position.SetZ(hit->GetPosition().Z());
    hit->SetPosition(new_hit_position);
    return 1;
  } else if ( M3!=99 && M3!=23 && M3!=24 && (M1==99 || M1==23 || M1==24) && (M2==99 || M2==23 || M2==24)){
    new_hit_position.SetX(hit->GetPosition().X()-fDeltaMisa[M3][0]);
    new_hit_position.SetY(hit->GetPosition().Y()-fDeltaMisa[M3][1]); 
    new_hit_position.SetZ(hit->GetPosition().Z());
    hit->SetPosition(new_hit_position);
    return 1;
  } else if ( M1==M2 && M1!=M3 && M3!=99 && M3!=23 && M3!=24 && M1!=99 && M1!=23 && M1!=24){
    Double_t pendM1M3 = (Point1_Y-Point3_Y)/(Point1_X-Point3_X);
    Double_t costM1M3 = (Point3_Y*Point1_X-Point3_X*Point1_Y)/(Point1_X-Point3_X);
    Double_t pendM2M3 = (Point2_Y-Point3_Y)/(Point2_X-Point3_X);
    Double_t costM2M3 = (Point3_Y*Point2_X-Point3_X*Point2_Y)/(Point2_X-Point3_X);
    Double_t bourder13_X, bourder13_Y;
    if (ac[M1][M3]==0.) {
      bourder13_X=bcost[M1][M3];
      bourder13_Y=pendM1M3*bcost[M1][M3]+costM1M3;
    } else {
      bourder13_X= (costM1M3-bcost[M1][M3])/(ac[M1][M3]-pendM1M3);
      bourder13_Y= (costM1M3*ac[M1][M3]-bcost[M1][M3]*pendM1M3)/(ac[M1][M3]-pendM1M3);
    }
    Double_t bourder23_X, bourder23_Y;
    if (ac[M2][M3]==0.) {
      bourder23_X=bcost[M2][M3];
      bourder23_Y=pendM2M3*bcost[M2][M3]+costM2M3;
    } else {
      bourder23_X= (costM2M3-bcost[M2][M3])/(ac[M2][M3]-pendM2M3);
      bourder23_Y= (costM2M3*ac[M2][M3]-bcost[M2][M3]*pendM2M3)/(ac[M2][M3]-pendM2M3);
    }
    Double_t area1 = fabs(Point1_X*Point2_Y+Point2_X*bourder13_Y+bourder13_X*bourder23_Y+bourder23_X*Point1_Y-Point2_X*Point1_Y-bourder13_X*Point2_Y-bourder23_X*bourder13_Y-Point1_X*bourder23_Y)/2.;
    if (area1>2000.) area1=1.;
    Double_t misPoint1_X = Point1_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
    Double_t misPoint1_Y = Point1_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
    Double_t misPoint2_X = Point2_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
    Double_t misPoint2_Y = Point2_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
    Double_t misPoint3_X = Point3_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
    Double_t misPoint3_Y = Point3_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
    Double_t mispendM1M3 = (misPoint1_Y-misPoint3_Y)/(misPoint1_X-misPoint3_X);
    Double_t miscostM1M3 = (misPoint3_Y*misPoint1_X-misPoint3_X*misPoint1_Y)/(misPoint1_X-misPoint3_X);
    Double_t mispendM2M3 = (misPoint2_Y-misPoint3_Y)/(misPoint2_X-misPoint3_X);
    Double_t miscostM2M3 = (misPoint3_Y*misPoint2_X-misPoint3_X*misPoint2_Y)/(misPoint2_X-misPoint3_X);
    Double_t misbourder13_X, misbourder13_Y;
    if (ac[M1][M3]==0.) {
      misbourder13_X=bcost[M1][M3];
      misbourder13_Y=mispendM1M3*bcost[M1][M3]+miscostM1M3;
    } else {
      misbourder13_X= (miscostM1M3-bcost[M1][M3])/(ac[M1][M3]-mispendM1M3);
      misbourder13_Y= (miscostM1M3*ac[M1][M3]-bcost[M1][M3]*mispendM1M3)/(ac[M1][M3]-mispendM1M3);
    }
    Double_t misbourder23_X, misbourder23_Y;
    if (ac[M2][M3]==0.) {
      misbourder23_X=bcost[M2][M3];
      misbourder23_Y=mispendM2M3*bcost[M2][M3]+miscostM2M3;
    } else {
      misbourder23_X= (miscostM2M3-bcost[M2][M3])/(ac[M2][M3]-mispendM2M3);
      misbourder23_Y= (miscostM2M3*ac[M2][M3]-bcost[M2][M3]*mispendM2M3)/(ac[M2][M3]-mispendM2M3);
    }
    Double_t area2=fabs(misPoint3_X*misbourder13_Y+misbourder13_X*misbourder23_Y+misbourder23_X*misPoint3_Y-misbourder13_X*misPoint3_Y-misbourder23_X*misbourder13_Y-misPoint3_X*misbourder23_Y)/2.;
    if (area2>2000.) area2=1.;
    correction = ComputeCorrection(M1,M3,area1,area2);
  } else if (M1==M3 && M1!=M2 && M2!=99 && M2!=23 && M2!=24  && M1!=99 && M1!=23 && M1!=24){
    Double_t pendM1M2 = (Point1_Y-Point2_Y)/(Point1_X-Point2_X);
    Double_t costM1M2 = (Point2_Y*Point1_X-Point2_X*Point1_Y)/(Point1_X-Point2_X);
    Double_t pendM2M3 = (Point2_Y-Point3_Y)/(Point2_X-Point3_X);
    Double_t costM2M3 = (Point3_Y*Point2_X-Point3_X*Point2_Y)/(Point2_X-Point3_X);
    Double_t bourder12_X, bourder12_Y;
    if (ac[M1][M2]==0.) {
      bourder12_X=bcost[M1][M2];
      bourder12_Y=pendM1M2*bcost[M1][M2]+costM1M2;
    } else {
      bourder12_X= (costM1M2-bcost[M1][M2])/(ac[M1][M2]-pendM1M2);
      bourder12_Y= (costM1M2*ac[M1][M2]-bcost[M1][M2]*pendM1M2)/(ac[M1][M2]-pendM1M2);
    }
    Double_t bourder23_X, bourder23_Y;
    if (ac[M2][M3]==0.) {
      bourder23_X=bcost[M2][M3];
      bourder23_Y=pendM2M3*bcost[M2][M3]+costM2M3;
    } else {
      bourder23_X= (costM2M3-bcost[M2][M3])/(ac[M2][M3]-pendM2M3);
      bourder23_Y= (costM2M3*ac[M2][M3]-bcost[M2][M3]*pendM2M3)/(ac[M2][M3]-pendM2M3);
    }
    Double_t area1=fabs(Point1_X*Point3_Y+Point3_X*bourder12_Y+bourder12_X*bourder23_Y+bourder23_X*Point1_Y-Point3_X*Point1_Y-bourder12_X*Point3_Y-bourder23_X*bourder12_Y-Point1_X*bourder23_Y)/2.;
    if (area1>2000.) area1=1.;
    Double_t misPoint1_X=Point1_X+fDeltaMisa[M1][0]-fDeltaMisa[M2][0];
    Double_t misPoint1_Y=Point1_Y+fDeltaMisa[M1][1]-fDeltaMisa[M2][1];
    Double_t misPoint2_X=Point2_X+fDeltaMisa[M1][0]-fDeltaMisa[M2][0];
    Double_t misPoint2_Y=Point2_Y+fDeltaMisa[M1][1]-fDeltaMisa[M2][1];
    Double_t misPoint3_X=Point3_X+fDeltaMisa[M1][0]-fDeltaMisa[M2][0];
    Double_t misPoint3_Y=Point3_Y+fDeltaMisa[M1][1]-fDeltaMisa[M2][1];
    Double_t mispendM1M2=(misPoint1_Y-misPoint2_Y)/(misPoint1_X-misPoint2_X);
    Double_t miscostM1M2=(misPoint2_Y*misPoint1_X-misPoint2_X*misPoint1_Y)/(misPoint1_X-misPoint2_X);
    Double_t mispendM2M3=(misPoint2_Y-misPoint3_Y)/(misPoint2_X-misPoint3_X);
    Double_t miscostM2M3=(misPoint3_Y*misPoint2_X-misPoint3_X*misPoint2_Y)/(misPoint2_X-misPoint3_X);
    Double_t misbourder12_X, misbourder12_Y;
    if (ac[M1][M2]==0.) {
      misbourder12_X=bcost[M1][M2];
      misbourder12_Y=mispendM1M2*bcost[M1][M2]+miscostM1M2;
    } else {
      misbourder12_X= (miscostM1M2-bcost[M1][M2])/(ac[M1][M2]-mispendM1M2);
      misbourder12_Y= (miscostM1M2*ac[M1][M2]-bcost[M1][M2]*mispendM1M2)/(ac[M1][M2]-mispendM1M2);
    }
    Double_t misbourder23_X, misbourder23_Y;
    if (ac[M2][M3]==0.) {
      misbourder23_X=bcost[M2][M3];
      misbourder23_Y=mispendM2M3*bcost[M2][M3]+miscostM2M3;
    } else {
      misbourder23_X= (miscostM2M3-bcost[M2][M3])/(ac[M2][M3]-mispendM2M3);
      misbourder23_Y= (miscostM2M3*ac[M2][M3]-bcost[M2][M3]*mispendM2M3)/(ac[M2][M3]-mispendM2M3);
    }
    Double_t area2=fabs(misPoint2_X*misbourder12_Y+misbourder12_X*misbourder23_Y+misbourder23_X*misPoint2_Y-misbourder12_X*misPoint2_Y-misbourder23_X*misbourder12_Y-misPoint2_X*misbourder23_Y)/2.;
    if (area2>2000.) area2=1.;
    correction = ComputeCorrection(M1,M2,area1,area2);
  } else if (M2==M3 && M1!=M2  && M1!=99 && M1!=23 && M1!=24 && M2!=99 && M2!=23 && M2!=24){
    Double_t pendM1M3 = (Point1_Y-Point3_Y)/(Point1_X-Point3_X);
    Double_t costM1M3 = (Point3_Y*Point1_X-Point3_X*Point1_Y)/(Point1_X-Point3_X);
    Double_t pendM1M2 = (Point1_Y-Point2_Y)/(Point1_X-Point2_X);
    Double_t costM1M2 = (Point2_Y*Point1_X-Point2_X*Point1_Y)/(Point1_X-Point2_X);
    Double_t bourder13_X, bourder13_Y;
    if (ac[M1][M3]==0.) {
      bourder13_X=bcost[M1][M3];
      bourder13_Y=pendM1M3*bcost[M1][M3]+costM1M3;
    } else {
      bourder13_X= (costM1M3-bcost[M1][M3])/(ac[M1][M3]-pendM1M3);
      bourder13_Y= (costM1M3*ac[M1][M3]-bcost[M1][M3]*pendM1M3)/(ac[M1][M3]-pendM1M3);
    }
    Double_t bourder12_X, bourder12_Y;
    if (ac[M1][M2]==0.) {
      bourder12_X=bcost[M1][M2];
      bourder12_Y=pendM1M2*bcost[M1][M2]+costM1M2;
    } else {
      bourder12_X= (costM1M2-bcost[M1][M2])/(ac[M1][M2]-pendM1M2);
      bourder12_Y= (costM1M2*ac[M1][M2]-bcost[M1][M2]*pendM1M2)/(ac[M1][M2]-pendM1M2);
    }
    Double_t area1=fabs(Point1_X*bourder12_Y+bourder12_X*bourder13_Y+bourder13_X*Point1_Y-bourder12_X*Point1_Y-bourder13_X*bourder12_Y-Point1_X*bourder13_Y)/2.;
    if (area1>2000.) area1=1.;
    Double_t misPoint1_X = Point1_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
    Double_t misPoint1_Y = Point1_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
    Double_t misPoint2_X = Point2_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
    Double_t misPoint2_Y = Point2_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
    Double_t misPoint3_X = Point3_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
    Double_t misPoint3_Y = Point3_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
    Double_t mispendM1M3 = (misPoint1_Y-misPoint3_Y)/(misPoint1_X-misPoint3_X);
    Double_t miscostM1M3 = (misPoint3_Y*misPoint1_X-misPoint3_X*misPoint1_Y)/(misPoint1_X-misPoint3_X);
    Double_t mispendM1M2 = (misPoint1_Y-misPoint2_Y)/(misPoint1_X-misPoint2_X);
    Double_t miscostM1M2 = (misPoint2_Y*misPoint1_X-misPoint2_X*misPoint1_Y)/(misPoint1_X-misPoint2_X);
    Double_t misbourder13_X, misbourder13_Y;
    if (ac[M1][M3]==0.) {
      misbourder13_X=bcost[M1][M3];
      misbourder13_Y=mispendM1M3*bcost[M1][M3]+miscostM1M3;
    } else {
      misbourder13_X= (miscostM1M3-bcost[M1][M3])/(ac[M1][M3]-mispendM1M3);
      misbourder13_Y= (miscostM1M3*ac[M1][M3]-bcost[M1][M3]*mispendM1M3)/(ac[M1][M3]-mispendM1M3);
    }
    Double_t misbourder12_X, misbourder12_Y;
    if (ac[M1][M2]==0.) {
      misbourder12_X=bcost[M1][M2];
      misbourder12_Y=mispendM1M2*bcost[M1][M2]+miscostM1M2;
    } else {
      misbourder12_X= (miscostM1M2-bcost[M1][M2])/(ac[M1][M2]-mispendM1M2);
      misbourder12_Y= (miscostM1M2*ac[M1][M2]-bcost[M1][M2]*mispendM1M2)/(ac[M1][M2]-mispendM1M2);
    }
    Double_t area2=fabs(misPoint3_X*misPoint2_Y+misPoint2_X*misbourder12_Y+misbourder12_X*misbourder13_Y+misbourder13_X*misPoint3_Y-misPoint2_X*misPoint3_Y-misbourder12_X*misPoint2_Y-misbourder13_X*misbourder12_Y-misPoint3_X*misbourder13_Y)/2.;
    if (area2>2000.) area2=1.;
    correction = ComputeCorrection(M1,M2,area1,area2);
  } else if (M1!=M2 && M1!=M3 && M2!=M3  && M2!=99 && M2!= 23 && M2!=24 && M3!=99 && M3!=23 && M3!=24 && (M1==99 || M1==23 || M1==24)){
    Double_t misPoint2_X = Point2_X-fDeltaMisa[M2][0];
    Double_t misPoint2_Y = Point2_Y-fDeltaMisa[M2][1];
    Double_t misPoint3_X = Point3_X-fDeltaMisa[M2][0];
    Double_t misPoint3_Y = Point3_Y-fDeltaMisa[M2][1];
    Double_t mispendM2M3 = (misPoint2_Y-misPoint3_Y)/(misPoint2_X-misPoint3_X);
    Double_t miscostM2M3 = (misPoint3_Y*misPoint2_X-misPoint3_X*misPoint2_Y)/(misPoint2_X-misPoint3_X);
    Double_t misbourder23_X, misbourder23_Y;
    if (ac[M2][M3]==0.) {
      misbourder23_X=bcost[M2][M3];
      misbourder23_Y=mispendM2M3*bcost[M2][M3]+miscostM2M3;
    } else {
      misbourder23_X= (miscostM2M3-bcost[M2][M3])/(ac[M2][M3]-mispendM2M3);
      misbourder23_Y= (miscostM2M3*ac[M2][M3]-bcost[M2][M3]*mispendM2M3)/(ac[M2][M3]-mispendM2M3);
    }
    Double_t area1=fabs(Point1_X*misbourder23_Y+misbourder23_X*misPoint2_Y+misPoint2_X*Point1_Y-misbourder23_X*Point1_Y-misPoint2_X*misbourder23_Y-Point1_X*misPoint2_Y)/2.;
    if (area1>2000.) area1=1.;
    Double_t mismisPoint2_X=Point2_X-fDeltaMisa[M3][0];
    Double_t mismisPoint2_Y=Point2_Y-fDeltaMisa[M3][1];
    Double_t mismisPoint3_X=Point3_X-fDeltaMisa[M3][0];
    Double_t mismisPoint3_Y=Point3_Y-fDeltaMisa[M3][1];
    Double_t mismispendM2M3=(mismisPoint2_Y-mismisPoint3_Y)/(mismisPoint2_X-mismisPoint3_X);
    Double_t mismiscostM2M3=(mismisPoint3_Y*mismisPoint2_X-mismisPoint3_X*mismisPoint2_Y)/(mismisPoint2_X-mismisPoint3_X);
    Double_t mismisbourder23_X, mismisbourder23_Y;
    if (ac[M2][M3]==0.) {
      mismisbourder23_X=bcost[M2][M3];
      mismisbourder23_Y=mismispendM2M3*bcost[M2][M3]+mismiscostM2M3;
    } else {
      mismisbourder23_X= (mismiscostM2M3-bcost[M2][M3])/(ac[M2][M3]-mismispendM2M3);
      mismisbourder23_Y= (mismiscostM2M3*ac[M2][M3]-bcost[M2][M3]*mismispendM2M3)/(ac[M2][M3]-mismispendM2M3);
    }
    Double_t area2=fabs(Point1_X*mismisbourder23_Y+mismisbourder23_X*mismisPoint3_Y+mismisPoint3_X*Point1_Y-mismisbourder23_X*Point1_Y-mismisPoint3_X*mismisbourder23_Y-Point1_X*mismisPoint3_Y)/2.;
    if (area2>2000.) area2=1.;
  } else if (M1!=M2 && M1!=M3 && M2!=M3  && M1!=99 && M1!= 23 && M1!=24 && M3!=99 && M3!=23 && M3!=24 && (M2==99 || M2==23 || M2==24)){
    if (fakemirror[M1][M3]==0 || M2!=23 || M2!=24){
      Double_t pendM1M3=(Point1_Y-Point3_Y)/(Point1_X-Point3_X);
      Double_t costM1M3=(Point3_Y*Point1_X-Point3_X*Point1_Y)/(Point1_X-Point3_X);
      Double_t bourder13_X, bourder13_Y;
      if (ac[M1][M3]==0.) {
        bourder13_X=bcost[M1][M3];
        bourder13_Y=pendM1M3*bcost[M1][M3]+costM1M3;
      } else {
        bourder13_X= (costM1M3-bcost[M1][M3])/(ac[M1][M3]-pendM1M3);
        bourder13_Y= (costM1M3*ac[M1][M3]-bcost[M1][M3]*pendM1M3)/(ac[M1][M3]-pendM1M3);
      }
      Double_t area1=fabs(Point1_X*bourder13_Y+bourder13_X*Point2_Y+Point2_X*Point1_Y-bourder13_X*Point1_Y-Point2_X*bourder13_Y-Point1_X*Point2_Y)/2.;
      if (area1>2000.) area1=1.;
      Double_t misPoint1_X = Point1_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
      Double_t misPoint1_Y = Point1_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
      Double_t misPoint2_X = Point2_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
      Double_t misPoint2_Y = Point2_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
      Double_t misPoint3_X = Point3_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
      Double_t misPoint3_Y = Point3_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
      Double_t mispendM1M3=(misPoint1_Y-misPoint3_Y)/(misPoint1_X-misPoint3_X);
      Double_t miscostM1M3=(misPoint3_Y*misPoint1_X-misPoint3_X*misPoint1_Y)/(misPoint1_X-misPoint3_X);
      Double_t misbourder13_X, misbourder13_Y;
      if (ac[M1][M3]==0.) {
        misbourder13_X=bcost[M1][M3];
        misbourder13_Y=mispendM1M3*bcost[M1][M3]+miscostM1M3;
      } else {
        misbourder13_X= (miscostM1M3-bcost[M1][M3])/(ac[M1][M3]-mispendM1M3);
        misbourder13_Y= (miscostM1M3*ac[M1][M3]-bcost[M1][M3]*mispendM1M3)/(ac[M1][M3]-mispendM1M3);
      }
      Double_t area2=fabs(misPoint3_X*misbourder13_Y+misbourder13_X*misPoint2_Y+misPoint2_X*misPoint3_Y-misbourder13_X*misPoint3_Y-misPoint2_X*misbourder13_Y-misPoint3_X*misPoint2_Y)/2.;
      if (area2>2000.) area2=1.;
    } else {
      if (M1<M3) {
        vertex_X=v[0][M1][M3][0];
        vertex_Y=v[0][M1][M3][1];
      } else {
        vertex_X=v[0][M3][M1][0];
        vertex_Y=v[0][M3][M1][1];
      }
      Double_t pendM1M3=(Point1_Y-Point3_Y)/(Point1_X-Point3_X);
      Double_t costM1M3=(Point3_Y*Point1_X-Point3_X*Point1_Y)/(Point1_X-Point3_X);
      Double_t pendM1M2=(Point1_Y-Point2_Y)/(Point1_X-Point2_X);
      Double_t costM1M2=(Point2_Y*Point1_X-Point2_X*Point1_Y)/(Point1_X-Point2_X);
      Double_t bourder13_X, bourder13_Y;
      if (ac[M1][M3]==0.) {
        bourder13_X=bcost[M1][M3];
        bourder13_Y=pendM1M3*bcost[M1][M3]+costM1M3;
      } else {
        bourder13_X= (costM1M3-bcost[M1][M3])/(ac[M1][M3]-pendM1M3);
        bourder13_Y= (costM1M3*ac[M1][M3]-bcost[M1][M3]*pendM1M3)/(ac[M1][M3]-pendM1M3);
      }
      Double_t bourder12_X, bourder12_Y;
      if (ac[M1][fakemirror[M1][M3]]==0.) {
        bourder12_X=bcost[M1][fakemirror[M1][M3]];
        bourder12_Y=pendM1M2*bcost[M1][fakemirror[M1][M3]]+costM1M2;
      } else {
        bourder12_X= (costM1M2-bcost[M1][fakemirror[M1][M3]])/(ac[M1][fakemirror[M1][M3]]-pendM1M2);
        bourder12_Y= (costM1M2*ac[M1][fakemirror[M1][M3]]-bcost[M1][fakemirror[M1][M3]]*pendM1M2)/(ac[M1][fakemirror[M1][M3]]-pendM1M2);
      }
      Double_t area1=fabs(Point1_X*vertex_Y+vertex_X*bourder13_Y+bourder13_X*bourder12_Y+bourder12_X*Point1_Y-vertex_X*Point1_Y-bourder13_X*vertex_Y-bourder12_X*bourder13_Y-Point1_X*bourder12_Y)/2.;
      if (area1>2000.) area1=1.;
      Double_t misPoint1_X=Point1_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
      Double_t misPoint1_Y=Point1_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
      Double_t misPoint2_X=Point2_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
      Double_t misPoint2_Y=Point2_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
      Double_t misPoint3_X=Point3_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
      Double_t misPoint3_Y=Point3_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
      Double_t mispendM1M3=(misPoint1_Y-misPoint3_Y)/(misPoint1_X-misPoint3_X);
      Double_t miscostM1M3=(misPoint3_Y*misPoint1_X-misPoint3_X*misPoint1_Y)/(misPoint1_X-misPoint3_X);
      Double_t mispendM2M3=(misPoint2_Y-misPoint3_Y)/(misPoint2_X-misPoint3_X);
      Double_t miscostM2M3=(misPoint3_Y*misPoint2_X-misPoint3_X*misPoint2_Y)/(misPoint2_X-misPoint3_X);
      Double_t misbourder13_X, misbourder13_Y;
      if (ac[M1][M3]==0.) {
        misbourder13_X=bcost[M1][M3];
        misbourder13_Y=mispendM1M3*bcost[M1][M3]+miscostM1M3;
      } else {
        misbourder13_X= (miscostM1M3-bcost[M1][M3])/(ac[M1][M3]-mispendM1M3);
        misbourder13_Y= (miscostM1M3*ac[M1][M3]-bcost[M1][M3]*mispendM1M3)/(ac[M1][M3]-mispendM1M3);
      }
      Double_t misbourder23_X, misbourder23_Y;
      if (ac[M3][fakemirror[M3][M1]]==0.) {
        misbourder23_X=bcost[M3][fakemirror[M3][M1]];
        misbourder23_Y=mispendM2M3*bcost[M3][fakemirror[M3][M1]]+miscostM2M3;
      } else {
        misbourder23_X= (miscostM2M3-bcost[M3][fakemirror[M3][M1]])/(ac[M3][fakemirror[M3][M1]]-mispendM2M3);
        misbourder23_Y= (miscostM2M3*ac[M3][fakemirror[M3][M1]]-bcost[M3][fakemirror[M3][M1]]*mispendM2M3)/(ac[M3][fakemirror[M3][M1]]-mispendM2M3);
      }
      Double_t area2=fabs(misPoint3_X*vertex_Y+vertex_X*misbourder13_Y+misbourder13_X*misbourder23_Y+misbourder23_X*misPoint3_Y-vertex_X*misPoint3_Y-misbourder13_X*vertex_Y-misbourder23_X*misbourder13_Y-misPoint3_X*misbourder23_Y)/2.;
      if (area2>2000.) area2=1.;
    }
  } else if (M1!=M2 && M1!=M3 && M2!=M3 && M1!=99 && M1!=23 && M1!=24 && M2!=99 && M2!=23 && M2!=24 && (M3==99 || M3==23 || M3==24)){
    if (fakemirror[M1][M2]==0 || M3!=23 || M3!=24){
      Double_t pendM1M2=(Point1_Y-Point2_Y)/(Point1_X-Point2_X);
      Double_t costM1M2=(Point2_Y*Point1_X-Point2_X*Point1_Y)/(Point1_X-Point2_X);
      Double_t bourder12_X, bourder12_Y;
      if (ac[M1][M2]==0.) {
        bourder12_X=bcost[M1][M2];
        bourder12_Y=pendM1M2*bcost[M1][M2]+costM1M2;
      } else {
        bourder12_X= (costM1M2-bcost[M1][M2])/(ac[M1][M2]-pendM1M2);
        bourder12_Y= (costM1M2*ac[M1][M2]-bcost[M1][M2]*pendM1M2)/(ac[M1][M2]-pendM1M2);
      }
      Double_t area1=fabs(Point1_X*bourder12_Y+bourder12_X*Point3_Y+Point3_X*Point1_Y-bourder12_X*Point1_Y-Point3_X*bourder12_Y-Point1_X*Point3_Y)/2.;
      if (area1>2000.) area1=1.;
      Double_t misPoint1_X=Point1_X+fDeltaMisa[M1][0]-fDeltaMisa[M2][0];
      Double_t misPoint1_Y=Point1_Y+fDeltaMisa[M1][1]-fDeltaMisa[M2][1];
      Double_t misPoint2_X=Point2_X+fDeltaMisa[M1][0]-fDeltaMisa[M2][0];
      Double_t misPoint2_Y=Point2_Y+fDeltaMisa[M1][1]-fDeltaMisa[M2][1];
      Double_t misPoint3_X=Point3_X+fDeltaMisa[M1][0]-fDeltaMisa[M2][0];
      Double_t misPoint3_Y=Point3_Y+fDeltaMisa[M1][1]-fDeltaMisa[M2][1];
      Double_t mispendM1M2=(misPoint1_Y-misPoint2_Y)/(misPoint1_X-misPoint2_X);
      Double_t miscostM1M2=(misPoint2_Y*misPoint1_X-misPoint2_X*misPoint1_Y)/(misPoint1_X-misPoint2_X);
      Double_t misbourder12_X, misbourder12_Y;
      if (ac[M1][M2]==0.) {
        misbourder12_X=bcost[M1][M2];
        misbourder12_Y=mispendM1M2*bcost[M1][M2]+miscostM1M2;
      } else {
        misbourder12_X= (miscostM1M2-bcost[M1][M2])/(ac[M1][M2]-mispendM1M2);
        misbourder12_Y= (miscostM1M2*ac[M1][M2]-bcost[M1][M2]*mispendM1M2)/(ac[M1][M2]-mispendM1M2);
      }
      Double_t area2=fabs(misPoint2_X*misbourder12_Y+misbourder12_X*misPoint3_Y+misPoint3_X*misPoint2_Y-misbourder12_X*misPoint2_Y-misPoint3_X*misbourder12_Y-misPoint2_X*misPoint3_Y)/2.;
      if (area2>2000.) area2=1.;
    } else {
      if (M1<M2) {
        vertex_X=v[0][M1][M2][0];
        vertex_Y=v[0][M1][M2][1];
      } else {
        vertex_X=v[0][M2][M1][0];
        vertex_Y=v[0][M2][M1][1];
      }
      Double_t pendM1M3=(Point1_Y-Point3_Y)/(Point1_X-Point3_X);
      Double_t costM1M3=(Point3_Y*Point1_X-Point3_X*Point1_Y)/(Point1_X-Point3_X);
      Double_t pendM1M2=(Point1_Y-Point2_Y)/(Point1_X-Point2_X);
      Double_t costM1M2=(Point2_Y*Point1_X-Point2_X*Point1_Y)/(Point1_X-Point2_X);
      Double_t bourder13_X, bourder13_Y;
      if (ac[M1][fakemirror[M1][M2]]==0.) {
        bourder13_X=bcost[M1][fakemirror[M1][M2]];
        bourder13_Y=pendM1M3*bcost[M1][fakemirror[M1][M2]]+costM1M3;
      } else {
        bourder13_X= (costM1M3-bcost[M1][fakemirror[M1][M2]])/(ac[M1][fakemirror[M1][M2]]-pendM1M3);
        bourder13_Y= (costM1M3*ac[M1][fakemirror[M1][M2]]-bcost[M1][fakemirror[M1][M2]]*pendM1M3)/(ac[M1][fakemirror[M1][M2]]-pendM1M3);
      }
      Double_t bourder12_X, bourder12_Y;
      if (ac[M1][M2]==0.) {
        bourder12_X=bcost[M1][M2];
        bourder12_Y=pendM1M2*bcost[M1][M2]+costM1M2;
      } else {
        bourder12_X= (costM1M2-bcost[M1][M2])/(ac[M1][M2]-pendM1M2);
        bourder12_Y= (costM1M2*ac[M1][M2]-bcost[M1][M2]*pendM1M2)/(ac[M1][M2]-pendM1M2);
      }
      Double_t area1=fabs(Point1_X*vertex_Y+vertex_X*bourder13_Y+bourder13_X*bourder12_Y+bourder12_X*Point1_Y-vertex_X*Point1_Y-bourder13_X*vertex_Y-bourder12_X*bourder13_Y-Point1_X*bourder12_Y)/2.;
      if (area1>2000.) area1=1.;
      Double_t misPoint1_X=Point1_X+fDeltaMisa[M1][0]-fDeltaMisa[M2][0];
      Double_t misPoint1_Y=Point1_Y+fDeltaMisa[M1][1]-fDeltaMisa[M2][1];
      Double_t misPoint2_X=Point2_X+fDeltaMisa[M1][0]-fDeltaMisa[M2][0];
      Double_t misPoint2_Y=Point2_Y+fDeltaMisa[M1][1]-fDeltaMisa[M2][1];
      Double_t misPoint3_X=Point3_X+fDeltaMisa[M1][0]-fDeltaMisa[M2][0];
      Double_t misPoint3_Y=Point3_Y+fDeltaMisa[M1][1]-fDeltaMisa[M2][1];
      Double_t mispendM1M2=(misPoint1_Y-misPoint2_Y)/(misPoint1_X-misPoint2_X);
      Double_t miscostM1M2=(misPoint2_Y*misPoint1_X-misPoint2_X*misPoint1_Y)/(misPoint1_X-misPoint2_X);
      Double_t mispendM2M3=(misPoint2_Y-misPoint3_Y)/(misPoint2_X-misPoint3_X);
      Double_t miscostM2M3=(misPoint3_Y*misPoint2_X-misPoint3_X*misPoint2_Y)/(misPoint2_X-misPoint3_X);
      Double_t misbourder12_X, misbourder12_Y;
      if (ac[M1][M2]==0.) {
        misbourder12_X=bcost[M1][M2];
        misbourder12_Y=mispendM1M2*bcost[M1][M2]+miscostM1M2;
      } else {
        misbourder12_X= (miscostM1M2-bcost[M1][M2])/(ac[M1][M2]-mispendM1M2);
        misbourder12_Y= (miscostM1M2*ac[M1][M2]-bcost[M1][M2]*mispendM1M2)/(ac[M1][M2]-mispendM1M2);
      }
      Double_t misbourder23_X, misbourder23_Y;
      if (ac[M2][fakemirror[M2][M1]]==0.) {
        misbourder23_X=bcost[M2][fakemirror[M2][M1]];
        misbourder23_Y=mispendM2M3*bcost[M2][fakemirror[M2][M1]]+miscostM2M3;
      } else {
        misbourder23_X= (miscostM2M3-bcost[M2][fakemirror[M2][M1]])/(ac[M2][fakemirror[M2][M1]]-mispendM2M3);
        misbourder23_Y= (miscostM2M3*ac[M2][fakemirror[M2][M1]]-bcost[M2][fakemirror[M2][M1]]*mispendM2M3)/(ac[M2][fakemirror[M2][M1]]-mispendM2M3);
      }
      Double_t area2=fabs(misPoint2_X*vertex_Y+vertex_X*misbourder12_Y+misbourder12_X*misbourder23_Y+misbourder23_X*misPoint2_Y-vertex_X*misPoint2_Y-misbourder12_X*vertex_Y-misbourder23_X*misbourder12_Y-misPoint2_X*misbourder23_Y)/2.;
      if (area2>2000.) area2=1.;
    }
  } else if (M1!=M2 && M1!=M3 && M2!=M3  && M1!=99 && M1!=23 && M1!=24 && M2!=99 && M2!=23 && M2!=24 && M3!=99 && M3!=23 && M3!=24){
    if (M1<M2 && M2<M3 && M1<M3) {
      vertex_X=v[M1][M2][M3][0];
      vertex_Y=v[M1][M2][M3][1];
    } else if (M1<M2 && M2>M3 && M1<M3) {
      vertex_X=v[M1][M3][M2][0];
      vertex_Y=v[M1][M3][M2][1];
    } else if (M1>M2 && M2<M3 && M1>M3) {
      vertex_X=v[M2][M3][M1][0];
      vertex_Y=v[M2][M3][M1][1];
    } else if (M1>M2 && M2<M3 && M1<M3) {
      vertex_X=v[M2][M1][M3][0];
      vertex_Y=v[M2][M1][M3][1];
    } else if (M1>M2 && M2>M3 && M1>M3) {
      vertex_X=v[M3][M2][M1][0];
      vertex_Y=v[M3][M2][M1][1];
    } else if (M1<M2 && M2>M3 && M1>M3) {
      vertex_X=v[M3][M1][M2][0];
      vertex_Y=v[M3][M1][M2][1];
    }
    Double_t pendM1M3=(Point1_Y-Point3_Y)/(Point1_X-Point3_X);
    Double_t costM1M3=(Point3_Y*Point1_X-Point3_X*Point1_Y)/(Point1_X-Point3_X);
    Double_t pendM1M2=(Point1_Y-Point2_Y)/(Point1_X-Point2_X);
    Double_t costM1M2=(Point2_Y*Point1_X-Point2_X*Point1_Y)/(Point1_X-Point2_X);
    Double_t bourder13_X, bourder13_Y;
    if (ac[M1][M3]==0.) {
      bourder13_X=bcost[M1][M3];
      bourder13_Y=pendM1M3*bcost[M1][M3]+costM1M3;
    } else {
      bourder13_X= (costM1M3-bcost[M1][M3])/(ac[M1][M3]-pendM1M3);
      bourder13_Y= (costM1M3*ac[M1][M3]-bcost[M1][M3]*pendM1M3)/(ac[M1][M3]-pendM1M3);
    }
    Double_t bourder12_X, bourder12_Y;
    if (ac[M1][M2]==0.) {
      bourder12_X=bcost[M1][M2];
      bourder12_Y=pendM1M2*bcost[M1][M2]+costM1M2;
    } else {
      bourder12_X= (costM1M2-bcost[M1][M2])/(ac[M1][M2]-pendM1M2);
      bourder12_Y= (costM1M2*ac[M1][M2]-bcost[M1][M2]*pendM1M2)/(ac[M1][M2]-pendM1M2);
    }
    Double_t area1=fabs(Point1_X*vertex_Y+vertex_X*bourder13_Y+bourder13_X*bourder12_Y+bourder12_X*Point1_Y-vertex_X*Point1_Y-bourder13_X*vertex_Y-bourder12_X*bourder13_Y-Point1_X*bourder12_Y)/2.;
    if (area1>2000.) area1=1.;
    Double_t misPoint1_X=Point1_X+fDeltaMisa[M1][0]-fDeltaMisa[M2][0];
    Double_t misPoint1_Y=Point1_Y+fDeltaMisa[M1][1]-fDeltaMisa[M2][1];
    Double_t misPoint2_X=Point2_X+fDeltaMisa[M1][0]-fDeltaMisa[M2][0];
    Double_t misPoint2_Y=Point2_Y+fDeltaMisa[M1][1]-fDeltaMisa[M2][1];
    Double_t misPoint3_X=Point3_X+fDeltaMisa[M1][0]-fDeltaMisa[M2][0];
    Double_t misPoint3_Y=Point3_Y+fDeltaMisa[M1][1]-fDeltaMisa[M2][1];
    Double_t mispendM1M2=(misPoint1_Y-misPoint2_Y)/(misPoint1_X-misPoint2_X);
    Double_t miscostM1M2=(misPoint2_Y*misPoint1_X-misPoint2_X*misPoint1_Y)/(misPoint1_X-misPoint2_X);
    Double_t mispendM2M3=(misPoint2_Y-misPoint3_Y)/(misPoint2_X-misPoint3_X);
    Double_t miscostM2M3=(misPoint3_Y*misPoint2_X-misPoint3_X*misPoint2_Y)/(misPoint2_X-misPoint3_X);
    Double_t misbourder12_X, misbourder12_Y;
    if (ac[M1][M2]==0.) {
      misbourder12_X=bcost[M1][M2];
      misbourder12_Y=mispendM1M2*bcost[M1][M2]+miscostM1M2;
    } else {
      misbourder12_X= (miscostM1M2-bcost[M1][M2])/(ac[M1][M2]-mispendM1M2);
      misbourder12_Y= (miscostM1M2*ac[M1][M2]-bcost[M1][M2]*mispendM1M2)/(ac[M1][M2]-mispendM1M2);
    }
    Double_t misbourder23_X, misbourder23_Y;
    if (ac[M2][M3]==0.) {
      misbourder23_X=bcost[M2][M3];
      misbourder23_Y=mispendM2M3*bcost[M2][M3]+miscostM2M3;
    } else {
      misbourder23_X= (miscostM2M3-bcost[M2][M3])/(ac[M2][M3]-mispendM2M3);
      misbourder23_Y= (miscostM2M3*ac[M2][M3]-bcost[M2][M3]*mispendM2M3)/(ac[M2][M3]-mispendM2M3);
    }
    Double_t area2=fabs(misPoint2_X*vertex_Y+vertex_X*misbourder12_Y+misbourder12_X*misbourder23_Y+misbourder23_X*misPoint2_Y-vertex_X*misPoint2_Y-misbourder12_X*vertex_Y-misbourder23_X*misbourder12_Y-misPoint2_X*misbourder23_Y)/2.;
    if (area2>2000.) area2=1.;
    Double_t mismisPoint1_X=Point1_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
    Double_t mismisPoint1_Y=Point1_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
    Double_t mismisPoint2_X=Point2_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
    Double_t mismisPoint2_Y=Point2_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
    Double_t mismisPoint3_X=Point3_X+fDeltaMisa[M1][0]-fDeltaMisa[M3][0];
    Double_t mismisPoint3_Y=Point3_Y+fDeltaMisa[M1][1]-fDeltaMisa[M3][1];
    Double_t mismispendM1M3=(mismisPoint1_Y-mismisPoint3_Y)/(mismisPoint1_X-mismisPoint3_X);
    Double_t mismiscostM1M3=(mismisPoint3_Y*mismisPoint1_X-mismisPoint3_X*mismisPoint1_Y)/(mismisPoint1_X-mismisPoint3_X);
    Double_t mismispendM2M3=(mismisPoint2_Y-mismisPoint3_Y)/(mismisPoint2_X-mismisPoint3_X);
    Double_t mismiscostM2M3=(mismisPoint3_Y*mismisPoint2_X-mismisPoint3_X*mismisPoint2_Y)/(mismisPoint2_X-mismisPoint3_X);
    Double_t mismisbourder13_X, mismisbourder13_Y;
    if (ac[M1][M3]==0.) {
      mismisbourder13_X=bcost[M1][M3];
      mismisbourder13_Y=mismispendM1M3*bcost[M1][M3]+mismiscostM1M3;
    } else {
      mismisbourder13_X= (mismiscostM1M3-bcost[M1][M3])/(ac[M1][M3]-mismispendM1M3);
      mismisbourder13_Y= (mismiscostM1M3*ac[M1][M3]-bcost[M1][M3]*mismispendM1M3)/(ac[M1][M3]-mismispendM1M3);
    }
    Double_t mismisbourder23_X, mismisbourder23_Y;
    if (ac[M3][M2]==0.) {
      mismisbourder23_X=bcost[M3][M2];
      mismisbourder23_Y=mismispendM2M3*bcost[M3][M2]+mismiscostM2M3;
    } else {
      mismisbourder23_X= (mismiscostM2M3-bcost[M3][M2])/(ac[M3][M2]-mismispendM2M3);
      mismisbourder23_Y= (mismiscostM2M3*ac[M3][M2]-bcost[M3][M2]*mismispendM2M3)/(ac[M3][M2]-mismispendM2M3);
    }
    Double_t area3=fabs(mismisPoint3_X*vertex_Y+vertex_X*mismisbourder13_Y+mismisbourder13_X*mismisbourder23_Y+mismisbourder23_X*mismisPoint3_Y-vertex_X*mismisPoint3_Y-mismisbourder13_X*vertex_Y-mismisbourder23_X*mismisbourder13_Y-mismisPoint3_X*mismisbourder23_Y)/2.;
    if (area3>2000.) area3=1.;
    new_hit_position.SetX(hit->GetPosition().X()-(fDeltaMisa[M1][0]*pow(area1,2)+fDeltaMisa[M2][0]*pow(area2,2)+fDeltaMisa[M3][0]*pow(area3,2))/(pow(area1,2)+pow(area2,2)+pow(area3,2)));
    new_hit_position.SetY(hit->GetPosition().Y()-(fDeltaMisa[M1][1]*pow(area1,2)+fDeltaMisa[M2][1]*pow(area2,2)+fDeltaMisa[M3][1]*pow(area3,2))/(pow(area1,2)+pow(area2,2)+pow(area3,2))); 
    new_hit_position.SetZ(hit->GetPosition().Z());
    hit->SetPosition(new_hit_position);
    return 1;
  } else {
    Int_t M1_up=tools->MirrorSurface(Point1_X, Point1_Y,10.,1);
    Int_t M1_dw=tools->MirrorSurface(Point1_X, Point1_Y,-10.,1);
    if (M1_up>=1 && M1_up<23 && M1_dw>=1 && M1_dw<23) {
       new_hit_position.SetX(hit->GetPosition().X()-(fDeltaMisa[M1_up][0]+fDeltaMisa[M1_dw][0])/2.);
       new_hit_position.SetY(hit->GetPosition().Y()-(fDeltaMisa[M1_up][1]+fDeltaMisa[M1_dw][1])/2.); 
    } else if (M1_up>=1 && M1_up<23){
       new_hit_position.SetX(hit->GetPosition().X()-fDeltaMisa[M1_up][0]);
       new_hit_position.SetY(hit->GetPosition().Y()-fDeltaMisa[M1_up][1]); 
    } else if (M1_dw>=1 && M1_dw<23){
       new_hit_position.SetX(hit->GetPosition().X()-fDeltaMisa[M1_dw][0]);
       new_hit_position.SetY(hit->GetPosition().Y()-fDeltaMisa[M1_dw][1]); 
    }
    new_hit_position.SetZ(hit->GetPosition().Z());


    //hit->SetPosition(new_hit_position);
    return 1;
  }

  // Apply the mirror correction
  hit->SetPosition(hit->GetPosition()-correction);
  hit->SetPosition(hit->GetPosition());
  return 1;
}

TVector3 RICHImprovedRing::ComputeCorrection(Int_t ma, Int_t mb, Double_t a1, Double_t a2) {
  TVector3 corr;
//  Double_t denom = a1+a2;
//  Double_t f1x = fDeltaMisa[ma][0];
//  Double_t f2x = fDeltaMisa[mb][0];
//  Double_t f1y = fDeltaMisa[ma][1];
//  Double_t f2y = fDeltaMisa[mb][1];
//  corr.SetX((f1x*a1+f2x*a2)/denom);
//  corr.SetY((f1y*a1+f2y*a2)/denom);
////  corr.SetX((f1x));
////  corr.SetY((f1y));
//  corr.SetZ(0);

//  Double_t denom = a1*a1+a2*a2;
//  Double_t f1x = fDeltaMisa[ma][0];
//  Double_t f2x = fDeltaMisa[mb][0];
//  Double_t f1y = fDeltaMisa[ma][1];
//  Double_t f2y = fDeltaMisa[mb][1];
//  corr.SetX((f1x*a1*a1+f2x*a2*a2)/denom);
//  corr.SetY((f1y*a1*a1+f2y*a2*a2)/denom);
//  corr.SetZ(0);

  // No weight applied
  Double_t f1x = fDeltaMisa[ma][0];
  Double_t f2x = fDeltaMisa[mb][0];
  Double_t f1y = fDeltaMisa[ma][1];
  Double_t f2y = fDeltaMisa[mb][1];
  corr.SetX(f1x);
  corr.SetY(f1y);
  corr.SetZ(0);
  return corr;
}

void RICHImprovedRing::InitMirrorPositions()
{
  p[0][0][0]=0.;
  p[0][0][1]=0.;
  p[0][1][0]=0.;
  p[0][1][1]=0.;
  p[0][2][0]=0.;
  p[0][2][1]=0.;
  p[0][3][0]=0.;
  p[0][3][1]=0.;
  p[0][4][0]=0.;
  p[0][4][1]=0.;
  p[0][5][0]=0.;
  p[0][5][1]=0.;

  p[9][0][0]=-604.9822;
  p[9][0][1]=1407.273;
  p[9][1][0]=-302.1713;
  p[9][1][1]=1232.108;
  p[9][2][0]=-302.5643;
  p[9][2][1]=880.618;
  p[9][3][0]=-605.4903;
  p[9][3][1]=705.981;
  p[9][4][0]=-907.3783;
  p[9][4][1]=880.507;
  p[9][5][0]=-908.4193;
  p[9][5][1]=1232.219;

  p[17][0][0]=1.140884;
  p[17][0][1]=1410.055;
  p[17][1][0]=304.2186;
  p[17][1][1]=1234.862;
  p[17][2][0]=304.1706;
  p[17][2][1]=883.6355;
  p[17][3][0]=1.159304;
  p[17][3][1]=708.5905;
  p[17][4][0]=-301.5964;
  p[17][4][1]=883.6275;
  p[17][5][0]=-301.7034;
  p[17][5][1]=1234.87;

  p[4][0][0]=616.0545;
  p[4][0][1]=1411.085;
  p[4][1][0]=919.1667;
  p[4][1][1]=1236.041;
  p[4][2][0]=919.0697;
  p[4][2][1]=884.8148;
  p[4][3][0]=616.1178;
  p[4][3][1]=709.8158;
  p[4][4][0]=313.3247;
  p[4][4][1]=884.8058;
  p[4][5][0]=313.2847;
  p[4][5][1]=1236.05;

  p[10][0][0]=-906.7862;
  p[10][0][1]=879.7446;
  p[10][1][0]=-603.9165;
  p[10][1][1]=704.7246;
  p[10][2][0]=-603.9505;
  p[10][2][1]=353.3726;
  p[10][3][0]=-906.9273;
  p[10][3][1]=178.4316;
  p[10][4][0]=-1209.824;
  p[10][4][1]=353.3636;
  p[10][5][0]=-1209.788;
  p[10][5][1]=704.7336;

  p[20][0][0]=-299.2179;
  p[20][0][1]=879.3759;
  p[20][1][0]=3.730599;
  p[20][1][1]=704.5339;
  p[20][2][0]=3.581612;
  p[20][2][1]=353.1239;
  p[20][3][0]=-299.2848;
  p[20][3][1]=178.0999;
  p[20][4][0]=-602.1384;
  p[20][4][1]=353.2229;
  p[20][5][0]=-602.0394;
  p[20][5][1]=704.4349;

  p[12][0][0]=313.3375;
  p[12][0][1]=883.5286;
  p[12][1][0]=616.2212;
  p[12][1][1]=708.5726;
  p[12][2][0]=616.0572;
  p[12][2][1]=357.1526;
  p[12][3][0]=313.071;
  p[12][3][1]=182.1756;
  p[12][4][0]=10.27421;
  p[12][4][1]=357.1916;
  p[12][5][0]=10.37623;
  p[12][5][1]=708.5336;

  p[21][0][0]=921.0824;
  p[21][0][1]=883.6236;
  p[21][1][0]=1224.16;
  p[21][1][1]=708.5726;
  p[21][2][0]=1224.072;
  p[21][2][1]=357.3046;
  p[21][3][0]=921.1059;
  p[21][3][1]=182.3506;
  p[21][4][0]=618.2943;
  p[21][4][1]=357.3456;
  p[21][5][0]=618.2583;
  p[21][5][1]=708.5316;

  p[6][0][0]=-1208.882;
  p[6][0][1]=351.8385;
  p[6][1][0]=-905.842;
  p[6][1][1]=176.7515;
  p[6][2][0]=-906.043;
  p[6][2][1]=-174.5805;
  p[6][3][0]=-1209.044;
  p[6][3][1]=-349.4885;
  p[6][4][0]=-1511.81;
  p[6][4][1]=-174.5345;
  p[6][5][0]=-1511.817;
  p[6][5][1]=176.7055;

  p[5][0][0]=-601.9675;
  p[5][0][1]=352.9382;
  p[5][1][0]=-298.9482;
  p[5][1][1]=177.7952;
  p[5][2][0]=-298.9742;
  p[5][2][1]=-173.5028;
  p[5][3][0]=-602.0141;
  p[5][3][1]=-348.5198;
  p[5][4][0]=-904.9082;
  p[5][4][1]=-173.5058;
  p[5][5][0]=-904.8583;
  p[5][5][1]=177.7982;

  p[3][0][0]=617.0106;
  p[3][0][1]=353.4884;
  p[3][1][0]=919.8736;
  p[3][1][1]=178.6014;
  p[3][2][0]=919.7276;
  p[3][2][1]=-172.8446;
  p[3][3][0]=616.5578;
  p[3][3][1]=-347.7606;
  p[3][4][0]=313.9646;
  p[3][4][1]=-172.6736;
  p[3][5][0]=314.0006;
  p[3][5][1]=178.4304;

  p[1][0][0]=1224.873;
  p[1][0][1]=352.6879;
  p[1][1][0]=1527.851;
  p[1][1][1]=177.4779;
  p[1][2][0]=1527.732;
  p[1][2][1]=-173.6701;
  p[1][3][0]=1224.303;
  p[1][3][1]=-348.801;
  p[1][4][0]=921.3879;
  p[1][4][1]=-173.6471;
  p[1][5][0]=921.3839;
  p[1][5][1]=177.4549;

  p[11][0][0]=-904.0635;
  p[11][0][1]=-177.051;
  p[11][1][0]=-600.9554;
  p[11][1][1]=-352.213;
  p[11][2][0]=-601.0044;
  p[11][2][1]=-703.431;
  p[11][3][0]=-904.0333;
  p[11][3][1]=-878.454;
  p[11][4][0]=-1206.766;
  p[11][4][1]=-703.424;
  p[11][5][0]=-1206.844;
  p[11][5][1]=-352.22;

  p[14][0][0]=-297.6943;
  p[14][0][1]=-175.6188;
  p[14][1][0]=5.315837;
  p[14][1][1]=-350.7058;
  p[14][2][0]=5.157848;
  p[14][2][1]=-702.1018;
  p[14][3][0]=-297.7042;
  p[14][3][1]=-877.0388;
  p[14][4][0]=-600.7572;
  p[14][4][1]=-702.0688;
  p[14][5][0]=-600.6882;
  p[14][5][1]=-350.7388;

  p[13][0][0]=314.3028;
  p[13][0][1]=-175.2481;
  p[13][1][0]=617.264;
  p[13][1][1]=-350.2241;
  p[13][2][0]=617.0899;
  p[13][2][1]=-701.700;
  p[13][3][0]=314.1854;
  p[13][3][1]=-876.6871;
  p[13][4][0]=11.19794;
  p[13][4][1]=-701.6591;
  p[13][5][0]=11.25895;
  p[13][5][1]=-350.2651;

  p[16][0][0]=922.5141;
  p[16][0][1]=-175.0949;
  p[16][1][0]=1225.809;
  p[16][1][1]=-350.3329;
  p[16][2][0]=1225.91;
  p[16][2][1]=-701.1989;
  p[16][3][0]=923.0109;
  p[16][3][1]=-876.4109;
  p[16][4][0]=619.8021;
  p[16][4][1]=-701.2689;
  p[16][5][0]=619.6551;
  p[16][5][1]=-350.2629;

  p[15][0][0]=620.0246;
  p[15][0][1]=-703.531;
  p[15][1][0]=923.1095;
  p[15][1][1]=-878.518;
  p[15][2][0]=923.1674;
  p[15][2][1]=-1230.048;
  p[15][3][0]=620.1311;
  p[15][3][1]=-1405.198;
  p[15][4][0]=317.2905;
  p[15][4][1]=-1230.037;
  p[15][5][0]=317.3815;
  p[15][5][1]=-878.529;

  p[22][0][0]=14.00862;
  p[22][0][1]=-704.6452;
  p[22][1][0]=316.9342;
  p[22][1][1]=-879.6272;
  p[22][2][0]=316.9592;
  p[22][2][1]=-1230.817;
  p[22][3][0]=14.05039;
  p[22][3][1]=-1405.803;
  p[22][4][0]=-288.6738;
  p[22][4][1]=-1230.762;
  p[22][5][0]=-288.6788;
  p[22][5][1]=-879.6822;

  p[8][0][0]=-599.8333;
  p[8][0][1]=-704.3596;
  p[8][1][0]=-297.1782;
  p[8][1][1]=-879.2995;
  p[8][2][0]=-297.3792;
  p[8][2][1]=-1230.788;
  p[8][3][0]=-600.1864;
  p[8][3][1]=-1405.685;
  p[8][4][0]=-903.0912;
  p[8][4][1]=-1230.608;
  p[8][5][0]=-902.8402;
  p[8][5][1]=-879.4796;

  p[23][0][0]=9.536897;
  p[23][0][1]=352.1638;
  p[23][1][0]=311.489;
  p[23][1][1]=177.4148;
  p[23][2][0]=311.276;
  p[23][2][1]=-173.6332;
  p[23][3][0]=9.081099;
  p[23][3][1]=-348.3822;

  p[24][0][0]=2.469554;
  p[24][0][1]=352.7277;
  p[24][1][0]=-299.931;
  p[24][1][1]=177.9067;
  p[24][2][0]=-299.835;
  p[24][2][1]=-173.3314;
  p[24][3][0]=2.710454;
  p[24][3][1]=-348.1524;

  v[0][9][17][0]=(p[9][1][0]+p[17][5][0])/2.;
  v[0][9][17][1]=(p[9][1][1]+p[17][5][1])/2.;

  v[0][4][17][0]=(p[17][1][0]+p[4][5][0])/2.;
  v[0][4][17][1]=(p[17][1][1]+p[4][5][1])/2.;

  v[0][9][10][0]=(p[9][4][0]+p[10][0][0])/2.;
  v[0][9][10][1]=(p[9][4][1]+p[10][0][1])/2.;

  v[9][17][20][0]=(p[9][2][0]+p[17][4][0]+p[20][0][0])/3.;
  v[9][17][20][1]=(p[9][2][1]+p[17][4][1]+p[20][0][1])/3.;

  v[4][12][17][0]=(p[17][2][0]+p[4][4][0]+p[12][0][0])/3.;
  v[4][12][17][1]=(p[17][2][1]+p[4][4][1]+p[12][0][1])/3.;

  v[0][4][21][0]=(p[4][2][0]+p[0][4][0]+p[21][0][0])/2.;
  v[0][4][21][1]=(p[4][2][1]+p[0][4][1]+p[21][0][1])/2.;

  v[9][10][20][0]=(p[9][3][0]+p[10][1][0]+p[20][5][0])/3.;
  v[9][10][20][1]=(p[9][3][1]+p[10][1][1]+p[20][5][1])/3.;

  v[12][17][20][0]=(p[17][3][0]+p[20][1][0]+p[12][5][0])/3.;
  v[12][17][20][1]=(p[17][3][1]+p[20][1][1]+p[12][5][1])/3.;

  v[4][12][21][0]=(p[4][3][0]+p[12][1][0]+p[21][5][0])/3.;
  v[4][12][21][1]=(p[4][3][1]+p[12][1][1]+p[21][5][1])/3.;

  v[0][6][10][0]=(p[0][2][0]+p[10][4][0]+p[6][0][0])/2.;
  v[0][6][10][1]=(p[0][2][1]+p[10][4][1]+p[6][0][1])/2.;

  v[5][10][20][0]=(p[10][2][0]+p[20][4][0]+p[5][0][0])/3.;
  v[5][10][20][1]=(p[10][2][1]+p[20][4][1]+p[5][0][1])/3.;

  v[0][12][20][0]=(p[20][2][0]+p[12][4][0]+p[0][0][0])/2.;
  v[0][12][20][1]=(p[20][2][1]+p[12][4][1]+p[0][0][1])/2.;

  v[3][12][21][0]=(p[12][2][0]+p[21][4][0]+p[3][0][0])/3.;
  v[3][12][21][1]=(p[12][2][1]+p[21][4][1]+p[3][0][1])/3.;

  v[0][1][21][0]=(p[21][2][0]+p[0][4][0]+p[1][0][0])/2.;
  v[0][1][21][1]=(p[21][2][1]+p[0][4][1]+p[1][0][1])/2.;

  v[5][6][10][0]=(p[10][3][0]+p[6][1][0]+p[5][5][0])/3.;
  v[5][6][10][1]=(p[10][3][1]+p[6][1][1]+p[5][5][1])/3.;

  v[0][5][20][0]=(p[20][3][0]+p[5][1][0]+p[0][5][0])/2.;
  v[0][5][20][1]=(p[20][3][1]+p[5][1][1]+p[0][5][1])/2.;

  v[0][3][12][0]=(p[12][3][0]+p[0][1][0]+p[3][5][0])/2.;
  v[0][3][12][1]=(p[12][3][1]+p[0][1][1]+p[3][5][1])/2.;

  v[1][3][21][0]=(p[21][3][0]+p[3][1][0]+p[1][5][0])/3.;
  v[1][3][21][1]=(p[21][3][1]+p[3][1][1]+p[1][5][1])/3.;

  v[5][6][11][0]=(p[6][2][0]+p[5][4][0]+p[11][0][0])/3.;
  v[5][6][11][1]=(p[6][2][1]+p[5][4][1]+p[11][0][1])/3.;

  v[0][5][14][0]=(p[5][2][0]+p[0][4][0]+p[14][0][0])/2.;
  v[0][5][14][1]=(p[5][2][1]+p[0][4][1]+p[14][0][1])/2.;

  v[0][3][13][0]=(p[0][2][0]+p[3][4][0]+p[13][0][0])/2.;
  v[0][3][13][1]=(p[0][2][1]+p[3][4][1]+p[13][0][1])/2.;

  v[1][3][16][0]=(p[3][2][0]+p[1][4][0]+p[16][0][0])/3.;
  v[1][3][16][1]=(p[3][2][1]+p[1][4][1]+p[16][0][1])/3.;

  v[0][6][11][0]=(p[6][3][0]+p[0][1][0]+p[11][5][0])/2.;
  v[0][6][11][1]=(p[6][3][1]+p[0][1][1]+p[11][5][1])/2.;

  v[5][11][14][0]=(p[5][3][0]+p[11][1][0]+p[14][5][0])/3.;
  v[5][11][14][1]=(p[5][3][1]+p[11][1][1]+p[14][5][1])/3.;

  v[0][13][14][0]=(p[0][3][0]+p[14][1][0]+p[13][5][0])/2.;
  v[0][13][14][1]=(p[0][3][1]+p[14][1][1]+p[13][5][1])/2.;

  v[3][13][16][0]=(p[3][3][0]+p[13][1][0]+p[16][5][0])/3.;
  v[3][13][16][1]=(p[3][3][1]+p[13][1][1]+p[16][5][1])/3.;

  v[0][1][16][0]=(p[1][3][0]+p[16][1][0]+p[0][5][0])/2.;
  v[0][1][16][1]=(p[1][3][1]+p[16][1][1]+p[0][5][1])/2.;

  v[8][11][14][0]=(p[11][2][0]+p[14][4][0]+p[8][0][0])/3.;
  v[8][11][14][1]=(p[11][2][1]+p[14][4][1]+p[8][0][1])/3.;

  v[13][14][22][0]=(p[14][2][0]+p[13][4][0]+p[22][0][0])/3.;
  v[13][14][22][1]=(p[14][2][1]+p[13][4][1]+p[22][0][1])/3.;

  v[13][15][16][0]=(p[13][2][0]+p[16][4][0]+p[15][0][0])/3.;
  v[13][15][16][1]=(p[13][2][1]+p[16][4][1]+p[15][0][1])/3.;

  v[0][8][11][0]=(p[11][3][0]+p[0][1][0]+p[8][5][0])/2.;
  v[0][8][11][1]=(p[11][3][1]+p[0][1][1]+p[8][5][1])/2.;

  v[8][14][22][0]=(p[14][3][0]+p[8][1][0]+p[22][5][0])/3.;
  v[8][14][22][1]=(p[14][3][1]+p[8][1][1]+p[22][5][1])/3.;

  v[13][15][22][0]=(p[13][3][0]+p[22][1][0]+p[15][5][0])/3.;
  v[13][15][22][1]=(p[13][3][1]+p[22][1][1]+p[15][5][1])/3.;

  v[0][15][16][0]=(p[16][3][0]+p[15][1][0]+p[0][5][0])/2.;
  v[0][15][16][1]=(p[16][3][1]+p[15][1][1]+p[0][5][1])/2.;

  v[0][8][22][0]=(p[8][2][0]+p[22][4][0]+p[0][0][0])/2.;
  v[0][8][22][1]=(p[8][2][1]+p[22][4][1]+p[0][0][1])/2.;

  v[0][15][22][0]=(p[22][2][0]+p[15][4][0]+p[0][0][0])/2.;
  v[0][15][22][1]=(p[22][2][1]+p[15][4][1]+p[0][0][1])/2.;

  fakemirror[9][17]=25;
  fakemirror[9][10]=27;
  fakemirror[17][9]=26;
  fakemirror[17][4]=25;
  fakemirror[4][17]=27;
  fakemirror[4][21]=26;
  fakemirror[10][9]=26;
  fakemirror[10][6]=25;
  fakemirror[21][4]=25;
  fakemirror[21][1]=26;
  fakemirror[6][10]=27;
  fakemirror[6][11]=25;
  fakemirror[1][21]=25;
  fakemirror[1][16]=27;
  fakemirror[11][6]=26;
  fakemirror[11][8]=25;
  fakemirror[16][1]=25;
  fakemirror[16][15]=26;
  fakemirror[8][11]=27;
  fakemirror[8][22]=25;
  fakemirror[22][8]=26;
  fakemirror[22][15]=25;
  fakemirror[15][16]=25;
  fakemirror[15][22]=26;
  fakemirror[12][20]=25;
  fakemirror[12][3]=25;
  fakemirror[3][12]=25;
  fakemirror[3][13]=25;
  fakemirror[13][3]=25;
  fakemirror[13][14]=25;
  fakemirror[14][13]=25;
  fakemirror[14][5]=25;
  fakemirror[5][14]=25;
  fakemirror[5][20]=25;
  fakemirror[20][12]=25;
  fakemirror[20][5]=25;

  ac[9][25]=((p[9][0][1]+p[0][4][1])-(p[9][1][1]+p[0][3][1]))/((p[9][0][0]+p[0][4][0])-(p[9][1][0]+p[0][3][0]));
  bcost[9][25]=((p[9][1][1]+p[0][3][1])*(p[9][0][0]+p[0][4][0])-(p[9][1][0]+p[0][3][0])*(p[9][0][1]+p[0][4][1]))/(2.*((p[9][0][0]+p[0][4][0])-(p[9][1][0]+p[0][3][0])));
  ac[9][17]=0.;
  bcost[9][17]=(p[9][1][0]+p[17][5][0])/4.+(p[9][2][0]+p[17][4][0])/4.;
  ac[9][20]=((p[9][2][1]+p[20][0][1])/2.-(p[9][3][1]+p[20][5][1])/2.)/((p[9][2][0]+p[20][0][0])/2.-(p[9][3][0]+p[20][5][0])/2.);
  bcost[9][20]=((p[9][3][1]+p[20][5][1])*(p[9][2][0]+p[20][0][0])/4.-(p[9][3][0]+p[20][5][0])*(p[9][2][1]+p[20][0][1])/4.)/((p[9][2][0]+p[20][0][0])/2.-(p[9][3][0]+p[20][5][0])/2.);
  ac[9][10]=((p[9][3][1]+p[10][1][1])/2.-(p[9][4][1]+p[10][0][1])/2.)/((p[9][3][0]+p[10][1][0])/2.-(p[9][4][0]+p[10][0][0])/2.);
  bcost[9][10]=((p[9][4][1]+p[10][0][1])*(p[9][3][0]+p[10][1][0])/4.-(p[9][4][0]+p[10][0][0])*(p[9][3][1]+p[10][1][1])/4.)/((p[9][3][0]+p[10][1][0])/2.-(p[9][4][0]+p[10][0][0])/2.);
  ac[9][26]=0.;
  bcost[9][26]=(p[9][4][0]+p[0][2][0])/4.+(p[9][5][0]+p[0][1][0])/4.;
  ac[9][27]=((p[9][5][1]+p[0][3][1])-(p[9][0][1]+p[0][2][1]))/((p[9][5][0]+p[0][3][0])-(p[9][0][0]+p[0][2][0]));
  bcost[9][27]=((p[9][0][1]+p[0][2][1])*(p[9][5][0]+p[0][3][0])-(p[9][0][0]+p[0][2][0])*(p[9][5][1]+p[0][3][1]))/(2.*((p[9][5][0]+p[0][3][0])-(p[9][0][0]+p[0][2][0])));

  ac[17][25]=((p[17][0][1]+p[0][4][1])-(p[17][1][1]+p[0][3][1]))/((p[17][0][0]+p[0][4][0])-(p[17][1][0]+p[0][3][0]));
  bcost[17][25]=((p[17][1][1]+p[0][3][1])*(p[17][0][0]+p[0][4][0])-(p[17][1][0]+p[0][3][0])*(p[17][0][1]+p[0][4][1]))/(2.*((p[17][0][0]+p[0][4][0])-(p[17][1][0]+p[0][3][0])));
  ac[17][4]=0.;
  bcost[17][4]=(p[17][1][0]+p[4][5][0])/4.+(p[17][2][0]+p[4][4][0])/4.;
  ac[17][12]=((p[17][2][1]+p[12][0][1])/2.-(p[17][3][1]+p[12][5][1])/2.)/((p[17][2][0]+p[12][0][0])/2.-(p[17][3][0]+p[12][5][0])/2.);
  bcost[17][12]=((p[17][3][1]+p[12][5][1])*(p[17][2][0]+p[12][0][0])/4.-(p[17][3][0]+p[12][5][0])*(p[17][2][1]+p[12][0][1])/4.)/((p[17][2][0]+p[12][0][0])/2.-(p[17][3][0]+p[12][5][0])/2.);
  ac[17][20]=((p[17][3][1]+p[20][1][1])/2.-(p[17][4][1]+p[20][0][1])/2.)/((p[17][3][0]+p[20][1][0])/2.-(p[17][4][0]+p[20][0][0])/2.);
  bcost[17][20]=((p[17][4][1]+p[20][0][1])*(p[17][3][0]+p[20][1][0])/4.-(p[17][4][0]+p[20][0][0])*(p[17][3][1]+p[20][1][1])/4.)/((p[17][3][0]+p[20][1][0])/2.-(p[17][4][0]+p[20][0][0])/2.);
  ac[17][9]=0.;
  bcost[17][9]=(p[17][4][0]+p[9][2][0])/4.+(p[17][5][0]+p[9][1][0])/4.;
  ac[17][26]=((p[17][5][1]+p[0][3][1])-(p[17][0][1]+p[0][2][1]))/((p[17][5][0]+p[0][3][0])-(p[17][0][0]+p[0][2][0]));
  bcost[17][26]=((p[17][0][1]+p[0][2][1])*(p[17][5][0]+p[0][3][0])/4.-(p[17][0][0]+p[0][2][0])*(p[17][5][1]+p[0][3][1])/4.)/((p[17][5][0]+p[0][3][0])/2.-(p[17][0][0]+p[0][2][0])/2.);

  ac[4][25]=((p[4][0][1]+p[0][4][1])-(p[4][1][1]+p[0][3][1]))/((p[4][0][0]+p[0][4][0])-(p[4][1][0]+p[0][3][0]));
  bcost[4][25]=((p[4][1][1]+p[0][3][1])*(p[4][0][0]+p[0][4][0])-(p[4][1][0]+p[0][3][0])*(p[4][0][1]+p[0][4][1]))/(2.*((p[4][0][0]+p[0][4][0])-(p[4][1][0]+p[0][3][0])));
  ac[4][26]=0.;
  bcost[4][26]=(p[4][1][0]+p[0][5][0])/4.+(p[4][2][0]+p[0][4][0])/4.;
  ac[4][21]=((p[4][2][1]+p[21][0][1])/2.-(p[4][3][1]+p[21][5][1])/2.)/((p[4][2][0]+p[21][0][0])/2.-(p[4][3][0]+p[21][5][0])/2.);
  bcost[4][21]=((p[4][3][1]+p[21][5][1])*(p[4][2][0]+p[21][0][0])/4.-(p[4][3][0]+p[21][5][0])*(p[4][2][1]+p[21][0][1])/4.)/((p[4][2][0]+p[21][0][0])/2.-(p[4][3][0]+p[21][5][0])/2.);
  ac[4][12]=((p[4][3][1]+p[12][1][1])/2.-(p[4][4][1]+p[12][0][1])/2.)/((p[4][3][0]+p[12][1][0])/2.-(p[4][4][0]+p[12][0][0])/2.);
  bcost[4][12]=((p[4][4][1]+p[12][0][1])*(p[4][3][0]+p[12][1][0])/4.-(p[4][4][0]+p[12][0][0])*(p[4][3][1]+p[12][1][1])/4.)/((p[4][3][0]+p[12][1][0])/2.-(p[4][4][0]+p[12][0][0])/2.);
  ac[4][17]=0.;
  bcost[4][17]=(p[4][4][0]+p[17][2][0])/4.+(p[4][5][0]+p[17][1][0])/4.;
  ac[4][27]=((p[4][5][1]+p[0][3][1])-(p[4][0][1]+p[0][2][1]))/((p[4][5][0]+p[0][3][0])-(p[4][0][0]+p[0][2][0]));
  bcost[4][27]=((p[4][0][1]+p[0][2][1])*(p[4][5][0]+p[0][3][0])-(p[4][0][0]+p[0][2][0])*(p[4][5][1]+p[0][3][1]))/(2.*((p[4][5][0]+p[0][3][0])-(p[4][0][0]+p[0][2][0])));

  ac[10][9]=((p[10][0][1]+p[9][4][1])/2.-(p[10][1][1]+p[9][3][1])/2.)/((p[10][0][0]+p[9][4][0])/2.-(p[10][1][0]+p[9][3][0])/2.);
  bcost[10][9]=((p[10][1][1]+p[9][3][1])*(p[10][0][0]+p[9][4][0])/4.-(p[10][1][0]+p[9][3][0])*(p[10][0][1]+p[9][4][1])/4.)/((p[10][0][0]+p[9][4][0])/2.-(p[10][1][0]+p[9][3][0])/2.);
  ac[10][20]=0.;
  bcost[10][20]=(p[10][1][0]+p[20][5][0])/4.+(p[10][2][0]+p[20][4][0])/4.;
  ac[10][5]=((p[10][2][1]+p[5][0][1])/2.-(p[10][3][1]+p[5][5][1])/2.)/((p[10][2][0]+p[5][0][0])/2.-(p[10][3][0]+p[5][5][0])/2.);
  bcost[10][5]=((p[10][3][1]+p[5][5][1])*(p[10][2][0]+p[5][0][0])/4.-(p[10][3][0]+p[5][5][0])*(p[10][2][1]+p[5][0][1])/4.)/((p[10][2][0]+p[5][0][0])/2.-(p[10][3][0]+p[5][5][0])/2.);
  ac[10][6]=((p[10][3][1]+p[6][1][1])/2.-(p[10][4][1]+p[6][0][1])/2.)/((p[10][3][0]+p[6][1][0])/2.-(p[10][4][0]+p[6][0][0])/2.);
  bcost[10][6]=((p[10][4][1]+p[6][0][1])*(p[10][3][0]+p[6][1][0])/4.-(p[10][4][0]+p[6][0][0])*(p[10][3][1]+p[6][1][1])/4.)/((p[10][3][0]+p[6][1][0])/2.-(p[10][4][0]+p[6][0][0])/2.);
  ac[10][25]=0.;
  bcost[10][25]=(p[10][4][0]+p[0][2][0])/4.+(p[10][5][0]+p[0][1][0])/4.;
  ac[10][26]=((p[10][5][1]+p[0][3][1])-(p[10][0][1]+p[0][2][1]))/((p[10][5][0]+p[0][3][0])-(p[10][0][0]+p[0][2][0]));
  bcost[10][26]=((p[10][0][1]+p[0][2][1])*(p[10][5][0]+p[0][3][0])-(p[10][0][0]+p[0][2][0])*(p[10][5][1]+p[0][3][1]))/(2.*((p[10][5][0]+p[0][3][0])-(p[10][0][0]+p[0][2][0])));

  ac[20][17]=((p[20][0][1]+p[17][4][1])/2.-(p[20][1][1]+p[17][3][1])/2.)/((p[20][0][0]+p[17][4][0])/2.-(p[20][1][0]+p[17][3][0])/2.);
  bcost[20][17]=((p[20][1][1]+p[17][3][1])*(p[20][0][0]+p[17][4][0])/4.-(p[20][1][0]+p[17][3][0])*(p[20][0][1]+p[17][4][1])/4.)/((p[20][0][0]+p[17][4][0])/2.-(p[20][1][0]+p[17][3][0])/2.);
  ac[20][12]=0.;
  bcost[20][12]=(p[20][1][0]+p[12][5][0])/4.+(p[20][2][0]+p[12][4][0])/4.;
  ac[20][25]=((p[20][2][1]+p[0][0][1])-(p[20][3][1]+p[0][5][1]))/((p[20][2][0]+p[0][0][0])-(p[20][3][0]+p[0][5][0]));
  bcost[20][25]=((p[20][3][1]+p[0][5][1])*(p[20][2][0]+p[0][0][0])-(p[20][3][0]+p[0][5][0])*(p[20][2][1]+p[0][0][1]))/(2.*((p[20][2][0]+p[0][0][0])-(p[20][3][0]+p[0][5][0])));
  ac[20][5]=((p[20][3][1]+p[5][1][1])/2.-(p[20][4][1]+p[5][0][1])/2.)/((p[20][3][0]+p[5][1][0])/2.-(p[20][4][0]+p[5][0][0])/2.);
  bcost[20][5]=((p[20][4][1]+p[5][0][1])*(p[20][3][0]+p[5][1][0])/4.-(p[20][4][0]+p[5][0][0])*(p[20][3][1]+p[5][1][1])/4.)/((p[20][3][0]+p[5][1][0])/2.-(p[20][4][0]+p[5][0][0])/2.);
  ac[20][10]=0.;
  bcost[20][10]=(p[20][4][0]+p[10][2][0])/4.+(p[20][5][0]+p[10][1][0])/4.;
  ac[20][9]=((p[20][5][1]+p[9][3][1])/2.-(p[20][0][1]+p[9][2][1])/2.)/((p[20][5][0]+p[9][3][0])/2.-(p[20][0][0]+p[9][2][0])/2.);
  bcost[20][9]=((p[20][0][1]+p[9][2][1])*(p[20][5][0]+p[9][3][0])/4.-(p[20][0][0]+p[9][2][0])*(p[20][5][1]+p[9][3][1])/4.)/((p[20][5][0]+p[9][3][0])/2.-(p[20][0][0]+p[9][2][0])/2.);

  ac[12][4]=((p[12][0][1]+p[4][4][1])/2.-(p[12][1][1]+p[4][3][1])/2.)/((p[12][0][0]+p[4][4][0])/2.-(p[12][1][0]+p[4][3][0])/2.);
  bcost[12][4]=((p[12][1][1]+p[4][3][1])*(p[12][0][0]+p[4][4][0])/4.-(p[12][1][0]+p[4][3][0])*(p[12][0][1]+p[4][4][1])/4.)/((p[12][0][0]+p[4][4][0])/2.-(p[12][1][0]+p[4][3][0])/2.);
  ac[12][21]=0.;
  bcost[12][21]=(p[12][1][0]+p[21][5][0])/4.+(p[12][2][0]+p[21][4][0])/4.;
  ac[12][3]=((p[12][2][1]+p[3][0][1])/2.-(p[12][3][1]+p[3][5][1])/2.)/((p[12][2][0]+p[3][0][0])/2.-(p[12][3][0]+p[3][5][0])/2.);
  bcost[12][3]=((p[12][3][1]+p[3][5][1])*(p[12][2][0]+p[3][0][0])/4.-(p[12][3][0]+p[3][5][0])*(p[12][2][1]+p[3][0][1])/4.)/((p[12][2][0]+p[3][0][0])/2.-(p[12][3][0]+p[3][5][0])/2.);
  ac[12][25]=((p[12][3][1]+p[0][1][1])-(p[12][4][1]+p[0][0][1]))/((p[12][3][0]+p[0][1][0])-(p[12][4][0]+p[0][0][0]));
  bcost[12][25]=((p[12][4][1]+p[0][0][1])*(p[12][3][0]+p[0][1][0])-(p[12][4][0]+p[0][0][0])*(p[12][3][1]+p[0][1][1]))/(2.*((p[12][3][0]+p[0][1][0])-(p[12][4][0]+p[0][0][0])));
  ac[12][20]=0.;
  bcost[12][20]=(p[12][4][0]+p[20][2][0])/4.+(p[12][5][0]+p[20][1][0])/4.;
  ac[12][17]=((p[12][5][1]+p[17][3][1])/2.-(p[12][0][1]+p[17][2][1])/2.)/((p[12][5][0]+p[17][3][0])/2.-(p[12][0][0]+p[17][2][0])/2.);
  bcost[12][17]=((p[12][0][1]+p[17][2][1])*(p[12][5][0]+p[17][3][0])/4.-(p[12][0][0]+p[17][2][0])*(p[12][5][1]+p[17][3][1])/4.)/((p[12][5][0]+p[17][3][0])/2.-(p[12][0][0]+p[17][2][0])/2.);

  ac[21][25]=((p[21][0][1]+p[0][4][1])-(p[21][1][1]+p[0][3][1]))/((p[21][0][0]+p[0][4][0])-(p[21][1][0]+p[0][3][0]));
  bcost[21][25]=((p[21][1][1]+p[0][3][1])*(p[21][0][0]+p[0][4][0])-(p[21][1][0]+p[0][3][0])*(p[21][0][1]+p[0][4][1]))/(2.*((p[21][0][0]+p[0][4][0])-(p[21][1][0]+p[0][3][0])));
  ac[21][26]=0.;
  bcost[21][26]=(p[21][1][0]+p[0][5][0])/4.+(p[21][2][0]+p[0][4][0])/4.;
  ac[21][1]=((p[21][2][1]+p[1][0][1])/2.-(p[21][3][1]+p[1][5][1])/2.)/((p[21][2][0]+p[1][0][0])/2.-(p[21][3][0]+p[1][5][0])/2.);
  bcost[21][1]=((p[21][3][1]+p[1][5][1])*(p[21][2][0]+p[1][0][0])/4.-(p[21][3][0]+p[1][5][0])*(p[21][2][1]+p[1][0][1])/4.)/((p[21][2][0]+p[1][0][0])/2.-(p[21][3][0]+p[1][5][0])/2.);
  ac[21][3]=((p[21][3][1]+p[3][1][1])/2.-(p[21][4][1]+p[3][0][1])/2.)/((p[21][3][0]+p[3][1][0])/2.-(p[21][4][0]+p[3][0][0])/2.);
  bcost[21][3]=((p[21][4][1]+p[3][0][1])*(p[21][3][0]+p[3][1][0])/4.-(p[21][4][0]+p[3][0][0])*(p[21][3][1]+p[3][1][1])/4.)/((p[21][3][0]+p[3][1][0])/2.-(p[21][4][0]+p[3][0][0])/2.);
  ac[21][12]=0.;
  bcost[21][12]=(p[21][4][0]+p[12][2][0])/4.+(p[21][5][0]+p[12][1][0])/4.;
  ac[21][4]=((p[21][5][1]+p[4][3][1])/2.-(p[21][0][1]+p[4][2][1])/2.)/((p[21][5][0]+p[4][3][0])/2.-(p[21][0][0]+p[4][2][0])/2.);
  bcost[21][4]=((p[21][0][1]+p[4][2][1])*(p[21][5][0]+p[4][3][0])/4.-(p[21][0][0]+p[4][2][0])*(p[21][5][1]+p[4][3][1])/4.)/((p[21][5][0]+p[4][3][0])/2.-(p[21][0][0]+p[4][2][0])/2.);

  ac[6][10]=((p[6][0][1]+p[10][4][1])/2.-(p[6][1][1]+p[10][3][1])/2.)/((p[6][0][0]+p[10][4][0])/2.-(p[6][1][0]+p[10][3][0])/2.);
  bcost[6][10]=((p[6][1][1]+p[10][3][1])*(p[6][0][0]+p[10][4][0])/4.-(p[6][1][0]+p[10][3][0])*(p[6][0][1]+p[10][4][1])/4.)/((p[6][0][0]+p[10][4][0])/2.-(p[6][1][0]+p[10][3][0])/2.);
  ac[6][5]=0.;
  bcost[6][5]=(p[6][1][0]+p[5][5][0])/4.+(p[6][2][0]+p[5][4][0])/4.;
  ac[6][11]=((p[6][2][1]+p[11][0][1])/2.-(p[6][3][1]+p[11][5][1])/2.)/((p[6][2][0]+p[11][0][0])/2.-(p[6][3][0]+p[11][5][0])/2.);
  bcost[6][11]=((p[6][3][1]+p[11][5][1])*(p[6][2][0]+p[11][0][0])/4.-(p[6][3][0]+p[11][5][0])*(p[6][2][1]+p[11][0][1])/4.)/((p[6][2][0]+p[11][0][0])/2.-(p[6][3][0]+p[11][5][0])/2.);
  ac[6][25]=((p[6][3][1]+p[0][1][1])-(p[6][4][1]+p[0][0][1]))/((p[6][3][0]+p[0][1][0])-(p[6][4][0]+p[0][0][0]));
  bcost[6][25]=((p[6][4][1]+p[0][0][1])*(p[6][3][0]+p[0][1][0])-(p[6][4][0]+p[0][0][0])*(p[6][3][1]+p[0][1][1]))/(2.*((p[6][3][0]+p[0][1][0])-(p[6][4][0]+p[0][0][0])));
  ac[6][26]=0.;
  bcost[6][26]=(p[6][4][0]+p[0][2][0])/4.+(p[6][5][0]+p[0][1][0])/4.;
  ac[6][27]=((p[6][5][1]+p[0][3][1])-(p[6][0][1]+p[0][2][1]))/((p[6][5][0]+p[0][3][0])-(p[6][0][0]+p[0][2][0]));
  bcost[6][27]=((p[6][0][1]+p[0][2][1])*(p[6][5][0]+p[0][3][0])-(p[6][0][0]+p[0][2][0])*(p[6][5][1]+p[0][3][1]))/(2.*((p[6][5][0]+p[0][3][0])-(p[6][0][0]+p[0][2][0])));

  ac[5][20]=((p[5][0][1]+p[20][4][1])/2.-(p[5][1][1]+p[20][3][1])/2.)/((p[5][0][0]+p[20][4][0])/2.-(p[5][1][0]+p[20][3][0])/2.);
  bcost[5][20]=((p[5][1][1]+p[20][3][1])*(p[5][0][0]+p[20][4][0])/4.-(p[5][1][0]+p[20][3][0])*(p[5][0][1]+p[20][4][1])/4.)/((p[5][0][0]+p[20][4][0])/2.-(p[5][1][0]+p[20][3][0])/2.);
  ac[5][25]=0.;
  bcost[5][25]=(p[5][1][0]+p[0][5][0])/4.+(p[5][2][0]+p[0][4][0])/4.;
  ac[5][14]=((p[5][2][1]+p[14][0][1])/2.-(p[5][3][1]+p[14][5][1])/2.)/((p[5][2][0]+p[14][0][0])/2.-(p[5][3][0]+p[14][5][0])/2.);
  bcost[5][14]=((p[5][3][1]+p[14][5][1])*(p[5][2][0]+p[14][0][0])/4.-(p[5][3][0]+p[14][5][0])*(p[5][2][1]+p[14][0][1])/4.)/((p[5][2][0]+p[14][0][0])/2.-(p[5][3][0]+p[14][5][0])/2.);
  ac[5][11]=((p[5][3][1]+p[11][1][1])/2.-(p[5][4][1]+p[11][0][1])/2.)/((p[5][3][0]+p[11][1][0])/2.-(p[5][4][0]+p[11][0][0])/2.);
  bcost[5][11]=((p[5][4][1]+p[11][0][1])*(p[5][3][0]+p[11][1][0])/4.-(p[5][4][0]+p[11][0][0])*(p[5][3][1]+p[11][1][1])/4.)/((p[5][3][0]+p[11][1][0])/2.-(p[5][4][0]+p[11][0][0])/2.);
  ac[5][6]=0.;
  bcost[5][6]=(p[5][4][0]+p[6][2][0])/4.+(p[5][5][0]+p[6][1][0])/4.;
  ac[5][10]=((p[5][5][1]+p[10][3][1])/2.-(p[5][0][1]+p[10][2][1])/2.)/((p[5][5][0]+p[10][3][0])/2.-(p[5][0][0]+p[10][2][0])/2.);
  bcost[5][10]=((p[5][0][1]+p[10][2][1])*(p[5][5][0]+p[10][3][0])/4.-(p[5][0][0]+p[10][2][0])*(p[5][5][1]+p[10][3][1])/4.)/((p[5][5][0]+p[10][3][0])/2.-(p[5][0][0]+p[10][2][0])/2.);

  ac[3][21]=((p[3][0][1]+p[21][4][1])/2.-(p[3][1][1]+p[21][3][1])/2.)/((p[3][0][0]+p[21][4][0])/2.-(p[3][1][0]+p[21][3][0])/2.);
  bcost[3][21]=((p[3][1][1]+p[21][3][1])*(p[3][0][0]+p[21][4][0])/4.-(p[3][1][0]+p[21][3][0])*(p[3][0][1]+p[21][4][1])/4.)/((p[3][0][0]+p[21][4][0])/2.-(p[3][1][0]+p[21][3][0])/2.);
  ac[3][1]=0.;
  bcost[3][1]=(p[3][1][0]+p[1][5][0])/4.+(p[3][2][0]+p[1][4][0])/4.;
  ac[3][16]=((p[3][2][1]+p[16][0][1])/2.-(p[3][3][1]+p[16][5][1])/2.)/((p[3][2][0]+p[16][0][0])/2.-(p[3][3][0]+p[16][5][0])/2.);
  bcost[3][16]=((p[3][3][1]+p[16][5][1])*(p[3][2][0]+p[16][0][0])/4.-(p[3][3][0]+p[16][5][0])*(p[3][2][1]+p[16][0][1])/4.)/((p[3][2][0]+p[16][0][0])/2.-(p[3][3][0]+p[16][5][0])/2.);
  ac[3][13]=((p[3][3][1]+p[13][1][1])/2.-(p[3][4][1]+p[13][0][1])/2.)/((p[3][3][0]+p[13][1][0])/2.-(p[3][4][0]+p[13][0][0])/2.);
  bcost[3][13]=((p[3][4][1]+p[13][0][1])*(p[3][3][0]+p[13][1][0])/4.-(p[3][4][0]+p[13][0][0])*(p[3][3][1]+p[13][1][1])/4.)/((p[3][3][0]+p[13][1][0])/2.-(p[3][4][0]+p[13][0][0])/2.);
  ac[3][25]=0.;
  bcost[3][25]=(p[3][4][0]+p[0][2][0])/4.+(p[3][5][0]+p[0][1][0])/4.;
  ac[3][12]=((p[3][5][1]+p[12][3][1])/2.-(p[3][0][1]+p[12][2][1])/2.)/((p[3][5][0]+p[12][3][0])/2.-(p[3][0][0]+p[12][2][0])/2.);
  bcost[3][12]=((p[3][0][1]+p[12][2][1])*(p[3][5][0]+p[12][3][0])/4.-(p[3][0][0]+p[12][2][0])*(p[3][5][1]+p[12][3][1])/4.)/((p[3][5][0]+p[12][3][0])/2.-(p[3][0][0]+p[12][2][0])/2.);

  ac[1][25]=((p[1][0][1]+p[0][4][1])-(p[1][1][1]+p[0][3][1]))/((p[1][0][0]+p[0][4][0])-(p[1][1][0]+p[0][3][0]));
  bcost[1][25]=((p[1][1][1]+p[0][3][1])*(p[1][0][0]+p[0][4][0])-(p[1][1][0]+p[0][3][0])*(p[1][0][1]+p[0][4][1]))/(2.*((p[1][0][0]+p[0][4][0])-(p[1][1][0]+p[0][3][0])));
  ac[1][26]=0.;
  bcost[1][26]=(p[1][1][0]+p[0][5][0])/4.+(p[1][2][0]+p[0][4][0])/4.;
  ac[1][27]=((p[1][2][1]+p[0][0][1])-(p[1][3][1]+p[0][5][1]))/((p[1][2][0]+p[0][0][0])-(p[1][3][0]+p[0][5][0]));
  bcost[1][27]=((p[1][3][1]+p[0][5][1])*(p[1][2][0]+p[0][0][0])-(p[1][3][0]+p[0][5][0])*(p[1][2][1]+p[0][0][1]))/(2.*((p[1][2][0]+p[0][0][0])-(p[1][3][0]+p[0][5][0])));
  ac[1][16]=((p[1][3][1]+p[16][1][1])/2.-(p[1][4][1]+p[16][0][1])/2.)/((p[1][3][0]+p[16][1][0])/2.-(p[1][4][0]+p[16][0][0])/2.);
  bcost[1][16]=((p[1][4][1]+p[16][0][1])*(p[1][3][0]+p[16][1][0])/4.-(p[1][4][0]+p[16][0][0])*(p[1][3][1]+p[16][1][1])/4.)/((p[1][3][0]+p[16][1][0])/2.-(p[1][4][0]+p[16][0][0])/2.);
  ac[1][3]=0.;
  bcost[1][3]=(p[1][4][0]+p[3][2][0])/4.+(p[1][5][0]+p[3][1][0])/4.;
  ac[1][21]=((p[1][5][1]+p[21][3][1])/2.-(p[1][0][1]+p[21][2][1])/2.)/((p[1][5][0]+p[21][3][0])/2.-(p[1][0][0]+p[21][2][0])/2.);
  bcost[1][21]=((p[1][0][1]+p[21][2][1])*(p[1][5][0]+p[21][3][0])/4.-(p[1][0][0]+p[21][2][0])*(p[1][5][1]+p[21][3][1])/4.)/((p[1][5][0]+p[21][3][0])/2.-(p[1][0][0]+p[21][2][0])/2.);

  ac[11][5]=((p[11][0][1]+p[5][4][1])/2.-(p[11][1][1]+p[5][3][1])/2.)/((p[11][0][0]+p[5][4][0])/2.-(p[11][1][0]+p[5][3][0])/2.);
  bcost[11][5]=((p[11][1][1]+p[5][3][1])*(p[11][0][0]+p[5][4][0])/4.-(p[11][1][0]+p[5][3][0])*(p[11][0][1]+p[5][4][1])/4.)/((p[11][0][0]+p[5][4][0])/2.-(p[11][1][0]+p[5][3][0])/2.);
  ac[11][14]=0.;
  bcost[11][14]=(p[11][1][0]+p[14][5][0])/4.+(p[11][2][0]+p[14][4][0])/4.;
  ac[11][8]=((p[11][2][1]+p[8][0][1])/2.-(p[11][3][1]+p[8][5][1])/2.)/((p[11][2][0]+p[8][0][0])/2.-(p[11][3][0]+p[8][5][0])/2.);
  bcost[11][8]=((p[11][3][1]+p[8][5][1])*(p[11][2][0]+p[8][0][0])/4.-(p[11][3][0]+p[8][5][0])*(p[11][2][1]+p[8][0][1])/4.)/((p[11][2][0]+p[8][0][0])/2.-(p[11][3][0]+p[8][5][0])/2.);
  ac[11][25]=((p[11][3][1]+p[0][1][1])-(p[11][4][1]+p[0][0][1]))/((p[11][3][0]+p[0][1][0])-(p[11][4][0]+p[0][0][0]));
  bcost[11][25]=((p[11][4][1]+p[0][0][1])*(p[11][3][0]+p[0][1][0])-(p[11][4][0]+p[0][0][0])*(p[11][3][1]+p[0][1][1]))/(2.*((p[11][3][0]+p[0][1][0])-(p[11][4][0]+p[0][0][0])));
  ac[11][26]=0.;
  bcost[11][26]=(p[11][4][0]+p[0][2][0])/4.+(p[11][5][0]+p[0][1][0])/4.;
  ac[11][6]=((p[11][5][1]+p[6][3][1])/2.-(p[11][0][1]+p[6][2][1])/2.)/((p[11][5][0]+p[6][3][0])/2.-(p[11][0][0]+p[6][2][0])/2.);
  bcost[11][6]=((p[11][0][1]+p[6][2][1])*(p[11][5][0]+p[6][3][0])/4.-(p[11][0][0]+p[6][2][0])*(p[11][5][1]+p[6][3][1])/4.)/((p[11][5][0]+p[6][3][0])/2.-(p[11][0][0]+p[6][2][0])/2.);

  ac[14][25]=((p[14][0][1]+p[0][4][1])-(p[14][1][1]+p[0][3][1]))/((p[14][0][0]+p[0][4][0])-(p[14][1][0]+p[0][3][0]));
  bcost[14][25]=((p[14][1][1]+p[0][3][1])*(p[14][0][0]+p[0][4][0])-(p[14][1][0]+p[0][3][0])*(p[14][0][1]+p[0][4][1]))/(2.*((p[14][0][0]+p[0][4][0])-(p[14][1][0]+p[0][3][0])));
  ac[14][13]=0.;
  bcost[14][13]=(p[14][1][0]+p[13][5][0])/4.+(p[14][2][0]+p[13][4][0])/4.;
  ac[14][22]=((p[14][2][1]+p[22][0][1])/2.-(p[14][3][1]+p[22][5][1])/2.)/((p[14][2][0]+p[22][0][0])/2.-(p[14][3][0]+p[22][5][0])/2.);
  bcost[14][22]=((p[14][3][1]+p[22][5][1])*(p[14][2][0]+p[22][0][0])/4.-(p[14][3][0]+p[22][5][0])*(p[14][2][1]+p[22][0][1])/4.)/((p[14][2][0]+p[22][0][0])/2.-(p[14][3][0]+p[22][5][0])/2.);
  ac[14][8]=((p[14][3][1]+p[8][1][1])/2.-(p[14][4][1]+p[8][0][1])/2.)/((p[14][3][0]+p[8][1][0])/2.-(p[14][4][0]+p[8][0][0])/2.);
  bcost[14][8]=((p[14][4][1]+p[8][0][1])*(p[14][3][0]+p[8][1][0])/4.-(p[14][4][0]+p[8][0][0])*(p[14][3][1]+p[8][1][1])/4.)/((p[14][3][0]+p[8][1][0])/2.-(p[14][4][0]+p[8][0][0])/2.);
  ac[14][11]=0.;
  bcost[14][11]=(p[14][4][0]+p[11][2][0])/4.+(p[14][5][0]+p[11][1][0])/4.;
  ac[14][5]=((p[14][5][1]+p[5][3][1])/2.-(p[14][0][1]+p[5][2][1])/2.)/((p[14][5][0]+p[5][3][0])/2.-(p[14][0][0]+p[5][2][0])/2.);
  bcost[14][5]=((p[14][0][1]+p[6][2][1])*(p[14][5][0]+p[5][3][0])/4.-(p[14][0][0]+p[5][2][0])*(p[14][5][1]+p[5][3][1])/4.)/((p[14][5][0]+p[5][3][0])/2.-(p[14][0][0]+p[5][2][0])/2.);

  ac[13][3]=((p[13][0][1]+p[3][4][1])/2.-(p[13][1][1]+p[3][3][1])/2.)/((p[13][0][0]+p[3][4][0])/2.-(p[13][1][0]+p[3][3][0])/2.);
  bcost[13][3]=((p[13][1][1]+p[3][3][1])*(p[13][0][0]+p[3][4][0])/4.-(p[13][1][0]+p[3][3][0])*(p[13][0][1]+p[3][4][1])/4.)/((p[13][0][0]+p[3][4][0])/2.-(p[13][1][0]+p[3][3][0])/2.);
  ac[13][16]=0.;
  bcost[13][16]=(p[13][1][0]+p[16][5][0])/4.+(p[13][2][0]+p[16][4][0])/4.;
  ac[13][15]=((p[13][2][1]+p[15][0][1])/2.-(p[13][3][1]+p[15][5][1])/2.)/((p[13][2][0]+p[15][0][0])/2.-(p[13][3][0]+p[15][5][0])/2.);
  bcost[13][15]=((p[13][3][1]+p[15][5][1])*(p[13][2][0]+p[15][0][0])/4.-(p[13][3][0]+p[15][5][0])*(p[13][2][1]+p[15][0][1])/4.)/((p[13][2][0]+p[15][0][0])/2.-(p[13][3][0]+p[15][5][0])/2.);
  ac[13][22]=((p[13][3][1]+p[22][1][1])/2.-(p[13][4][1]+p[22][0][1])/2.)/((p[13][3][0]+p[22][1][0])/2.-(p[13][4][0]+p[22][0][0])/2.);
  bcost[13][22]=((p[13][4][1]+p[22][0][1])*(p[13][3][0]+p[22][1][0])/4.-(p[13][4][0]+p[22][0][0])*(p[13][3][1]+p[22][1][1])/4.)/((p[13][3][0]+p[22][1][0])/2.-(p[13][4][0]+p[22][0][0])/2.);
  ac[13][14]=0.;
  bcost[13][14]=(p[13][4][0]+p[14][2][0])/4.+(p[13][5][0]+p[14][1][0])/4.;
  ac[13][25]=((p[13][5][1]+p[0][3][1])-(p[13][0][1]+p[0][2][1]))/((p[13][5][0]+p[0][3][0])-(p[13][0][0]+p[0][2][0]));
  bcost[13][25]=((p[13][0][1]+p[0][2][1])*(p[13][5][0]+p[0][3][0])-(p[13][0][0]+p[0][2][0])*(p[13][5][1]+p[0][3][1]))/(2.*((p[13][5][0]+p[0][3][0])-(p[13][0][0]+p[0][2][0])));

  ac[16][1]=((p[16][0][1]+p[1][4][1])/2.-(p[16][1][1]+p[1][3][1])/2.)/((p[16][0][0]+p[1][4][0])/2.-(p[16][1][0]+p[1][3][0])/2.);
  bcost[16][1]=((p[16][1][1]+p[1][3][1])*(p[16][0][0]+p[1][4][0])/4.-(p[16][1][0]+p[1][3][0])*(p[16][0][1]+p[1][4][1])/4.)/((p[16][0][0]+p[1][4][0])/2.-(p[16][1][0]+p[1][3][0])/2.);
  ac[16][25]=0.;
  bcost[16][25]=(p[16][1][0]+p[0][5][0])/4.+(p[16][2][0]+p[0][4][0])/4.;
  ac[16][26]=((p[16][2][1]+p[0][0][1])-(p[16][3][1]+p[0][5][1]))/((p[16][2][0]+p[0][0][0])-(p[16][3][0]+p[0][5][0]));
  bcost[16][26]=((p[16][3][1]+p[0][5][1])*(p[16][2][0]+p[0][0][0])-(p[16][3][0]+p[0][5][0])*(p[16][2][1]+p[0][0][1]))/(2.*((p[16][2][0]+p[0][0][0])-(p[16][3][0]+p[0][5][0])));
  ac[16][15]=((p[16][3][1]+p[15][1][1])/2.-(p[16][4][1]+p[15][0][1])/2.)/((p[16][3][0]+p[15][1][0])/2.-(p[16][4][0]+p[15][0][0])/2.);
  bcost[16][15]=((p[16][4][1]+p[15][0][1])*(p[16][3][0]+p[15][1][0])/4.-(p[16][4][0]+p[15][0][0])*(p[16][3][1]+p[15][1][1])/4.)/((p[16][3][0]+p[15][1][0])/2.-(p[16][4][0]+p[15][0][0])/2.);
  ac[16][13]=0.;
  bcost[16][13]=(p[16][4][0]+p[13][2][0])/4.+(p[16][5][0]+p[13][1][0])/4.;
  ac[16][3]=((p[16][5][1]+p[3][3][1])/2.-(p[16][0][1]+p[3][2][1])/2.)/((p[16][5][0]+p[3][3][0])/2.-(p[16][0][0]+p[3][2][0])/2.);
  bcost[16][3]=((p[16][0][1]+p[3][2][1])*(p[16][5][0]+p[3][3][0])/4.-(p[16][0][0]+p[3][2][0])*(p[16][5][1]+p[3][3][1])/4.)/((p[16][5][0]+p[3][3][0])/2.-(p[16][0][0]+p[3][2][0])/2.);

  ac[8][14]=((p[8][0][1]+p[14][4][1])/2.-(p[8][1][1]+p[14][3][1])/2.)/((p[8][0][0]+p[14][4][0])/2.-(p[8][1][0]+p[14][3][0])/2.);
  bcost[8][14]=((p[8][1][1]+p[14][3][1])*(p[8][0][0]+p[14][4][0])/4.-(p[8][1][0]+p[14][3][0])*(p[8][0][1]+p[14][4][1])/4.)/((p[8][0][0]+p[14][4][0])/2.-(p[8][1][0]+p[14][3][0])/2.);
  ac[8][22]=0.;
  bcost[8][22]=(p[8][1][0]+p[22][5][0])/4.+(p[8][2][0]+p[22][4][0])/4.;
  ac[8][25]=((p[8][2][1]+p[0][0][1])-(p[8][3][1]+p[0][5][1]))/((p[8][2][0]+p[0][0][0])-(p[8][3][0]+p[0][5][0]));
  bcost[8][25]=((p[8][3][1]+p[0][5][1])*(p[8][2][0]+p[0][0][0])-(p[8][3][0]+p[0][5][0])*(p[8][2][1]+p[0][0][1]))/(2.*((p[8][2][0]+p[0][0][0])-(p[8][3][0]+p[0][5][0])));
  ac[8][26]=((p[8][3][1]+p[0][1][1])-(p[8][4][1]+p[0][0][1]))/((p[8][3][0]+p[0][1][0])-(p[8][4][0]+p[0][0][0]));
  bcost[8][26]=((p[8][4][1]+p[0][0][1])*(p[8][3][0]+p[0][1][0])-(p[8][4][0]+p[0][0][0])*(p[8][3][1]+p[0][1][1]))/(2.*((p[8][3][0]+p[0][1][0])-(p[8][4][0]+p[0][0][0])));
  ac[8][27]=0.;
  bcost[8][27]=(p[8][4][0]+p[0][2][0])/4.+(p[8][5][0]+p[0][1][0])/4.;
  ac[8][11]=((p[8][5][1]+p[11][3][1])/2.-(p[8][0][1]+p[11][2][1])/2.)/((p[8][5][0]+p[11][3][0])/2.-(p[8][0][0]+p[11][2][0])/2.);
  bcost[8][11]=((p[8][0][1]+p[11][2][1])*(p[8][5][0]+p[11][3][0])/4.-(p[8][0][0]+p[11][2][0])*(p[8][5][1]+p[11][3][1])/4.)/((p[8][5][0]+p[11][3][0])/2.-(p[8][0][0]+p[11][2][0])/2.);

  ac[22][13]=((p[22][0][1]+p[13][4][1])/2.-(p[22][1][1]+p[13][3][1])/2.)/((p[22][0][0]+p[13][4][0])/2.-(p[22][1][0]+p[13][3][0])/2.);
  bcost[22][13]=((p[22][1][1]+p[13][3][1])*(p[22][0][0]+p[13][4][0])/4.-(p[22][1][0]+p[13][3][0])*(p[22][0][1]+p[13][4][1])/4.)/((p[22][0][0]+p[13][4][0])/2.-(p[22][1][0]+p[13][3][0])/2.);
  ac[22][15]=0.;
  bcost[22][15]=(p[22][1][0]+p[15][5][0])/4.+(p[22][2][0]+p[15][4][0])/4.;
  ac[22][25]=((p[22][2][1]+p[0][0][1])-(p[22][3][1]+p[0][5][1]))/((p[22][2][0]+p[0][0][0])-(p[22][3][0]+p[0][5][0]));
  bcost[22][25]=((p[22][3][1]+p[0][5][1])*(p[22][2][0]+p[0][0][0])-(p[22][3][0]+p[0][5][0])*(p[22][2][1]+p[0][0][1]))/(2.*((p[22][2][0]+p[0][0][0])-(p[22][3][0]+p[0][5][0])));
  ac[22][26]=((p[22][3][1]+p[0][1][1])-(p[22][4][1]+p[0][0][1]))/((p[22][3][0]+p[0][1][0])-(p[22][4][0]+p[0][0][0]));
  bcost[22][26]=((p[22][4][1]+p[0][0][1])*(p[22][3][0]+p[0][1][0])-(p[22][4][0]+p[0][0][0])*(p[22][3][1]+p[0][1][1]))/(2.*((p[22][3][0]+p[0][1][0])-(p[22][4][0]+p[0][0][0])));
  ac[22][8]=0.;
  bcost[22][8]=(p[22][4][0]+p[8][2][0])/4.+(p[22][5][0]+p[8][1][0])/4.;
  ac[22][14]=((p[22][5][1]+p[14][3][1])/2.-(p[22][0][1]+p[14][2][1])/2.)/((p[22][5][0]+p[14][3][0])/2.-(p[22][0][0]+p[14][2][0])/2.);
  bcost[22][14]=((p[22][0][1]+p[14][2][1])*(p[22][5][0]+p[14][3][0])/4.-(p[22][0][0]+p[14][2][0])*(p[22][5][1]+p[14][3][1])/4.)/((p[22][5][0]+p[14][3][0])/2.-(p[22][0][0]+p[14][2][0])/2.);

  ac[15][16]=((p[15][0][1]+p[16][4][1])/2.-(p[15][1][1]+p[16][3][1])/2.)/((p[15][0][0]+p[16][4][0])/2.-(p[15][1][0]+p[16][3][0])/2.);
  bcost[15][16]=((p[15][1][1]+p[16][3][1])*(p[15][0][0]+p[16][4][0])/4.-(p[15][1][0]+p[16][3][0])*(p[15][0][1]+p[16][4][1])/4.)/((p[15][0][0]+p[16][4][0])/2.-(p[15][1][0]+p[16][3][0])/2.);
  ac[15][25]=0.;
  bcost[15][25]=(p[15][1][0]+p[0][5][0])/4.+(p[15][2][0]+p[0][4][0])/4.;
  ac[15][26]=((p[15][2][1]+p[0][0][1])-(p[15][3][1]+p[0][5][1]))/((p[15][2][0]+p[0][0][0])-(p[15][3][0]+p[0][5][0]));
  bcost[15][26]=((p[15][3][1]+p[0][5][1])*(p[15][2][0]+p[0][0][0])-(p[15][3][0]+p[0][5][0])*(p[15][2][1]+p[0][0][1]))/(2.*((p[15][2][0]+p[0][0][0])-(p[15][3][0]+p[0][5][0])));
  ac[15][27]=((p[15][3][1]+p[0][1][1])-(p[15][4][1]+p[0][0][1]))/((p[15][3][0]+p[0][1][0])-(p[15][4][0]+p[0][0][0]));
  bcost[15][27]=((p[15][4][1]+p[0][0][1])*(p[15][3][0]+p[0][1][0])-(p[15][4][0]+p[0][0][0])*(p[15][3][1]+p[0][1][1]))/(2.*((p[15][3][0]+p[0][1][0])-(p[15][4][0]+p[0][0][0])));
  ac[15][22]=0.;
  bcost[15][22]=(p[15][4][0]+p[22][2][0])/4.+(p[15][5][0]+p[22][1][0])/4.;
  ac[15][13]=((p[15][5][1]+p[13][3][1])/2.-(p[15][0][1]+p[13][2][1])/2.)/((p[15][5][0]+p[13][3][0])/2.-(p[15][0][0]+p[13][2][0])/2.);
  bcost[15][13]=((p[15][0][1]+p[13][2][1])*(p[15][5][0]+p[13][3][0])/4.-(p[15][0][0]+p[13][2][0])*(p[15][5][1]+p[13][3][1])/4.)/((p[15][5][0]+p[13][3][0])/2.-(p[15][0][0]+p[13][2][0])/2.);

  // Mirror offsets
//  fDeltaMisa[1][0]=0.;
//  fDeltaMisa[1][1]=0.;
//  fDeltaMisa[3][0]=-0.44;
//  fDeltaMisa[3][1]=-18.4;
//  fDeltaMisa[4][0]=-18.;
//  fDeltaMisa[4][1]=-17.;
//  fDeltaMisa[5][0]=-4.;
//  fDeltaMisa[5][1]=3.14;
//  fDeltaMisa[6][0]=-0.44;
//  fDeltaMisa[6][1]=-5.63;
//  fDeltaMisa[8][0]=-12.0;
//  fDeltaMisa[8][1]=-2.5;
//  fDeltaMisa[9][0]=-9.0;
//  fDeltaMisa[9][1]=-4.9;
//  fDeltaMisa[10][0]=-8.69;
//  fDeltaMisa[10][1]=-1.92;
//  fDeltaMisa[11][0]=-16.09;
//  fDeltaMisa[11][1]=-19.97;
//  fDeltaMisa[12][0]=-6.89;
//  fDeltaMisa[12][1]=-13.0;
//  fDeltaMisa[13][0]=19.06;
//  fDeltaMisa[13][1]=-22.84;
//  fDeltaMisa[14][0]=6.1;
//  fDeltaMisa[14][1]=8.14;
//  fDeltaMisa[15][0]=-9.;
//  fDeltaMisa[15][1]=-18.;
//  fDeltaMisa[16][0]=5.6;
//  fDeltaMisa[16][1]=-14.0;
//  fDeltaMisa[17][0]=-7.07;
//  fDeltaMisa[17][1]=-11.6;
//  fDeltaMisa[20][0]=-6.99;
//  fDeltaMisa[20][1]=13.76;
//  fDeltaMisa[21][0]=-4.0;
//  fDeltaMisa[21][1]=-5.4;
//  fDeltaMisa[22][0]=4.61;
//  fDeltaMisa[22][1]=-4.4;
//  fDeltaMisa[23][0]=0.;
//  fDeltaMisa[23][1]=0.;
//  fDeltaMisa[24][0]=0.;
//  fDeltaMisa[24][1]=0.;
////  fDeltaMisa[2][0]=0.;
////  fDeltaMisa[2][1]=0.;
//  fDeltaMisa[2][0]=-31.;
//  fDeltaMisa[2][1]=-76.;
//  fDeltaMisa[7][0]=0.;
//  fDeltaMisa[7][1]=0.;
//  fDeltaMisa[18][0]=0.;
//  fDeltaMisa[18][1]=0.;
//  fDeltaMisa[19][0]=0.;
//  fDeltaMisa[19][1]=0.;
//  fDeltaMisa[99][0]=0.;
//  fDeltaMisa[99][1]=0.;

//////  fDeltaMisa[1][0]=0.;
//////  fDeltaMisa[1][1]=0.;
//////
//////  fDeltaMisa[3][0]=-6.9;
//////  fDeltaMisa[3][1]=-13.2;
//////
//////  fDeltaMisa[4][0]=-18.;
//////  fDeltaMisa[4][1]=-17.;
//////
//////  fDeltaMisa[5][0]=8.3;
//////  fDeltaMisa[5][1]=9.3;
//////
//////  fDeltaMisa[6][0]=-20.2;
//////  fDeltaMisa[6][1]=14.7;
//////
//////  fDeltaMisa[8][0]=-12.0;
//////  fDeltaMisa[8][1]=-2.5;
//////
//////  fDeltaMisa[9][0]=-9.0;
//////  fDeltaMisa[9][1]=-4.9;
//////
//////  fDeltaMisa[10][0]=8.8;
//////  fDeltaMisa[10][1]=10.5;
//////
//////  fDeltaMisa[11][0]=7.1;
//////  fDeltaMisa[11][1]=12.9;
//////
//////  fDeltaMisa[12][0]=-6.8;
//////  fDeltaMisa[12][1]=-11.4;
//////
//////  fDeltaMisa[13][0]=-6.8;
//////  fDeltaMisa[13][1]=-11.2;
//////
//////  fDeltaMisa[14][0]=7.8;
//////  fDeltaMisa[14][1]=11.1;
//////
//////  fDeltaMisa[15][0]=-9.;
//////  fDeltaMisa[15][1]=-18.;
//////
//////  fDeltaMisa[16][0]=5.6;
//////  fDeltaMisa[16][1]=-14.0;
//////
//////  fDeltaMisa[17][0]=10.1;
//////  fDeltaMisa[17][1]=14.3;
//////
//////  fDeltaMisa[20][0]=7.6;
//////  fDeltaMisa[20][1]=10.4;
//////
//////  fDeltaMisa[21][0]=-4.0;
//////  fDeltaMisa[21][1]=-5.4;
//////
//////  fDeltaMisa[22][0]=3.9;
//////  fDeltaMisa[22][1]=-6.8;
//////
//////  fDeltaMisa[23][0]=0.;
//////  fDeltaMisa[23][1]=0.;
//////  fDeltaMisa[24][0]=0.;
//////  fDeltaMisa[24][1]=0.;
////////  fDeltaMisa[2][0]=0.;
////////  fDeltaMisa[2][1]=0.;
//////  fDeltaMisa[2][0]=-31.;
//////  fDeltaMisa[2][1]=-76.;
//////  fDeltaMisa[7][0]=0.;
//////  fDeltaMisa[7][1]=0.;
//////  fDeltaMisa[18][0]=0.;
//////  fDeltaMisa[18][1]=0.;
//////  fDeltaMisa[19][0]=0.;
//////  fDeltaMisa[19][1]=0.;
//////  fDeltaMisa[99][0]=0.;
//////  fDeltaMisa[99][1]=0.;
}

/********************************************//**
 * Mirror focus position corrections
 ************************************************/
TVector2 RICHImprovedRing::GetChPosAngCorr(Int_t iCh)
{
  TVector2 Beta;
  if (iCh < 1952/2 || (1952 <= iCh && iCh < 1952+(243/2))) Beta = TVector2(127,0);
  else if ((1952/2 <= iCh && iCh < 1952) || (1952+(243/2) <= iCh && iCh <1952+243 )) Beta = TVector2(177,0);
  TVector2 Correction = Beta;// the correction due to the mirror rotation is already multiplied by the fFocalLength and it is given in mm;
  return Correction;
}
