#include "AnalysisTools_G.hh"
//#include "BTubeField.hh"
#include "TRecoSpectrometerCandidate.hh"
#include "TRecoLKrCandidate.hh"
#include "TMath.h"

//#include "Particle.hh"
//#include "BlueTubeTracker.hh"

AnalysisTools_G* AnalysisTools_G::fInstance = 0;

AnalysisTools_G::AnalysisTools_G() {
//  fBTubeField = new BTubeField();
//  fTracker = new BlueTubeTracker();
  InitializeMirrorCoordinate();
  fIsClusterCorrected = 0;
}

AnalysisTools_G *AnalysisTools_G::GetInstance() {
  if (fInstance==0) {fInstance = new AnalysisTools_G();}
  return fInstance;
}

/////////////////////////
// Single track vertex //
/////////////////////////
TVector3 AnalysisTools_G::SingleTrackVertex(TVector3 *vect1, TVector3 *vect2, TVector3 *vectb, TVector3 *vectc, Double_t *cda) {
  TVector3 v1 = *vect1;
  TVector3 v2 = *vect2;
  TVector3 pos1 = *vectb;
  TVector3 pos2 = *vectc;
  TVector3 r12 = pos1-pos2;
  Double_t v1xv2 = v1.Dot(v2);
  Double_t det   = pow(v1xv2,2)-v1.Mag2()*v2.Mag2();
  if (!det) return TVector3(-9999,-9999,-9999);
  Double_t t1 = (v2.Mag2()*r12.Dot(v1)-v1.Dot(v2)*r12.Dot(v2))/det;
  Double_t t2 = (v1.Dot(v2)*r12.Dot(v1)-v1.Mag2()*r12.Dot(v2))/det;
  TVector3 q1 = pos1+t1*v1;
  TVector3 q2 = pos2+t2*v2;
  TVector3 vertex = 0.5*(q1+q2);
  *cda = (q1-q2).Mag();
  return vertex;
}

///////////////////////////////////////////
// Determination of a multi track vertex //
// Simple analytical derivation          //
// nTracks: number of tracks             //
// ptracks: 4-momentum before magnet     //
// postracks: position before magnet     //
// cda: output variable                  //              
///////////////////////////////////////////
TVector3 AnalysisTools_G::MultiTrackVertexSimple(Int_t nTracks, TLorentzVector *ptracks, TVector3 *postracks, Double_t *cda) {
  TVector3 avPosition(0,0,0);
  TVector3 avSlope(0,0,0);
  TVector3 avSlope2(0,0,0);
  TVector3 avMixed(0,0,0);

  // Compute Z as the position of minimum apporach between tracks
  Double_t z0 = 0;
  for (Int_t j=0; j<nTracks; j++) {
    TVector3 position = postracks[j];
    TLorentzVector momentum = ptracks[j]; 
    avPosition += position;
    TVector3 ddz = momentum.Vect()*(1./momentum.Vect().Z());
    avSlope += ddz;
    avSlope2 += TVector3(ddz.X()*ddz.X(),ddz.Y()*ddz.Y(),ddz.Z()*ddz.Z());
    avMixed += TVector3(position.X()*ddz.X(),position.Y()*ddz.Y(),position.Z()*ddz.Z());
    z0 = position.Z();
  }
  avPosition = (1./nTracks)*avPosition;
  avSlope = (1./nTracks)*avSlope;
  avSlope2 = (1./nTracks)*avSlope2;
  avMixed = (1./nTracks)*avMixed;
  Double_t num = nTracks*(avMixed.X()+avMixed.Y())-nTracks*(avPosition.X()*avSlope.X()+avPosition.Y()*avSlope.Y());
  Double_t den = nTracks*(avSlope2.X()+avSlope2.Y())-nTracks*(avSlope.X()*avSlope.X()+avSlope.Y()*avSlope.Y());
  Double_t zvertex = z0-num/den;

  // Compute the trasnverse position and the cda
  TVector3 avPosVtx(0,0,0);
  TVector3 avPosVtx2(0,0,0);
  for (Int_t j=0; j<nTracks; j++) {
    TVector3 position = postracks[j];
    TLorentzVector momentum = ptracks[j]; 
    TVector3 posvtx = position+momentum.Vect()*(1./momentum.Vect().Z())*(zvertex-position.Z());
    avPosVtx += posvtx;
    avPosVtx2 += TVector3(posvtx.X()*posvtx.X(),posvtx.Y()*posvtx.Y(),posvtx.Z()*posvtx.Z());
  }
  avPosVtx = (1./nTracks)*avPosVtx;
  avPosVtx2 = (1./nTracks)*avPosVtx2;
  *cda = sqrt(avPosVtx2.X()+avPosVtx2.Y()-avPosVtx.X()*avPosVtx.X()-avPosVtx.Y()*avPosVtx.Y());

  return TVector3(avPosVtx.X(),avPosVtx.Y(),zvertex);
}

//////////////////////////////////////
// Compute 4-Momentum of a particle //
// patrack: 0 |momentum|            //
//          1 theta XZ              //
//          2 theta YZ              //
//          3 mass                  //
//////////////////////////////////////
TLorentzVector AnalysisTools_G::Get4Momentum(Double_t *partrack) {
  TLorentzVector pmom;
  Double_t thetaX = partrack[1];
  Double_t thetaY = partrack[2];
  Double_t pmag = partrack[0];
  Double_t pmomz = pmag/sqrt(1.+thetaX*thetaX+thetaY*thetaY);
  Double_t pmomx = pmomz*thetaX;
  Double_t pmomy = pmomz*thetaY;
  pmom.SetXYZM(pmomx/1000,pmomy/1000,pmomz/1000,partrack[3]);
  return pmom;
}

///////////////////////////////////
// Get track position at plane Z //
///////////////////////////////////
TVector3 AnalysisTools_G::GetPositionAtZ(TRecoSpectrometerCandidate *cand, Double_t zpos) {
  Double_t posx;
  Double_t posy;
  if (zpos<196345) { // front face of MNP33
    posx = cand->GetPositionBeforeMagnet().X()+cand->GetSlopeXBeforeMagnet()*(zpos-cand->GetPositionBeforeMagnet().Z());
    posy = cand->GetPositionBeforeMagnet().Y()+cand->GetSlopeYBeforeMagnet()*(zpos-cand->GetPositionBeforeMagnet().Z());
  } else {
    posx = cand->GetPositionAfterMagnet().X()+cand->GetSlopeXAfterMagnet()*(zpos-cand->GetPositionAfterMagnet().Z());
    posy = cand->GetPositionAfterMagnet().Y()+cand->GetSlopeYAfterMagnet()*(zpos-cand->GetPositionAfterMagnet().Z());
  }
  return TVector3(posx,posy,zpos);
}

///////////////////////
// Mirror reflection //
///////////////////////
void AnalysisTools_G::MirrorVector(TVector3 *vec, TVector3 *axis) {
  TVector3 S = ((vec->Dot(*axis))/(axis->Mag2()))*(*axis);
  TVector3 d = S-*vec;
  TVector3 ret = S+d;
  *vec = ret.Unit();
}

////////////////////////
// Find the mirror ID //
////////////////////////
Int_t AnalysisTools_G::MirrorSurface(Double_t xin, Double_t yin, Double_t SafeFactor, Bool_t flag) {
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

//////////////////////////////
// LKr Clusters Corrections //
//////////////////////////////
void AnalysisTools_G::ClusterCorrections(Int_t runNumber, TRecoLKrCandidate *thisCluster) {
//  if (fIsClusterCorrected) return;
  Double_t fEScale = 1.;
  if (runNumber<1600) fEScale = 1.02;
  else fEScale = 1.03;

  // ZS non linearity
  Double_t ue = thisCluster->GetClusterEnergy()/1000.;
  Double_t ce = ue;

  if (thisCluster->GetNCells()>9) {
    if (runNumber<1600) {
      if (ue<25.) ce = ue/(0.724761+0.0220857*ue-0.000873029*ue*ue+0.0000126505*ue*ue*ue); 
      if (ue>=25. && ue<65.) ce = ue/(0.928916625+0.00262869*(ue-25.)-0.000063694*(ue-25.)*(ue-25.)+0.00000068890*(ue-25.)*(ue-25.)*(ue-25.));
      if (ue>65.) ce = ue/0.976243425;
    } else {
       if (ue<22) ce = ue/(0.7666+0.0573489*log(ue));
       if (ue>=22 && ue<65) ce = ue/(0.828962+0.0369797*log(ue));
       if (ue>=65) ce = ue/(0.828962+0.0369797*log(65));
    }
  }

  thisCluster->SetClusterEnergy(ce*1000.);
  thisCluster->SetClusterEnergy(thisCluster->GetClusterEnergy()*fEScale); 

  // Energy loss in the hole
  Double_t energy = thisCluster->GetClusterEnergy()/1000.;
  Double_t radius = sqrt(thisCluster->GetClusterX()/10.*thisCluster->GetClusterX()/10.+thisCluster->GetClusterY()/10.*thisCluster->GetClusterY()/10.);
  if (radius>=14 && radius<=18.5) energy = energy/(0.97249+0.0014692*radius)*0.9999;
  Double_t rcorr2 = 0;
  if (radius>=14 && radius<18) rcorr2 = 0.0042-3.7e-4*radius; 
  if (radius>=18 && radius<20) rcorr2 = -0.00211; 
  if (radius>=20 && radius<22) rcorr2 = -0.01694-7.769e-4*radius; 
  thisCluster->SetClusterEnergy(energy*(1-rcorr2)*1000.);

  // Energy loss at low energy
  energy = thisCluster->GetClusterEnergy()/1000.;
  if (energy<15) thisCluster->SetClusterEnergy(1000*(energy+0.015)/(15+0.015)*15*0.9999);

  // Alignment correction
  thisCluster->SetClusterX(thisCluster->GetClusterX()+10*0.13);
  thisCluster->SetClusterY(thisCluster->GetClusterY()+10*0.08);

  // Global scale factor
  thisCluster->SetClusterSeedEnergy(thisCluster->GetClusterSeedEnergy()*fEScale); 
  thisCluster->SetCluster77Energy(thisCluster->GetCluster77Energy()*fEScale); 
 
  fIsClusterCorrected = 1;
}

////////////////////////////
// Get LAV station from Z //
////////////////////////////
Int_t AnalysisTools_G::GetLAVStation(Double_t z) {
  if (z>120000 && z<123000) return 0;
  if (z>128000 && z<130000) return 1;
  if (z>136000 && z<139000) return 2;
  if (z>142000 && z<146000) return 3;
  if (z>150000 && z<154000) return 4;
  if (z>164000 && z<168000) return 5;
  if (z>170000 && z<176000) return 6;
  if (z>178000 && z<182000) return 7;
  if (z>192000 && z<196000) return 8;
  if (z>202000 && z<206000) return 9;
  if (z>214000 && z<222000) return 10;
  if (z>236000 && z<240000) return 11;
  return -1;
}

//////////////////////////////////////////////
// Track propagation. It takes into account //
// the magnetic field of the spectrometer   //
//////////////////////////////////////////////
TVector3 AnalysisTools_G::Propagate(Int_t charge, TVector3 *momentum, TVector3 *position, Double_t *thetaxafter, Double_t fZEnd)
{
  Double_t fEC = TMath::C()* 1.e-9 * 1.e-4 * 1.e-2;
  TVector3 fB(0.,0.6928*10000,0.);
  TVector3 fPosExp;
  Double_t dMag = 0.1*1300;
  Double_t zMag = 0.1*196345;
//  Double_t fStartX = 0.1*position->X();
//  Double_t fStartY = 0.1*position->Y();
//  Double_t fStartZ = 0.1*position->Z();
  Double_t fStartX = position->X();
  Double_t fStartY = position->Y();
  Double_t fStartZ = position->Z();
  Double_t fPartP = momentum->Mag();
  Double_t fPartThetaX = momentum->X()/momentum->Z(); 
  Double_t fPartThetaY = momentum->Y()/momentum->Z(); 
//  fZEnd *= 0.1;
  Int_t fPartQ = charge;

  // fZEnd before magnet
  if ((fZEnd<=zMag&&fStartZ<=zMag) || (fZEnd>zMag&&fStartZ>zMag)) {
    fPosExp.SetX(fStartX+fPartThetaX*(fZEnd-fStartZ));
    fPosExp.SetY(fStartY+fPartThetaY*(fZEnd-fStartZ));
    fPosExp.SetZ(fZEnd);
//    return fPosExp*10;
    return fPosExp;
  }

  // fZEnd after MNP33
  fPosExp.SetX(fStartX+fPartThetaX*(zMag-fStartZ));
  fPosExp.SetY(fStartY+fPartThetaY*(zMag-fStartZ));
  fPosExp.SetZ(zMag);
  TVector3 fP;
  fP.SetZ(fPartP/sqrt(1.+fPartThetaX*fPartThetaX+fPartThetaY*fPartThetaY));
  fP.SetX(fP.Z()*fPartThetaX);
  fP.SetY(fP.Z()*fPartThetaY);
  Int_t qb = fB.Y()>0 ? 1 : -1;
  Double_t rho = (fP.Cross(fB)).Mag()/(fPartQ*fEC*fB.Mag2());
  Double_t delta = dMag/rho;
  Double_t sint = sin(atan(fPartThetaX));
  Double_t cost = cos(atan(fPartThetaX));
  Double_t dx = qb*rho*(-cost+sqrt(1-(delta-qb*sint)*(delta-qb*sint)));
  fPosExp.SetX(fPosExp.X()+dx);
  fPosExp.SetY(fPosExp.Y()+fPartThetaY*dMag);
  fPosExp.SetZ(fPosExp.Z()+dMag);
  Double_t fThetaXAfter = -qb*(delta-qb*sint)/sqrt(1.-(delta-qb*sint)*(delta-qb*sint));
  fPosExp.SetX(fPosExp.X()+fThetaXAfter*(fZEnd-fPosExp.Z()));
  fPosExp.SetY(fPosExp.Y()+fPartThetaY*(fZEnd-fPosExp.Z()));
  fPosExp.SetZ(fZEnd);
  *thetaxafter = fThetaXAfter;
 
//  return fPosExp*10;  
  return fPosExp;  
}

//void AnalysisTools_G::BTubeCorrection(TVector3 *vertex, Double_t *dxdz, Double_t *dydz, Double_t pmom) {
//  fBTubeField->GetCorrection(vertex,dxdz,dydz,pmom);
//} 

Bool_t AnalysisTools_G::LKrAcceptance(Double_t x, Double_t y, Double_t rmin, Double_t lmax) {
  Double_t r2 = x*x+y*y;
  Double_t ax = fabs(x);
  Double_t ay = fabs(y);
  if(r2<=rmin*rmin) return false;
  if(ax>=lmax) return false;
  if(ay>=lmax) return false;
  if((ax+ay)>=sqrt(2.)*lmax) return false;
  return true;
}

void AnalysisTools_G::InitializeMirrorCoordinate() {
  fMirrorPos[0][0][0]=0.;
  fMirrorPos[0][0][1]=0.;
  fMirrorPos[0][1][0]=0.;
  fMirrorPos[0][1][1]=0.;
  fMirrorPos[0][2][0]=0.;
  fMirrorPos[0][2][1]=0.;
  fMirrorPos[0][3][0]=0.;
  fMirrorPos[0][3][1]=0.;
  fMirrorPos[0][4][0]=0.;
  fMirrorPos[0][4][1]=0.;
  fMirrorPos[0][5][0]=0.;
  fMirrorPos[0][5][1]=0.;
  fMirrorPos[9][0][0]=-604.9822;
  fMirrorPos[9][0][1]=1407.273;
  fMirrorPos[9][1][0]=-302.1713;
  fMirrorPos[9][1][1]=1232.108;
  fMirrorPos[9][2][0]=-302.5643;
  fMirrorPos[9][2][1]=880.618;
  fMirrorPos[9][3][0]=-605.4903;
  fMirrorPos[9][3][1]=705.981;
  fMirrorPos[9][4][0]=-907.3783;
  fMirrorPos[9][4][1]=880.507;
  fMirrorPos[9][5][0]=-908.4193;
  fMirrorPos[9][5][1]=1232.219;
  fMirrorPos[17][0][0]=1.140884;
  fMirrorPos[17][0][1]=1410.055;
  fMirrorPos[17][1][0]=304.2186;
  fMirrorPos[17][1][1]=1234.862;
  fMirrorPos[17][2][0]=304.1706;
  fMirrorPos[17][2][1]=883.6355;
  fMirrorPos[17][3][0]=1.159304;
  fMirrorPos[17][3][1]=708.5905;
  fMirrorPos[17][4][0]=-301.5964;
  fMirrorPos[17][4][1]=883.6275;
  fMirrorPos[17][5][0]=-301.7034;
  fMirrorPos[17][5][1]=1234.87;
  fMirrorPos[4][0][0]=616.0545;
  fMirrorPos[4][0][1]=1411.085;
  fMirrorPos[4][1][0]=919.1667;
  fMirrorPos[4][1][1]=1236.041;
  fMirrorPos[4][2][0]=919.0697;
  fMirrorPos[4][2][1]=884.8148;
  fMirrorPos[4][3][0]=616.1178;
  fMirrorPos[4][3][1]=709.8158;
  fMirrorPos[4][4][0]=313.3247;
  fMirrorPos[4][4][1]=884.8058;
  fMirrorPos[4][5][0]=313.2847;
  fMirrorPos[4][5][1]=1236.05;
  fMirrorPos[10][0][0]=-906.7862;
  fMirrorPos[10][0][1]=879.7446;
  fMirrorPos[10][1][0]=-603.9165;
  fMirrorPos[10][1][1]=704.7246;
  fMirrorPos[10][2][0]=-603.9505;
  fMirrorPos[10][2][1]=353.3726;
  fMirrorPos[10][3][0]=-906.9273;
  fMirrorPos[10][3][1]=178.4316;
  fMirrorPos[10][4][0]=-1209.824;
  fMirrorPos[10][4][1]=353.3636;
  fMirrorPos[10][5][0]=-1209.788;
  fMirrorPos[10][5][1]=704.7336;
  fMirrorPos[20][0][0]=-299.2179;
  fMirrorPos[20][0][1]=879.3759;
  fMirrorPos[20][1][0]=3.730599;
  fMirrorPos[20][1][1]=704.5339;
  fMirrorPos[20][2][0]=3.581612;
  fMirrorPos[20][2][1]=353.1239;
  fMirrorPos[20][3][0]=-299.2848;
  fMirrorPos[20][3][1]=178.0999;
  fMirrorPos[20][4][0]=-602.1384;
  fMirrorPos[20][4][1]=353.2229;
  fMirrorPos[20][5][0]=-602.0394;
  fMirrorPos[20][5][1]=704.4349;
  fMirrorPos[12][0][0]=313.3375;
  fMirrorPos[12][0][1]=883.5286;
  fMirrorPos[12][1][0]=616.2212;
  fMirrorPos[12][1][1]=708.5726;
  fMirrorPos[12][2][0]=616.0572;
  fMirrorPos[12][2][1]=357.1526;
  fMirrorPos[12][3][0]=313.071;
  fMirrorPos[12][3][1]=182.1756;
  fMirrorPos[12][4][0]=10.27421;
  fMirrorPos[12][4][1]=357.1916;
  fMirrorPos[12][5][0]=10.37623;
  fMirrorPos[12][5][1]=708.5336;
  fMirrorPos[21][0][0]=921.0824;
  fMirrorPos[21][0][1]=883.6236;
  fMirrorPos[21][1][0]=1224.16;
  fMirrorPos[21][1][1]=708.5726;
  fMirrorPos[21][2][0]=1224.072;
  fMirrorPos[21][2][1]=357.3046;
  fMirrorPos[21][3][0]=921.1059;
  fMirrorPos[21][3][1]=182.3506;
  fMirrorPos[21][4][0]=618.2943;
  fMirrorPos[21][4][1]=357.3456;
  fMirrorPos[21][5][0]=618.2583;
  fMirrorPos[21][5][1]=708.5316;
  fMirrorPos[6][0][0]=-1208.882;
  fMirrorPos[6][0][1]=351.8385;
  fMirrorPos[6][1][0]=-905.842;
  fMirrorPos[6][1][1]=176.7515;
  fMirrorPos[6][2][0]=-906.043;
  fMirrorPos[6][2][1]=-174.5805;
  fMirrorPos[6][3][0]=-1209.044;
  fMirrorPos[6][3][1]=-349.4885;
  fMirrorPos[6][4][0]=-1511.81;
  fMirrorPos[6][4][1]=-174.5345;
  fMirrorPos[6][5][0]=-1511.817;
  fMirrorPos[6][5][1]=176.7055;
  fMirrorPos[5][0][0]=-601.9675;
  fMirrorPos[5][0][1]=352.9382;
  fMirrorPos[5][1][0]=-298.9482;
  fMirrorPos[5][1][1]=177.7952;
  fMirrorPos[5][2][0]=-298.9742;
  fMirrorPos[5][2][1]=-173.5028;
  fMirrorPos[5][3][0]=-602.0141;
  fMirrorPos[5][3][1]=-348.5198;
  fMirrorPos[5][4][0]=-904.9082;
  fMirrorPos[5][4][1]=-173.5058;
  fMirrorPos[5][5][0]=-904.8583;
  fMirrorPos[5][5][1]=177.7982;
  fMirrorPos[3][0][0]=617.0106;
  fMirrorPos[3][0][1]=353.4884;
  fMirrorPos[3][1][0]=919.8736;
  fMirrorPos[3][1][1]=178.6014;
  fMirrorPos[3][2][0]=919.7276;
  fMirrorPos[3][2][1]=-172.8446;
  fMirrorPos[3][3][0]=616.5578;
  fMirrorPos[3][3][1]=-347.7606;
  fMirrorPos[3][4][0]=313.9646;
  fMirrorPos[3][4][1]=-172.6736;
  fMirrorPos[3][5][0]=314.0006;
  fMirrorPos[3][5][1]=178.4304;
  fMirrorPos[1][0][0]=1224.873;
  fMirrorPos[1][0][1]=352.6879;
  fMirrorPos[1][1][0]=1527.851;
  fMirrorPos[1][1][1]=177.4779;
  fMirrorPos[1][2][0]=1527.732;
  fMirrorPos[1][2][1]=-173.6701;
  fMirrorPos[1][3][0]=1224.303;
  fMirrorPos[1][3][1]=-348.801;
  fMirrorPos[1][4][0]=921.3879;
  fMirrorPos[1][4][1]=-173.6471;
  fMirrorPos[1][5][0]=921.3839;
  fMirrorPos[1][5][1]=177.4549;
  fMirrorPos[11][0][0]=-904.0635;
  fMirrorPos[11][0][1]=-177.051;
  fMirrorPos[11][1][0]=-600.9554;
  fMirrorPos[11][1][1]=-352.213;
  fMirrorPos[11][2][0]=-601.0044;
  fMirrorPos[11][2][1]=-703.431;
  fMirrorPos[11][3][0]=-904.0333;
  fMirrorPos[11][3][1]=-878.454;
  fMirrorPos[11][4][0]=-1206.766;
  fMirrorPos[11][4][1]=-703.424;
  fMirrorPos[11][5][0]=-1206.844;
  fMirrorPos[11][5][1]=-352.22;
  fMirrorPos[14][0][0]=-297.6943;
  fMirrorPos[14][0][1]=-175.6188;
  fMirrorPos[14][1][0]=5.315837;
  fMirrorPos[14][1][1]=-350.7058;
  fMirrorPos[14][2][0]=5.157848;
  fMirrorPos[14][2][1]=-702.1018;
  fMirrorPos[14][3][0]=-297.7042;
  fMirrorPos[14][3][1]=-877.0388;
  fMirrorPos[14][4][0]=-600.7572;
  fMirrorPos[14][4][1]=-702.0688;
  fMirrorPos[14][5][0]=-600.6882;
  fMirrorPos[14][5][1]=-350.7388;
  fMirrorPos[13][0][0]=314.3028;
  fMirrorPos[13][0][1]=-175.2481;
  fMirrorPos[13][1][0]=617.264;
  fMirrorPos[13][1][1]=-350.2241;
  fMirrorPos[13][2][0]=617.0899;
  fMirrorPos[13][2][1]=-701.700;
  fMirrorPos[13][3][0]=314.1854;
  fMirrorPos[13][3][1]=-876.6871;
  fMirrorPos[13][4][0]=11.19794;
  fMirrorPos[13][4][1]=-701.6591;
  fMirrorPos[13][5][0]=11.25895;
  fMirrorPos[13][5][1]=-350.2651;
  fMirrorPos[16][0][0]=922.5141;
  fMirrorPos[16][0][1]=-175.0949;
  fMirrorPos[16][1][0]=1225.809;
  fMirrorPos[16][1][1]=-350.3329;
  fMirrorPos[16][2][0]=1225.91;
  fMirrorPos[16][2][1]=-701.1989;
  fMirrorPos[16][3][0]=923.0109;
  fMirrorPos[16][3][1]=-876.4109;
  fMirrorPos[16][4][0]=619.8021;
  fMirrorPos[16][4][1]=-701.2689;
  fMirrorPos[16][5][0]=619.6551;
  fMirrorPos[16][5][1]=-350.2629;
  fMirrorPos[15][0][0]=620.0246;
  fMirrorPos[15][0][1]=-703.531;
  fMirrorPos[15][1][0]=923.1095;
  fMirrorPos[15][1][1]=-878.518;
  fMirrorPos[15][2][0]=923.1674;
  fMirrorPos[15][2][1]=-1230.048;
  fMirrorPos[15][3][0]=620.1311;
  fMirrorPos[15][3][1]=-1405.198;
  fMirrorPos[15][4][0]=317.2905;
  fMirrorPos[15][4][1]=-1230.037;
  fMirrorPos[15][5][0]=317.3815;
  fMirrorPos[15][5][1]=-878.529;
  fMirrorPos[22][0][0]=14.00862;
  fMirrorPos[22][0][1]=-704.6452;
  fMirrorPos[22][1][0]=316.9342;
  fMirrorPos[22][1][1]=-879.6272;
  fMirrorPos[22][2][0]=316.9592;
  fMirrorPos[22][2][1]=-1230.817;
  fMirrorPos[22][3][0]=14.05039;
  fMirrorPos[22][3][1]=-1405.803;
  fMirrorPos[22][4][0]=-288.6738;
  fMirrorPos[22][4][1]=-1230.762;
  fMirrorPos[22][5][0]=-288.6788;
  fMirrorPos[22][5][1]=-879.6822;
  fMirrorPos[8][0][0]=-599.8333;
  fMirrorPos[8][0][1]=-704.3596;
  fMirrorPos[8][1][0]=-297.1782;
  fMirrorPos[8][1][1]=-879.2995;
  fMirrorPos[8][2][0]=-297.3792;
  fMirrorPos[8][2][1]=-1230.788;
  fMirrorPos[8][3][0]=-600.1864;
  fMirrorPos[8][3][1]=-1405.685;
  fMirrorPos[8][4][0]=-903.0912;
  fMirrorPos[8][4][1]=-1230.608;
  fMirrorPos[8][5][0]=-902.8402;
  fMirrorPos[8][5][1]=-879.4796;
  fMirrorPos[23][0][0]=9.536897;
  fMirrorPos[23][0][1]=352.1638;
  fMirrorPos[23][1][0]=311.489;
  fMirrorPos[23][1][1]=177.4148;
  fMirrorPos[23][2][0]=311.276;
  fMirrorPos[23][2][1]=-173.6332;
  fMirrorPos[23][3][0]=9.081099;
  fMirrorPos[23][3][1]=-348.3822;
  fMirrorPos[24][0][0]=2.469554;
  fMirrorPos[24][0][1]=352.7277;
  fMirrorPos[24][1][0]=-299.931;
  fMirrorPos[24][1][1]=177.9067;
  fMirrorPos[24][2][0]=-299.835;
  fMirrorPos[24][2][1]=-173.3314;
  fMirrorPos[24][3][0]=2.710454;
  fMirrorPos[24][3][1]=-348.1524;
  fMirrorPos[2][0][0]=0;
  fMirrorPos[2][0][1]=350.;
  fMirrorPos[2][1][0]=315.;
  fMirrorPos[2][1][1]=175.;
  fMirrorPos[2][2][0]=315.;
  fMirrorPos[2][2][1]=-175.;
  fMirrorPos[2][3][0]=0;
  fMirrorPos[2][3][1]=-350;
  fMirrorPos[2][4][0]=-315;
  fMirrorPos[2][4][1]=-175.;
  fMirrorPos[2][5][0]=-315;
  fMirrorPos[2][5][1]=175.;
}

/////////////////////////////////////
// Trigger Mask                    //
// 1 RICH                          //
// 2 RICHx!MUV3 or RICHx!MUV3x!LAV //
// 3 RICHx!MUV3x!LAV               //
// 4 RICHx2 MUV3                   //
/////////////////////////////////////
Int_t TriggerMask(Int_t triggerMask) {
  if (triggerMask&1) return 1;
  if (triggerMask&2||triggerMask&4) return 2;
  if (triggerMask&4) return 3;
  if (triggerMask&3) return 4;
}
/*
///////////////////////////
// Blue field correction //
///////////////////////////
TVector3 AnalysisTools_G::BlueFieldCorrection(Particle *track, Double_t zvertex, Double_t mass) {
  Double_t zmax = 183311.;
  Double_t zmin = 101800.;
  Double_t gtom = 1000.;
  TVector3 initPos = track->GetPosition();
  if (initPos.Z()>zmax) {
    initPos.SetX(track->GetPosition().X()+(track->GetMomentum().X()/track->GetMomentum().Z())*(zmax-track->GetPosition().Z()));    
    initPos.SetY(track->GetPosition().Y()+(track->GetMomentum().Y()/track->GetMomentum().Z())*(zmax-track->GetPosition().Z()));    
    initPos.SetZ(zmax);
  }
  fTracker->SetCharge(1); 
  fTracker->SetInitialPosition(initPos); // Unit: mm
  fTracker->SetInitialMomentum(track->GetMomentum().X()*gtom,track->GetMomentum().Y()*gtom,track->GetMomentum().Z()*gtom); // Unit: MeV
  Double_t zEndTracker = zvertex;
  if (zvertex<zmin) zEndTracker = zmin;
  if (zvertex>=zmax) zEndTracker = zmax;
  fTracker->SetZFinal(zEndTracker); // mm
  fTracker->TrackParticle(); // 
  TVector3 vertexCorr = fTracker->GetFinalPosition(); // New vertex
  TVector3 correctedMomentum = fTracker->GetFinalMomentum();
  TLorentzVector momentumCorr;
  momentumCorr.SetXYZM(correctedMomentum.X()/gtom,correctedMomentum.Y()/gtom,correctedMomentum.Z()/gtom,mass); // New 4-momentum
  if (zvertex<zmin) vertexCorr = vertexCorr+correctedMomentum*(1./correctedMomentum.Z())*(zvertex-vertexCorr.Z());
  track->SetMomentum(momentumCorr); 
  return vertexCorr;
} 

TVector3 AnalysisTools_G::BlueFieldCorrection(TVector3 *pmom, TVector3 initPos, Int_t charge, Double_t zvertex) {
  Double_t zmax = 183311.;
  Double_t zmin = 101800.;
  Double_t gtom = 1000.;
  if (initPos.Z()>zmax) {
    initPos.SetX(initPos.X()+(pmom->X()/pmom->Z())*(zmax-initPos.Z()));    
    initPos.SetY(initPos.Y()+(pmom->Y()/pmom->Z())*(zmax-initPos.Z()));    
    initPos.SetZ(zmax);
  }
  if (initPos.Z()<zmin) {
    initPos.SetX(initPos.X()+(pmom->X()/pmom->Z())*(zmin-initPos.Z()));    
    initPos.SetY(initPos.Y()+(pmom->Y()/pmom->Z())*(zmin-initPos.Z()));    
    initPos.SetZ(zmin);
  }
  fTracker->SetCharge(charge); 
  fTracker->SetInitialPosition(initPos); // Unit: mm
  fTracker->SetInitialMomentum(pmom->X()*gtom,pmom->Y()*gtom,pmom->Z()*gtom); // Unit: MeV
  Double_t zEndTracker = zvertex;
  if (zvertex<zmin) zEndTracker = zmin;
  if (zvertex>=zmax) zEndTracker = zmax;
  fTracker->SetZFinal(zEndTracker); // mm
  fTracker->TrackParticle(); // 
  TVector3 vertexCorr = fTracker->GetFinalPosition(); // New vertex
  TVector3 correctedMomentum = fTracker->GetFinalMomentum();
  pmom->SetXYZ(correctedMomentum.X()/gtom,correctedMomentum.Y()/gtom,correctedMomentum.Z()/gtom); // New momentum
  if (zvertex<zmin) vertexCorr = vertexCorr+correctedMomentum*(1./correctedMomentum.Z())*(zvertex-vertexCorr.Z());
  return vertexCorr;
} 
*/
