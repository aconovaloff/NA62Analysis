#include "AnalysisTools.hh"
//#include "BlueTubeTracker.hh"
/*TVector3 AnalysisTools::*/TVector3 prop(TVector3 InitPos, TVector3 InitP, double PropPosZ)
{
                                TVector3 EndPos;
                                double MagStartPos = 197000;    //in mm, actually center of magnet (1.3 m thickness)//
                                double kick = -270; // -270; //in MeV//
//                              TVector3 InitP = fMCSimple["pi+"][0]->GetInitialMomentum();
//                              TVector3 InitPos = fMCSimple["pi+"][0]->GetProdPos().Vect();
//                              double PropPosZ = LKrStartPos*10;                                               

                                double propDist = PropPosZ - InitPos.Z();
//                              TVector3 PiPlusMom(InitP.X(), InitP.Y(), InitP.Z()); 
                                TVector3 PartMomDir = InitP.Unit();
              //                          if (propDist < PropPosZ - MagStartPos){
              				  if (PropPosZ < MagStartPos){
                                                double propVectorMag = propDist/PartMomDir.Z();
                                                TVector3 propVector = PartMomDir*propVectorMag;
                                                EndPos = propVector + InitPos;
                                                                                    }
                                        else                                        {
                                                double propDist2 = (MagStartPos - InitPos.Z());
                                                double propVectorMag2 = propDist2/PartMomDir.Z();
                                                TVector3 propVector2 = PartMomDir*propVectorMag2;
                                                InitP.SetX(InitP.X() + kick);
                                                double propDist = PropPosZ - MagStartPos;
                                                TVector3 PartMomDir = InitP.Unit();
                                                double propVectorMag = propDist/PartMomDir.Z();
                                                TVector3 propVector = PartMomDir*propVectorMag;
                                                EndPos = propVector + propVector2 + InitPos;
                                                                                   }
                                return EndPos;
}
TVector3 prop_neg(TVector3 InitPos, TVector3 InitP, double PropPosZ)
{
                                TVector3 EndPos;
                                double MagStartPos = 197000;    //in mm, actually center of magnet (1.3 m thickness)//
                                double kick = 270; // -270; //in MeV//
                                double propDist = PropPosZ - InitPos.Z();
                                TVector3 PartMomDir = InitP.Unit();
                                          if (PropPosZ < MagStartPos){
                                                double propVectorMag = propDist/PartMomDir.Z();
                                                TVector3 propVector = PartMomDir*propVectorMag;
                                                EndPos = propVector + InitPos;
                                                                                    }
                                        else                                        {
                                                double propDist2 = (MagStartPos - InitPos.Z());
                                                double propVectorMag2 = propDist2/PartMomDir.Z();
                                                TVector3 propVector2 = PartMomDir*propVectorMag2;
                                                InitP.SetX(InitP.X() + kick);
                                                double propDist = PropPosZ - MagStartPos;
                                                TVector3 PartMomDir = InitP.Unit();
                                                double propVectorMag = propDist/PartMomDir.Z();
                                                TVector3 propVector = PartMomDir*propVectorMag;
                                                EndPos = propVector + propVector2 + InitPos;
                                                                                   }
                                return EndPos;
}
/*
TVector3 prop_neg(TVector3 InitPos, TVector3 InitP, double PropPosZ)
{
                                TVector3 EndPos;
                                double MagStartPos = 197000;    //in mm, actually center of magnet (1.3 m thickness)//
                                double kick = 270; // -270; //in MeV//
                                double propDist = PropPosZ - InitPos.Z(); 
                                TVector3 PartMomDir = InitP.Unit();
                                        if (propDist < PropPosZ - MagStartPos){
                                                double propVectorMag = propDist/PartMomDir.Z();
                                                TVector3 propVector = PartMomDir*propVectorMag;
                                                EndPos = propVector + InitPos;
                                                                                    }
                                        else                                        {
                                                double propDist2 = (MagStartPos - InitPos.Z());
                                                double propVectorMag2 = propDist2/PartMomDir.Z();
                                                TVector3 propVector2 = PartMomDir*propVectorMag2;
                                                InitP.SetX(InitP.X() + kick);
                                                double propDist = PropPosZ - MagStartPos;
                                                TVector3 PartMomDir = InitP.Unit();
                                                double propVectorMag = propDist/PartMomDir.Z();
                                                TVector3 propVector = PartMomDir*propVectorMag;
                                                EndPos = propVector + propVector2 + InitPos;
                                                                                   }
                                return EndPos;
}
*/
Double_t /*AnalysisTools::*/LKrZeroSup(Double_t cells, Double_t energy)
{

   Double_t ue = energy;
  Double_t ce;
  if(cells > 5) {
        if (ue<25.) ce = ue/(0.724761+0.0220857*ue-0.000873029*ue*ue+0.0000126505*ue*ue*ue);
        if (ue>=25. && ue<65.) ce = ue/(0.928916625+0.00262869*(ue-25.)-0.000063694*(ue-25.)*(ue-25.)+0.00000068890*(ue-25.)*(ue-25.)*(ue-25.));
        if (ue>65.) ce = ue/0.976243425;
                }
  else ce = ue;
  return ce;
}

TVector3 /*AnalysisTools::*/Vertex(TVector3 vect1, TVector3 vect2, TVector3 vectb, TVector3 vectc)
{
  TVector3 v1 = vect1;
  TVector3 v2 = vect2;
  TVector3 pos1 = vectb;
  TVector3 pos2 = vectc;
  TVector3 r12 = pos1-pos2;
  Double_t v1xv2 = v1.Dot(v2);
  Double_t det   = pow(v1xv2,2)-v1.Mag2()*v2.Mag2();
  if (!det) return TVector3(-9999,-9999,-9999);
  Double_t t1 = (v2.Mag2()*r12.Dot(v1)-v1.Dot(v2)*r12.Dot(v2))/det;
  Double_t t2 = (v1.Dot(v2)*r12.Dot(v1)-v1.Mag2()*r12.Dot(v2))/det;
  TVector3 q1 = pos1+t1*v1;
  TVector3 q2 = pos2+t2*v2;
  TVector3 vertex = 0.5*(q1+q2);
  Double_t cda = (q1-q2).Mag();
  return vertex;
}

Double_t /*AnalysisTools::*/cda(TVector3 vect1, TVector3 vect2, TVector3 vectb, TVector3 vectc)
{
  TVector3 v1 = vect1;
  TVector3 v2 = vect2;
  TVector3 pos1 = vectb;
  TVector3 pos2 = vectc;
  TVector3 r12 = pos1-pos2;
  Double_t v1xv2 = v1.Dot(v2);
  Double_t det   = pow(v1xv2,2)-v1.Mag2()*v2.Mag2();
  if (!det) return -99999;
  Double_t t1 = (v2.Mag2()*r12.Dot(v1)-v1.Dot(v2)*r12.Dot(v2))/det;
  Double_t t2 = (v1.Dot(v2)*r12.Dot(v1)-v1.Mag2()*r12.Dot(v2))/det;
  TVector3 q1 = pos1+t1*v1;
  TVector3 q2 = pos2+t2*v2;
  TVector3 vertex = 0.5*(q1+q2);
  Double_t cda = (q1-q2).Mag();
  return cda;
}

TVector3 /*AnalysisTools::*/Propagate(Int_t charge, TVector3 *momentum, TVector3 *position, Double_t fZEnd)
{

  Double_t fEC = TMath::C()* 1.e-9 * 1.e-4 * 1.e-2;
  TVector3  fB;
  fB.SetXYZ(0.,0.6928*10000,0.);

  TVector3 fPosExp;
  Double_t dMag = 0.1*1300;
  Double_t zMag = 0.1*196345;
  Double_t fStartX = 0.1*position->X();
  Double_t fStartY = 0.1*position->Y();
  Double_t fStartZ = 0.1*position->Z();
  Double_t fPartP = momentum->Mag();
  Double_t fPartThetaX = momentum->X()/momentum->Z();
  Double_t fPartThetaY = momentum->Y()/momentum->Z();
  fZEnd *= 0.1;
  Int_t fPartQ = 1;

  // fZEnd before magnet
  if (fZEnd<=zMag)
  {
    fPosExp.SetX(fStartX+fPartThetaX*(fZEnd-fStartZ));
    fPosExp.SetY(fStartY+fPartThetaY*(fZEnd-fStartZ));
    fPosExp.SetZ(fZEnd);
    return fPosExp*10;
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

  return fPosExp*10;
}

                                                                                                                             
TVector3 backprop(TVector3 InitPos, TVector3 InitP, double PropPosZ)
{
                                TVector3 EndPos;
                                double propDist =  InitPos.Z() - PropPosZ;
                                TVector3 PartMomDir = InitP.Unit();

                                                double propVectorMag = propDist/PartMomDir.Z();
                                                TVector3 propVector = PartMomDir*propVectorMag;
                                                EndPos = InitPos - propVector;
                                                                                  
                                return EndPos;
}
/*
TVector3 BlueTubeProp(BlueTubeTracker *ftracker, TVector3 K_pos, TVector3 DecayVertex, TVector3 initPos, TVector3 PionMomentum){
 Double_t zvertex = DecayVertex.Z();
 Double_t zmax = 183311.;
  Double_t zmin = 101800.;
  Double_t gtom = 1000.;
  TVector3 initPos = PionPosition;
  if (initPos.Z()>zmax) {
    initPos.SetX(PionPosition.X()+(PionMomentum.X()/PionMomentum.Z())*(zmax-PionPosition.Z()));
    initPos.SetY(PionPosition.Y()+(PionMomentum.Y()/PionMomentum.Z())*(zmax-PionPosition.Z()));
    initPos.SetZ(zmax);
  }
  fTracker->SetCharge(1);
  fTracker->SetInitialPosition(initPos); // Unit: mm
  fTracker->SetInitialMomentum(PionMomentum.X()*gtom,PionMomentum.Y()*gtom,PionMomentum.Z()*gtom); // Unit: MeV
  Double_t zEndTracker = zvertex;
  if (zvertex<zmin) zEndTracker = zmin;
  if (zvertex>=zmax) zEndTracker = zmax;
  fTracker->SetZFinal(zEndTracker); // mm
  fTracker->TrackParticle(); // 
  TVector3 vertexCorr = fTracker->GetFinalPosition(); // New vertex
  TVector3 correctedMomentum = fTracker->GetFinalMomentum();
  TLorentzVector pionMomentumCorr;
pionMomentumCorr.SetXYZM(correctedMomentum.X()/gtom,correctedMomentum.Y()/gtom,correctedMomentum.Z()/gtom,0.13957018); // New pion momentum
  if (zvertex<zmin) vertexCorr = vertexCorr+correctedMomentum*(1./correctedMomentum.Z())*(zvertex-vertexCorr.Z());
TVector3 K_pos= Propagate(1,&KaonMomentum,&KaonPosition, Vert.Z());
TVector2 cda_new_vect(K_pos.X()-vertexCorr.X(), K_pos.Y()-vertexCorr.Y());
*/

/*
TVector3 InitPos = fInitPos;
TVector3 InitP = fInitP;
Double_t PropPosZ = fPropPosZ; 
BlueTubeTracker *tracker = new BlueTubeTracker();
tracker->SetCharge(+1);
tracker->SetInitialPosition(InitPos.X(), InitPos.Y(), InitPos.Z());  // Unit: mm
tracker->SetInitialMomentum(InitP.X(), InitP.Y(),  InitP.Z());  // Unit: MeV/c
tracker->SetZFinal(PropPosZ);                         // Unit: mm
tracker->TrackParticle();
TVector3 NewMomentum = tracker->GetFinalMomentum(); // Unit: MeV/c
TVector3 NewPosition = tracker->GetFinalPosition(); // Unit: mm
return NewPosition;
 
}
*/
