#include "MUV2Analysis.hh"
#include "AnalysisTools.hh"
#include "TRecoMUV2Event.hh"
#include "Event.hh"

MUV2Analysis::MUV2Analysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}

void MUV2Analysis::Set(TRecoMUV2Event* fevent, Double_t freftime, Double_t ftLow, Double_t ftHi, bool fMCflag ){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
MCflag = fMCflag;
TRecoMUV2Candidate* MUV2Cand;
for(int i=0; i<event->GetNCandidates(); i++){
  MUV2Cand = (TRecoMUV2Candidate*)event->GetCandidate(i);
  Double_t Tdiff = reftime - MUV2Cand->GetTime();
  fUserMethods->FillHisto("MUV2_Set_reft", Tdiff);
}
}

Double_t MUV2Analysis::TotalEnergy(){
Double_t TotalEnergy = 0;
TRecoMUV2Candidate *MUV2Cand;
for(int i=0; i<event->GetNCandidates(); i++) {
 MUV2Cand = (TRecoMUV2Candidate*)event->GetCandidate(i);
 Double_t Tdiff = reftime - MUV2Cand->GetTime();
 if(!MCflag){if(Tdiff<tLow||Tdiff>tHi)continue;}
 TotalEnergy = TotalEnergy + MUV2Cand->GetEnergy(); 
}
return TotalEnergy;

}

Double_t MUV2Analysis::TotalEnergyHits(){
Double_t TotalEnergy = 0;
TClonesArray& Hits = (*(event->GetHits()));
 for (Int_t i=0; i<event->GetNHits(); i++) {
   TRecoMUV2Hit *hit = (TRecoMUV2Hit*)Hits[i];
   hit = (TRecoMUV2Hit*)event->GetHit(i);
   Double_t Tdiff = reftime - hit->GetTime();
   if(!MCflag){if(Tdiff<tLow||Tdiff>tHi)continue;}
   TotalEnergy = TotalEnergy + hit->GetEnergy();
 }
return TotalEnergy;
}

Int_t MUV2Analysis::MatchedTrack(TVector3 Pos){
                MatchedClusterE = 0;
                TRecoMUV2Candidate *MUV2Cand;
                TVector3 Pion_Prop_MUV2=Pos;
                int MUV2 = 0;
                Int_t MUV2_minD;
                MUV2ID = -1;
                TVector2 Pion_Prop_MUV2_2D(Pion_Prop_MUV2.X(), Pion_Prop_MUV2.Y());
                double distanceMUV2 = 99999;
                for(int i=0; i<event->GetNCandidates(); i++)        {
               MUV2Cand = (TRecoMUV2Candidate*)event->GetCandidate(i);
                double Tdiff_MUV2 = reftime - MUV2Cand->GetTime();
               if(!MCflag){if(Tdiff_MUV2 < tLow || Tdiff_MUV2 > tHi) continue;}
                TVector2 HitPos(MUV2Cand->GetPosition().X() , MUV2Cand->GetPosition().Y());
                TVector2 distMUV2 = Pion_Prop_MUV2_2D - HitPos;
                if (distMUV2.Mod() < distanceMUV2)      {
                        distanceMUV2 = distMUV2.Mod();
                        if (distanceMUV2<600)   {
                        MUV2++;
                        MUV2_minD = i;
                                                }
                                                        }
                                                                   }
	if(distanceMUV2<99999)fUserMethods->FillHisto("MUV2_MatchedTrack_MinD", distanceMUV2);
                if(MUV2>0){
                        MUV2Cand = (TRecoMUV2Candidate*)event->GetCandidate(MUV2_minD);
                        MatchedClusterE=MUV2Cand->GetEnergy();
//                        fUserMethods->FillHisto("MUV2_MatchedTrack_MinD", distanceMUV2);
                        MUV2ID = MUV2_minD;
			Distance = distanceMUV2; 
                        return 1;
 }
                else return 0;
}


Int_t MUV2Analysis::MuonID(){ 
int MUV2_N = 0;
TRecoMUV2Candidate *MUV2Cand;
for(int i=0; i<event->GetNCandidates(); i++) {
 MUV2Cand = (TRecoMUV2Candidate*)event->GetCandidate(i);
 Double_t Tdiff = reftime - MUV2Cand->GetTime();
 if(!MCflag){if(Tdiff<tLow||Tdiff>tHi)continue;}
 if(MUV2Cand->GetEnergy()<1.5)MUV2_N++;
}
if (MUV2_N>1)return 0;
else return 1;

}

void MUV2Analysis::BookHistos_MuonID(){
}
void MUV2Analysis::BookHistos_MatchedTrack(){
fUserMethods->BookHisto(new TH1F("MUV2_MatchedTrack_MinD","Min distance MUV2 track", 1500, 0,1500));
}

void MUV2Analysis::BookHistos_Set(){
fUserMethods->BookHisto(new TH1I("MUV2_Set_reft", "reference time", 1200, -50, 50));
}

void MUV2Analysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
}
