#include "MUV1Analysis.hh"
#include "AnalysisTools.hh"
#include "TRecoMUV1Event.hh"
#include "Event.hh"

MUV1Analysis::MUV1Analysis(NA62Analysis::Core::BaseAnalysis *ba){
fUserMethods = new UserMethods(ba);
}

void MUV1Analysis::Set(TRecoMUV1Event* fevent, double freftime, double ftLow, double ftHi, bool fMCflag ){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
MCflag = fMCflag;
TRecoMUV1Candidate* MUV1Cand;
for(int i=0; i<event->GetNCandidates(); i++){
  MUV1Cand = (TRecoMUV1Candidate*)event->GetCandidate(i);
  Double_t Tdiff = reftime - MUV1Cand->GetTime();
  fUserMethods->FillHisto("MUV1_Set_reft", Tdiff);
}
}

Int_t MUV1Analysis::MuonID(){ 
int MUV1_N = 0;
TRecoMUV1Candidate *MUV1Cand;
for(int i=0; i<event->GetNCandidates(); i++) {
 MUV1Cand = (TRecoMUV1Candidate*)event->GetCandidate(i);
 double Tdiff = reftime - MUV1Cand->GetTime();
// fUserMethods->FillHisto("MUV1_MuonID_reft", Tdiff); 
if (!MCflag){ if(Tdiff<tLow||Tdiff>tHi)continue;}
 if(MUV1Cand->GetEnergy()<1.5)MUV1_N++;
}
if (MUV1_N>1)return 0;
else return 1;

}

Double_t MUV1Analysis::TotalEnergyHits(){
Double_t TotalEnergy = 0;
TClonesArray& Hits = (*(event->GetHits()));
 for (Int_t i=0; i<event->GetNHits(); i++) {
   TRecoMUV1Hit *hit = (TRecoMUV1Hit*)Hits[i];
   hit = (TRecoMUV1Hit*)event->GetHit(i);
   Double_t Tdiff = reftime - hit->GetTime();
   if(!MCflag){if(Tdiff<tLow||Tdiff>tHi)continue;}
   TotalEnergy = TotalEnergy + hit->GetEnergy();
 }
return TotalEnergy;
}


Int_t MUV1Analysis::MatchedTrack(TVector3 Pos){
		MatchedClusterE = 0; 
		TRecoMUV1Candidate *MUV1Cand;
                TVector3 Pion_Prop_MUV1=Pos;
                int MUV1 = 0;
                Int_t MUV1_minD;
		MUV1ID = -1;
                TVector2 Pion_Prop_MUV1_2D(Pion_Prop_MUV1.X(), Pion_Prop_MUV1.Y());
                double distanceMUV1 = 99999;
                for(int i=0; i<event->GetNCandidates(); i++)        {
               MUV1Cand = (TRecoMUV1Candidate*)event->GetCandidate(i);
                double Tdiff_MUV1 = reftime - MUV1Cand->GetTime();
               if(!MCflag){if(Tdiff_MUV1 < tLow || Tdiff_MUV1 > tHi) continue;}
                TVector2 HitPos(MUV1Cand->GetPosition().X() , MUV1Cand->GetPosition().Y());
                TVector2 distMUV1 = Pion_Prop_MUV1_2D - HitPos;
                if (distMUV1.Mod() < distanceMUV1)      {
                        distanceMUV1 = distMUV1.Mod();
                        if (distanceMUV1<600)   {
                        MUV1++;
                        MUV1_minD = i;
                                                }
                                                        }
                                                                   }
if(distanceMUV1<99999)fUserMethods->FillHisto("MUV1_MatchedTrack_MinD", distanceMUV1);
		if(MUV1>0){
			MUV1Cand = (TRecoMUV1Candidate*)event->GetCandidate(MUV1_minD);
			MatchedClusterE=MUV1Cand->GetEnergy(); 
//			fUserMethods->FillHisto("MUV1_MatchedTrack_MinD", distanceMUV1); 
			MUV1ID = MUV1_minD; 
			Distance = distanceMUV1; 
			return 1;
 }
		else return 0;
}

void MUV1Analysis::BookHistos_MatchedTrack(){
fUserMethods->BookHisto(new TH1F("MUV1_MatchedTrack_MinD","Min distance MUV1 track", 1500, 0,1500)); 
}
void MUV1Analysis::BookHistos_MuonID(){
//fUserMethods->BookHisto(new TH1F("MUV1_MuonID_reft","reference time", 200, -100, 100));
}
void MUV1Analysis::BookHistos_Set(){
fUserMethods->BookHisto(new TH1F("MUV1_Set_reft","reference time", 800, -100, 100));
}

void MUV1Analysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
}
