#include "LKrAnalysis.hh"
#include "AnalysisTools.hh"
#include "TRecoLKrEvent.hh"
#include "Event.hh"
#include "UserMethods.hh"
using namespace TMath;
LKrAnalysis::LKrAnalysis(NA62Analysis::Core::BaseAnalysis *ba ){
fUserMethods = new UserMethods(ba); 
}

void LKrAnalysis::Set(TRecoLKrEvent* fevent, double freftime, double ftLow, double ftHi, bool fMCflag ){
event = fevent;
reftime = freftime;
tLow = ftLow;
tHi = ftHi;
fMCflag = MCflag;
clusterE = 0;
TRecoLKrCandidate* LKrCand;
MatchedPion_unmatched=-1;  
MatchedMuon_unmatched=-1;
for(int i=0; i<event->GetNCandidates(); i++){
  LKrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
  Double_t Tdiff = reftime - LKrCand->GetTime();
  fUserMethods->FillHisto("LKr_Set_reft", Tdiff);
}
}

Int_t LKrAnalysis::MuonID(TVector3 Pos){
				TRecoLKrCandidate *LKrCand;
                                TVector3 Pion_Prop_LKr2 = Pos;
                                double mindist = 99999;
                                double LKrStartPos = 240413;
                                int m = 0;
                                int unmatched = 0;
                                int minDCand;  				
				double EoverP=-1;
                                TVector2 Pion_Prop_LKr2D(Pion_Prop_LKr2.X(), Pion_Prop_LKr2.Y());
                                for(int i=0; i<event->GetNCandidates(); i++) {
                                LKrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
				double Tdiff_LKr = reftime-LKrCand->GetTime();
 				if(!MCflag){if(Tdiff_LKr < tLow || Tdiff_LKr > tHi)continue;}
                                TVector2 ClusterPos(LKrCand->GetClusterX(), LKrCand->GetClusterY());
                                TVector2 dist = Pion_Prop_LKr2D - ClusterPos;
                                Double_t ce = LKrCand->GetClusterEnergy();
 //                             Double_t cells = LKrCand->GetNCells();
   //                           Double_t ce = LKrZeroSup(cells, ue);
					if (dist.Mod() > 150) unmatched++;
                                        if (dist.Mod() > 150) continue;
                                        if (dist.Mod() < mindist)       {
                                                mindist = dist.Mod();
                                                minDCand = i;
                                                m++;
                                                                        }

                                                                                }
                                Double_t ce;
//				if (m==0||unmatched>0) return 0;
				if(unmatched>0) return 0; 
                                if(m>0) {
                                LKrCand = (TRecoLKrCandidate*)event->GetCandidate(minDCand);
//                                Double_t ue = LKrCand->GetClusterEnergy();
//                                Double_t cells = LKrCand->GetNCells();
                                clusterE = LKrCand->GetClusterEnergy();
				if (clusterE > 3) return 0; 
				fUserMethods->FillHisto("LKr_MuonID_clusterE", LKrCand->GetClusterEnergy()*1000);
                                fUserMethods->FillHisto("LKr_MuonID_minD", mindist);  
                                        }
//                                if (ce>1) return 0;
//				else return 1;
				return 1;
}
void LKrAnalysis::BookHistos_MuonID(){
fUserMethods->BookHisto(new TH1I("LKr_MuonID_clusterE", "Matched Cluster Energy", 50, 0, 50000));
fUserMethods->BookHisto(new TH1I("LKr_MuonID_minD", "Matched Cluster Distance", 150, 0, 150));
}

Int_t LKrAnalysis::SinglePionCluster(TVector3 PionMomentum, TVector3 Pos, bool bg_flag){
				RMS = 0;
				SinglePionCluster_unmatched=0;
				clusterE = 0; 
				SinglePionCluster_minID = -1;  
				Double_t P = PionMomentum.Mag(); 
                                TRecoLKrCandidate *LKrCand;
                                TVector3 Pion_Prop_LKr2 = Pos;
                                mindist = 99999;
                                double LKrStartPos = 240413;
                                int m = 0;
                                int unmatched = 0;
                                int minDCand;
                                double EoverP=-1;
                                TVector2 Pion_Prop_LKr2D(Pion_Prop_LKr2.X(), Pion_Prop_LKr2.Y());
                                for(int i=0; i<event->GetNCandidates(); i++) {
                                LKrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
			 	double Tdiff_LKr = reftime-LKrCand->GetTime();
// 				if(Tdiff_LKr < tLow || Tdiff_LKr > tHi)continue;
                                TVector2 ClusterPos(LKrCand->GetClusterX(), LKrCand->GetClusterY());
                                TVector2 dist = Pion_Prop_LKr2D - ClusterPos;
                                Double_t ce = LKrCand->GetClusterEnergy();
 				 if (dist.Mod() > 150)   {
                                                if(ce>=2&&((Tdiff_LKr>=-15&&Tdiff_LKr<=10)||(Tdiff_LKr>=22&&Tdiff_LKr<=28)))SinglePionCluster_unmatched++;
                                                if(ce<=2&&(Tdiff_LKr>=-5&&Tdiff_LKr<=5))SinglePionCluster_unmatched++;
                                                         }
                                if(!MCflag){if(Tdiff_LKr < tLow || Tdiff_LKr > tHi)continue;}
//                                        if (dist.Mod() > 150) SinglePionCluster_unmatched++;
                                        if (dist.Mod() < mindist)       {
                                                mindist = dist.Mod();
                                                minDCand = i;
                                                if (dist.Mod()<150)m++;
                                                                        }

                                                                                }
				fUserMethods->FillHisto("LKr_SinglePionCluster_MinD", mindist);
				
                                if(m==1) {
                                LKrCand = (TRecoLKrCandidate*)event->GetCandidate(minDCand);
				SinglePionCluster_minID = minDCand; 
                                Double_t ce = LKrCand->GetClusterEnergy();
				clusterE = ce;
 				EoverP = ce/P;
 fUserMethods->FillHisto("LKr_SinglePionCluster_EtoP", EoverP);
				RMS = sqrt(LKrCand->GetClusterRMSY()*LKrCand->GetClusterRMSY()+LKrCand->GetClusterRMSX()*LKrCand->GetClusterRMSX());
                                        }
                                if (m==1&&EoverP<=0.8)return 1;
				else if (m==1&&EoverP>0.8)return 2;  
				else if(m>1)return 3; 
//if (!bg_flag){				if(unmatched>0) return 0;}
//				if(m>1)return 0; 
                                else return 0;
}

Int_t LKrAnalysis::SinglePhoton(TVector3 PionMomentum, TVector3 Pos, TVector3 DecayVert, TLorentzVector KaonFourMomentum, TLorentzVector PionFourMomentum){
				Double_t P = PionMomentum.Mag(); 
				Int_t photonID= -1; 
                                TRecoLKrCandidate *LKrCand;
                                TVector3 Pion_Prop_LKr2 = Pos;
                                double mindist = 99999;
                                double LKrStartPos = 240413;
                                int m = 0;
                                int unmatched = 0;
                                int minDCand;
                                double EoverP=-1;
                                TVector2 Pion_Prop_LKr2D(Pion_Prop_LKr2.X(), Pion_Prop_LKr2.Y());
                                for(int i=0; i<event->GetNCandidates(); i++) {
                                LKrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
			 	double Tdiff_LKr = reftime-LKrCand->GetTime();
 				if(Tdiff_LKr < tLow || Tdiff_LKr > tHi)continue;
                                TVector2 ClusterPos(LKrCand->GetClusterX(), LKrCand->GetClusterY());
                                TVector2 dist = Pion_Prop_LKr2D - ClusterPos;
                                Double_t ce = LKrCand->GetClusterEnergy();
                                        if (dist.Mod() > 150&&ce>3){ 
						unmatched++;
						photonID = i;
								   }
                                        if (dist.Mod() < mindist)       {
                                                mindist = dist.Mod();
                                                minDCand = i;
                                                if (dist.Mod()<150)m++;
                                                                        }

                                                                                }
//				fUserMethods->FillHisto("LKr_SinglePionCluster_MinD", mindist);
				if(m!=0||unmatched!=1)return 0;
				LKrCand = (TRecoLKrCandidate*)event->GetCandidate(photonID);
				Double_t ce = LKrCand->GetClusterEnergy();
                                TVector3 PosPhoton((LKrCand->GetClusterX()), LKrCand->GetClusterY(), LKrStartPos);
                                TVector3 Pos2 = PosPhoton - DecayVert;
                                TVector3 Pdir = Pos2.Unit();
                                double Px = Pdir.X()*ce;
                                double Py = Pdir.Y()*ce;
                                double Pz = Pdir.Z()*ce;
                                TLorentzVector gFourMom(Px, Py, Pz, ce);
				TLorentzVector p2 = KaonFourMomentum - PionFourMomentum - gFourMom;
//				TVector3 K_pos= Propagate(1,&KaonMomentum,&KaonPosition, Vert.Z());
				TVector3 p2three = p2.Vect();
				TVector3 p2pos = Propagate(0, &p2three, &DecayVert, LKrStartPos); 
cout<<"gothere";
				if (sqrt(p2pos.X()*p2pos.X()+p2pos.Y()*p2pos.Y())<1200)return 0; //missing LKr? don't know exact dimension, 1130 is detector acceptance
cout<<"gothere2";				
				TLorentzVector npion = p2 + gFourMom;
				Double_t miss = p2.E() - p2three.Mag(); 
				TVector3 npionvect = npion.Vect();
//		Double_t miss = npion.E()-npion.Mag(); 
				fUserMethods->FillHisto("EPmiss", miss); 
				if (miss>0.0004 || miss< -0.0004)return 0; 
				return 1; 

//                                if(m>0) {
//                                LKrCand = (TRecoLKrCandidate*)event->GetCandidate(minDCand);
//                                Double_t ce = LKrCand->GetClusterEnergy();
//				clusterE = ce;
// 				EoverP = ce/P;
// fUserMethods->FillHisto("LKr_SinglePionCluster_EtoP", EoverP);
//                                        }
//                                if (m>0&&EoverP>0.8)return 0; 
//				if(unmatched>0) return 0;
//                                else return 1;
}

void LKrAnalysis::BookHistos_SinglePhoton(){
fUserMethods->BookHisto(new TH1I("EPmiss", "E miss minus P miss", 100, -0.002, 0.002));
}

void LKrAnalysis::BookHistos_SinglePionCluster(){
fUserMethods->BookHisto(new TH1I("LKr_SinglePionCluster_EtoP", "E/p", 120, 0, 1.2));
//fUserMethods->BookHisto(new TH1I("clusterE", "cluster e", 100, 0, 100));
fUserMethods->BookHisto(new TH1I("LKr_SinglePionCluster_MinD", " Distance of minimum distance track to closest matched LKr Cluster", 200, 0, 200));
}

void LKrAnalysis::SaveAllPlots(){
fUserMethods->SaveAllPlots();
}

Int_t LKrAnalysis::OneNeutralPion(TLorentzVector Pi4, Int_t pflag){
          TLorentzVector PiPlus4Mom = Pi4; 
Int_t ParticleFlag = pflag; //0 for general, 1 for muons, 2 for electrons
	 nPion_flag = false;
	 TRecoLKrCandidate *lkrCand;
	 TRecoLKrCandidate *lkrCand2;
	double EnergyRatioMin = .28;
	double EnergyCellRatioMin = 1.66; 
	double dGammaMin = 20;
	double dDeadCellMin = 20; 
	double LKrOutRadius = 1130;
	double LKrInRadius = 150; 
        double LKrStartPos = 240413;
	double pi0Mass = .1349766;
	double CHODStartPos = 239009;
	double TrimStartPos = 101800;
	double trim = .0012; 
//Horizaont LKr position correction
double Xcorr = 0; //add
	double Ecorr = 1; //multiply
        Double_t tCedarhi =0.5;//
        Double_t tCedarlow =-0.7;//
        Double_t tCedarmid =-0.13; //
        double LAV_tLow = -3;
        double LAV_tHi = 4;
        double MUV3_tLow = -2;
        double MUV3_tHi =1;
        Double_t tCHANTILow = -6;
        Double_t tCHANTIHi = 2;
        Double_t tMUV1Low = -14;
        Double_t tMUV1Hi = 6;
        Double_t tMUV2Low = -5;
        Double_t tMUV2Hi = 5;
        Double_t tRICHLow = -1;
        Double_t tRICHHi = 0.5;
vector<int> CandID;
	vector<int> CandID2;

	for(int i=0; i<event->GetNCandidates(); i++){
				lkrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
				Double_t ce = lkrCand->GetClusterEnergy();
//      				Double_t cells = lkrCand->GetNCells();
//				Double_t ce = LKrZeroSup(cells, ue); 
//detector acceptance
if(sqrt((lkrCand->GetClusterX()+Xcorr)*(lkrCand->GetClusterX()+Xcorr)+lkrCand->GetClusterY()*lkrCand->GetClusterY())>LKrOutRadius || sqrt((lkrCand->GetClusterX()+Xcorr)*(lkrCand->GetClusterX()+Xcorr)+lkrCand->GetClusterY()*lkrCand->GetClusterY())<LKrInRadius)continue;				
					fUserMethods->FillHisto("LKr_OneNeutralPion_CellEnergy", ce*Ecorr, lkrCand->GetNCells());
					if(lkrCand->GetNCells()/ce > EnergyCellRatioMin) continue;
					 fUserMethods->FillHisto("LKr_OneNeutralPion_CellEnergy_Post", ce*Ecorr, lkrCand->GetNCells());
 					if(lkrCand->GetClusterDDeadCell() < dDeadCellMin) continue;
					if(ce < 3) continue;
					CandID.push_back(i);
						      	} //LKr for loop
	vector<TLorentzVector> gamma4;
	vector<TVector2> gammaPos;
	double Pi0Z;
	double decayZ;
//find decay vertex
bool n = false;
	vector<int> pairID;
//	TVector3 DecayVert;  
	if (CandID.size()<2)return 0; 
		for (int ii = 0; ii < CandID.size() - 1; ii++){
		lkrCand = (TRecoLKrCandidate*)event->GetCandidate(CandID[ii]);
	//	if ((lkrCand->GetClusterX()+Xcorr)>-43.5&&(lkrCand->GetClusterX()+Xcorr)<-39.5&&lkrCand->GetClusterY()>-39.6&&lkrCand->GetClusterY()<-37.5) continue;
		 	for (int jj = ii + 1; jj < CandID.size(); jj++){
				lkrCand2 = (TRecoLKrCandidate*)event->GetCandidate(CandID[jj]);
//				if ((lkrCand2->GetClusterX()+1)>-43.5&&(lkrCand2->GetClusterX()+1)<-39.5&&lkrCand2->GetClusterY()>-39.6&&lkrCand2->GetClusterY()<-37.5) continue;			
				double dgamma = sqrt(((lkrCand->GetClusterX()+Xcorr)-(lkrCand2->GetClusterX()))*((lkrCand->GetClusterX()+Xcorr)-(lkrCand2->GetClusterX())) + (lkrCand->GetClusterY()-lkrCand2->GetClusterY())*(lkrCand->GetClusterY()-lkrCand2->GetClusterY()));
				if (dgamma < dGammaMin) continue; 
	
				Double_t ce = lkrCand->GetClusterEnergy();//changed from ue to accomodate new software
 //     				Double_t cells = lkrCand->GetNCells();
//				Double_t ce = LKrZeroSup(cells, ue); 


				Double_t ce2 = lkrCand2->GetClusterEnergy();//changed from ue to accomodate new software
//      				Double_t cells2 = lkrCand2->GetNCells();
//				Double_t ce2 = LKrZeroSup(cells2, ue2); 

				Pi0Z = sqrt(ce*Ecorr*ce2*Ecorr*dgamma*dgamma)/pi0Mass;		
				decayZ = LKrStartPos - Pi0Z;
				if(decayZ > 180000) continue; 
				if(decayZ < 105000) continue;
					if (n) return 0;
					n = true; 
					pairID.push_back(CandID[ii]);
					pairID.push_back(CandID[jj]); 		   			 
						
 								          }
								}
//calculate both photon's momentum
if(!n) return 0;
	double decayX = Tan(trim)*(decayZ-TrimStartPos);	
				double decayY = 0;
				DecayVert.SetX(decayX);
				DecayVert.SetY(decayY);
				DecayVert.SetZ(decayZ);
				for(int j=0; j<2; j++){	
				lkrCand = (TRecoLKrCandidate*)event->GetCandidate(pairID[j]);				
				Double_t ce = lkrCand->GetClusterEnergy();
				TVector3 Pos((lkrCand->GetClusterX()+Xcorr), lkrCand->GetClusterY(), LKrStartPos); 
				TVector3 Pos2 = Pos - DecayVert;
				TVector3 Pdir = Pos2.Unit();
				double Px = Pdir.X()*ce*Ecorr;
				double Py = Pdir.Y()*ce*Ecorr;
				double Pz = Pdir.Z()*ce*Ecorr;
				TLorentzVector gFourMom(Px, Py, Pz, ce*Ecorr);
				gamma4.push_back(gFourMom);
				TVector2 GammaPos((lkrCand->GetClusterX()+Xcorr), lkrCand->GetClusterY());
				gammaPos.push_back(GammaPos); 
							 }
//use paired photon 4 momentum to get Pi+ missing mass
				 Pi04Mom = gamma4[0] + gamma4[1];
                        TLorentzVector KMom(Sin(trim)*75, 0, Cos(trim)*75, sqrt(75*75+.493677*.493677));
                        	PiPlus4Mom_LKr = KMom - Pi04Mom;
				PiPlus_mmsq = PiPlus4Mom_LKr.M2(); 
				PiPlus3Mom_LKr = PiPlus4Mom_LKr.Vect();
				nPion_flag = true; 
				double smallLKr;
				int minD_id;
				Int_t LKr = 0; 
				TVector3 InitP; 
				InitP.SetX(PiPlus4Mom.X());
				InitP.SetY(PiPlus4Mom.Y());
				InitP.SetZ(PiPlus4Mom.Z()); 
				TVector3 PiPlusLKr;
				Double_t distanceLKr = 99999;
				PiPlusLKr = prop(DecayVert, InitP*1000, LKrStartPos);	
 				TVector2 PiPlusLKr2D(PiPlusLKr.X(), PiPlusLKr.Y());
				for(int i=0; i<event->GetNCandidates(); i++){
					if (i == pairID[0] || i == pairID[1]) continue;
 
				lkrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
				TVector2 ClusterPos((lkrCand->GetClusterX()+Xcorr), lkrCand->GetClusterY());
				TVector2 dist = PiPlusLKr2D - ClusterPos; 
				TVector2 sepGamma1 = PiPlusLKr2D - gammaPos[0];
					if (sepGamma1.Mod() < 200)continue;
				TVector2 sepGamma2 = PiPlusLKr2D - gammaPos[1];
					if (sepGamma2.Mod() < 200)continue;
                                Double_t ce = lkrCand->GetClusterEnergy();
//                                Double_t cells = lkrCand->GetNCells();
 //                               Double_t ce = LKrZeroSup(cells, ue);
if(ParticleFlag==1){			if(ce>1) continue; }
if(ParticleFlag==2){
					if(lkrCand->GetNCells()/ce > EnergyCellRatioMin) continue;
 					if(lkrCand->GetClusterDDeadCell() < dDeadCellMin) continue;
					if(ce < 3) continue;
}
					if (dist.Mod()<distanceLKr)	{
						distanceLKr = dist.Mod();
						minD_id = i;
 						if (dist.Mod() < 150) LKr++;
									}
										}
                                if(LKr>0)lkrCand = (TRecoLKrCandidate*)event->GetCandidate(minD_id);
				if(LKr>0)clusterE = lkrCand->GetClusterEnergy(); 
				if (LKr == 0)return 0;
				if (LKr > 0)return 1;
}

void LKrAnalysis::BookHistos_OneNeutralPion(){
                TH2I* LKr_OneNeutralPion_CellEnergy = new TH2I("LKr_OneNeutralPion_CellEnergy", "Cell to Energy Ratio (Pi+ missing LKr acceptance)", 100, 0, 100, 110, 0, 110);
                LKr_OneNeutralPion_CellEnergy->SetOption("COLZ");
                LKr_OneNeutralPion_CellEnergy->Draw();
                fUserMethods->BookHisto(LKr_OneNeutralPion_CellEnergy);

                TH2I* LKr_OneNeutralPion_CellEnergy_Post = new TH2I("LKr_OneNeutralPion_CellEnergy_Post", "Cell to Energy Ratio (Pi+ missing LKr acceptance)", 100, 0, 100, 110, 0, 110);
                LKr_OneNeutralPion_CellEnergy_Post->SetOption("COLZ");
                LKr_OneNeutralPion_CellEnergy_Post->Draw();
                fUserMethods->BookHisto(LKr_OneNeutralPion_CellEnergy_Post);
}
void LKrAnalysis::BookHistos_OneNeutralPion_Straws(){
                TH2I* LKr_OneNeutralPion_Straws_CellEnergy = new TH2I("LKr_OneNeutralPion_Straws_CellEnergy", "Cell to Energy Ratio (Pi+ missing LKr acceptance)", 100, 0, 100, 110, 0, 110);
                LKr_OneNeutralPion_Straws_CellEnergy->SetOption("COLZ");
                LKr_OneNeutralPion_Straws_CellEnergy->Draw();               
               fUserMethods->BookHisto(LKr_OneNeutralPion_Straws_CellEnergy);
}

void LKrAnalysis::BookHistos_Set(){
fUserMethods->BookHisto(new TH1F("LKr_Set_reft", "Reference Time", 400, -50, 50));
}

void LKrAnalysis::Test(TRecoLKrEvent*& LKrEvent){
for (Int_t iClus=0; iClus<LKrEvent->GetNCandidates(); iClus++) {
TRecoLKrCandidate* Lcand = (TRecoLKrCandidate*)LKrEvent->GetCandidate(iClus);
Lcand->SetClusterEnergy(15); 
}
}

Int_t LKrAnalysis::MatchedCluster(TVector3 Pos){ 
				clusterE = 0;
                                TRecoLKrCandidate *LKrCand;
                                TVector3 Pion_Prop_LKr2 = Pos;
                                double mindist = 99999;
                                double LKrStartPos = 240413;
                                int m = 0;
                                int minDCand;
				MatchedCluster_minID = -1; 
                                TVector2 Pion_Prop_LKr2D(Pion_Prop_LKr2.X(), Pion_Prop_LKr2.Y());
                                for(int i=0; i<event->GetNCandidates(); i++) {
                                LKrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
                                double Tdiff_LKr = reftime-LKrCand->GetTime();
                                TVector2 ClusterPos(LKrCand->GetClusterX(), LKrCand->GetClusterY());
                                TVector2 dist = Pion_Prop_LKr2D - ClusterPos;
                                Double_t ce = LKrCand->GetClusterEnergy();
                                        if (dist.Mod() < mindist)       {
                                                mindist = dist.Mod();
                                                minDCand = i;
                                                if (dist.Mod()<150)m++;
                                                                        }

                                                                                }
				if(mindist<99999)fUserMethods->FillHisto("LKr_MatchedCluster_MinD", mindist);
				if (m>0){
                                LKrCand = (TRecoLKrCandidate*)event->GetCandidate(minDCand);
                                Double_t ce = LKrCand->GetClusterEnergy();
                                clusterE = ce;
				fUserMethods->FillHisto("clusterE", clusterE);
				MatchedCluster_minID = minDCand; 
				MatchedCluster_RMS = sqrt(LKrCand->GetClusterRMSY()*LKrCand->GetClusterRMSY()+LKrCand->GetClusterRMSX()*LKrCand->GetClusterRMSX()); 
					}		
                                if(m==0)return 0; 
				if(m>0)return 1; 
}

Int_t LKrAnalysis::MatchedPion(TVector3 PionMomentum, TVector3 Pos){
				RMS = 0;
				Double_t P = PionMomentum.Mag();
				MatchedPion_unmatched=0; 
				bool MatchedPion=false; 
                                TRecoLKrCandidate *LKrCand;
                                TVector3 Pion_Prop_LKr2 = Pos;
                                double mindist = 99999;
                                double LKrStartPos = 240413;
				double EoverP = -1; 
                                int m = 0;
                                MatchedPion_minID = -1;
                                TVector2 Pion_Prop_LKr2D(Pion_Prop_LKr2.X(), Pion_Prop_LKr2.Y());
                                for(int i=0; i<event->GetNCandidates(); i++) {
                                LKrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
				double Tdiff_LKr = reftime-LKrCand->GetTime();
//			if(!MCflag){if(Tdiff_LKr < -15 || Tdiff_LKr > 28)continue;} //per G's recommendation
  //                      if(!MCflag){if(Tdiff_LKr > 10 && Tdiff_LKr < 22)continue;} //per G's recommendation 

//                   if(!MCflag){if(Tdiff_LKr < tLow || Tdiff_LKr > tHi)continue;}
                                TVector2 ClusterPos(LKrCand->GetClusterX(), LKrCand->GetClusterY());
                                TVector2 dist = Pion_Prop_LKr2D - ClusterPos;
                                Double_t ce = LKrCand->GetClusterEnergy();
//				if (dist.Mod() > 150) MatchedPion_unmatched++; 
					 if (dist.Mod() > 150)   {
                                                if(ce>=2&&((Tdiff_LKr>=-15&&Tdiff_LKr<=10)||(Tdiff_LKr>=22&&Tdiff_LKr<=28)))MatchedPion_unmatched++;
                                                if(ce<=2&&(Tdiff_LKr>=-5&&Tdiff_LKr<=5))MatchedPion_unmatched++;
                                                                 }
				if(!MCflag){if(Tdiff_LKr < tLow || Tdiff_LKr > tHi)continue;}
                                        if (dist.Mod() < mindist)       {
                                                mindist = dist.Mod();
                                                MatchedPion_minID = i;
                                                if (dist.Mod()<150)m++;
                                                                        }

                                                                                }
//                                if(mindist<99999)fUserMethods->FillHisto("LKr_MatchedPion_MinD", mindist);
                                if (m>0){
//                              	if(m==1){
                                LKrCand = (TRecoLKrCandidate*)event->GetCandidate(MatchedPion_minID);
                                Double_t ce = LKrCand->GetClusterEnergy();
                                clusterE = ce;
				EoverP = ce/P;
				RMS = sqrt(LKrCand->GetClusterRMSY()*LKrCand->GetClusterRMSY()+LKrCand->GetClusterRMSX()*LKrCand->GetClusterRMSX());
				if (EoverP<0.8)MatchedPion=true; 
  //                          fUserMethods->FillHisto("clusterE", clusterE);
                                        }
                                if(MatchedPion)return 1;
                                else return 0;
}
Int_t LKrAnalysis::MatchedMuon(TVector3 Pos){
				RMS = 0; 
                                MatchedMuon_unmatched=0;
                                bool MatchedMuon=false;
                                TRecoLKrCandidate *LKrCand;
                                TVector3 Pion_Prop_LKr2 = Pos;
                                double mindist = 99999;
                                double LKrStartPos = 240413;
                                int m = 0;
				MatchedMuon_minID = -1;
                                int minDCand;
                                TVector2 Pion_Prop_LKr2D(Pion_Prop_LKr2.X(), Pion_Prop_LKr2.Y());
                                for(int i=0; i<event->GetNCandidates(); i++) {
                                LKrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
				double Tdiff_LKr = reftime-LKrCand->GetTime(); 
//                     if(!MCflag){if(Tdiff_LKr < tLow || Tdiff_LKr > tHi)continue;}
                                TVector2 ClusterPos(LKrCand->GetClusterX(), LKrCand->GetClusterY());
                                TVector2 dist = Pion_Prop_LKr2D - ClusterPos;
                                Double_t ce = LKrCand->GetClusterEnergy();
//					if (dist.Mod() > 150)MatchedMuon_unmatched++; 
						if (dist.Mod() > 150)   {
                                                if(ce>=2&&((Tdiff_LKr>=-15&&Tdiff_LKr<=10)||(Tdiff_LKr>=22&&Tdiff_LKr<=28)))MatchedMuon_unmatched++;
                                                if(ce<=2&&(Tdiff_LKr>=-5&&Tdiff_LKr<=5))MatchedMuon_unmatched++;
                                                                 }
		if(!MCflag){if(Tdiff_LKr < tLow || Tdiff_LKr > tHi)continue;}
                                        if (dist.Mod() < mindist)       {
                                                mindist = dist.Mod();
                                                minDCand = i;
                                                if (dist.Mod()<150)m++;
                                                                        }

                                                                                }
                                if (m==1){
                                LKrCand = (TRecoLKrCandidate*)event->GetCandidate(minDCand);
                                Double_t ce = LKrCand->GetClusterEnergy();
                                clusterE = ce;
				RMS = sqrt(LKrCand->GetClusterRMSY()*LKrCand->GetClusterRMSY()+LKrCand->GetClusterRMSX()*LKrCand->GetClusterRMSX()); 
                              //  if (ce<3) put back in when doing non calo rej stuf
				MatchedMuon=true;
				MatchedMuon_minID = minDCand;
                                        }
                                if(MatchedMuon)return 1;
                                else return 0;
}


void LKrAnalysis::BookHistos_MatchedCluster(){
fUserMethods->BookHisto(new TH1F("LKr_MatchedCluster_MinD", "Minimum Cluster Distance", 500, 0, 500));
fUserMethods->BookHisto(new TH1I("clusterE", "cluster e", 100, 0, 100));
}


Int_t LKrAnalysis::TwoNeutralPions(TLorentzVector Pi4, TVector3 fDecayVertexPi){
	TVector3 DecayVertexPi = fDecayVertexPi; 
          TLorentzVector PiPlus4Mom = Pi4; 
	 nPion_flag = false;
	 TRecoLKrCandidate *lkrCand;
	 TRecoLKrCandidate *lkrCand2;
	double EnergyRatioMin = .28;
	double EnergyCellRatioMin = 1.66; 
	double dGammaMin = 20;
	double dDeadCellMin = 2; 
	double LKrOutRadius = 1130;
	double LKrInRadius = 150; 
        double LKrStartPos = 240413;
	double pi0Mass = .1349766;
	double CHODStartPos = 239009;
	double TrimStartPos = 101800;
	double trim = .0012; 
//Horizaont LKr position correction
double Xcorr = 0; //add
	double Ecorr = 1; //multiply
        Double_t tCedarhi =0.5;//
        Double_t tCedarlow =-0.7;//
        Double_t tCedarmid =-0.13; //
        double LAV_tLow = -3;
        double LAV_tHi = 4;
        double MUV3_tLow = -2;
        double MUV3_tHi =1;
        Double_t tCHANTILow = -6;
        Double_t tCHANTIHi = 2;
        Double_t tMUV1Low = -14;
        Double_t tMUV1Hi = 6;
        Double_t tMUV2Low = -5;
        Double_t tMUV2Hi = 5;
        Double_t tRICHLow = -1;
        Double_t tRICHHi = 0.5;
vector<int> CandID;
	vector<int> CandID2;

	for(int i=0; i<event->GetNCandidates(); i++){
				lkrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
				Double_t ce = lkrCand->GetClusterEnergy();
//      				Double_t cells = lkrCand->GetNCells();
//				Double_t ce = LKrZeroSup(cells, ue); 
//detector acceptance
if(sqrt((lkrCand->GetClusterX()+Xcorr)*(lkrCand->GetClusterX()+Xcorr)+lkrCand->GetClusterY()*lkrCand->GetClusterY())>LKrOutRadius || sqrt((lkrCand->GetClusterX()+Xcorr)*(lkrCand->GetClusterX()+Xcorr)+lkrCand->GetClusterY()*lkrCand->GetClusterY())<LKrInRadius)continue;				
					fUserMethods->FillHisto("LKr_TwoNeutralPions_CellEnergy", ce*Ecorr, lkrCand->GetNCells());
					if(lkrCand->GetNCells()/ce > EnergyCellRatioMin) continue;
					 fUserMethods->FillHisto("LKr_TwoNeutralPions_CellEnergy_Post", ce*Ecorr, lkrCand->GetNCells());
 					if(lkrCand->GetClusterDDeadCell() < dDeadCellMin) continue;
					if(ce < 3) continue;
					CandID.push_back(i);
						      	} //LKr for loop
	vector<TLorentzVector> gamma4;
vector<TVector3> DecayVertex;
	vector<TVector2> gammaPos;
	double Pi0Z;
	double decayZ;
vector<pair<int,int>> PhotonPairs;
vector<Double_t> Zcoord;
//find decay vertex
	vector<int> pairID;
//	TVector3 DecayVert; 
	fUserMethods->FillHisto("CandID", CandID.size());  
	if (CandID.size()<4)return 0;  
cout<<"GotHere"<<endl;
fUserMethods->FillHisto("n", 1);  
		for (int ii = 0; ii < CandID.size() - 1; ii++){
		lkrCand = (TRecoLKrCandidate*)event->GetCandidate(CandID[ii]);
	//	if ((lkrCand->GetClusterX()+Xcorr)>-43.5&&(lkrCand->GetClusterX()+Xcorr)<-39.5&&lkrCand->GetClusterY()>-39.6&&lkrCand->GetClusterY()<-37.5) continue;
		 	for (int jj = ii + 1; jj < CandID.size(); jj++){
				lkrCand2 = (TRecoLKrCandidate*)event->GetCandidate(CandID[jj]);
//				if ((lkrCand2->GetClusterX()+1)>-43.5&&(lkrCand2->GetClusterX()+1)<-39.5&&lkrCand2->GetClusterY()>-39.6&&lkrCand2->GetClusterY()<-37.5) continue;			
				double dgamma = sqrt(((lkrCand->GetClusterX()+Xcorr)-(lkrCand2->GetClusterX()))*((lkrCand->GetClusterX()+Xcorr)-(lkrCand2->GetClusterX())) + (lkrCand->GetClusterY()-lkrCand2->GetClusterY())*(lkrCand->GetClusterY()-lkrCand2->GetClusterY()));
				if (dgamma < dGammaMin) continue; 
	
				Double_t ce = lkrCand->GetClusterEnergy();//changed from ue to accomodate new software
 //     				Double_t cells = lkrCand->GetNCells();
//				Double_t ce = LKrZeroSup(cells, ue); 


				Double_t ce2 = lkrCand2->GetClusterEnergy();//changed from ue to accomodate new software
//      				Double_t cells2 = lkrCand2->GetNCells();
//				Double_t ce2 = LKrZeroSup(cells2, ue2); 

				Pi0Z = sqrt(ce*Ecorr*ce2*Ecorr*dgamma*dgamma)/pi0Mass;		
				decayZ = LKrStartPos - Pi0Z;
				if(decayZ > 180000) continue; 
				if(decayZ < 105000) continue;
 				double decayX = Tan(trim)*(decayZ-TrimStartPos);
                                double decayY = 0;
                                DecayVert.SetX(decayX);
                                DecayVert.SetY(decayY);
                                DecayVert.SetZ(decayZ);
//				TVector3 difference = DecayVertexPi - DecayVert;
  //      			if (fabs(difference.Mag())>30)continue;
				PhotonPairs.push_back(std::make_pair(CandID[ii], CandID[jj]));
//		Zcoord.push_back(decayZ); 
				DecayVertex.push_back(DecayVert);
//					pairID.push_back(CandID[ii]);
//					pairID.push_back(CandID[jj]); 		   			 
						
 								          }
								}
fUserMethods->FillHisto("n", 2);
 
vector<pair<int,int>> BestPhotonPairs;
BestPhotonPairs.push_back(std::make_pair(-1, -1));
BestPhotonPairs.push_back(std::make_pair(-1, -1));
Int_t match = 0; 
Double_t dist = 999999;
if (PhotonPairs.size()<2)return 0;
DecayVertexAve.SetXYZ(0,0,0); 
for (int i = 0; i<DecayVertex.size(); i++){
	for (int j = i+1; j<DecayVertex.size(); j++){
	if (PhotonPairs[i].first==PhotonPairs[j].first)continue;
        if (PhotonPairs[i].second==PhotonPairs[j].second)continue;
        if (PhotonPairs[i].first==PhotonPairs[j].second)continue;
	if (PhotonPairs[i].second==PhotonPairs[j].first)continue; 
        if (PhotonPairs[j].first==PhotonPairs[i].first)continue;
        if (PhotonPairs[j].second==PhotonPairs[i].second)continue;
        if (PhotonPairs[j].first==PhotonPairs[i].second)continue;
        if (PhotonPairs[j].second==PhotonPairs[i].first)continue;

	TVector3 diff =DecayVertex[i] - DecayVertex[j];
	fUserMethods->FillHisto("x", diff.X());
	fUserMethods->FillHisto("z", diff.Z());	 
	if (diff.Mag()<dist) {
		dist = diff.Mag(); 
		BestPhotonPairs[0].first=PhotonPairs[i].first, 
		BestPhotonPairs[0].second=PhotonPairs[i].second;
                BestPhotonPairs[1].first=PhotonPairs[j].first,
                BestPhotonPairs[1].second=PhotonPairs[j].second;
		DecayVertexAve = 0.5*(DecayVertex[i] + DecayVertex[j]); 
		match++;
}
}		
}
if (match==0)return 0; 
fUserMethods->FillHisto("dist", dist); 
if (dist>8000)return 0;
if (PhotonPairs[0].first==-1)return 0;
        if (PhotonPairs[0].second==-1)return 0;
        if (PhotonPairs[1].first==-1)return 0;
        if (PhotonPairs[1].second==-1)return 0;

fUserMethods->FillHisto("Pairs", BestPhotonPairs[0].first); 
fUserMethods->FillHisto("Pairs", BestPhotonPairs[0].second);
fUserMethods->FillHisto("Pairs", BestPhotonPairs[1].first);
fUserMethods->FillHisto("Pairs",BestPhotonPairs[1].second); 
/*
vector<pair<int,int>> BestPhotonPairs;
Double_t zz = 99999;
Int_t minzID;
Int_t minz = 0;
for (int i = 0; i<Zcoord.size(); i++){
if (Zcoord[i]<zz){
zz == Zcoord[i];
minzID = i;
minz++;
}
}
if (minz ==0)return 0;
BestPhotonPairs.push_back(std::make_pair(PhotonPairs[minzID].first, PhotonPairs[minzID].second));

Double_t zz_2 = 99999;
Int_t minzID_2;
Int_t minz_2 = 0;
for (int i = 0; i<Zcoord.size(); i++){
if(i=minzID)continue; 
if (Zcoord[i]<zz_2){
zz_2 == Zcoord[i];
minzID_2 = i;
minz_2++;
}
}
if (minz_2==0)return 0; 
*/
//BestPhotonPairs.push_back(std::make_pair(PhotonPairs[minzID_2].first, PhotonPairs[minzID_2].second));
//this loop may best be moved above previous two loops, not sure
/*
Int_t repeat = 0;
fUserMethods->FillHisto("Size", BestPhotonPairs.size()); 
for (int k=0;k<BestPhotonPairs.size(); k++){
    for(int j=k+1; j <BestPhotonPairs.size(); j++){
	if(BestPhotonPairs[k].first ==BestPhotonPairs[j].second)repeat++;
}
}
fUserMethods->FillHisto("repeat", repeat); 
if (repeat > 0)return 0; 
*/
fUserMethods->FillHisto("n", 3);
vector<pair<TLorentzVector,TLorentzVector>> PhotonMomentum;
vector<pair<TVector2,TVector2>> PhotonPosition;

                                for(int j=0; j<2; j++){
                                lkrCand = (TRecoLKrCandidate*)event->GetCandidate(BestPhotonPairs[j].first);
                                Double_t ce = lkrCand->GetClusterEnergy();
                                TVector3 Pos((lkrCand->GetClusterX()+Xcorr), lkrCand->GetClusterY(), LKrStartPos);
                                TVector3 Pos2 = Pos - DecayVert;
                                TVector3 Pdir = Pos2.Unit();
                                double Px = Pdir.X()*ce*Ecorr;
                                double Py = Pdir.Y()*ce*Ecorr;
                                double Pz = Pdir.Z()*ce*Ecorr;
                                TLorentzVector gFourMom(Px, Py, Pz, ce*Ecorr);
//                                gamma4.push_back(gFourMom);
                                TVector2 GammaPos((lkrCand->GetClusterX()+Xcorr), lkrCand->GetClusterY());
//                                gammaPos.push_back(GammaPos);

                                lkrCand = (TRecoLKrCandidate*)event->GetCandidate(BestPhotonPairs[j].second);
                                ce = lkrCand->GetClusterEnergy();
                                TVector3 Pos_2((lkrCand->GetClusterX()+Xcorr), lkrCand->GetClusterY(), LKrStartPos);
                                TVector3 Pos2_2 = Pos - DecayVert;
                                TVector3 Pdir_2 = Pos2.Unit();
                                Px = Pdir_2.X()*ce*Ecorr;
                                Py = Pdir_2.Y()*ce*Ecorr;
                                Pz = Pdir_2.Z()*ce*Ecorr;
                                TLorentzVector gFourMom2(Px, Py, Pz, ce*Ecorr);
//                                gamma4.push_back(gFourMom);
                                TVector2 GammaPos2((lkrCand->GetClusterX()+Xcorr), lkrCand->GetClusterY());
//                                gammaPos.push_back(GammaPos);
				PhotonMomentum.push_back(std::make_pair(gFourMom, gFourMom2));
				PhotonPosition.push_back(std::make_pair(GammaPos, GammaPos2));
                                                         }
//good
				 Pi04Mom = PhotonMomentum[0].first + PhotonMomentum[0].second;
                                 TLorentzVector Pi04Mom2 = PhotonMomentum[1].first + PhotonMomentum[1].second;
                        TLorentzVector KMom(Sin(trim)*75, 0, Cos(trim)*75, sqrt(75*75+.493677*.493677));
                        	PiPlus4Mom_LKr = KMom - Pi04Mom - Pi04Mom2;
				PiPlus_mmsq = PiPlus4Mom_LKr.M2(); 
				PiPlus3Mom_LKr = PiPlus4Mom_LKr.Vect();
				nPion_flag = true; 
				double smallLKr;
				int minD_id;
				Int_t LKr = 0; 
				TVector3 InitP; 
				InitP.SetX(PiPlus4Mom.X());
				InitP.SetY(PiPlus4Mom.Y());
				InitP.SetZ(PiPlus4Mom.Z()); 
				TVector3 PiPlusLKr;
				Double_t distanceLKr = 99999;
				PiPlusLKr = prop(DecayVert, InitP*1000, LKrStartPos);	
 				TVector2 PiPlusLKr2D(PiPlusLKr.X(), PiPlusLKr.Y());
//good
return 0; 
//seg fault after this, need to find 
				for(int i=0; i<event->GetNCandidates(); i++){
					if (i == BestPhotonPairs[0].first || i == BestPhotonPairs[0].second ||i == BestPhotonPairs[1].first ||i == BestPhotonPairs[1].second ) continue;
 
				lkrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
				TVector2 ClusterPos((lkrCand->GetClusterX()+Xcorr), lkrCand->GetClusterY());
				TVector2 dist = PiPlusLKr2D - ClusterPos; 
				TVector2 sepGamma1 = PiPlusLKr2D - PhotonPosition[0].first;
					if (sepGamma1.Mod() < 200)continue;
				TVector2 sepGamma2 = PiPlusLKr2D - PhotonPosition[0].second;
					if (sepGamma2.Mod() < 200)continue;
                                TVector2 sepGamma3 = PiPlusLKr2D - PhotonPosition[1].first;
                                        if (sepGamma3.Mod() < 200)continue;
                                TVector2 sepGamma4 = PiPlusLKr2D - PhotonPosition[1].second;
                                        if (sepGamma4.Mod() < 200)continue;
 
                                Double_t ce = lkrCand->GetClusterEnergy();
					if (dist.Mod()<distanceLKr)	{
						distanceLKr = dist.Mod();
						minD_id = i;
 						if (dist.Mod() < 150) LKr++;
									}																			}
				if (LKr==0)return 0;
				if (LKr>0)return 1; 
}

Int_t LKrAnalysis::TrackCellMatching(TVector3 *posatlkr) {
  TRecoLKrEvent* fLKrEvent = event;
  Double_t trackTime = reftime; 
  fLKrCellCandidate = new TRecoLKrCandidate;  
  Double_t xpart = posatlkr->X();
  Double_t ypart = posatlkr->Y();
//Int_t minid = -1;

 // Look for cells associated to the track, if no cluster candidate matches the track
  Double_t eclus = 0.;
  Double_t tclus = 0.;
  Double_t xclus = 0.;
  Double_t yclus = 0.;
  Int_t ncells = 0;
  Double_t emax = -999999.;
  TClonesArray& Hits = (*(fLKrEvent->GetHits()));
  for (Int_t jHit=0; jHit<fLKrEvent->GetNHits(); jHit++) {
    TRecoLKrHit *hit = (TRecoLKrHit*)Hits[jHit];
    Int_t ixcell = (Int_t)hit->GetXCellID();
    Int_t iycell = (Int_t)hit->GetYCellID();
    Double_t xcell = hit->GetPosition().X();
    Double_t ycell = hit->GetPosition().Y();
    Double_t ecell = hit->GetEnergy()/1000.;
    if (ecell<=-999999.) continue;
//  Double_t t0cell = (ixcell>=0 && iycell>=0) ? fLKrCellT0[ixcell][iycell] : 0;
    Double_t tcell = hit->GetTime()+3.2;
    Double_t dist = sqrt((xcell-xpart)*(xcell-xpart)+(ycell-ypart)*(ycell-ypart));
    Double_t dtime = tcell-trackTime;
    if (dist<150. && fabs(dtime)<40) { // Select a region around the extrapolated track position at LKr
      eclus += ecell;
      if (emax<ecell) {
        emax = ecell;
        tclus = tcell+1.3;
      }
      xclus += xcell*ecell;
      yclus += ycell*ecell;
      ncells++;
//      minid = 1;
    }
  }
  if (!ncells || emax<0.04) { // no cells found
    tclus = -99999.;
    eclus = 0;
    xclus = -99999.;
    yclus = -99999.;
 //   minid = -1;
//    return 0;
  }
  
  xclus /= eclus;
  yclus /= eclus;
  Double_t ue = eclus;
  if (ncells>9) {
    if (ue<22) eclus = ue/(0.7666+0.0573489*log(ue));
    if (ue>=22 && ue<65) eclus = ue/(0.828962+0.0369797*log(ue));
    if (ue>=65) eclus = ue/(0.828962+0.0369797*log(65));
  } 
  eclus *= 1.03;
  eclus = eclus<15 ? (eclus+0.015)/(15+0.015)*15*0.9999 : eclus; 
  Double_t mindist = sqrt((xclus-xpart)*(xclus-xpart)+(yclus-ypart)*(yclus-ypart));
//  fLKrCellCandidate->SetDiscriminant(mindist);
//  fLKrCellCandidate->SetDeltaTime(tclus-trackTime);
  fLKrCellCandidate->SetClusterEnergy(eclus);
//  fLKrCellCandidate->SetSeedEnergy(emax);
fLKrCellCandidate->SetTime(tclus);
  fLKrCellCandidate->SetClusterX(xclus);
  fLKrCellCandidate->SetClusterY(yclus);
  fLKrCellCandidate->SetNCells(ncells);
/*
  if (fFlag==0) {
    fUserMethods->FillHisto("LKrAnalysis_mindistvst_cell",tclus-trackTime,mindist);
    fUserMethods->FillHisto("LKrAnalysis_mindistvstchod_cell",tclus-chodTime,mindist);
    fUserMethods->FillHisto("LKrAnalysis_minxy_cell",xclus-xpart,yclus-ypart);
    fUserMethods->FillHisto("LKrAnalysis_energyvsdist_cell",mindist,eclus);
  }
  if (fFlag==2) {
    fUserMethods->FillHisto("LKrAnalysis_2gammamindistvst_cell",tclus-trackTime,mindist);
    fUserMethods->FillHisto("LKrAnalysis_2gammamindistvstchod_cell",tclus-chodTime,mindist);
    fUserMethods->FillHisto("LKrAnalysis_2gammaminxy_cell",xclus-xpart,yclus-ypart);
    fUserMethods->FillHisto("LKrAnalysis_2gammaenergyvsdist_cell",mindist,eclus);
  }
*/
  if (!ncells || emax<0.04)return 0;
  else{ 
	clusterE = fLKrCellCandidate->GetClusterEnergy(); 	
	return 1;
 }
}

void LKrAnalysis::BookHistos_TwoNeutralPions(){
                TH2I* LKr_TwoNeutralPions_CellEnergy = new TH2I("LKr_TwoNeutralPions_CellEnergy", "Cell to Energy Ratio (Pi+ missing LKr acceptance)", 100, 0, 100, 110, 0, 110);
                LKr_TwoNeutralPions_CellEnergy->SetOption("COLZ");
                LKr_TwoNeutralPions_CellEnergy->Draw();
                fUserMethods->BookHisto(LKr_TwoNeutralPions_CellEnergy);

                TH2I* LKr_TwoNeutralPions_CellEnergy_Post = new TH2I("LKr_TwoNeutralPions_CellEnergy_Post", "Cell to Energy Ratio (Pi+ missing LKr acceptance)", 100, 0, 100, 110, 0, 110);
                LKr_TwoNeutralPions_CellEnergy_Post->SetOption("COLZ");
                LKr_TwoNeutralPions_CellEnergy_Post->Draw();
                fUserMethods->BookHisto(LKr_TwoNeutralPions_CellEnergy_Post);
fUserMethods->BookHisto(new TH1F("n", "hmm", 10, 0, 10));
fUserMethods->BookHisto(new TH1F("dist", "hmm", 180, 0, 180000));
fUserMethods->BookHisto(new TH1F("CandID", "# of photons", 10, 0, 10));
fUserMethods->BookHisto(new TH1F("x", "hmm", 360, -180000, 180000));
fUserMethods->BookHisto(new TH1F("z", "hmm", 360, -180000, 180000));
fUserMethods->BookHisto(new TH1F("repeat", "hmm", 50, 0, 50));
fUserMethods->BookHisto(new TH1F("Size", "hmm", 10, 0, 10));
fUserMethods->BookHisto(new TH1F("Pairs", "hmm", 20, -10, 10));

}

Int_t LKrAnalysis::OneNeutralPion_Straws(TVector3 Pi3, TVector3 PiPos, TLorentzVector K4, Int_t pflag){
	Int_t ParticleFlag = pflag; //0 general, 1 muons, 2 electrons
	OneNeutralPion_Straws_minID=-1;
	OneNeutralPion_Straws_RMS = 0;
	TLorentzVector KMom = K4;
//          TLorentzVector PiPlus4Mom = Pi4; 
	 TRecoLKrCandidate *lkrCand;
	 TRecoLKrCandidate *lkrCand2;
	double EnergyRatioMin = .28;
	double EnergyCellRatioMin = 1.66; 
	double dGammaMin = 20;
	double dDeadCellMin = 20; 
	double LKrOutRadius = 1130;
	double LKrInRadius = 150; 
        double LKrStartPos = 240413;
	double pi0Mass = .1349766;
	double CHODStartPos = 239009;
	double TrimStartPos = 101800;
	double trim = .0012; 
//Horizaont LKr position correction
double Xcorr = 0; //add
	double Ecorr = 1; //multiply
        Double_t tCedarhi =0.5;//
        Double_t tCedarlow =-0.7;//
        Double_t tCedarmid =-0.13; //
        double LAV_tLow = -3;
        double LAV_tHi = 4;
        double MUV3_tLow = -2;
        double MUV3_tHi =1;
        Double_t tCHANTILow = -6;
        Double_t tCHANTIHi = 2;
        Double_t tMUV1Low = -14;
        Double_t tMUV1Hi = 6;
        Double_t tMUV2Low = -5;
        Double_t tMUV2Hi = 5;
        Double_t tRICHLow = -1;
        Double_t tRICHHi = 0.5;
vector<int> CandID;
	vector<int> CandID2;

	for(int i=0; i<event->GetNCandidates(); i++){
				lkrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
				Double_t ce = lkrCand->GetClusterEnergy();
//      				Double_t cells = lkrCand->GetNCells();
//				Double_t ce = LKrZeroSup(cells, ue); 
//detector acceptance
if(sqrt((lkrCand->GetClusterX()+Xcorr)*(lkrCand->GetClusterX()+Xcorr)+lkrCand->GetClusterY()*lkrCand->GetClusterY())>LKrOutRadius || sqrt((lkrCand->GetClusterX()+Xcorr)*(lkrCand->GetClusterX()+Xcorr)+lkrCand->GetClusterY()*lkrCand->GetClusterY())<LKrInRadius)continue;				
					fUserMethods->FillHisto("LKr_OneNeutralPion_Straws_CellEnergy", ce*Ecorr, lkrCand->GetNCells());
					if(lkrCand->GetNCells()/ce > EnergyCellRatioMin) continue;
//					 fUserMethods->FillHisto("LKr_OneNeutralPion_CellEnergy_Post", ce*Ecorr, lkrCand->GetNCells());
 					if(lkrCand->GetClusterDDeadCell() < dDeadCellMin) continue;
					if(ce < 3) continue;
					CandID.push_back(i);
						      	} //LKr for loop
	vector<TLorentzVector> gamma4;
	vector<TVector2> gammaPos;
	double Pi0Z;
	double decayZ;
//find decay vertex
bool n = false;
	vector<int> pairID;
//	TVector3 DecayVert;  
	if (CandID.size()<2)return 0; 
		for (int ii = 0; ii < CandID.size() - 1; ii++){
		lkrCand = (TRecoLKrCandidate*)event->GetCandidate(CandID[ii]);
	//	if ((lkrCand->GetClusterX()+Xcorr)>-43.5&&(lkrCand->GetClusterX()+Xcorr)<-39.5&&lkrCand->GetClusterY()>-39.6&&lkrCand->GetClusterY()<-37.5) continue;
		 	for (int jj = ii + 1; jj < CandID.size(); jj++){
				lkrCand2 = (TRecoLKrCandidate*)event->GetCandidate(CandID[jj]);
//				if ((lkrCand2->GetClusterX()+1)>-43.5&&(lkrCand2->GetClusterX()+1)<-39.5&&lkrCand2->GetClusterY()>-39.6&&lkrCand2->GetClusterY()<-37.5) continue;			
				double dgamma = sqrt(((lkrCand->GetClusterX()+Xcorr)-(lkrCand2->GetClusterX()))*((lkrCand->GetClusterX()+Xcorr)-(lkrCand2->GetClusterX())) + (lkrCand->GetClusterY()-lkrCand2->GetClusterY())*(lkrCand->GetClusterY()-lkrCand2->GetClusterY()));
				if (dgamma < dGammaMin) continue; 
	
				Double_t ce = lkrCand->GetClusterEnergy();//changed from ue to accomodate new software
 //     				Double_t cells = lkrCand->GetNCells();
//				Double_t ce = LKrZeroSup(cells, ue); 


				Double_t ce2 = lkrCand2->GetClusterEnergy();//changed from ue to accomodate new software
//      				Double_t cells2 = lkrCand2->GetNCells();
//				Double_t ce2 = LKrZeroSup(cells2, ue2); 

				Pi0Z = sqrt(ce*Ecorr*ce2*Ecorr*dgamma*dgamma)/pi0Mass;		
				decayZ = LKrStartPos - Pi0Z;
				if(decayZ > 180000) continue; 
				if(decayZ < 105000) continue;
					if (n) return 0;
					n = true; 
					pairID.push_back(CandID[ii]);
					pairID.push_back(CandID[jj]); 		   			 
						
 								          }
								}
//calculate both photon's momentum
if(!n) return 0;
	double decayX = Tan(trim)*(decayZ-TrimStartPos);	
				double decayY = 0;
				DecayVert.SetX(decayX);
				DecayVert.SetY(decayY);
				DecayVert.SetZ(decayZ);
				for(int j=0; j<2; j++){	
				lkrCand = (TRecoLKrCandidate*)event->GetCandidate(pairID[j]);				
				Double_t ce = lkrCand->GetClusterEnergy();
				TVector3 Pos((lkrCand->GetClusterX()+Xcorr), lkrCand->GetClusterY(), LKrStartPos); 
				TVector3 Pos2 = Pos - DecayVert;
				TVector3 Pdir = Pos2.Unit();
				double Px = Pdir.X()*ce*Ecorr;
				double Py = Pdir.Y()*ce*Ecorr;
				double Pz = Pdir.Z()*ce*Ecorr;
				TLorentzVector gFourMom(Px, Py, Pz, ce*Ecorr);
				gamma4.push_back(gFourMom);
				TVector2 GammaPos((lkrCand->GetClusterX()+Xcorr), lkrCand->GetClusterY());
				gammaPos.push_back(GammaPos); 
							 }
//use paired photon 4 momentum to get Pi+ missing mass
				 Pi04Mom = gamma4[0] + gamma4[1];
//                        TLorentzVector KMom(Sin(trim)*75, 0, Cos(trim)*75, sqrt(75*75+.493677*.493677));
                        	PiPlus4Mom_LKr = KMom - Pi04Mom;
				PiPlus_mmsq = PiPlus4Mom_LKr.M2(); 
				PiPlus3Mom_LKr = PiPlus4Mom_LKr.Vect();
				nPion_flag = true; 
				double smallLKr;
//		OneNeutralPionStraws_minID=-1
				Int_t LKr = 0; 
				TVector3 InitP = Pi3;  
				TVector3 PionPosition = PiPos;
				TVector3 PiPlusLKr;
				Double_t distanceLKr = 150;
				PiPlusLKr = prop(PionPosition, InitP*1000, LKrStartPos);	
 				TVector2 PiPlusLKr2D(PiPlusLKr.X(), PiPlusLKr.Y());
				for(int i=0; i<event->GetNCandidates(); i++){
					if (i == pairID[0] || i == pairID[1]) continue;
 
				lkrCand = (TRecoLKrCandidate*)event->GetCandidate(i);
				TVector2 ClusterPos((lkrCand->GetClusterX()+Xcorr), lkrCand->GetClusterY());
				TVector2 dist = PiPlusLKr2D - ClusterPos; 
				TVector2 sepGamma1 = PiPlusLKr2D - gammaPos[0];
					if (sepGamma1.Mod() < 200)continue;
				TVector2 sepGamma2 = PiPlusLKr2D - gammaPos[1];
					if (sepGamma2.Mod() < 200)continue;
                                Double_t ce = lkrCand->GetClusterEnergy();
if(ParticleFlag==1){                    if(ce>1) continue; }
if(ParticleFlag==2){
                                        if(lkrCand->GetNCells()/ce > EnergyCellRatioMin) continue;
                                        if(lkrCand->GetClusterDDeadCell() < dDeadCellMin) continue;
                                        if(ce < 3) continue;
}

					if (dist.Mod()<distanceLKr)	{
						distanceLKr = dist.Mod();
						OneNeutralPion_Straws_minID = i;
						LKr++;
//					if (dist.Mod() < 150) LKr++;
									}
										}
                                if(LKr>0){
				lkrCand = (TRecoLKrCandidate*)event->GetCandidate(OneNeutralPion_Straws_minID);
				OneNeutralPion_Straws_RMS = sqrt(lkrCand->GetClusterRMSY()*lkrCand->GetClusterRMSY()+lkrCand->GetClusterRMSX()*lkrCand->GetClusterRMSX());
				clusterE = lkrCand->GetClusterEnergy(); 
				return 1;
 }
				else return 0;
}
