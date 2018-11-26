#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "RunCalculations.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "PreAnalyzers.hh"
#include "CalorimeterAnalysis.hh"
#include "PhotonAnalysis.hh"
#include "TrackAnalysis.hh"
#include "StrawAnalysis.hh"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;
RunCalculations::RunCalculations(Core::BaseAnalysis *ba) : Analyzer(ba, "RunCalculations")
{
	/// \MemberDescr
/// \EndMemberDescr
RequestTree("Spectrometer",new TRecoSpectrometerEvent);
        RequestTree("CHOD",new TRecoCHODEvent);
        RequestTree("Cedar",new TRecoCedarEvent);
        RequestTree("RICH",new TRecoRICHEvent);
        RequestTree("LKr",new TRecoLKrEvent);
         RequestTree("CHANTI",new TRecoCHANTIEvent);
        RequestTree("LAV",new TRecoLAVEvent);
         RequestTree("MUV1",new TRecoMUV1Event);
         RequestTree("MUV2",new TRecoMUV2Event);
        RequestTree("MUV3",new TRecoMUV3Event);
        RequestTree("GigaTracker",new TRecoGigaTrackerEvent);
        RequestTree("MUV0",new TRecoMUV0Event);
        RequestTree("SAC",new TRecoSACEvent);
        RequestTree("IRC",new TRecoIRCEvent);
 RequestTree("SAV",new TRecoSAVEvent);
RequestL0Data();
Tracka = new TrackAnalysis(ba);  
Strawa = new StrawAnalysis(ba); 
}

void RunCalculations::InitOutput(){

}

void RunCalculations::InitHist(){
	BookHisto(new TH1F("mismatch", "Mismatched bins in MC", 11, 0, 11));
	BookHisto(new TH1F("mmsq_SREvents", "Events in Signal Region", 250, -.35, .15));
	BookHisto(new TH1F("NEvents", "Events per burst", 30, 0, 30000));
 BookHisto(new TH1F("EventNumber", "Events per burst", 200, 0, 200));
	BookHisto(new TH1F("init_p", "Initial momentum in MC", 100, 0, 100));
 BookHisto(new TH1I("Retarded", "momentum", 100, 0, 100));
for(int i=0; i<4; i++){
  acc_mmsq[i].Form("acc_mmsq_%d", i);
  BookHisto(new TH1I(acc_mmsq[i].Data(), "mmsq of accepted events", 250, -.35, .15));
}
for(int i=0; i<4; i++){
 mc_pi[i].Form("mc_pi_%d", i);
  BookHisto(new TH1I(mc_pi[i].Data(), "momentum with no cuts", 100, 0, 100));

  mc_p[i].Form("mc_p_%d", i);
  BookHisto(new TH1I(mc_p[i].Data(), "initial momentum from mc", 100, 0, 100));
 mc_ps[i].Form("mc_ps_%d", i);
  BookHisto(new TH1I(mc_ps[i].Data(), "momentum from straws", 100, 0, 100));
 mc_chambers[i].Form("mc_chambers_%d", i);
  BookHisto(new TH1I(mc_chambers[i].Data(), "fourth chamber",5, 0, 5));
mc_chambersp[i].Form("mc_chambersp_%d", i);
  BookHisto(new TH1I(mc_chambersp[i].Data(), "fourth chamber",5, 0, 5));


}
BookHisto(new TH1F("events", "Event counter for MC events and GTK type", 125, 0, 125));
TH2F* pos = new TH2F("pos", "position at ch3", 300, -1500, 1500, 300, -1500, 1500);
pos->SetOption("COLZ");
pos->Draw();
BookHisto(pos);
	Tracka->One->BookHistos();
	Tracka->Calo->BookHistos();  
	Tracka->Photon->BookHistos(); 
	Tracka->BookHistos(); 

}

void RunCalculations::DefineMCSimple(){
 int kID = fMCSimple.AddParticle(0, 321); //ask for beam Kaon

        fMCSimple.AddParticle(kID, 211); //ask for positive pion from initial kaon decay
 //       int pi0ID = fMCSimple.AddParticle(kID, 111); //ask for positive pion from initial kaon decay
  //      fMCSimple.AddParticle(pi0ID, 22); //ask for positive pion from initial kaon decay
  //      fMCSimple.AddParticle(pi0ID, 22); //ask for positive pion from initial kaon decay

//	 fMCSimple.AddParticle(kID, -13); //ask for positive pion from initial kaon decay


}

void RunCalculations::StartOfRunUser(){
TypeCounter=0;
burst_number=0;
	/// \MemberDescr
	/// This method is called at the beginning of the processing (corresponding to a start of run in the normal NA62 data taking)\n
	/// Do here your start of run processing if any
	/// \EndMemberDescr
}

void RunCalculations::StartOfBurstUser(){

Tracka->Photon->LKra_G->StartBurst(2015,/*GetRawHeader()->GetRunID()*/3809);
NEvents=0;
burst_number++;
cout<<"Burst Number+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "<<burst_number<<endl;

	/// \MemberDescr
	/// This method is called when a new file is opened in the ROOT TChain (corresponding to a start/end of burst in the normal NA62 data taking) + at the beginning of the first file\n
	/// Do here your start/end of burst processing if any
	/// \EndMemberDescr
}

void RunCalculations::Process(int iEvent){
NEvents++;
//MC bullshit
//cout<<"start"<<endl;
/*
Double_t InitialMomentum=0;
 TRecoSpectrometerEvent *SpectrometerEvent = (TRecoSpectrometerEvent*)GetEvent("Spectrometer");
if(SpectrometerEvent->GetNCandidates() == 0)return;
TRecoSpectrometerCandidate *StrawCand;
StrawCand = (TRecoSpectrometerCandidate*)SpectrometerEvent->GetCandidate(0);
*/
//cout<<StrawCand->GetChamberId(3)<<endl;
Bool_t MCflag=false;
Bool_t get_particle=false;

if(get_particle==true){
 Double_t InitialMomentum=0;
 TRecoSpectrometerEvent *SpectrometerEvent = (TRecoSpectrometerEvent*)GetEvent("Spectrometer");
if(SpectrometerEvent->GetNCandidates() == 0)return;
TRecoSpectrometerCandidate *StrawCand;
StrawCand = (TRecoSpectrometerCandidate*)SpectrometerEvent->GetCandidate(0);

        bool withMC = true;
        if(fMCSimple.fStatus == MCSimple::kMissing){
//              Event* MCTruthEvent = GetMCEvent();
//              for(int i=0; i<MCTruthEvent->GetNKineParts(); i++){
//                      FillHisto("pdgID", ((KinePart*)MCTruthEvent->GetKineParts()->At(i))->GetParticleName(), 1);
//              }
                withMC = false;
        }
        if(fMCSimple.fStatus == MCSimple::kEmpty) withMC = false;
if(withMC==false)return; 
 
 Int_t particle = 1; //1 for kpinn/pipi0, 0 for mu nu
//        Double_t InitialMomentum;//= new TVector3();
//	Double_t InitialMomentum_mu;
	Double_t ProdPos;
//	Double_t ProdPos_mu;
//cout<<"with mc"<<endl; 
if (particle==1){
//	TVector3 *ProdPos_pi = new TVector3();
//Double_t ProdPos_pi;
	ProdPos = fMCSimple["pi+"][0]->GetProdPos().Z();
//        TVector3 *InitialMomentum_pi = new TVector3();
        InitialMomentum = fMCSimple["pi+"][0]->GetInitialMomentum().Mag()/1000;
                Int_t bin_i=10;
                        if (InitialMomentum>=15&&InitialMomentum<20)bin_i=0;
                        else if (InitialMomentum>=20&&InitialMomentum<25)bin_i=1;
                        else if (InitialMomentum>=25&&InitialMomentum<30)bin_i=2;
                        else if (InitialMomentum>=30&&InitialMomentum<35)bin_i=3;
                	FillHisto(mc_pi[bin_i].Data(), InitialMomentum);
}
else{
//	TVector3 *ProdPos_mu = new TVector3();
        ProdPos = fMCSimple["mu+"][0]->GetProdPos().Z();
//        TVector3 *InitialMomentum_mu = new TVector3();
        InitialMomentum = fMCSimple["mu+"][0]->GetInitialMomentum().Mag()/1000;
}
//cout<<"prod_pos"<<endl; 
	if(ProdPos < 105000||ProdPos > 165000)return; 
 	FillHisto("init_p", InitialMomentum);  
Int_t bin_p=10;
if (InitialMomentum>=15&&InitialMomentum<20)bin_p=0;
else if (InitialMomentum>=20&&InitialMomentum<25)bin_p=1;
else if (InitialMomentum>=25&&InitialMomentum<30)bin_p=2;
else if (InitialMomentum>=30&&InitialMomentum<35)bin_p=3;
//else bin=5;
for(int i=0; i<4; i++){
if(i==bin_p)FillHisto(mc_p[i].Data(), InitialMomentum);
}
for(int i=0; i<4; i++){
if (bin_p==10)continue; 
Int_t missed = -2;
if(i!=StrawCand->GetChamberId(i))missed=i;
if(missed!=-2)FillHisto(mc_chambers[bin_p].Data(), i);
}
//if(i==bin_p/*&&StrawCand->GetChamberId(3)==3*/)FillHisto(mc_chambers[i].Data(), StrawCand->GetNChambers());
//if(i==bin_p&&missed!=-2)FillHisto(mc_chambers[bin_p].Data(), StrawCand->GetChamberId(missed));
if(StrawCand->GetChamberId(3)!=3){
TVector3 PiPosStart = StrawCand->GetPositionAfterMagnet();
Double_t Posx = StrawCand->GetSlopeXAfterMagnet()*(218885-PiPosStart.Z())+PiPosStart.X();
Double_t Posy = StrawCand->GetSlopeYAfterMagnet()*(218885-PiPosStart.Z())+PiPosStart.Y();
FillHisto("pos", Posx, Posy);
}

	if (bin_p !=10){
//	TRecoSpectrometerEvent *SpectrometerEvent = (TRecoSpectrometerEvent*)GetEvent("Spectrometer");
	Strawa->Set(SpectrometerEvent,1,1,1, MCflag);
		if (Strawa->SingleTrack()==1){
		Int_t bin_s=10;
		TVector3 pmom = Strawa->GetPionMomentum();
			if (pmom.Mag()>=15&&pmom.Mag()<20)bin_s=0;
			else if (pmom.Mag()>=20&&pmom.Mag()<25)bin_s=1;
			else if (pmom.Mag()>=25&&pmom.Mag()<30)bin_s=2;
			else if (pmom.Mag()>=30&&pmom.Mag()<35)bin_s=3;
			if(bin_s==2)FillHisto("Retarded",pmom.Mag());
			 if(bin_p==2)cout<<"wtf????"<<pmom.Mag()<<endl;
cout<<"bin#:  "<<bin_p<<endl;
	//
//if(pmom.Mag()>35||pmom.Mag()<15)return; 
if(pmom.Mag()<35&&pmom.Mag()>15){if (bin_p==bin_s)FillHisto(mc_ps[bin_p].Data(), pmom.Mag());};
//if(bin_p==bin_s/*&&StrawCand->GetChamberId(3)==3*/)FillHisto(mc_chambersp[bin_p].Data(), StrawCand->GetNChambers());
if(bin_p==bin_s&&StrawCand->GetChamberId(2)==2)FillHisto(mc_chambersp[bin_p].Data(), StrawCand->GetChamberId(2));
// 		if (bin_s!=bin_p)FillHisto("mismatch", bin_s);
 //		else FillHisto("mismatch", 9); 
 }
 }
}
// TRecoSpectrometerEvent *SpectrometerEvent = (TRecoSpectrometerEvent*)GetEvent("Spectrometer");
//cout<<"anything"<<endl;

        TRecoMUV3Event *MUV3Event = (TRecoMUV3Event*)GetEvent("MUV3");
        TRecoMUV1Event *MUV1Event = (TRecoMUV1Event*)GetEvent("MUV1");
        TRecoMUV2Event *MUV2Event = (TRecoMUV2Event*)GetEvent("MUV2");
        TRecoLKrEvent *LKrEvent = (TRecoLKrEvent*)GetEvent("LKr");
        TRecoCHANTIEvent *CHANTIEvent = (TRecoCHANTIEvent*)GetEvent("CHANTI");
        TRecoLAVEvent *LAVEvent = (TRecoLAVEvent*)GetEvent("LAV");
        TRecoRICHEvent *RICHEvent = (TRecoRICHEvent*)GetEvent("RICH");
        TRecoSpectrometerEvent *SpectrometerEvent = (TRecoSpectrometerEvent*)GetEvent("Spectrometer");
        TRecoGigaTrackerEvent *GigaTrackerEvent = (TRecoGigaTrackerEvent*)GetEvent("GigaTracker");
        TRecoCHODEvent *CHODEvent = (TRecoCHODEvent*)GetEvent("CHOD");
        TRecoCedarEvent *CedarEvent = (TRecoCedarEvent*)GetEvent("Cedar");
        TRecoMUV0Event *MUV0Event = (TRecoMUV0Event*)GetEvent("MUV0");
        TRecoSACEvent *SACEvent = (TRecoSACEvent*)GetEvent("SAC");
        TRecoSAVEvent *SAVEvent = (TRecoSAVEvent*)GetEvent("SAV");
        TRecoIRCEvent *IRCEvent = (TRecoIRCEvent*)GetEvent("IRC");
	NA62Analysis::UserMethods::OutputState isEvtValid;
//	 NA62Analysis::UserMethods::OutputState isTypeValid;
//	TRecoGigaTrackerEvent* GigaTrackerEvent = (TRecoGigaTrackerEvent*)GetOutput("GigaTrackerEvtReco.GigaTrackerEvent", isEvtValid);

//	if (isEvtValid==kOInvalid)return;
//TypeCounter++; 
//if(TypeCounter>3461)return;//GTK123 line 26 of list 
//if(TypeCounter>4263)return;//GTK12 line 86 of list
//if(TypeCounter>8929)return; //GTK13 line 38 of list
//if(TypeCounter>234)return; //GTK23 line 58 of list (these are all relative
//Int_t events = 23;
//FillHisto("events", events);

//Int_t EventNumber = (Int_t*)GetOutput("GigaTrackerEvtReco2.EventNumberw", isTypeValid);
//if(isTypeValid)FillHisto("EventNumberw", EventNumberw); 

//Bool_t MCflag = false; 

PreAnalyzers* Pre = new PreAnalyzers();
Pre->LKrMToG(LKrEvent);
L0TPData *L0 = GetL0Data();
Int_t RunNumber = 3809;//GetRawHeader()->GetRunID(); 
Tracka->Set(CedarEvent, GigaTrackerEvent, CHANTIEvent, SpectrometerEvent, CHODEvent, LKrEvent, RICHEvent, MUV1Event, MUV2Event, MUV3Event, LAVEvent, IRCEvent, SACEvent, SAVEvent, L0, MCflag, RunNumber); 
//Tracka->One->Set(SpectrometerEvent, CHODEvent, CedarEvent, GigaTrackerEvent, CHANTIEvent, LKrEvent, RICHEvent,  L0, MCflag, 3813); 
if(Tracka->One->OneTrackSelection()==1){
Tracka->Calo->Set(MUV1Event, MUV2Event, MUV3Event, Tracka->One->GetLKrCand(), Tracka->One->GetReft(), Tracka->One->Getmmsq(), Tracka->One->GetPionPosition(), Tracka->One->GetPionMomentum(), MCflag);
Tracka->Photon->Set(LKrEvent, LAVEvent, IRCEvent, SACEvent, SAVEvent, Tracka->One->GetLKrCand(), Tracka->One->GetReft(), Tracka->One->Getmmsq(), Tracka->One->LKra->GetSinglePionCluster_unmatched(), Tracka->One->GetLKrPosition(), MCflag);
Double_t mmsq = Tracka->One->Getmmsq(); 
Int_t bin = 10;

//unblock for acceptance calculations
//if (InitialMomentum<15)bin=0;
/*
 if (InitialMomentum>=15&&InitialMomentum<20)bin=0;
else if (InitialMomentum>=20&&InitialMomentum<25)bin=1;
else if (InitialMomentum>=25&&InitialMomentum<30)bin=2;
else if (InitialMomentum>=30&&InitialMomentum<35)bin=3;
//else bin=5;
for(int i=0; i<4; i++){
if(i==bin&&((mmsq>0&&mmsq<0.01)||(mmsq>0.026&&mmsq<0.068)))FillHisto(acc_mmsq[i].Data(), Tracka->One->Getmmsq());
}
*/
//Tracka->PhotonRejectionFactor();
//if (Tracka->PhotonRejectionFactor()==1)FillHisto("mmsq_SREvents",Tracka->One->Getmmsq());
//Tracka->CalorimeterRejectionFactor();
//Tracka->RICHRejectionFactor();
//Tracka->PionKineFactor();
//Tracka->MuonKineFactor();
//Tracka->MuonAccidentalPhotonFactor(); 
Tracka->RICHAnalysisPion();
//Tracka->CalorimeterPlots();
//Tracka->RICHAnalysisMuon();
//Tracka->ep(); 
//Tracka->MUV3Efficiency(); 
//Tracka->Kmu2Normalization(); 
/*
Double_t reftime = Tracka->One->GetReft(); 
TRecoSAVCandidate* SAVCand;
for(int i=0; i<SAVEvent->GetNCandidates(); i++){
  SAVCand = (TRecoSAVCandidate*)SAVEvent->GetCandidate(i);
  Double_t Tdiff = reftime - SAVCand->GetTime();
  FillHisto("SAV_reft", Tdiff);
}
*/

}
/*
else if (Tracka->One->BackgroundSelection()==1){
Tracka->Calo->Set(MUV1Event, MUV2Event, MUV3Event, Tracka->One->GetLKrCand(), Tracka->One->GetReft(), Tracka->One->Getmmsq(), Tracka->One->GetPionPosition(), Tracka->One->GetPionMomentum(), MCflag);
Tracka->Photon->Set(LKrEvent, LAVEvent, IRCEvent, SACEvent, SAVEvent, Tracka->One->GetLKrCand(), Tracka->One->GetReft(), Tracka->One->Getmmsq(), Tracka->One->LKra->GetSinglePionCluster_unmatched(), Tracka->One->GetLKrPosition(), MCflag);
Tracka->BeamAccidentalPhotonFactor(); 
}
*/
//else return;
 
}

void RunCalculations::PostProcess(){
	/// \MemberDescr
	/// This function is called after an event has been processed by all analyzers. It could be used to free some memory allocated
	/// during the Process.
	/// \EndMemberDescr

}

void RunCalculations::EndOfBurstUser(){
	/// \MemberDescr
	/// This method is called when a new file is opened in the ROOT TChain (corresponding to a start/end of burst in the normal NA62 data taking) + at the end of the last file\n
	/// Do here your start/end of burst processing if any
	/// \EndMemberDescr
FillHisto("NEvents", NEvents); 
}

void RunCalculations::EndOfRunUser(){
	/// \MemberDescr
	/// This method is called at the end of the processing (corresponding to a end of run in the normal NA62 data taking)\n
	/// Do here your end of run processing if any\n
	/// \n
	/// You can also use this space to save plots in your output file.\n
	/// If you want to save them all, just call\n
	/// \code
	 	SaveAllPlots();
	Tracka->One->SaveAllPlots(); 
//	Tracka->Calo->SaveAllPlots();
//	Tracka->Photon->SaveAllPlots(); 
	Tracka->SaveAllPlots();
//	Tracka->CaloPlot->SaveAllPlots();  
	/// \endcode
	/// Or you can just save the ones you want with\n
	/// \code
	/// 	histogram->Write()\n
	///		fHisto.Get...("histoname")->Write();
	/// \endcode
	/// \n
	/// To run over a set of histograms you can use Iterators (HistoHandler::IteratorTH1,
	/// HistoHandler::IteratorTH2, HistoHandler::IteratorTGraph). You can use it to run over
	/// all the histograms or only a subset of histogram by using one of the two forms of
	/// GetIterator...  (replace ... by TH1, TH2 or TGraph)\n
	/// \code
	/// 	GetIterator...()
	/// \endcode
	/// will get an Iterator running over all histograms of this type while
	/// \code
	/// 	GetIterator...("baseName")
	/// \endcode
	/// will get an Iterator running only over the histograms of this type whose name starts
	/// with baseName.\n
	/// For more details and examples on how to use the Iterator after getting it, please refer
	/// to the HistoHandler::Iterator documentation.\n
	/// Although this is described here, Iterators can be used anywhere after the
	/// histograms have been booked.
	/// \EndMemberDescr

}

void RunCalculations::DrawPlot(){
	/// \MemberDescr
	/// This method is called at the end of processing to draw plots when the -g option is used.\n
	/// If you want to draw all the plots, just call\n
	/// \code
	 	DrawAllPlots();
	/// \endcode
	/// Or get the pointer to the histogram with\n
	/// \code
	/// 	fHisto.GetTH1("histoName");// for TH1
	/// 	fHisto.GetTH2("histoName");// for TH2
	/// 	fHisto.GetGraph("graphName");// for TGraph and TGraphAsymmErrors
	///     fHisto.GetHisto("histoName");// for TH1 or TH2 (returns a TH1 pointer)
	/// \endcode
	/// and manipulate it as usual (TCanvas, Draw, ...)\n
	/// \EndMemberDescr
}

RunCalculations::~RunCalculations(){
	/// \MemberDescr
	/// Destructor of the Analyzer. If you allocated any memory for class
	/// members, delete them here.
	/// \EndMemberDescr
}
