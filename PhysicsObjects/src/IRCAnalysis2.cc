#include "IRCAnalysis2.hh"
//#include "AnalysisTools.hh"
//#include "Parameters.hh"

IRCAnalysis2::IRCAnalysis2(NA62Analysis::Core::BaseAnalysis *ba) { 
//  fPar = Parameters::GetInstance(); 
  
  // IRC priority definition
  fIRCPriorityMask[0] = 0; // --SNH--
  fIRCPriorityMask[1] = 4; // LL __ __ __
  fIRCPriorityMask[2] = 5; // __ LH __ __
   fIRCPriorityMask[3] =10; // LL LH __ __
   fIRCPriorityMask[4] = 2; // __ __ TH __
    fIRCPriorityMask[5] = 7; // LL __ TH __
   fIRCPriorityMask[6] =11; // __ LH TH __
  fIRCPriorityMask[7] =13; // LL LH TH __
  fIRCPriorityMask[8] = 1; // __ __ __ TL
  fIRCPriorityMask[9] =12; // LL __ __ TL
  fIRCPriorityMask[10]= 6; // __ LH __ TL
  fIRCPriorityMask[11]=14; // LL LH __ TL
  fIRCPriorityMask[12]= 3; // __ __ TH TL
  fIRCPriorityMask[13]= 8; // LL __ TH TL
  fIRCPriorityMask[14]= 9; // __ LH TH TL
  fIRCPriorityMask[15]=15; // LL LH TH TL

 // fTools = AnalysisTools::GetInstance();
  fUserMethods = new UserMethods(ba);

  // Histo booking
  fUserMethods->BookHisto(new TH2F("IRCAnalysis_hit_chtime","",10,0,10,400,-50,50));
  fUserMethods->BookHisto(new TH2F("IRCAnalysis_hit_tottime","",400,0,200,400,-50,50));
  fUserMethods->BookHisto(new TH2F("IRCAnalysis_hit_besttime","",400,0,200,400,-50,50));

  // Photon candidates
  for (Int_t kk=0; kk<20; kk++) fPhotonCandidate[kk] = new PhotonVetoCandidate();
}

IRCAnalysis2::~IRCAnalysis2()
{ }

void IRCAnalysis2::StartBurst(Int_t burstid) {
}

void IRCAnalysis2::Clear() {
  for (Int_t kk=0; kk<20; kk++) fPhotonCandidate[kk]->Clear(); 
}

Int_t IRCAnalysis2::MakeCandidate(Double_t reftime, TRecoIRCEvent* event) {

  TClonesArray& Hits = (*(event->GetHits()));

  // At least a good hit within 5 ns from the reference time
  Int_t bestIRCHitType = 0;
  for (Int_t jHit=0; jHit<event->GetNHits(); jHit++) {
    TRecoIRCHit *hit = (TRecoIRCHit*)Hits[jHit];
    Int_t chid = hit->GetChannelID();
    Int_t edge = hit->GetEdgeMask();
    Double_t dtime = hit->GetTime()-reftime;
    if (fIRCPriorityMask[edge]>bestIRCHitType && fabs(dtime)<5) bestIRCHitType = fIRCPriorityMask[edge];
  }   
  if (bestIRCHitType<9) return 0;

  Int_t nirc = 0;
  Double_t etot = 0;
  Double_t maxtot = 0;
  Double_t maxtime = 999999.;
  for (Int_t jHit=0; jHit<event->GetNHits(); jHit++) {
    TRecoIRCHit *hit = (TRecoIRCHit*)Hits[jHit];
    Int_t chid = hit->GetChannelID();
    Int_t edge = hit->GetEdgeMask();
    Double_t dtime = hit->GetTime()-reftime;
    Double_t tot = (hit->GetTrailingEdgeLow() && hit->GetLeadingEdgeLow()) ? hit->GetTrailingEdgeLow()-hit->GetLeadingEdgeLow() : 0;
    Double_t tslew = 0;
    if (tot>2&&tot<40) tslew = 6.38-0.303*tot+0.003578*tot*tot;
    if (tot>=40 && tot<60) tslew = 6.38-0.303*40+0.003578*40*40;
    if (tot>2&&tot<15) tslew += 1.2;
    dtime -= tslew;
    fUserMethods->FillHisto("IRCAnalysis_hit_chtime",chid,dtime);
    fUserMethods->FillHisto("IRCAnalysis_hit_tottime",tot,dtime);
    etot += tot;
    if (tot>=maxtot) {
      maxtot = tot;
      maxtime = dtime;
    }
    if (nirc<20) {
      fPhotonCandidate[nirc]->SetTime(dtime+reftime);
      fPhotonCandidate[nirc]->SetID(jHit);
    }
    nirc++;
  }
  if (nirc) fUserMethods->FillHisto("IRCAnalysis_hit_besttime",etot,maxtime);

  if (nirc && fabs(maxtime)<5) return 1; 
  return 0;
}
