// ---------------------------------------------------------------
// History:
//
// Created by Tommaso Spadaro (tommaso.spadaro@lnf.infn.it) 2016-02-26
//
// ---------------------------------------------------------------

#include <iostream>
#include "SAVMatchingMC.hh"

using namespace std;

/// \class SAVMatchingMC
/// \Brief
/// Provides the output of various algorithms for time-matching SAV information with respect to a given reference detector. Only TDC-based data is used.
/// \EndBrief
/// 
/// \Detailed
/// \n
/// Provided a reference time is given in input, the class will search for IRC, SAC hits in time, withing a given window (default +-5ns).
/// Among all of the reconstructed hits in time, only specific edge-combination types are considered.
/// Edge combination ("edge mask") is ranked from the most-ill defined (isolated trailings for low or high threshold),
/// to the most complete (both edges of low and high threshold). 
/// A "priority" of the edgemasks is therefore build, defined as follows:<br>
///  SAVPriorityMask = 0; // --SNH--    <br>
///  SAVPriorityMask = 1; // __ __ __ TL<br>
///  SAVPriorityMask = 2; // __ __ TH __<br>
///  SAVPriorityMask = 3; // __ __ TH TL<br>
///  SAVPriorityMask = 4; // LL __ __ __<br>
///  SAVPriorityMask = 5; // __ LH __ __ <br>
///  SAVPriorityMask = 6; // __ LH __ TL<br>
///  SAVPriorityMask = 7; // LL __ TH __<br>
///  SAVPriorityMask = 8; // LL __ TH TL<br>
///  SAVPriorityMask = 9; // __ LH TH TL <-- to be selected <br>
///  SAVPriorityMask =10; // LL LH __ __ <-- to be selected <br>
///  SAVPriorityMask =11; // __ LH TH __ <-- to be selected <br>
///  SAVPriorityMask =12; // LL __ __ TL <-- to be selected <br>
///  SAVPriorityMask =13; // LL LH TH __ <-- to be selected <br>
///  SAVPriorityMask =14; // LL LH __ TL <-- to be selected <br>
///  SAVPriorityMask =15; // LL LH TH TL <-- to be selected <br>
///  \n
/// The integer SAVHasTimeMatching method returns a word with bit 0(1) set if at least one IRC(SAC) hit is found in time with priority >=9.
/// The list of the blocks in time with the reference detector (of all priorities) can be obtained as output via the GetIndexOfIRCMatchedBlocks and GetIndexOfSACMatchedBlocks methods.
/// The highest priority available ("bestEdgeMask") can be output via the GetBestEdgeMaskOfIRCMatchedBlocks and GetBestEdgeMaskOfSACMatchedBlocks methods.
/// The time difference with the smallest absolute value can be output via the GetBestTimeOfIRCMatchedBlocks and GetBestTimeOfSACMatchedBlocks methods. 
/// \n Example use-case code:
/// \code
/// // Init part of an Analyzer
/// fSAVMatchingMC = new SAVMatchingMC(); // fSAVMatchingMC is a private SAVMatchingMC* variable
/// ...
/// // Event process part of an Analyzer:
/// fSAVMatchingMC->SetReferenceTime(CHODTime);
/// Bool_t matched = fSAVMatchingMC->SAVHasTimeMatching(IRCEvent, SACEvent);
/// if (fVerboseFlag) {	
///   cout << "Did any of SAV match? " << matched << endl;
///   fSAVMatchingMC->Print();
/// }
/// \endcode
/// \author Tommaso Spadaro (tommaso.spadaro@lnf.infn.it)
/// \EndDetailed

SAVMatchingMC::SAVMatchingMC() {
  Clear();
  fRefTime = -999;

  fSAVPriorityMask[0] = 0; // --SNH--
  fSAVPriorityMask[1] = 4; // LL __ __ __ 
  fSAVPriorityMask[2] = 5; // __ LH __ __ 
  fSAVPriorityMask[3] =10; // LL LH __ __ <-- to be selected
  fSAVPriorityMask[4] = 2; // __ __ TH __
  fSAVPriorityMask[5] = 7; // LL __ TH __
  fSAVPriorityMask[6] =11; // __ LH TH __ <-- to be selected
  fSAVPriorityMask[7] =13; // LL LH TH __ <-- to be selected
  fSAVPriorityMask[8] = 1; // __ __ __ TL
  fSAVPriorityMask[9] =12; // LL __ __ TL <-- to be selected
  fSAVPriorityMask[10]= 6; // __ LH __ TL
  fSAVPriorityMask[11]=14; // LL LH __ TL <-- to be selected
  fSAVPriorityMask[12]= 3; // __ __ TH TL
  fSAVPriorityMask[13]= 8; // LL __ TH TL
  fSAVPriorityMask[14]= 9; // __ LH TH TL <-- to be selected
  fSAVPriorityMask[15]=15; // LL LH TH TL <-- to be selected

  fLowThrSAVTimeCut[0] = 5.; // [ns] default for IRC
  fHighThrSAVTimeCut[0] = 5.; // [ns] default for IRC
  fLowThrSAVTimeCut[1] = 5.; // [ns] default for SAC
  fHighThrSAVTimeCut[1] = 5.; // [ns] default for SAC
}

Int_t SAVMatchingMC::SAVHasTimeMatching(TRecoIRCEvent* IRCEvent, TRecoSACEvent* SACEvent) {

  for (Int_t i=0; i<2; i++) {
    fNMatched[i] = 0;  
    fBestEdgeMask[i] = 0; 
    fDeltaTBest[i] = 999.;
  }

  TRecoVEvent* fSAVEvent[2] = {IRCEvent, SACEvent};
  Int_t outputFlag = 0;

  for (Int_t iDet = 0; iDet<2; iDet++) { // loop for IRC, SAC matching

    TClonesArray& hitArray = (* (fSAVEvent[iDet]->GetHits()));    

    for (Int_t i=0; i< fSAVEvent[iDet]->GetNHits(); i++) {
      Int_t edgeMask;
      if (iDet==0) edgeMask = ((TRecoIRCHit*) hitArray[i])->GetEdgeMask();
      else edgeMask = ((TRecoSACHit*) hitArray[i])->GetEdgeMask();

      TRecoVHit* hit = (TRecoVHit*) hitArray[i];
      int chid = hit->GetChannelID();
      double dt = hit->GetTime() - fRefTime;
    
//      if ( (edgeMask & 1 && edgeMask & 2)/* && TMath::Abs(dt)>fHighThrSAVTimeCut[iDet]*/ ) continue; 
//      if (!(edgeMask & 1 && edgeMask & 2)/* && TMath::Abs(dt)>fLowThrSAVTimeCut[iDet]*/ ) continue;

      bool noisymatch = kFALSE;
      for (int is = 0; is< (int) fNoisyChannels[iDet].size(); is++){
	if (fNoisyChannels[iDet].at(is) == hit->GetChannelID()) {
	  noisymatch = kTRUE;
	  break;
	}
      }
      if (noisymatch) continue;    

      if (fSAVPriorityMask[edgeMask] > fBestEdgeMask[iDet]) fBestEdgeMask[iDet] = fSAVPriorityMask[edgeMask];

      if (fSAVPriorityMask[edgeMask] >= 9) {
	if (fabs(dt) < fabs(fDeltaTBest[iDet])){
	  fDeltaTBest[iDet] = dt;
	}

	if (fNMatched[iDet]<100) {
	  fIndexMatched[iDet][fNMatched[iDet]] = i;
	  fChannelIDMatched[iDet][fNMatched[iDet]] = chid;
	  fNMatched[iDet]++;
	}
      }
    }
    if (fBestEdgeMask[iDet]>=9) outputFlag |= (1>>iDet); 
  }

  return outputFlag;
}


void SAVMatchingMC::Clear() {
  for (Int_t i=0; i<2; i++) {
    fNMatched[i] = 0;  
    fBestEdgeMask[i] = 0; 
    fDeltaTBest[i] = 999; 
    fNoisyChannels[i].clear();  
  }
  fRefTime = -1000000; 
}

void SAVMatchingMC::Print() {
  cout << "[SAVMatchingMC]: Printing---" << endl;
  cout << "Have matched " << fNMatched[0] << " IRC blocks, the best of which has edgemask priority (ranging from 1 to 15) = " << fBestEdgeMask[0] << endl;
  cout << "Have matched " << fNMatched[1] << " SAC blocks, the best of which has edgemask priority (ranging from 1 to 15) = " << fBestEdgeMask[1] << endl;

  if (!(fSAVEvent[0]) || !(fSAVEvent[1])) {
    cout << "Wrong SAVEvents provided " << fSAVEvent[0] << " " << fSAVEvent[1] << endl;
    return;
  }
  TString DetS[2] = {"IRC","SAC"};
  for (Int_t iDet = 0; iDet<2; iDet++) {
    TClonesArray& hitArray = (* (fSAVEvent[iDet]->GetHits()));  

    for (Int_t i=0; i<fNMatched[iDet]; i++) {
      if (fIndexMatched[iDet][i]<0 || fIndexMatched[iDet][i]>= fSAVEvent[iDet]->GetNHits()) {
	cout << " SAVMatchingMC Print: wrong index stored or wrong SAVEvent provided " << i << " " << fIndexMatched[iDet][i] << " " << fSAVEvent[iDet]->GetNHits() << DetS[iDet].Data()<<endl;
	return;
      }
      // TRecoVHit* hit = (TRecoVHit*) hitArray[fIndexMatched[iDet][i]];
      Int_t edgeMask;
      if (iDet==0) edgeMask = ((TRecoIRCHit*) hitArray[i])->GetEdgeMask();
      else edgeMask = ((TRecoSACHit*) hitArray[i])->GetEdgeMask();

      cout << " Matched block " << i << " chid " << fChannelIDMatched[iDet][i] << " edge " << edgeMask << " priority " << fSAVPriorityMask[edgeMask] << DetS[iDet].Data()<< endl;
    }
  }
  cout << "[SAVMatchingMC]: Printing end" << endl;

}
