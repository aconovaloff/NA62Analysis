// ---------------------------------------------------------------
// History:
//
// Created by Tommaso Spadaro (tommaso.spadaro@lnf.infn.it) 2015-12-13
//
// ---------------------------------------------------------------

#include <iostream>
#include "LAVMatchingMC.hh"

using namespace std;

/// \class LAVMatchingMC
/// \Brief
/// Provides the output of various algorithms for time-matching LAV information with respect to a given reference detector
/// \EndBrief
/// 
/// \Detailed
/// \n
/// Provided a reference time is given in input, the class will search for LAV blocks hits in time, withing a given window (default +-5ns).
/// By default, certain blocks are masked as noisy when allocating an object of type LAVMatchingMC.
/// Among all of the reconstructed hits in time, only specific edge-combination types are considered.
/// Edge combination ("edge mask") is ranked from the most-ill defined (isolated trailings for low or high threshold),
/// to the most complete (both edges of low and high threshold). 
/// A "priority" of the edgemasks is therefore build, defined as follows:<br>
///  LAVPriorityMask = 0; // --SNH--     <br>
///  LAVPriorityMask = 1; // __ __ __ TL <br>
///  LAVPriorityMask = 2; // __ __ TH __ <br>
///  LAVPriorityMask = 3; // __ __ TH TL <br>
///  LAVPriorityMask = 4; // LL __ __ __ <br>
///  LAVPriorityMask = 5; // __ LH __ __ <br>
///  LAVPriorityMask = 6; // __ LH __ TL <br>
///  LAVPriorityMask = 7; // LL __ TH __ <br>
///  LAVPriorityMask = 8; // LL __ TH TL <br>
///  LAVPriorityMask = 9; // __ LH TH TL <-- to be selected <br>
///  LAVPriorityMask =10; // LL LH __ __ <-- to be selected <br>
///  LAVPriorityMask =11; // __ LH TH __ <-- to be selected <br>
///  LAVPriorityMask =12; // LL __ __ TL <-- to be selected <br>
///  LAVPriorityMask =13; // LL LH TH __ <-- to be selected <br>
///  LAVPriorityMask =14; // LL LH __ TL <-- to be selected <br>
///  LAVPriorityMask =15; // LL LH TH TL <-- to be selected <br>
///  \n
/// The first L/T letters in letter pairs above mean leading/trailing edges; the second L/H mean low/high thresholds.
/// The boolean LAVHasTimeMatching method returns true if at least one block is found in time with priority >=9.
/// The list of the blocks in time with the reference detector (of all priorities) can be obtained as output via the GetIndexOfMatchedBlocks method.
/// The highest priority available ("bestEdgeMask") can be output via the GetBestEdgeMaskOfMatchedBlocks method.
/// The time difference with the smallest absolute value can be output via the GetBestTimeOfMatchedBlocks. 
/// For further details, see: https://indico.cern.ch/event/489423/contribution/1/attachments/1217640/1778693/LAVMatchingMCClass_pinunu_wg_26_1_2016-1.pdf
/// \n Example use-case code:
/// \code
/// // Init part of an Analyzer
/// fLAVMatchingMC = new LAVMatchingMC(); // fLAVMatchingMC is a private LAVMatchingMC* variable
/// ...
/// // Event process part of an Analyzer:
/// fLAVMatchingMC->SetReferenceTime(CHODTime);
/// Bool_t matched = fLAVMatchingMC->LAVHasTimeMatching(LAVEvent);
/// if (fVerboseFlag) {	
///   cout << "Did LAV match? " << matched << endl;
///   fLAVMatchingMC->Print();
/// }
/// \endcode
/// \author Tommaso Spadaro (tommaso.spadaro@lnf.infn.it)
/// \EndDetailed

LAVMatchingMC::LAVMatchingMC() {
  Clear();
  fRefTime = -999;

  fLAVPriorityMask[0] = 0; // --SNH--
  fLAVPriorityMask[1] = 4; // LL __ __ __ 
  fLAVPriorityMask[2] = 5; // __ LH __ __ 
  fLAVPriorityMask[3] =10; // LL LH __ __ <-- to be selected
  fLAVPriorityMask[4] = 2; // __ __ TH __
  fLAVPriorityMask[5] = 7; // LL __ TH __
  fLAVPriorityMask[6] =11; // __ LH TH __ <-- to be selected
  fLAVPriorityMask[7] =13; // LL LH TH __ <-- to be selected
  fLAVPriorityMask[8] = 1; // __ __ __ TL
  fLAVPriorityMask[9] =12; // LL __ __ TL <-- to be selected
  fLAVPriorityMask[10]= 6; // __ LH __ TL
  fLAVPriorityMask[11]=14; // LL LH __ TL <-- to be selected
  fLAVPriorityMask[12]= 3; // __ __ TH TL
  fLAVPriorityMask[13]= 8; // LL __ TH TL
  fLAVPriorityMask[14]= 9; // __ LH TH TL <-- to be selected
  fLAVPriorityMask[15]=15; // LL LH TH TL <-- to be selected

  fLowThrLAVTimeCut  = 5.; // ns
  fHighThrLAVTimeCut = 5.; // ns

  // Default for noisy channels for runs 3809,4068/9,4141 - to be checked for other runs
  SetNoisyChannel(84083);  //  227 of LAV8
  SetNoisyChannel(90083);  //   35 of LAV9
  SetNoisyChannel(103040); //  196 of LAV10
  SetNoisyChannel(121101); //  105 of LAV12
  SetNoisyChannel(123083); //  227 of LAV12
  for (Int_t i=0; i<12; i++) fMaskedStations[i] = kFALSE;
}

Bool_t LAVMatchingMC::ChannelIsNoisy(Int_t id) {
  vector<Int_t>::iterator i = find(fNoisyChannels.begin(), fNoisyChannels.end(), id);
  return (i!=fNoisyChannels.end());
}

Bool_t LAVMatchingMC::LAVHasTimeMatching(TRecoLAVEvent* LAVEvent) {
  fLAVEvent = LAVEvent;
  TClonesArray& hitArray = (*(LAVEvent->GetHits()));

  fNMatched = 0;  
  fNMatchedLowThr = 0; 
  fNMatchedHighThr = 0;
  fBestEdgeMask = 0; 
  fDeltaTBest = -99999.;

  for (Int_t i=0; i<LAVEvent->GetNHits(); i++) {

    TRecoLAVHit* hit = (TRecoLAVHit*) hitArray[i];
    int station = hit->GetLAVID()-1;
//cout<<1<<endl;
    if (fMaskedStations[station]) continue;
//cout<<2<<endl;
    Int_t chid = hit->GetChannelID();
    Double_t dt = hit->GetTime() - fRefTime;
//cout<<hit->GetEdgeMask()<<endl;
//    if ( (hit->GetEdgeMask() & 1 && hit->GetEdgeMask() & 2)/* && TMath::Abs(dt)>fHighThrLAVTimeCut*/ ) continue; 
//    if (!(hit->GetEdgeMask() & 1 && hit->GetEdgeMask() & 2)/* && TMath::Abs(dt)>fLowThrLAVTimeCut */ ) continue;
//cout<<3<<endl;
    Bool_t noisymatch = kFALSE;
    for (UInt_t is=0; is<fNoisyChannels.size(); is++){
      if (fNoisyChannels.at(is) == hit->GetChannelID()) {
	noisymatch = kTRUE;
	break;
      }
    }
    if (noisymatch) continue;
//cout<<4<<endl;
    if (hit->GetEdgeMask() & 1 && hit->GetEdgeMask() & 8) fNMatchedLowThr++;      
    if (hit->GetEdgeMask() & 2 && hit->GetEdgeMask() & 4) fNMatchedHighThr++;
    if (fLAVPriorityMask[hit->GetEdgeMask()] > fBestEdgeMask) fBestEdgeMask = fLAVPriorityMask[hit->GetEdgeMask()];

    if (fLAVPriorityMask[hit->GetEdgeMask()] >= 9) {
      if (fabs(dt) < fabs(fDeltaTBest)) {
	fDeltaTBest = dt;
      }
      if (fNMatched<100) {
	fIndexMatched[fNMatched] = i;
	fChannelIDMatched[fNMatched] = chid;
	fNMatched++;
      }
    }
  }
  return (fBestEdgeMask>=9);
}

void LAVMatchingMC::SetMaskedStation(Int_t val, Bool_t masked) {
  if (val<1 || val>12) { // station IDs are 1-12
    cout << "LAVMatchingMC::SetMaskedStation >> Error: input station ID out of range. Allowed range is [1,12]. Input given is " << val << endl;
    return;
  }
  fMaskedStations[val-1] = masked; // indices of fMaskedStations are 0-11
}

void LAVMatchingMC::MaskAllStations() {
  for (Int_t i=1; i<=12; i++) SetMaskedStation(i);
}
void LAVMatchingMC::UnmaskAllStations() {
  for (Int_t i=1; i<=12; i++) SetMaskedStation(i, kFALSE);
}

void LAVMatchingMC::Clear() {
  fNMatched = 0;
  fNMatchedLowThr = 0;
  fNMatchedHighThr = 0;
  fBestEdgeMask = 0;
  fRefTime = -1000000;
  fNoisyChannels.clear();
}

void LAVMatchingMC::Print() {
  cout << "[LAVMatchingMC]: Printing---" << endl;
  cout << "Have matched " << fNMatched << " blocks, the best of which has edgemask priority (ranging from 1 to 15) = " << fBestEdgeMask << endl;
  cout << "Have matched " << fNMatchedLowThr << " low-threshold blocks " << endl;
  cout << "Have matched " << fNMatchedHighThr << " high-threshold blocks " << endl;
  if (!(fLAVEvent)) {
    cout << "Wrong LAVEvent provided " << endl;
    return;
  }

  TClonesArray& hitArray = (*(fLAVEvent->GetHits()));

  for (Int_t i=0; i<TMath::Min(fNMatched,100); i++) {
    if (fIndexMatched[i]<0 || fIndexMatched[i]>= fLAVEvent->GetNHits()) {
      cout << " LAVMatchingMC Print: wrong index stored or wrong LAVEvent provided " << i << " " << fIndexMatched[i] << " " << fLAVEvent->GetNHits() << endl;
      return;
    }
    TRecoLAVHit* hit = (TRecoLAVHit*) hitArray[fIndexMatched[i]];
    cout << " Matched block " << i << " chid " << hit->GetChannelID() << " edge " << hit->GetEdgeMask() << " priority " << fLAVPriorityMask[hit->GetEdgeMask()] << endl;
  }
  cout << "[LAVMatchingMC]: Printing end" << endl;
}

void LAVMatchingMC::PrintMaskedStations() {
  cout << "[LAVMatchingMC] Masked stations are:";
  Int_t NumberOfMaskedStations = 0;
  for (Int_t i=0; i<12; i++) {
    if (fMaskedStations[i]) {
      cout << " " << i+1;
      NumberOfMaskedStations++;
    }
  }
  if (!NumberOfMaskedStations) cout << " none";
  cout << endl;
}
