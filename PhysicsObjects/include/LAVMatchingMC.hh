// ---------------------------------------------------------------
// History:
//
// Created by Tommaso Spadaro (tommaso.spadaro@lnf.infn.it) 2015-12-13
//
// ---------------------------------------------------------------

#ifndef LAVMatchingMC_HH
#define LAVMatchingMC_HH 1

#include "Persistency.hh"
#include <iostream>
using namespace std;

class LAVMatchingMC 
{
public:
  LAVMatchingMC();
  ~LAVMatchingMC() {}

  void SetTimeCuts(Double_t valLow, Double_t valHigh) { fLowThrLAVTimeCut = valLow; fHighThrLAVTimeCut = valHigh;} ///< Set time cuts w/wo low & high threshold leading edges (ns)
  void SetNoisyChannel(Int_t val)           { fNoisyChannels.push_back(val); } ///< Define channel val as noisy (taken from hit->GetChannelID())
  Bool_t ChannelIsNoisy(Int_t id); ///< Check is a channel is declared noisy or not
  void SetMaskedStation(Int_t val, Bool_t masked = kTRUE);  ///< Set station with ID val to be masked/unmasked when trying to match. Input should range between 1 and 12. Second argument: kTRUE to mask, kFALSE to unmask.
  void MaskAllStations();   ///< Mask all stations; useful as a first step for matching in LAV12 only
  void UnmaskAllStations(); ///< Unmask all stations
  void SetReferenceTime(Double_t val)       { fRefTime = val;                } ///< Set time of reference detector to be compared with hit time
  Bool_t LAVHasTimeMatching(TRecoLAVEvent*); ///< Perform time matching, return true if at least 1 block is found 

  Int_t   GetNumberOfMatchedBlocks() { return fNMatched;    }
  Int_t*  GetIndexOfMatchedBlocks()  { return fIndexMatched;} ///< Returns an array of indices of matched LAV hits in the LAV hits array
  Int_t*  GetChannelIDOfMatchedBlocks()  { return fChannelIDMatched;}
  Int_t   GetNumberOfMatchedBlocksLowThr() { return fNMatchedLowThr; }
  Int_t   GetNumberOfMatchedBlocksHighThr(){ return fNMatchedHighThr;}
  Int_t   GetBestEdgeMaskOfMatchedBlocks() { return fBestEdgeMask;}
  Double_t GetBestTimeOfMatchedBlocks(){return fDeltaTBest;}

  void Clear(); ///< Clear event-by-event counters
  void Print(); ///< Print relevant information about an event (to be run after LAVHasTimeMatching)
  void PrintMaskedStations(); ///< Print which stations are currently masked

private:
  Double_t fRefTime;             ///< Reference time [ns]
  Double_t fLowThrLAVTimeCut;    ///< Cut wrt reference time [ns]
  Double_t fHighThrLAVTimeCut;   ///< Cut wrt reference time [ns]

  TRecoLAVEvent* fLAVEvent;

  Int_t    fNMatched;            ///< Number of matched blocks 
  Int_t    fIndexMatched[100];   ///< Index of matched blocks (store them up to 100)
  Int_t    fChannelIDMatched[100];///< ChannelID of matched blocks (store them up to 100)
  Int_t    fNMatchedLowThr;      ///< Number of matched blocks low threshold
  Int_t    fNMatchedHighThr;     ///< Number of matched blocks high threshold
  Int_t    fBestEdgeMask;        ///< Best edgemask matched
  Double_t fDeltaTBest;          ///< Best time matching among blocks with chosen edgemasks
  vector<Int_t> fNoisyChannels;  ///< Vector of channel IDs labelled as noisy: they are ignored by the matchig algorithm
  Int_t fLAVPriorityMask[16];
  Bool_t fMaskedStations[12];    ///< Set i-th component of array to 1 if station ID i+1 has to be masked when trying to match
};
#endif
