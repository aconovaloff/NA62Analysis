#include "RICHCandidate.hh"

RICHCandidate::RICHCandidate()
{}

RICHCandidate::~RICHCandidate()
{}

void RICHCandidate::Clear() {
  fIsRICHCandidate = 0;
  fDiscriminant = 99999999.;
  fDX           = -9999999.;
  fDY           = -9999999.;
  fXcenter      = -99999.;
  fYcenter      = -99999.;
  fRadius       = 0.;
  fDeltaTime    = 99999999.;
  fChi2         = 99999999.;
  fNHits        = 0;
  fMirrorID     = 99;
}
