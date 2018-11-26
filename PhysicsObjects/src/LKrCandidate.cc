#include "LKrCandidate.hh"

LKrCandidate::LKrCandidate()
{}

LKrCandidate::~LKrCandidate()
{}

void LKrCandidate::Clear() {
  fIsLKrCandidate = 0;
  fDiscriminant = 99999999.;
  fDeltaTime = 99999999.;
  fX = -99999.;
  fY = -99999.;
  fEnergy = 0.;
  fSeedEnergy = 0.;
  fNCells = 0;
  fTime = -999999.;
}
