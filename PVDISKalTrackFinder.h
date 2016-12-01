#ifndef ROOT_PVDIS_KAL_TRACK_FINDER
#define ROOT_PVDIS_KAL_TRACK_FINDER
//c++
#include <iostream>
#include <map>
#include <cassert>
#include <vector>
//ROOT
#include "TClonesArray.h"
#include "TSeqCollection.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVector2.h"

//SoLIDTracking
#include "SoLIDTrack.h"
#include "SoLIDGEMHit.h"
#include "SoLIDUtility.h"
#include "SoLIDFieldMap.h"
#include "SoLIDGEMTracker.h"
#include "SoLIDGEMChamber.h"
#include "SoLIDECal.h"
#include "SoLKalMatrix.h"
#include "SoLKalFieldStepper.h"
#include "SoLKalTrackFinder.h"


using namespace std;

class PVDISKalTrackFinder : public SoLKalTrackFinder
{
public:
  PVDISKalTrackFinder() {;}
  PVDISKalTrackFinder(bool isMC);
  ~PVDISKalTrackFinder();
  
  void ProcessHits(TClonesArray* theTracks);
  void Clear( Option_t* opt="" );
  void SetGEMDetector(vector<SoLIDGEMTracker*> thetrackers);
  
  bool GetSeedEfficiency(int /*i*/) const { return fSeedEfficiency; }
  bool GetMCTrackEfficiency(int /*i*/) const { return fMcTrackEfficiency; }
  
protected:
  void FindDoubletSeed(Int_t planej, Int_t planek);
  void TrackFollow();
  void MergeSeed();
  void FindandAddVertex();
  void FinalSelection(TClonesArray* theTracks);
  
  SoLKalTrackSite & SiteInitWithSeed(DoubletSeed* thisSeed);
  Bool_t     ECCoarseCheck(SoLIDGEMHit* theHit, Int_t & index);
  double     CalDeltaPhi(const double & phi1, const double & phi2); 
  double     CalDeltaR(const double & r1, const double & r2);
  void       Rotate(Double_t& x, Double_t& y);
  Double_t   StraightLinePredict(const Double_t& x1, const Double_t& z1, const Double_t& x2, 
                                 const Double_t& z2, const Double_t& targetZ);
  int GetHitsInWindow(int plane, double x, double wx, double y, double wy, bool flag = false);
  SoLIDGEMHit* FindCloestHitInWindow(double &x, double &y);
  Bool_t CheckChargeAsy(SoLKalTrackSystem* theSystem);
  Double_t FindVertexZ(SoLKalTrackState* thisState);
  void CopyTrack(SoLIDTrack* soltrack, SoLKalTrackSystem* kaltrack);  
  
  
  bool fIsMC;
  bool fSeedEfficiency;
  bool fMcTrackEfficiency;
  Double_t fRefPhi;
  Double_t fRefSin;
  Double_t fRefCos;
  vector<SoLIDGEMHit*> fWindowHits;
  map< Int_t, vector<SoLIDGEMHit*> > fGoodHits;
  Int_t fNGoodTrack;
};
#endif
