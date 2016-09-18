#ifndef ROOT_SOL_KAL_TRACK_FINDER
#define ROOT_SOL_KAL_TRACK_FINDER
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


using namespace std;

class SoLKalTrackFinder
{
  public:
  SoLKalTrackFinder() {;}
  SoLKalTrackFinder(bool isMC, Int_t ntrackers);
  ~SoLKalTrackFinder();
  
  void SetGEMDetector(vector<SoLIDGEMTracker*> thetrackers);
  void SetECalDetector(SoLIDECal* theECal) { fECal = theECal; }
  void ProcessHits(TClonesArray* theTracks);
  
  void SetBPM(Double_t x, Double_t y);
  void CalCircle(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x3,
                 Double_t y3, Double_t* R,Double_t* Xc, Double_t* Yc);
                           
  void Clear( Option_t* opt="" );
  int  GetNSeeds() const { return fNSeeds; }
  bool GetSeedEfficiency(int i) const { return fSeedEfficiency[i]; }
  bool GetMCTrackEfficiency(int i) const { return fMcTrackEfficiency[i]; }

  
  
  protected:
  
  struct DoubletSeed{
    SeedType type;
    SoLIDGEMHit* hita;
    SoLIDGEMHit* hitb;
    Bool_t isActive;
    Double_t initMom;
    Double_t initTheta;
    Double_t initPhi;
    Double_t charge;
    ECType flag;
    
    DoubletSeed() {}
    DoubletSeed(SeedType t, SoLIDGEMHit* a, SoLIDGEMHit* b, Double_t mom,
                Double_t theta, Double_t phi, Double_t q, ECType f)
    :type(t), hita(a), hitb(b), isActive(kTRUE), initMom(mom), 
     initTheta(theta), initPhi(phi), charge(q), flag(f) {}
    
    void Deactive() { isActive = kFALSE; }
  };
  
  //Main analysis functions
  void FindDoubletSeed(Int_t planej, Int_t planek, ECType type = kFAEC);
  void MergeSeed();
  void TrackFollow();
  void CoarseCheckVertex();
  void FindandAddVertex();
  void FinalSelection(TClonesArray* theTracks);
  void CopyTrack(SoLIDTrack* soltrack, SoLKalTrackSystem* kaltrack);
  void ECalFinalMatch();
  
  
  //assistent functions
  SoLKalTrackSite & SiteInitWithSeed(DoubletSeed* thisSeed);
  Bool_t TriggerCheck(SoLIDGEMHit* theHit, ECType type);
  SoLIDGEMHit* FindCloestHitInWindow(double &x, double &y);
  double CalDeltaPhi(double phi1, double phi2); 
  double CalDeltaR(double r1, double r2);
  double PredictR(Int_t &plane, SoLIDGEMHit* hit1, SoLIDGEMHit* hit2);
  int GetHitsInWindow(int plane, double x, double wx, double y, double wy, bool flag = false);
  Double_t FindVertexZ(SoLKalTrackState* thisState);
  Bool_t CheckChargeAsy(SoLKalTrackSystem* theSystem);
  void GetHitChamberList(vector<Int_t> &theList, Int_t thisChamber, Int_t size);
  Int_t GetChamIDFromPos(Double_t &x, Double_t &y, Int_t TrackerID);
  Bool_t CalInitParForPair(SoLIDGEMHit* hita, SoLIDGEMHit* hitb, Double_t &charge, 
                           Double_t& mom, Double_t& theta, Double_t& phi, ECType& type);
  Int_t BinarySearchForR(TSeqCollection* array, Double_t &lowr);
    
  bool fIsMC;
  std::vector<SoLIDGEMTracker*> fGEMTracker;
  SoLIDECal *fECal;
  int fNSeeds;
  bool fSeedEfficiency[2];
  bool fMcTrackEfficiency[2];
  Int_t fNTrackers;
  Int_t fEventNum;
  vector<SoLIDGEMHit*> fWindowHits;
  SoLKalFieldStepper* fFieldStepper;
  TClonesArray*      fCoarseTracks;
  Double_t fTargetPlaneZ;
  Double_t fTargetCenter;
  Double_t fTargetLength;
  Double_t fBPMX;
  Double_t fBPMY;
  map< Int_t, vector<SoLIDGEMHit*> > fGoodHits;
  map< SeedType, vector<DoubletSeed> > fSeedPool;
  Int_t fNGoodTrack;
  Double_t fChi2PerNDFCut;
  vector<SoLIDCaloHit> * fCaloHits;  
};

#endif

