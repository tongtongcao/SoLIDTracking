//*************************************************//
//abstract base class for the track finder of      //
//various detector configuration of SoLID          //
//*************************************************//

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
#include "SoLIDGEMTracker.h"
#include "SoLIDECal.h"
#include "SoLIDUtility.h"
#include "SoLIDGEMHit.h"

using namespace std;

#define MAXWINDOWHIT 200
#define MAXNSEEDS 2000

class SoLKalFieldStepper;

class SoLKalTrackFinder 
{
public:
  SoLKalTrackFinder();
  virtual ~SoLKalTrackFinder() {;}
  
  virtual void SetGEMDetector(vector<SoLIDGEMTracker*> thetrackers);
  void SetECalDetector(SoLIDECal* theECal) { fECal = theECal; }
  void SetBPM(Double_t x, Double_t y);
  void SetTargetGeometry(Double_t& z, Double_t& center, Double_t& length);
  int  GetNSeeds() const { return fNSeeds; }
  
  //pure virtual function to be implimented in derived classes
#ifdef MCDATA
  virtual bool GetSeedEfficiency(int i) const = 0;
  virtual bool GetMCTrackEfficiency(int i) const = 0;
#endif
  
  virtual void Clear( Option_t* opt="" ) = 0;
  virtual void ProcessHits(TClonesArray* theTracks) = 0;
  
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
    DoubletSeed(SeedType& t, SoLIDGEMHit* a, SoLIDGEMHit* b, Double_t& mom,
                Double_t& theta, Double_t& phi, Double_t& q, ECType& f)
    :type(t), hita(a), hitb(b), isActive(kTRUE), initMom(mom), 
     initTheta(theta), initPhi(phi), charge(q), flag(f) {}
    
    void Deactive() { isActive = kFALSE; }
  };

  void CalCircle(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x3,
                 Double_t y3, Double_t* R,Double_t* Xc, Double_t* Yc);
  
  SoLKalFieldStepper*                  fFieldStepper;
  std::vector<SoLIDGEMTracker*>        fGEMTracker;
  SoLIDECal*                           fECal;
  TClonesArray*                        fCoarseTracks;
  Int_t                                fNTrackers;
  Int_t                                fNSeeds;
  Int_t                                fEventNum;
  Double_t                             fBPMX;
  Double_t                             fBPMY;
  Double_t                             fTargetPlaneZ;
  Double_t                             fTargetCenter;
  Double_t                             fTargetLength;
  Double_t                             fChi2PerNDFCut;
  vector<SoLIDCaloHit>*                fCaloHits;
  map< SeedType, vector<DoubletSeed> > fSeedPool;
  
  ClassDef(SoLKalTrackFinder,0)
};



#endif
