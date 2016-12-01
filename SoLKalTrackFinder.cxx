//SoLIDTracking
#include "SoLKalTrackFinder.h"
#include "SoLKalFieldStepper.h"
#include "SoLKalTrackSystem.h"
#include "TVector2.h"

#define MAXNTRACKS 1000
ClassImp(SoLKalTrackFinder)
SoLKalTrackFinder::SoLKalTrackFinder()
: fGEMTracker(nullptr), fECal(nullptr), fNTrackers(0),fNSeeds(0), fEventNum(0),
  fBPMX(0), fBPMY(0), fChi2PerNDFCut(30.)
{
  fFieldStepper = SoLKalFieldStepper::GetInstance();
  fCoarseTracks = new TClonesArray("SoLKalTrackSystem", MAXNTRACKS, kTRUE);
  vector<DoubletSeed> midBackSeed;
  midBackSeed.reserve(MAXNSEEDS);
  midBackSeed.clear();
  vector<DoubletSeed> frontMidSeed;
  frontMidSeed.reserve(MAXNSEEDS);
  frontMidSeed.clear();
  vector<DoubletSeed> frontBackSeed;
  frontBackSeed.reserve(MAXNSEEDS);
  frontBackSeed.clear();
  
  fSeedPool[kMidBack] = midBackSeed;
  fSeedPool[kFrontMid] = frontMidSeed;
  fSeedPool[kFrontBack] = frontBackSeed;
}
//__________________________________________________________________________
void SoLKalTrackFinder::SetGEMDetector(vector<SoLIDGEMTracker*> thetrackers)
{
  fGEMTracker = thetrackers;
  fNTrackers  = (Int_t)thetrackers.size();
}
//__________________________________________________________________________
void SoLKalTrackFinder::SetBPM(Double_t x, Double_t y)
{
  fBPMX = x + gRandom->Gaus(0., 3e-4);
  fBPMY = y + gRandom->Gaus(0., 3e-4);
}
//__________________________________________________________________________
void SoLKalTrackFinder::SetTargetGeometry(Double_t& z, Double_t& center, Double_t& length)
{
  fTargetPlaneZ = z;
  fTargetCenter = center;
  fTargetLength = length;
}
//____________________________________________________________________________
void SoLKalTrackFinder::CalCircle(Double_t x1,Double_t y1,Double_t x2,Double_t y2,Double_t x3,
                                         Double_t y3, Double_t* R,Double_t* Xc, Double_t* Yc){
  if (y1!=y2&& (x1-x2)*(y2-y3)!=(x2-x3)*(y1-y2)){
    *Xc = ((y2-y3)*(x1*x1-x2*x2+y1*y1-y2*y2)/2.-(y1-y2)*(x2*x2-x3*x3+y2*y2-y3*y3)/2.)/((x1-x2)*(y2-y3)-(x2-x3)*(y1-y2));
    *Yc = ((x1*x1-x2*x2+y1*y1-y2*y2)/2.-(*Xc)*(x1-x2))/(y1-y2);
    *R = sqrt((*Xc-x1)*(*Xc-x1)+(*Yc-y1)*(*Yc-y1));
    
  }else if (y2!=y3&&(x2-x3)*(y3-y1)!=(x3-x1)*(y2-y3)){
    *Xc = ((y3-y1)*(x2*x2-x3*x3+y2*y2-y3*y3)/2.-(y2-y3)*(x3*x3-x1*x1+y3*y3-y1*y1)/2.)/((x2-x3)*(y3-y1)-(x3-x1)*(y2-y3));
    *Yc = ((x2*x2-x3*x3+y2*y2-y3*y3)/2.-(*Xc)*(x2-x3))/(y2-y3);
    *R = sqrt((*Xc-x2)*(*Xc-x2)+(*Yc-y2)*(*Yc-y2));
   
  }else if (y3!=y1&&(x3-x1)*(y1-y2)!=(x1-x2)*(y3-y1)){
    *Xc = ((y1-y2)*(x3*x3-x1*x1+y3*y3-y1*y1)/2.-(y3-y1)*(x1*x1-x2*x2+y1*y1-y2*y2)/2.)/((x3-x1)*(y1-y2)-(x1-x2)*(y3-y1));
    *Yc = ((x3*x3-x1*x1+y3*y3-y1*y1)/2.-(*Xc)*(x3-x1))/(y3-y1);
    *R = sqrt((*Xc-x3)*(*Xc-x3)+(*Yc-y3)*(*Yc-y3));
    
  }else{
    *Xc = 0;
    *Yc = 0;
    *R = 0;
  }
}
