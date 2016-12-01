//c++
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
//SoLIDTracking
#include "SoLKalTrackSite.h"
ClassImp(SoLKalTrackSite)

SoLKalTrackSite::SoLKalTrackSite(Int_t m, Int_t p, Double_t chi2)
:TObjArray(2), fCurStatePtr(0), fM(m,1), fV(m,m),
 fH(m,p), fHt(p,m), fResVec(m,1), fR(m,m), fDeltaChi2(0.), fMaxDeltaChi2(chi2)
{
}
//_________________________________________________________________________
SoLKalTrackSite::SoLKalTrackSite(SoLIDGEMHit *ht, Int_t m, Int_t  p, Double_t chi2)
 : TObjArray(2), fCurStatePtr(0), fM(m,1), fV(m,m),
 fH(m,p), fHt(p,m), fResVec(m,1), fR(m,m), fDeltaChi2(0.), fMaxDeltaChi2(chi2),
 fGEMHit(ht) 
{
  fM(kIdxX0, 0) = ht->GetX(); 
  fM(kIdxY0, 0) = ht->GetY();
  Double_t phi = atan2(ht->GetY(), ht->GetX());
  Double_t dr = 5.e-4;
  Double_t drphi = 5.4e-5;
  
  Double_t dx = sqrt( pow( cos(phi)*dr, 2) + pow( sin(phi)*drphi, 2) );
  Double_t dy = sqrt( pow( sin(phi)*dr, 2) + pow( cos(phi)*drphi, 2) );
  fV(kIdxX0, kIdxX0) = pow(dx, 2);
  fV(kIdxY0, kIdxY0) = pow(dy, 2);
  fZ0 = ht->GetZ();
}
//_________________________________________________________________________
SoLKalTrackSite::~SoLKalTrackSite()
{
  SetOwner(kTRUE);
  Delete();
}
//_________________________________________________________________________
void SoLKalTrackSite::SetHitResolution(Double_t ex, Double_t ey)
{
  fV(kIdxX0, kIdxX0) = TMath::Power(ex, 2);
  fV(kIdxY0, kIdxY0) = TMath::Power(ey, 2);
}
//__________________________________________________________________________
void SoLKalTrackSite::SetMeasurement(Double_t x, Double_t y)
{
  fM(kIdxX0, 0) = x;
  fM(kIdxY0, 0) = y;
}
//__________________________________________________________________________
inline Int_t SoLKalTrackSite::CalcXexp(const SoLKalTrackState &a, TVector3 &xx) const
{
  xx.SetXYZ(a(0,0), a(1,0), fZ0);
  return 1;
}
//_________________________________________________________________________
Int_t SoLKalTrackSite::CalcExpectedMeasVec(const SoLKalTrackState &a, SoLKalMatrix &h)
{
   TVector3 xxv;
   if (!CalcXexp(a,xxv)) return 0;  // no hit
   h(0, 0) = xxv.X();
   h(1, 0) = xxv.Y();
   return 1;
}
//_________________________________________________________________________
Int_t SoLKalTrackSite::CalcMeasVecDerivative(const SoLKalTrackState&/* &a*/, SoLKalMatrix &H)
{
  H.Zero();
  H(0, 0) = 1.;
  H(1, 1) = 1.;
  return 1;
}
//__________________________________________________________________________
Bool_t SoLKalTrackSite::IsAccepted()
{
   if (fDeltaChi2 < fMaxDeltaChi2){
     return true;
   }else{
     return false;
   }
}
//___________________________________________________________________________
Bool_t SoLKalTrackSite::Filter()
{
   // prea and preC should be preset by SoLKalTrackState::Propagate()
   SoLKalTrackState &prea = GetState(SoLKalTrackSite::kPredicted);
   SoLKalMatrix h = fM;
   
   if (!CalcExpectedMeasVec(prea,h)) return kFALSE;
   
   SoLKalMatrix pull  = fM - h;
   
   SoLKalMatrix preC  = GetState(SoLKalTrackSite::kPredicted).GetCovMat();
   
   if (!CalcMeasVecDerivative(prea,fH)) return kFALSE;
   fHt = SoLKalMatrix(SoLKalMatrix::kTransposed, fH);

   // Calculate covariance matrix of residual

   SoLKalMatrix preR    = fV + fH * preC * fHt;
   
   SoLKalMatrix preRinv = SoLKalMatrix(SoLKalMatrix::kInverted, preR);
   if (preR.Determinant() == 0) {cout<<"1"<<endl; preR.Print(); fV.Print(); fH.Print(); preC.Print(); }
   // Calculate kalman gain matrix

   SoLKalMatrix preCinv = SoLKalMatrix(SoLKalMatrix::kInverted, preC);
   if (preC.Determinant() == 0) {cout<<"2"<<endl; preC.Print();}
   SoLKalMatrix G       = SoLKalMatrix(SoLKalMatrix::kInverted, fV);

   // Calculate filtered state vector

   SoLKalMatrix curCinv = preCinv + fHt * G * fH;
   SoLKalMatrix curC    = SoLKalMatrix(SoLKalMatrix::kInverted, curCinv);
   if (curC.Determinant() == 0) {cout<<"3"<<endl; curCinv.Print();}
   SoLKalMatrix K       = curC * fHt * G;

   SoLKalMatrix Kpull  = K * pull;
   SoLKalMatrix Kpullt = SoLKalMatrix(SoLKalMatrix::kTransposed,Kpull);
   SoLKalMatrix av     = prea + Kpull;
   SoLKalTrackState &a     = CreateState(av,curC,SoLKalTrackSite::kFiltered);
   SoLKalTrackState *aPtr  = &a;

   Add(aPtr);
   SetOwner();

   // Calculate chi2 increment

   fR      = fV - fH * curC *fHt;
   if (!CalcExpectedMeasVec(a,h)) return kFALSE;
   fResVec = fM - h;
   SoLKalMatrix curResVect = SoLKalMatrix(SoLKalMatrix::kTransposed, fResVec);
   fDeltaChi2 = (curResVect * G * fResVec + Kpullt * preCinv * Kpull)(0,0);

   if (IsAccepted()) return kTRUE;
   else              return kFALSE;
}
//______________________________________________________________________________
void SoLKalTrackSite::Smooth(SoLKalTrackSite &pre)
{
   if (&GetState(SoLKalTrackSite::kSmoothed)) return;

   SoLKalTrackState &cura  = GetState(SoLKalTrackSite::kFiltered);
   SoLKalTrackState &prea  = pre.GetState(SoLKalTrackSite::kPredicted);
   SoLKalTrackState &sprea = pre.GetState(SoLKalTrackSite::kSmoothed);

   SoLKalMatrix curC    = cura.GetCovMat();
   SoLKalMatrix curFt   = cura.GetPropMat("T");
   SoLKalMatrix preC    = prea.GetCovMat();
   SoLKalMatrix spreC   = sprea.GetCovMat();
   SoLKalMatrix preCinv = SoLKalMatrix(SoLKalMatrix::kInverted, preC);
   SoLKalMatrix curA    = curC * curFt * preCinv;
   SoLKalMatrix curAt   = SoLKalMatrix(SoLKalMatrix::kTransposed, curA);
   SoLKalMatrix scurC   = curC + curA * (spreC - preC) * curAt;

   SoLKalMatrix sv = cura + curA * (sprea - prea);
   Add(&CreateState(sv,scurC,SoLKalTrackSite::kSmoothed));
   SetOwner();

   // Update residual vector

   fR       = fV - fH * scurC *fHt;
   fResVec -= fH * (sv - cura);
   SoLKalMatrix curResVect = SoLKalMatrix(SoLKalMatrix::kTransposed, fResVec);
   SoLKalMatrix curRinv    = SoLKalMatrix(SoLKalMatrix::kInverted, fR);
   fDeltaChi2 = (curResVect * curRinv * fResVec)(0,0);
}
//______________________________________________________________________________
void SoLKalTrackSite::InvFilter()
{
   if (&GetState(SoLKalTrackSite::kInvFiltered)) return;

   SoLKalTrackState &sa = GetState(SoLKalTrackSite::kSmoothed);
   SoLKalMatrix pull = fResVec;

   SoLKalMatrix sC     = sa.GetCovMat();
   SoLKalMatrix sR     = fH * sC * fHt - fV;
   SoLKalMatrix sRinv  = SoLKalMatrix(SoLKalMatrix::kInverted, sR);
   SoLKalMatrix Kstar  = sC * fHt * sRinv;
   SoLKalMatrix svstar = sa + Kstar * pull;
   SoLKalMatrix sCinv  = SoLKalMatrix(SoLKalMatrix::kInverted, sC);
   SoLKalMatrix G      = SoLKalMatrix(SoLKalMatrix::kInverted, fV);
   SoLKalMatrix Cstar  = SoLKalMatrix(SoLKalMatrix::kInverted, sCinv + fHt * G * fH);
   Add(&CreateState(svstar,Cstar,SoLKalTrackSite::kInvFiltered));
   SetOwner();
}
//_______________________________________________________________________________
SoLKalMatrix SoLKalTrackSite::GetResVec (SoLKalTrackSite::EStType t)
{
   using namespace std;
   SoLKalTrackState &a  = GetState(t);
   SoLKalTrackState &sa = (&GetState(SoLKalTrackSite::kSmoothed) != 0
                    ? GetState(SoLKalTrackSite::kSmoothed)
                    : GetState(SoLKalTrackSite::kFiltered));
   if (!&a || !&sa) {
     cerr << ":::::: ERROR in SoLKalTrackSite::GetResVec(EStType) " << endl
          << " Invalid states requested"                      << endl
          << " &a = " << &a << " &sa = " << &sa               << endl
          << " Abort!"                                        << endl;
     ::abort();
   }
   if (&a == &sa) {
      return fResVec;
   } else {
      return (fResVec - fH * (a - sa));
   }
}
//__________________________________________________________________________
void SoLKalTrackSite::Add(TObject *obj)
{
   TObjArray::Add(obj);
   fCurStatePtr = static_cast<SoLKalTrackState*>(obj);
   fCurStatePtr->SetSitePtr(this);
}
//___________________________________________________________________________
SoLKalTrackState & SoLKalTrackSite::GetState(SoLKalTrackSite::EStType t)
{
   SoLKalTrackState *ap = 0;
   if (t >= 0 && t < GetEntries()) {
      ap = static_cast<SoLKalTrackState*>(UncheckedAt(t));
   }
   
   return *ap;
}       
//___________________________________________________________________________
SoLKalTrackState & SoLKalTrackSite::CreateState(const SoLKalMatrix &sv, Int_t type)
{
   SetOwner();
   return *(new SoLKalTrackState(sv,*this,type));
}
//____________________________________________________________________________
SoLKalTrackState & SoLKalTrackSite::CreateState(const SoLKalMatrix &sv,
                                                const SoLKalMatrix &c,
                                                Int_t       type)
{
   SetOwner();
   return *(new SoLKalTrackState(sv,c,*this,type));
}
//_____________________________________________________________________________
SoLIDGEMHit* SoLKalTrackSite::GetPredInfoHit()
{
  SoLKalTrackState &a  = GetState(kPredicted);
  
  fGEMHit->SetPredictHit( a(kIdxX0, 0), a(kIdxY0, 0), 
  sqrt(a.GetCovMat()(kIdxX0, kIdxX0)), 
  sqrt(a.GetCovMat()(kIdxY0, kIdxY0)) );
  a = (&GetState(SoLKalTrackSite::kSmoothed) != 0
      ? GetState(SoLKalTrackSite::kSmoothed)
      : GetState(SoLKalTrackSite::kFiltered));

  //a = GetState(SoLKalTrackSite::kFiltered);

  Double_t momentum = TMath::Abs(1./a(kIdxQP, 0)); 
  TVector3 tmp;
  tmp.SetZ( 1./(TMath::Sqrt(a(kIdxTX, 0)*a(kIdxTX, 0) + a(kIdxTY, 0)*a(kIdxTY, 0) + 1. )) );
  tmp.SetX(a(kIdxTX, 0) * tmp.Z());
  tmp.SetY(a(kIdxTY, 0) * tmp.Z());
  tmp = tmp.Unit();

  fGEMHit->SetMomentum(momentum*tmp.X(), momentum*tmp.Y(), momentum*tmp.Z());
 
  return fGEMHit;
}





























