//ROOT
#include "TMath.h"
//SoLIDTracking
#include "SoLIDGEMHit.h"

ClassImp(Hit)
ClassImp(SoLIDRawHit)
ClassImp(SoLIDGEMHit)
#ifdef MCDATA
ClassImp(SoLIDMCRawHit)
ClassImp(SoLIDMCGEMHit)
#endif
//-------------------------------  Hit  -----------------------------------//
//___________________________________________________________________________
 inline
 Int_t Hit::Compare( const TObject* obj ) const
 {
   // Used for sorting hits in a TSeqCollection (e.g. TList, TClonesArray).
   // A hit is "less than" another hit if its position is smaller.
   // Returns -1 if this is smaller than rhs, 0 if equal, +1 if greater.
   assert( dynamic_cast<const Hit*>(obj) );
   const Hit* rhs = static_cast<const Hit*>(obj);
   assert( fReadOut == rhs->fReadOut );
   if( fPos  < rhs->fPos )  return -1;
   if( fPos  > rhs->fPos )  return  1;
   return 0;
}

//___________________________________________________________________________
inline
Int_t Hit::Compare( const Hit* rhs, Double_t maxdist ) const
{
  // Determine if two hits are within maxdist of each other.
  // Returns -1 if this<rhs, 0 if overlap, +1 if this>rhs.
  if( fPos+maxdist < rhs->fPos )
    // this hit is "smaller than" (to the left of) rhs
    return -1;
  if( rhs->fPos+maxdist < fPos )
    // this hit is "larger than" (to the right of) rhs
    return +1;
  // The hits overlap within the maxdist tolerance
  return 0;
}
//____________________________________________________________________________
void Hit::Print( Option_t* opt ) const
{
  // Print hit info
  
     cout << "Hit: plane=" << (fReadOut ? fReadOut->GetName() : "??")
            << " pos="       << GetPos()<<endl;
                 //  << " z="         << GetZ()
                   //       << " res="       << GetResolution();
                            if( *opt != 'C' )
                                cout << endl;
}
//____________________________________________________________________________
//---------------------------------------------------------------------------//

//----------------------------- SoLIDRawHit ---------------------------------//
//_____________________________________________________________________________
void SoLIDRawHit::Print( Option_t* opt ) const
{
  // Print hit info

  Hit::Print("C");
  cout << " size="  << GetSize()
       << " type="  << GetType();
  if( *opt != 'C' )
    cout << endl;
}
//_____________________________________________________________________________
//---------------------------------------------------------------------------//

#ifdef MCDATA
//----------------------------- SoLIDMCRawHit -------------------------------//
//_____________________________________________________________________________
void SoLIDMCRawHit::Print( Option_t* ) const
{
  // Print hit info

  SoLIDRawHit::Print("C");
  MCPrint();
}
//______________________________________________________________________________
//-----------------------------------------------------------------------------//
#endif
//---------------------------- SoLIDGEMHit ------------------------------------//
//______________________________________________________________________________
SoLIDGEMHit::SoLIDGEMHit(Int_t chamberID, Int_t trackerID, 
Double_t r, Double_t phi, Double_t z, Hit* uhit, Hit* vhit)
:fIsUsed(kFALSE), fChamberID(chamberID), fTrackerID(trackerID), fR(r), fPhi(phi), fZ(z), 
fUHit(uhit), fVHit(vhit)
{
  fX = r*TMath::Cos(phi);
  fY = r*TMath::Sin(phi);
}
//______________________________________________________________________________
void SoLIDGEMHit::SetPredictHit(Double_t x, Double_t y, Double_t ex, Double_t ey)
{
  fPredictX = x;
  fPredictY = y;
  fPredicteX = ex;
  fPredicteY = ey;
}
//_______________________________________________________________________________
void SoLIDGEMHit::SetMomentum(Double_t px, Double_t py, Double_t pz)
{
  fPX = px;
  fPY = py;
  fPZ = pz;
}
//_______________________________________________________________________________
inline Int_t SoLIDGEMHit::Compare( const TObject* obj ) const
 {
   // Used for sorting hits in a TSeqCollection (e.g. TList, TClonesArray).
   // A hit is "less than" another hit if its radius (measure from z) is smaller
   // it might turn out that sorting in phi is better
   
   assert( dynamic_cast<const SoLIDGEMHit*>(obj) );
   const SoLIDGEMHit* rhs = static_cast<const SoLIDGEMHit*>(obj);
   if( fR  < rhs->GetR() )  return -1;
   if( fR  > rhs->GetR() )  return  1;
   return 0;
}
//_______________________________________________________________________________
//------------------------------------------------------------------------------//

#ifdef MCDATA
//---------------------------- SoLIDGEMHit ------------------------------------//
//______________________________________________________________________________
Int_t SoLIDMCGEMHit::IsSignalHit()
{
  Int_t mctracku = dynamic_cast<SoLIDMCRawHit*>(fUHit)->fMCTrack;
  Int_t mctrackv = dynamic_cast<SoLIDMCRawHit*>(fVHit)->fMCTrack;
  if (mctracku > 0 && mctracku == mctrackv) return mctracku;
  else return 0;
}
//_______________________________________________________________________________
#endif

























