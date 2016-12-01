#ifndef ROOT_SoLID_GEM_Chamber
#define ROOT_SoLID_GEM_Chamber
//c++
#include <vector>
//root
#include "TClonesArray.h"
#include "TVector2.h"
//Hall A analyzer
#include "THaSubDetector.h"
//SoLID Tracking
#include "SoLIDGEMReadOut.h"
#include "SoLIDGEMTracker.h"
#include "SoLIDUtility.h"

class SoLIDGEMTracker;
class SoLIDGEMReadOut;

class SoLIDGEMChamber : public THaSubDetector{
  public:
  SoLIDGEMChamber(Int_t ichamber, const char* name, const char* description = "",
                  THaDetectorBase* parent = 0 );
  SoLIDGEMChamber() : fHits(0) {}
  virtual ~SoLIDGEMChamber();

  virtual void      Clear( Option_t* opt="" );
  virtual Int_t     Decode( const THaEvData& );
  virtual EStatus   Init( const TDatime& date );
  virtual void      Print( Option_t* opt="" ) const;
  virtual void      PrintDataBase(Int_t level) const;
  virtual Int_t     Begin( THaRunBase* r=0 );
  virtual Int_t     End( THaRunBase* r=0 );

  void              SetChamberID(Int_t i) { fChamberID = i; }
  Int_t             GetChamberID() const { return fChamberID; }
  Double_t          GetPhi() const { return fPhi; }//in rad from -pi to pi
  Double_t          GetPhiCover() const { return fPhiCover; }
  Double_t          GetPhiOffset() const { return fPhiOffset; }
  TVector2          GetReference() const { return fReference; }
  const TVector3&   GetCenterPos() const { return fChamberCenter; }
  Double_t          GetZ() const { return fTrackerZ + fDz; }
  Int_t             GetNHits() const { return fHits->GetLast()+1; }
  TSeqCollection*   GetHits() const { return fHits; }
  Double_t          GetPhiInLab() { return fPhiInLab; }

  //TSeqCollection*   GetRawHits(Int_t i) const { return fGEMReadOut[i]->GetHits(); }


  SoLIDGEMReadOut * GetReadOut(Int_t i) const {
    if (i >=0 && i< fNReadOut) { return fGEMReadOut[i]; }
    else { return 0; }
  }

  Double_t GetUOccupancy() const { return fUOccupancy; }
  Double_t GetUHitOccupancy() const { return fUHitOccupancy; }
  Double_t GetVOccupancy() const { return fVOccupancy; }
  Double_t GetVHitOccupancy() const { return fVHitOccupancy; }

  protected:
  virtual Int_t     ReadDatabase( const TDatime& date );
  virtual Int_t     ReadGeometry( FILE* file, const TDatime& date,
                                  Bool_t required = kTRUE );
  virtual Int_t     DefineVariables( EMode mode = kDefine );
  Int_t             ProcessRawHits(TSeqCollection* uhits, TSeqCollection* vhits);
  void              UVtoCylinCoor(Double_t upos, Double_t vpos, Double_t* r, Double_t* phi);
  Bool_t            Contains(Double_t *r, Double_t *phi);
  Bool_t            CheckChargeAsymmetry(Double_t qu, Double_t qv);
  void              RotateToLab(Double_t *phi);

  Int_t             fChamberID;    //sector is the same as chamber, but chamber sounds nicer IMO
  Int_t             fParentTrackerID;
  //for database
  Int_t             fNReadOut;     //number of readout plane for in this chamber
  Int_t             fDo3DAmCorr;   //whether we do charge amplitude correlation cut
  Double_t          f3DAmCorrCut;  //cut on the charge amplitude correlation
  Double_t          fRMin;         //inner radius of the chamber
  Double_t          fRMax;         //outer radius of the chamber
  Double_t          fDz;           //shift in z with respect to the parent tracker
                                   //should be 0 for PVDIS, but not always for SIDIS due to overlap
  Double_t          fPhi;         //the phi angle of this chamber with respect to the phi angle of
                                  //the tracker system it belongs to
  Double_t          fPhiCover;    //the phi coverage of this chamber with respect to the reference
  Double_t          fPhiOffset;   //phi offset in addition to fPhi of this chamber
  Double_t          fPhiInLab;    //the total phi angle with respect to the lab system
  TVector2          fReference;   //the reference point of this chamber, see comment in ReadGeometry
  Double_t          fPhiMin;
  Double_t          fPhiMax;
  Double_t          fTrackerZ;
  TClonesArray*     fHits;        //collection for a pair of raw hits, checked by amplitude matching,
                                  // whether inside active area of the GEM detector, and so on
  Double_t          fUOccupancy;
  Double_t          fVOccupancy;
  Double_t          fUHitOccupancy;
  Double_t          fVHitOccupancy;

  TVector3          fChamberCenter;

  std::vector<SoLIDGEMReadOut*> fGEMReadOut; //The readout plane that the chamber has, usually 2
  ClassDef(SoLIDGEMChamber,0)
};
#endif
