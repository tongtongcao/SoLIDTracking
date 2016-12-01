#ifndef ROOT_SoLID_GEM_Tracker
#define ROOT_SoLID_GEM_Tracker
//c++
#include <vector>
//root
#include "TClonesArray.h"
//Hall A analyzer
#include "THaSubDetector.h"
//SoLID Tracking
#include "SoLIDGEMChamber.h"

class SoLIDGEMChamber;

class SoLIDGEMTracker : public THaSubDetector{
  public:
  SoLIDGEMTracker(Int_t iGEM, const char* name, const char* description = "",
           THaDetectorBase* parent = 0 );
  SoLIDGEMTracker(): fHits(0) {}
  virtual ~SoLIDGEMTracker();
  
  virtual void      Clear( Option_t* opt="" );
  virtual Int_t     Decode( const THaEvData& );
  virtual EStatus   Init( const TDatime& date );
  virtual void      Print( Option_t* opt="" ) const;
  virtual void      PrintDataBase(Int_t level) const;
  virtual Int_t     Begin( THaRunBase* r=0 );
  virtual Int_t     End( THaRunBase* r=0 );
  void              SetTrackerID(Int_t i) { fTrackerID = i; }
  Int_t             GetTrackerID() const { return fTrackerID; }
  Int_t             GetNChamber() const { return fNChamber; }
  Double_t&         GetZ()  { return fTrackerZ; }
  Int_t             CombineChamberHits();
  SoLIDGEMChamber * GetChamber(Int_t i) const {
    if (i >= 0 && i < fNChamber){ return fGEMChamber[i]; }
    else { return 0; }
  }
  Int_t             GetNHits() const { return fNHits; }
  
  Double_t GetUMeanOccu();
  Double_t GetUMeanHitOccu();
  Double_t GetVMeanOccu();
  Double_t GetVMeanHitOccu();
  
  Int_t             KillCrossTalk() const { return fKillCrossTalk; }
  Double_t          GetCrossTalkThres() const { return fCrossTalkThres; }
  Int_t             GetCrossStripApart() const { return fCrossStripApart; }
  
  protected:
  virtual Int_t     ReadDatabase( const TDatime& date );
  virtual Int_t     DefineVariables( EMode mode = kDefine );
  //for database
  Int_t             fNChamber;     //number of GEM sectors/chambers this tracker has
  Double_t          fTrackerZ;     //the z position of a tracker, indivial chamber may
                                   // has a offset with this z, like the case of SIDIS
  Int_t             fDoCombineHits;//whether combine all the hits from Chamber for ouput
  Int_t             fKillCrossTalk;
  Int_t             fCrossStripApart;
  Double_t          fCrossTalkThres;
  
  Int_t             fTrackerID;   //GEM ID, 0~5 for SIDIS, 0~4 for PVDIS
  TClonesArray*     fHits;   
  std::vector<SoLIDGEMChamber*> fGEMChamber;
  Int_t             fNHits;
  
  ClassDef(SoLIDGEMTracker,0)  // One Tracker plane coordinate direction
};
#endif

