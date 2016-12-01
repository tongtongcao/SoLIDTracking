#ifndef ROOT_SoLID_GEM_ReadOut
#define ROOT_SoLID_GEM_ReadOut
//c++
#include <vector>
//ROOT
#include "TClonesArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TH1.h"
//Hall A analyzer
#include "THaAnalysisObject.h"
#include "THaSubDetector.h"
#include "SoLIDGEMHit.h"
#define ALL(c) (c).begin(), (c).end()

using namespace std;
namespace Podd {
  class MCHitInfo;
}

struct StripData_t {
  Float_t adcraw;
  Float_t adc;
  Float_t time;
  Bool_t  pass;
  StripData_t() {}
  StripData_t( Float_t _raw, Float_t _adc, Float_t _time, Bool_t _pass )
    : adcraw(_raw), adc(_adc), time(_time), pass(_pass) {}
};

class SoLIDGEMReadOut : public THaSubDetector{
  public:
  typedef std::vector<Int_t>    Vint_t;
  typedef std::vector<Float_t>  Vflt_t;
  typedef std::vector<Bool_t>   Vbool_t;
  
  SoLIDGEMReadOut(Int_t ireadout, const char* name, const char* description = "",
           THaDetectorBase* parent = 0 );
  SoLIDGEMReadOut() : fADCraw(0), fADC(0), fHitTime(0), fADCcor(0), fGoodHit(0),
                      fHits(0) {}
  virtual ~SoLIDGEMReadOut();
  
  virtual void      Clear( Option_t* opt="" );
  virtual Int_t     Decode( const THaEvData& );
  virtual EStatus   Init( const TDatime& date );
  virtual void      Print( Option_t* opt="" ) const;
  virtual void      PrintDataBase(Int_t level) const;
  virtual Int_t     Begin( THaRunBase* r=0 );
  virtual Int_t     End( THaRunBase* r=0 );
  
  void              SetReadOutID(Int_t i) { fReadOutID = i; }
  Int_t             GetReadOutID() const { return fReadOutID; }
  Double_t          GetStripAngle() const { return fStripAngle; }
  Double_t          GetSinStripAngle() const { return fSinStripAngle; }
  Double_t          GetCosStripAngle() const { return fCosStripAngle; }
  Double_t          GetStart() const { return fStartPos+fCoorOffset; }
  TSeqCollection*   GetHits()        const { return fHits; }
  Int_t             GetNhits()       const { return fHits->GetLast()+1; }
  Double_t          GetResolution()  const { return fResolution; }
  Int_t             GetNSigStrips()  const { return fSigStrips.size(); }
  Double_t          GetPitch() const { return fStripPitch; }
  Double_t          GetOccupancy() const { return fOccupancy; }
  Double_t          GetHitOccupancy() const { return fHitOcc; }
  Bool_t            IsNoisyEvent() {
    if ( (fHits->GetLast()+1) >= fMaxHits) return kTRUE;
    else return kFALSE;
  }
  
  protected:
  
  virtual Int_t     ReadDatabase( const TDatime& date );
  virtual Int_t     DefineVariables( EMode mode = kDefine );
  StripData_t       ChargeDep( const vector<Float_t>& amp );
  StripData_t       ChipChargeDep( const vector<Float_t>& amp);
  Int_t             MapChannel( Int_t idx ) const;
  void              AddStrip( Int_t istrip );
  void              UpdateOffset();
  void              KillCrossTalk();
  
  Int_t             fReadOutID;       //0 means u and 1 means v   
  //for database
  Double_t          fStripAngle;      //strip angle of the readout strip with respect
                                      // to the symmetric axis of the chamber
  Double_t          fSinStripAngle;   //
  Double_t          fCosStripAngle;   //
  Double_t          fResolution;      // assumed resolution of the reconstructed hit
  Int_t             fMaxClusterSize;  // maximum cluster size we allow before splitting
  Int_t             fADCMin;          // ADC cut
  Double_t          fSplitFrac;       // used in clustering for splitting a cluster
  Int_t             fMaxHits;         // Maximum number of hits we allow in reconstruction
  Int_t             fMaxTimeSample;   // The maximum number of time sample from APV25
  Double_t          fADCSigma;        //
  Int_t             fDoNoise;         //  
  Int_t             fCheckPulseShape; //
  Int_t             fDoHisto;         //
  Double_t          fStripPitch;      // Pitch of the readout strips
  Double_t          fDz;              //
  Double_t          fDPhi;            //
  Double_t          fStartPos;        // Starting position of the readout strips
  Int_t             fNStrip;          // total number of strips on this readout
  Vflt_t            fPed;             // pedestal for each channal
  Int_t             fDeconMode;       // 1 if the APV25 chip is in the deconvolution mode
  Int_t             fKillCrossTalk;   // 1 to apply cross talk signal elimination algorithm
  Double_t          fCrossTalkThres;  // Threshold for eliminating cross talk signal
  Int_t             fCrossStripApart; // the number of strip that the induced signal is away
                                      // from the main signal
  
  // Event data, hits etc.
  Float_t*          fADCraw;          // [fNelem] Integral of raw ADC samples
  Float_t*          fADC;             // [fNelem] Integral of deconvoluted ADC samples
  Float_t*          fHitTime;         // [fNelem] Leading-edge time of deconv signal (ns)
  Float_t*          fADCcor;          // [fNelem] fADC corrected for pedestal & noise
  Byte_t*           fGoodHit;         // [fNelem] Strip data passed pulse shape test
  Double_t          fDNoise;          // Event-by-event noise (avg below fMinAmpl)
  Vint_t            fSigStrips;       // Ordered strip numbers with signal (adccor > minampl)
  Vbool_t           fStripsSeen;      // Flags for duplicate strip number detection

  UInt_t            fNRawStrips;      // Statistics: strips with any data
  UInt_t            fNHitStrips;      // Statistics: strips > 0
  Double_t          fHitOcc;          // Statistics: hit occupancy fNhitStrips/fNelem
  Double_t          fOccupancy;       // Statistics: occupancy GetNsigStrips/fNelem
  Double_t          fCoorOffset;      // The offset needed to transport the hit coordinate
                                      // from chamber frame to lab frame
  TClonesArray*     fHits;            // Cluster data (groups of hits)


  // Hardware channel mapping
  enum EChanMapType { kOneToOne, kReverse, kGassiplexAdapter1,
                      kGassiplexAdapter2, kTable };

  EChanMapType      fMapType;         // Type of hardware channel mapping to use
  Vint_t            fChanMap;         // [fNelem] Optional hardware channel mapping
  TH1*              fADCMap;          // Histogram of strip numbers weighted by ADC
  
  
#ifdef MCDATA
  // Only set when analyzing simulation data
  Vint_t            fMCHitList;       // Elements in fMCHitInfo with data
  Podd::MCHitInfo*  fMCHitInfo;       // [fNStrip] MC truth data for fHits
#endif
  // Only for TESTCODE, but kept for binary compatibility
  TH1*          fHitMap;              // Histogram of active sensor numbers
    
  ClassDef(SoLIDGEMReadOut,0)  // One Tracker plane coordinate direction
};
#endif
