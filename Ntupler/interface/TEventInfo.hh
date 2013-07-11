#ifndef HIGGSANA_NTUPLER_TEVENTINFO_HH
#define HIGGSANA_NTUPLER_TEVENTINFO_HH

#include <TObject.h>

namespace mithep 
{
  class TEventInfo : public TObject
  {
    public:
      TEventInfo(){}
      ~TEventInfo(){}

      UInt_t runNum; 			    // run number in data
      UInt_t evtNum; 			    // event number in data
      UInt_t lumiSec;			    // lumi section      
      Float_t eventweight;
      UInt_t nPUEvents;                     // number of pileup events.
      UInt_t nPUMinusOne;
      UInt_t nPUPlusOne;
      Float_t RhoKt6PFJetsForIso25;    
      Float_t RhoKt6PFJetsCentralNeutral;
      Float_t RhoKt6PFJets;

      ULong_t triggerBits;		    // HLT trigger bits 
      UInt_t l1triggerBits;		    // L1 trigger bits 

      Bool_t hasGoodPV;                     // event has a good PV?
      UInt_t nPV0;                          // number of reconstructed primary vertices in event
      Float_t pvx, pvy, pvz;		    // primary vertex with the most associated tracks 
      Float_t bsx, bsy, bsz;		    // beamspot					  

      Float_t tcMEx, tcMEy, tcSumET;	    // track-corrected MET
      Float_t pfMEx, pfMEy, pfSumET;	    // particle flow MET
      Float_t pfTrackMEx, pfTrackMEy, pfTrackSumET; // particle flow track MET
      Float_t pfNeutralMEx, pfNeutralMEy, pfNeutralSumET; // particle flow neutral MET
      Float_t pfNeutralNoFwdMEx, pfNeutralNoFwdMEy, pfNeutralNoFwdSumET; // particle flow neutral MET, up to eta 2.5


    ClassDef(TEventInfo,1)
  };
}
#endif
