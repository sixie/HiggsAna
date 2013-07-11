#ifndef HIGGSANA_NTUPLER_TPHOTON_HH
#define HIGGSANA_NTUPLER_TPHOTON_HH

#include <TObject.h>

namespace mithep 
{
  class TPhoton : public TObject
  {
    public:
      TPhoton(){}
      ~TPhoton(){}

      Float_t et, eta, phi; 	              // kinematics
      Float_t scEt, scEta, scPhi;             // supercluster
      Float_t trkIso03Hollow, trkIso03Solid;  // track isolation
      Float_t emIso03;                        // ECAL-based isolation
      Float_t hadIso03;                       // HCAL-based isolation
      Float_t HoverE;		              // H/E
      Float_t R9;		              // ratio of energies in 3x3 to SC
      Float_t sigiEtaiEta;                    // eta-width of shower in number of crystals
      UInt_t  hltMatchBits;  	              // bits from matching with HLT primitives
      UInt_t  scID;                           // supercluster ID (for matching to electron superclusters)
      Bool_t  hasPixelSeed;                   // supercluster has pixel seed?
    
    ClassDef(TPhoton,1)
  };  
}
#endif
