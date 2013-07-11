#ifndef HIGGSANA_NTUPLER_TELECTRON_HH
#define HIGGSANA_NTUPLER_TELECTRON_HH

#include <TObject.h>

namespace mithep
{
  class TElectron : public TObject
  {
    public:
      TElectron(){}
      ~TElectron(){}
    
      Float_t pt, eta, phi, p;     // kinematics, p is bestTrack momentum
      Float_t pfPt, pfEta, pfPhi;  // kinematics      
      Float_t scEt, scEta, scPhi;  // supercluster
      Int_t q;                     // charge
      Int_t  isMCReal;        
      ULong_t hltMatchBits;         // bits for matching with HLT primitives
      UInt_t l1TriggerMatchBits;   // bits for matching with L1 Seeds
      Bool_t isEcalDriven;         // is ECAL seeded electron?
      Bool_t isTrackerDriven;      // is TrackerDriven electron
      Bool_t  isEB;                // is in barrel

      Float_t d0, d0Err, dz;         // impact parameter
      Float_t ip3d, ip3dSig;         //IP3D PV
      Int_t   isConv;                // is conversion? (vertexing method)
      Float_t nExpHitsInner;       // number of hits expected before first hit
      Float_t partnerDeltaCot;     // cot(theta) difference with conversion partner track       
      Float_t partnerDist;         // distance in x-y plane to nearest conversion partner track
      Float_t partnerRadius;       // radius of helix intersection with conversion partner track

      Int_t   nBrem;               // Number of Brems
      Float_t fBrem, EOverP;       // fBrem, EOverP
      Float_t pIn;                 // mode of the gsf track at vertex
      Float_t ESeedClusterOverPIn; // ESeedClusterOverPout     
      Float_t ESeedClusterOverPout;// ESeedClusterOverPout    
      Float_t EEleClusterOverPout; 
      Float_t EcalEnergy;
      Float_t EcalEnergyError;
      Float_t deltaEtaIn;          // eta difference between track (at vertex) and SC
      Float_t deltaPhiIn;          // phi difference between track (at vertex) and SC
      Float_t dEtaCalo;
      Float_t dPhiCalo;
      Float_t sigiEtaiEta;         // eta-width of shower in number of crystals
      Float_t sigiPhiiPhi;         // phi-width of shower in number of crystals
      Float_t CovIEtaIPhi;
      Float_t SCEtaWidth;
      Float_t SCPhiWidth;
      Float_t R9;
      Float_t PreShowerOverRaw;

      Float_t HoverE;              // H / E
      Float_t GsfTrackChi2OverNdof;
      Float_t KFTrackChi2OverNdof;
      Float_t KFTrackNHits;
      Float_t KFTrackNLayersWithMeasurement;
      Float_t SeedE1x5OverE5x5;

      Float_t likelihood;          // likelihood
      Float_t mva;                 // mva

      Float_t trkIso03;            // track isolation
      Float_t emIso03;             // ECAL-based isolation
      Float_t hadIso03;            // HCAL-based isolation
      Float_t trkIso04;            // track isolation
      Float_t emIso04;             // ECAL-based isolation
      Float_t hadIso04;            // HCAL-based isolation

    ClassDef(TElectron,1)
  };
}
#endif
