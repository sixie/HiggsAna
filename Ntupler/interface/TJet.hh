#ifndef HIGGSANA_NTUPLER_TJET_HH
#define HIGGSANA_NTUPLER_TJET_HH

#include <TObject.h>

namespace mithep 
{
  class TJet : public TObject
  {
    public:
      TJet(){}
      ~TJet(){}

      Float_t pt, eta, phi, mass;  // kinematics
      Double_t TrackCountingHighEffBJetTagsDisc; //btag discriminator
      ULong_t hltMatchBits;   // bits from matching with HLT primitives
      Bool_t PassBetaVertexAssociationCut;
      UInt_t NConstituents;
      Float_t NeutralHadronFraction;
      Float_t NeutralEMFraction;
      Float_t ChargedHadronFraction;
      Float_t ChargedEMFraction;
      UInt_t ChargedMultiplicity;      

      Float_t TrackCountingHighPurBJetTagsDisc;
      Float_t GhostTrackBJetTagsDisc;
      Float_t SoftElectronByPtBJetTagsDisc;
      Float_t SoftElectronByIP3dBJetTagsDisc;
      Float_t SoftMuonByPtBJetTagsDisc;
      Float_t SoftMuonByIP3dBJetTagsDisc;
      Float_t SoftMuonBJetTagsDisc;
      Float_t SimpleSecondaryVertexHighPurBJetTagsDisc;
      Float_t SimpleSecondaryVertexHighEffBJetTagsDisc;
      Float_t SimpleSecondaryVertexBJetTagsDisc;
      Float_t CombinedSecondaryVertexBJetTagsDisc;
      Float_t CombinedSecondaryVertexMVABJetTagsDisc;
      Float_t JetProbabilityBJetTagsDisc;
      Float_t JetBProbabilityBJetTagsDisc;

      Float_t JetArea;

      Float_t DzAvg;
      Float_t DzPtSqrWeightedAvg;
      
      Float_t rawPt;

    ClassDef(TJet,1)
  };
}
#endif
