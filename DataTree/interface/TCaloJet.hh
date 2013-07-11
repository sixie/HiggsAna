#ifndef HIGGSANA_NTUPLER_TCALOJET_HH
#define HIGGSANA_NTUPLER_TCALOJET_HH

#include <TObject.h>

namespace higgsana
{
  class TCaloJet : public TObject
  {
    public:
      TCaloJet(){}
      ~TCaloJet(){}

      Float_t pt, eta, phi, mass;  // kinematics
      Double_t TrackCountingHighEffBJetTagsDisc; //btag discriminator
      ULong_t hltMatchBits;   // bits from matching with HLT primitives
      Bool_t PassBetaVertexAssociationCut;
      UInt_t NConstituents;
      Float_t EmEnergyInEB;
      Float_t EmEnergyInEE;
      Float_t EmEnergyInHF;
      Float_t HadEnergyInHO;
      Float_t HadEnergyInHB;
      Float_t HadEnergyInHF;
      Float_t HadEnergyInHE;
      Float_t EnergyFractionH;
      Float_t EnergyFractionEm;

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
      Float_t rawPt;

    ClassDef(TCaloJet,1)
  };
}
#endif
