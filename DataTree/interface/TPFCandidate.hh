#ifndef HIGGSANA_NTUPLER_TPFCANDIDATE_HH
#define HIGGSANA_NTUPLER_TPFCANDIDATE_HH

#include <TObject.h>

namespace higgsana
{
  class TPFCandidate : public TObject
  {
    public:
      TPFCandidate(){}
      ~TPFCandidate(){}
    
      Float_t pt, eta, phi, e;     // kinematics      
      Int_t q;                     // charge    
      Double_t dz;                 // dz to the primary vertex
      Bool_t   IsPFNoPU;
      UInt_t   pfType;              
      UInt_t   matchedObjectType;   // matched to reco ele,mu,photons
      UInt_t   matchedObjectIndex;  // index of the matched object
      //ULong_t  pfFlags;             

    ClassDef(TPFCandidate,2)
  };
}
#endif
