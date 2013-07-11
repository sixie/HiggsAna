#ifndef HIGGSANA_NTUPLER_TGenParticle_HH
#define HIGGSANA_NTUPLER_TGenParticle_HH

#include <TObject.h>

namespace mithep 
{
  class TGenParticle : public TObject
  {
    public:
      TGenParticle(){}
      ~TGenParticle(){} 
  
      Float_t pt, eta, phi, mass;     // kinematics
      Int_t pdgid;                    
      Int_t status;                   
      Int_t motherPdgID;

      ClassDef(TGenParticle,1)
        };  
}
#endif
