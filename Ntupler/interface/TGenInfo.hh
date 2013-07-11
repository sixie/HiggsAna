#ifndef HIGGSANA_NTUPLER_TGENINFO_HH
#define HIGGSANA_NTUPLER_TGENINFO_HH

#include <TObject.h>

namespace mithep 
{
  // Generator level info data object
  class TGenInfo : public TObject
  {
    public:
      TGenInfo(){}
      ~TGenInfo(){}
      
      UInt_t nGenPart;		                      // number of generated particles (status 1) in event
      UInt_t nGenCh;		                      // number of generated charged particles (status 1) in event
      UInt_t npho;  		                      // number of FSR photons
      Int_t id_1, id_2;		                      // parton PDF ID
      Float_t x_1, x_2;		                      // parton momentum fraction
      Float_t pdf_1, pdf_2; 	                      // parton PDF value
      Float_t weight;		                      // event weight
      Float_t scalePdf;		                      // Q-scale in PDF evaluation
      Float_t scale;		                      // event energy scale
      Float_t vmass_1, vpt_1, vy_1, veta_1, vphi_1;   // vector boson 1 kinematics
      Float_t vmass_2, vpt_2, vy_2, veta_2, vphi_2;   // vector boson 2 kinematics
      Float_t mass, pt, y, phi;     	              // dilepton/lepton-neutrino kinematics
      Float_t pt_1, eta_1, phi_1;                     // lepton kinematics
      Float_t pt_2, eta_2, phi_2;                     // anti-lepton/neutrino kinematics
      Float_t phopt, phoeta, phophi;                  // leading photon kinematics
      Float_t decx, decy, decz;	                      // boson decay vertex	  
      Float_t met, metPhi, ptBosonSystem;
      Float_t nGenJets;                               //
      Float_t jetpt_1, jetpt_2,jetpt_3,jetpt_4 ;      //

    ClassDef(TGenInfo,1)
  };
}
#endif
