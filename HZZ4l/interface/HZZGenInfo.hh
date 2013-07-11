#ifndef HIGGSANA_HZZ4L_GENINFO
#define HIGGSANA_HZZ4L_GENINFO

#include <TObject.h>

// #define VARLIST_GENINFO "weight/F:genQid/I:pid_1/I:pid_2/I:x_1/F:x_2/F:id_1/I:id_2/I:vmass_a/F:vpt_a/F:veta_a/F:vphi_a/F:vmass_b/F:vpt_b/F:veta_b/F:vphi_b/F:id_1_a/I:id_2_a/I:id_1_b/I:id_2_b/I:pt_1_a/F:eta_1_a/F:phi_1_a/F:pt_2_a/F:eta_2_a/F:phi_2_a/F:pt_1_b/F:eta_1_b/F:phi_1_b/F:pt_2_b/F:eta_2_b/F:phi_2_b/F:pt_zz/F:eta_zz/F:phi_ZZ/F:m_zz/F"

class HZZGenInfo : public TObject 
{ 
  public:
    HZZGenInfo(){}
    ~HZZGenInfo(){}

  float weight;				// event weight
  int   genQid;                         // largest abs(id) of quark in the hard interaction
  int   pid_1, pid_2;			// parton ID
  float x_1, x_2;			// parton momentum fraction
  int   id_a, id_b;	    		// boson IDs
  float vmass_a, vpt_a, veta_a, vphi_a;	// boson A info
  float vmass_b, vpt_b, veta_b, vphi_b;	// boson B info
  int   id_1_a, id_2_a;			// lepton/quark IDs
  int   id_1_b, id_2_b;			// lepton/quark IDs
  float pt_1_a, eta_1_a, phi_1_a;       // lepton info
  float pt_2_a, eta_2_a, phi_2_a;  
  float pt_1_b, eta_1_b, phi_1_b;
  float pt_2_b, eta_2_b, phi_2_b;  

  float pt_zz, eta_zz, phi_zz, m_zz;

  ClassDef(HZZGenInfo,1)      
};

#endif
