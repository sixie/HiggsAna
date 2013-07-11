#ifndef HIGGSANA_HZZ4L_KINEMATICS
#define HIGGSANA_HZZ4L_KINEMATICS

#include <TObject.h>

// #define VARLIST_GENINFO "l1type/I:l2type/I:l3type/I:l4type/I:l1pt/F:l2pt/F:l3pt/F:l4pt/F:l1eta/F:l2eta/F:l3eta/F:l4eta/F:l1phi/F:l2phi/F:l3phi/F:l4phi/F:Z1pt/F:Z1pt/F:ZZpt/F:Z1eta/F:Z1eta/F:ZZeta/F:mZ1/F:mZ2/F:m4l/F:channel:I"

class HZZKinematics : public TObject 
{ 
  public:
    HZZKinematics(){}
    ~HZZKinematics(){}
    
    int l1type,l2type,l3type,l4type;
    float l1pt,l2pt,l3pt,l4pt;
    float l1eta,l2eta,l3eta,l4eta;
    float l1phi,l2phi,l3phi,l4phi;
    float Z1pt, Z2pt, ZZpt;
    float Z1eta, Z2eta, ZZeta;
    float mZ1, mZ2, m4l;
    int channel;

    ClassDef(HZZKinematics,1)      
};

#endif
