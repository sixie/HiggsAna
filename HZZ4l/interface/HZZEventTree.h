#ifndef HZZEventTree_H
#define HZZEventTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

class HZZEventTree {

  public:



    /// variables
    Float_t                 fWeight;
    UInt_t                  fRunNumber;
    UInt_t                  fLumiSectionNumber;
    UInt_t                  fEventNumber;
    Float_t                 fRho;
    UInt_t                  fNVtx;
    Float_t                 fMet;

    Int_t                   fGenLep1Type;
    Float_t                 fGenLep1Pt; 
    Float_t                 fGenLep1Eta; 
    Float_t                 fGenLep1Phi; 
    Int_t                   fGenLep2Type;
    Float_t                 fGenLep2Pt; 
    Float_t                 fGenLep2Eta; 
    Float_t                 fGenLep2Phi; 
    Int_t                   fGenLep3Type;
    Float_t                 fGenLep3Pt; 
    Float_t                 fGenLep3Eta; 
    Float_t                 fGenLep3Phi; 
    Int_t                   fGenLep4Type;
    Float_t                 fGenLep4Pt; 
    Float_t                 fGenLep4Eta; 
    Float_t                 fGenLep4Phi; 

    Int_t                   fLep1Type;
    Float_t                 fLep1Pt; 
    Float_t                 fLep1Eta; 
    Float_t                 fLep1Phi; 
    Int_t                   fLep2Type;
    Float_t                 fLep2Pt; 
    Float_t                 fLep2Eta; 
    Float_t                 fLep2Phi; 
    Int_t                   fLep3Type;
    Float_t                 fLep3Pt; 
    Float_t                 fLep3Eta; 
    Float_t                 fLep3Phi; 
    Int_t                   fLep4Type;
    Float_t                 fLep4Pt; 
    Float_t                 fLep4Eta; 
    Float_t                 fLep4Phi; 


    Int_t                   fChannel;
    Float_t                 fZ1Mass; 
    Float_t                 fZ1Pt; 
    Float_t                 fZ1Eta; 
    Float_t                 fZ1Rapidity; 
    Float_t                 fZ2Mass; 
    Float_t                 fZ2Pt; 
    Float_t                 fZ2Eta; 
    Float_t                 fZ2Rapidity; 
    Float_t                 fFourLeptonPt; 
    Float_t                 fFourLeptonMass; 
    Float_t                 fFourLeptonRapidity; 

    Bool_t                  fPassFullSelection;


  public:
    /// this is the main element
    TTree *tree_;
    TFile *f_;
    TDirectory *dir_;
  
    /// hold the names of variables to facilitate things (filled during Init)
    std::vector<std::string> variables_;

    /// default constructor  
    HZZEventTree()  {};
    /// default destructor
    ~HZZEventTree(){ 
      if (f_) f_->Close();  
    };
    
    /// initialize varibles and fill list of available variables
    void InitVariables() {

      fWeight               =  0.0;
      fRunNumber            =  0.0;
      fLumiSectionNumber    =  0.0;
      fEventNumber          =  0.0;
      fRho                  =  0.0;
      fNVtx                 =  0.0;
      fMet                  =  0.0;
      fGenLep1Type          =  0.0;
      fGenLep1Pt            =  0.0;
      fGenLep1Eta           =  0.0;
      fGenLep1Phi           =  0.0;
      fGenLep2Type          =  0.0;
      fGenLep2Pt            =  0.0;
      fGenLep2Eta           =  0.0;
      fGenLep2Phi           =  0.0;
      fGenLep3Type          =  0.0;
      fGenLep3Pt            =  0.0;
      fGenLep3Eta           =  0.0;
      fGenLep3Phi           =  0.0;
      fGenLep4Type          =  0.0;
      fGenLep4Pt            =  0.0;
      fGenLep4Eta           =  0.0;
      fGenLep4Phi           =  0.0;
      fLep1Type             =  0.0;
      fLep1Pt               =  0.0;
      fLep1Eta              =  0.0;
      fLep1Phi              =  0.0;
      fLep2Type             =  0.0;
      fLep2Pt               =  0.0;
      fLep2Eta              =  0.0;
      fLep2Phi              =  0.0;
      fLep3Type             =  0.0;
      fLep3Pt               =  0.0;
      fLep3Eta              =  0.0;
      fLep3Phi              =  0.0;
      fLep4Type             =  0.0;
      fLep4Pt               =  0.0;
      fLep4Eta              =  0.0;
      fLep4Phi              =  0.0;
      fChannel              =  0.0;
      fZ1Mass               =  0.0;
      fZ1Pt                 =  0.0;
      fZ1Eta                =  0.0;
      fZ1Rapidity           =  0.0;
      fZ2Mass               =  0.0;
      fZ2Pt                 =  0.0;
      fZ2Eta                =  0.0;
      fZ2Rapidity           =  0.0;
      fFourLeptonPt         =  0.0;
      fFourLeptonMass       =  0.0;
      fFourLeptonRapidity   =  0.0;
      fPassFullSelection    =  kFALSE;

    }
    
    /// load a HZZEventTree
    void LoadTree(const char* file){
      f_ = TFile::Open(file);
      assert(f_);
      tree_ = dynamic_cast<TTree*>(f_->Get("probe_tree"));
      if (!tree_) {
        dir_ = dynamic_cast<TDirectory*>(f_->Get("zz4lTree"));
        assert(dir_);
        tree_ = dynamic_cast<TTree*>(dir_->Get("probe_tree"));
      }
      assert(tree_);
    }
    
    /// create a HZZEventTree
    void CreateTree(){
      tree_ = new TTree("probe_tree","probe_tree");
      f_ = 0;

      //book the branches

      tree_->Branch("weight", &fWeight, "weight/F");
      tree_->Branch("run",&fRunNumber,"run/i");     
      tree_->Branch("lumi",&fLumiSectionNumber, "lumi/i");
      tree_->Branch("event",&fEventNumber, "event/i");     
      tree_->Branch("rho",&fRho,"rho/F");
      tree_->Branch("recoVertices",&fNVtx,"recoVertices/i");
      tree_->Branch("pfmet",&fMet,"pfmet/F");    

      tree_->Branch("genl1pdgId",&fGenLep1Type ,"genl1pdgId/I");         
      tree_->Branch("genl1pt",&fGenLep1Pt   ,"genl1pt/F");         
      tree_->Branch("genl1eta",&fGenLep1Eta ,"genl1eta/F");          
      tree_->Branch("genl1phi",&fGenLep1Phi  ,"genl1phi/F");         
      tree_->Branch("genl2pdgId",&fGenLep2Type ,"genl2pdgId/I");         
      tree_->Branch("genl2pt",&fGenLep2Pt   ,"genl2pt/F");         
      tree_->Branch("genl2eta",&fGenLep2Eta ,"genl2eta/F");          
      tree_->Branch("genl2phi",&fGenLep2Phi  ,"genl2phi/F");         
      tree_->Branch("genl3pdgId",&fGenLep3Type ,"genl3pdgId/I");         
      tree_->Branch("genl3pt",&fGenLep3Pt   ,"genl3pt/F");         
      tree_->Branch("genl3eta",&fGenLep3Eta ,"genl3eta/F");          
      tree_->Branch("genl3phi",&fGenLep3Phi  ,"genl3phi/F");         
      tree_->Branch("genl4pdgId",&fGenLep4Type ,"genl4pdgId/I");         
      tree_->Branch("genl4pt",&fGenLep4Pt   ,"genl4pt/F");         
      tree_->Branch("genl4eta",&fGenLep4Eta ,"genl4eta/F");          
      tree_->Branch("genl4phi",&fGenLep4Phi  ,"genl4phi/F");         

      tree_->Branch("l1pdgId",&fLep1Type,"l1pdgId/I");         
      tree_->Branch("l1pt",&fLep1Pt,"l1pt/F");           
      tree_->Branch("l1eta",&fLep1Eta,"l1eta/F");           
      tree_->Branch("l1phi",&fLep1Phi,"l1phi/F");           
      tree_->Branch("l2pdgId",&fLep2Type,"l2pdgId/I");         
      tree_->Branch("l2pt",&fLep2Pt,"l2pt/F");           
      tree_->Branch("l2eta",&fLep2Eta,"l2eta/F");           
      tree_->Branch("l2phi",&fLep2Phi,"l2phi/F");           
      tree_->Branch("l3pdgId",&fLep3Type,"l3pdgId/I");         
      tree_->Branch("l3pt",&fLep3Pt,"l3pt/F");           
      tree_->Branch("l3eta",&fLep3Eta,"l3eta/F");           
      tree_->Branch("l3phi",&fLep3Phi,"l3phi/F");           
      tree_->Branch("l4pdgId",&fLep4Type,"l4pdgId/I");         
      tree_->Branch("l4pt",&fLep4Pt,"l4pt/F");           
      tree_->Branch("l4eta",&fLep4Eta,"l4eta/F");           
      tree_->Branch("l4phi",&fLep4Phi,"l4phi/F");           

      tree_->Branch("channel",&fChannel,"channel/I");            
      tree_->Branch("z1mass",&fZ1Mass,"z1mass/F");       
      tree_->Branch("z1pt",&fZ1Pt,"z1pt/F");       
      tree_->Branch("z1eta",&fZ1Eta,"z1eta/F");       
      tree_->Branch("z1rap",&fZ1Rapidity,"z1rap/F");       
      tree_->Branch("z2mass",&fZ2Mass,"z2mass/F");       
      tree_->Branch("z2pt",&fZ2Pt,"z2pt/F");       
      tree_->Branch("z2eta",&fZ2Eta,"z2eta/F");       
      tree_->Branch("z2rap",&fZ2Rapidity,"z2rap/F");       
      tree_->Branch("pt",&fFourLeptonPt,"pt/F");       
      tree_->Branch("mass",&fFourLeptonMass,"mass/F");
      tree_->Branch("rap",&fFourLeptonRapidity,"rap/F");  
      tree_->Branch("passFullSelection",&fPassFullSelection,"rap/O");  

    } 

    // initialze a HZZEventTree
    void InitTree(){
      assert(tree_);
      // don't forget to set pointers to zero before you set address
      // or you will fully appreciate that "ROOT sucks" :)
      InitVariables();
      //Set branch address
      Int_t currentState = gErrorIgnoreLevel;
      // gErrorIgnoreLevel = kError;
      gErrorIgnoreLevel = kBreak;


      tree_->SetBranchAddress("weight", &fWeight);
      tree_->SetBranchAddress("run",&fRunNumber);
      tree_->SetBranchAddress("lumi",&fLumiSectionNumber);
      tree_->SetBranchAddress("event",&fEventNumber);     
      tree_->SetBranchAddress("rho",&fRho);
      tree_->SetBranchAddress("recoVertices",&fNVtx);
      tree_->SetBranchAddress("pfmet",&fMet);
      tree_->SetBranchAddress("genl1pdgId",&fGenLep1Type);
      tree_->SetBranchAddress("genl1pt",&fGenLep1Pt);
      tree_->SetBranchAddress("genl1eta",&fGenLep1Eta );
      tree_->SetBranchAddress("genl1phi",&fGenLep1Phi  );
      tree_->SetBranchAddress("genl2pdgId",&fGenLep2Type );
      tree_->SetBranchAddress("genl2pt",&fGenLep2Pt   );
      tree_->SetBranchAddress("genl2eta",&fGenLep2Eta );
      tree_->SetBranchAddress("genl2phi",&fGenLep2Phi  );
      tree_->SetBranchAddress("genl3pdgId",&fGenLep3Type );
      tree_->SetBranchAddress("genl3pt",&fGenLep3Pt   );
      tree_->SetBranchAddress("genl3eta",&fGenLep3Eta );
      tree_->SetBranchAddress("genl3phi",&fGenLep3Phi  );
      tree_->SetBranchAddress("genl4pdgId",&fGenLep4Type );
      tree_->SetBranchAddress("genl4pt",&fGenLep4Pt   );
      tree_->SetBranchAddress("genl4eta",&fGenLep4Eta );
      tree_->SetBranchAddress("genl4phi",&fGenLep4Phi  );

      tree_->SetBranchAddress("l1pdgId",&fLep1Type);
      tree_->SetBranchAddress("l1pt",&fLep1Pt);
      tree_->SetBranchAddress("l1eta",&fLep1Eta);
      tree_->SetBranchAddress("l1phi",&fLep1Phi);
      tree_->SetBranchAddress("l2pdgId",&fLep2Type);
      tree_->SetBranchAddress("l2pt",&fLep2Pt);
      tree_->SetBranchAddress("l2eta",&fLep2Eta);
      tree_->SetBranchAddress("l2phi",&fLep2Phi);
      tree_->SetBranchAddress("l3pdgId",&fLep3Type);
      tree_->SetBranchAddress("l3pt",&fLep3Pt);
      tree_->SetBranchAddress("l3eta",&fLep3Eta);
      tree_->SetBranchAddress("l3phi",&fLep3Phi);
      tree_->SetBranchAddress("l4pdgId",&fLep4Type);
      tree_->SetBranchAddress("l4pt",&fLep4Pt);
      tree_->SetBranchAddress("l4eta",&fLep4Eta);
      tree_->SetBranchAddress("l4phi",&fLep4Phi);

      tree_->SetBranchAddress("channel",&fChannel);
      tree_->SetBranchAddress("z1mass",&fZ1Mass);
      tree_->SetBranchAddress("z1pt",&fZ1Pt);
      tree_->SetBranchAddress("z1eta",&fZ1Eta);
      tree_->SetBranchAddress("z1rap",&fZ1Rapidity);
      tree_->SetBranchAddress("z2mass",&fZ2Mass);
      tree_->SetBranchAddress("z2pt",&fZ2Pt);
      tree_->SetBranchAddress("z2eta",&fZ2Eta);
      tree_->SetBranchAddress("z2rap",&fZ2Rapidity);
      tree_->SetBranchAddress("pt",&fFourLeptonPt);
      tree_->SetBranchAddress("mass",&fFourLeptonMass);
      tree_->SetBranchAddress("rap",&fFourLeptonRapidity);
      tree_->SetBranchAddress("passFullSelection",&fPassFullSelection);

      gErrorIgnoreLevel = currentState;
    }

}; 



#endif
