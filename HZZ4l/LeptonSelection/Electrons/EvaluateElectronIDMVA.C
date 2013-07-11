//root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet0LowPt.root","output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt","Likelihood,LikelihoodD,BDT,BDTG,MLPBNN")'
//root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","output/ElectronNtuple.Real.Subdet0LowPt.root","Subdet0LowPt","Likelihood,LikelihoodD,BDT,BDTG,MLPBNN")'





/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TBranch.h"
#include "TRandom.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"

#include "HiggsAna/CommonData/interface/ElectronTree.h"
// #include "HiggsAna/HZZ4l/LeptonSelection/Electrons/TMVAGui.C"
#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace std;
using namespace TMVA;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

//--------------------------------------------------------------------

void fillUnderOverFlow(TH1F *h1, float value, float weight = 1)
{
  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void EvaluateElectronIDMVA(
string inputFile      = "", 
string outputFile     = "",
TString label         = "ElectronIDMVA"
) {   
#ifdef __CINT__
  gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

  //--------------------------------------------------------------------
  // path to weights dir (this is where MVA training info is stored)
  // output root file will be stored at [path]/output
  //--------------------------------------------------------------------
  
  vector<string> samples;
  samples.push_back(inputFile);

  //--------------------------------------------------------------------------------
  // IMPORTANT: set the following variables to the same set used for MVA training!!!
  //--------------------------------------------------------------------------------

  std::map<std::string,int> mvaVar;

  mvaVar[ "fbrem" ] = 1;
  mvaVar[ "deta" ] = 1;
  mvaVar[ "dphi" ] = 1;
  mvaVar[ "see" ] = 1;
  mvaVar[ "etawidth" ] = 1;
  mvaVar[ "phiwidth" ] = 1;
  mvaVar[ "HoE" ] = 1;
  mvaVar[ "EoP" ] = 1;
  mvaVar[ "e1x5e5x5" ] = 1;
  mvaVar[ "EoPout" ] = 1;
  mvaVar[ "detacalo" ] = 1;
  mvaVar[ "kfchi2" ] = 1;
  mvaVar[ "kfhits" ] = 1;
  mvaVar[ "spp" ] = 1;
  mvaVar[ "IoEmIoP" ] = 1;
  mvaVar[ "R9" ] = 1;
  mvaVar[ "gsfchi2" ] = 1;
  mvaVar[ "PreShowerOverRaw" ] = 1; 

  //---------------------------------------------------------------
  // specifies the selection applied to events in the training
  //---------------------------------------------------------------

  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;

  Use["BDTCat_BDTG_V0"] = 1;
  Use["BDTCat_BDTG_V1"] = 1;
  Use["BDTCat_BDTG_V2"] = 0;
  Use["BDTCat_BDTG_V3"] = 0;
  Use["BDTCat_BDTG_V4"] = 0;
  Use["BDTCat_BDTG_V5"] = 0;

  // ---------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "==> Start TMVAClassificationApplication" << std::endl;

//   // Select methods (don't look at this code - not of interest)
//   if (myMethodList != "") {
//     for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

//     std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
//     for (UInt_t i=0; i<mlist.size(); i++) {
//       std::string regMethod(mlist[i]);

//       if (Use.find(regMethod) == Use.end()) {
//         std::cout << "Method \"" << regMethod 
//                   << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
//         for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
//           std::cout << it->first << " ";
//         }
//         std::cout << std::endl;
//         return;
//       }
//       Use[regMethod] = 1;
//     }
//   }

  // --------------------------------------------------------------------------------------------------

  const unsigned int nsamples = samples.size();
  
  for( unsigned int i = 0 ; i < nsamples ; ++i ){

    // --- Create the Reader object

    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    
    TMVA::Reader *reader2 = new TMVA::Reader( "!Color:!Silent" );    
    TMVA::Reader  *fTMVAReader[6];
    TMVA::Reader  *fTMVAReaderICHEP2012[6];
    for(UInt_t j=0; j<6; ++j) {
      fTMVAReaderICHEP2012[j] = new TMVA::Reader( "!Color:!Silent:Error" );  
    }
    for(UInt_t j=0; j<6; ++j) {
      fTMVAReader[j] = new TMVA::Reader( "!Color:!Silent:Error" );  
    }
    TMVA::Reader  *fTMVAReader3[6];
    for(UInt_t j=0; j<6; ++j) {
      fTMVAReader3[j] = new TMVA::Reader( "!Color:!Silent:Error" );  
    }

    Float_t varfbrem;
    Float_t vardeta;
    Float_t vardphi;
    Float_t varsee;
    Float_t varetawidth;
    Float_t varphiwidth;
    Float_t varHoE;
    Float_t varEoP;
    Float_t varEoPout;
    Float_t vardetacalo;
    Float_t varkfchi2;
    Float_t varkfhits;
    Float_t varspp;
    Float_t varIoEmIoP;
    Float_t varR9;
    Float_t vargsfchi2;
    Float_t varPreShowerOverRaw;
    Float_t varkfhitsall;
    Float_t vareleEoPout;
    Float_t varOneMinuse1x5e5x5;
    Float_t varIoEmIoP_ICHEP2012;
    Float_t vardetaabs;
    Float_t vardphiabs;
    Float_t vardetacaloabs;
    Float_t varPt, varEta, varEventNumberParity;
    Int_t varCiC;

    for(UInt_t j=0; j<6; ++j) {

      //ICHEP2012 weights
      fTMVAReaderICHEP2012[j]->AddVariable( "fbrem", &varfbrem );
      fTMVAReaderICHEP2012[j]->AddVariable( "kfchi2", &varkfchi2 );
      fTMVAReaderICHEP2012[j]->AddVariable( "kfhits", &varkfhits );
      fTMVAReaderICHEP2012[j]->AddVariable( "gsfchi2", &vargsfchi2 );
      fTMVAReaderICHEP2012[j]->AddVariable( "deta", &vardetaabs );
      fTMVAReaderICHEP2012[j]->AddVariable( "dphi", &vardphiabs );
      fTMVAReaderICHEP2012[j]->AddVariable( "detacalo", &vardetacaloabs );
      fTMVAReaderICHEP2012[j]->AddVariable( "see", &varsee );
      fTMVAReaderICHEP2012[j]->AddVariable( "spp", &varspp );
      fTMVAReaderICHEP2012[j]->AddVariable( "etawidth", &varetawidth );
      fTMVAReaderICHEP2012[j]->AddVariable( "phiwidth", &varphiwidth );
      fTMVAReaderICHEP2012[j]->AddVariable( "e1x5e5x5", &varOneMinuse1x5e5x5 );
      fTMVAReaderICHEP2012[j]->AddVariable( "R9", &varR9 );
      fTMVAReaderICHEP2012[j]->AddVariable( "HoE", &varHoE );
      fTMVAReaderICHEP2012[j]->AddVariable( "EoP", &varEoP );
      fTMVAReaderICHEP2012[j]->AddVariable( "IoEmIoP", &varIoEmIoP_ICHEP2012 );
      fTMVAReaderICHEP2012[j]->AddVariable( "eleEoPout", &vareleEoPout );
      if (j==2 || j==5) {
        fTMVAReaderICHEP2012[j]->AddVariable( "PreShowerOverRaw", &varPreShowerOverRaw );
      }
      fTMVAReaderICHEP2012[j]->AddSpectator( "eta",   &varEta );
      fTMVAReaderICHEP2012[j]->AddSpectator( "pt" ,   &varPt  );

      //Si's weights
      fTMVAReader3[j]->AddVariable( "fbrem", &varfbrem );
      fTMVAReader3[j]->AddVariable( "kfchi2", &varkfchi2 );
      fTMVAReader3[j]->AddVariable( "kfhits", &varkfhits );
      fTMVAReader3[j]->AddVariable( "gsfchi2", &vargsfchi2 );
      fTMVAReader3[j]->AddVariable( "deta", &vardeta );
      fTMVAReader3[j]->AddVariable( "dphi", &vardphi );
      fTMVAReader3[j]->AddVariable( "detacalo", &vardetacalo );
      fTMVAReader3[j]->AddVariable( "see", &varsee );
      fTMVAReader3[j]->AddVariable( "spp", &varspp );
      fTMVAReader3[j]->AddVariable( "etawidth", &varetawidth );
      fTMVAReader3[j]->AddVariable( "phiwidth", &varphiwidth );
      fTMVAReader3[j]->AddVariable( "OneMinusSeedE1x5OverE5x5", &varOneMinuse1x5e5x5 );
      fTMVAReader3[j]->AddVariable( "R9", &varR9 );
      fTMVAReader3[j]->AddVariable( "HoE", &varHoE );
      fTMVAReader3[j]->AddVariable( "EoP", &varEoP );
      fTMVAReader3[j]->AddVariable( "IoEmIoP", &varIoEmIoP );
//       fTMVAReader3[j]->AddVariable( "eleEoPout", &vareleEoPout );
      fTMVAReader3[j]->AddVariable( "EEleoPout", &vareleEoPout );
      if (j==2 || j==5) {
        fTMVAReader3[j]->AddVariable( "PreShowerOverRaw", &varPreShowerOverRaw );
      }
      fTMVAReader3[j]->AddSpectator( "eta",   &varEta );
      fTMVAReader3[j]->AddSpectator( "pt" ,   &varPt  );


    }

    // --- Book the MVA methods

    //--------------------------------------------------------------------------------------
    // tell evaluateMVA_smurf_hww where to find the weights dir, which contains the trained MVA's. 
    // In this example, the weights dir is located at [path]/[dir]
    // and the output root file is written to [path]/[output]
    //--------------------------------------------------------------------------------------

    TString dir    = "weights/";
    TString outdir = "output/";
    TString prefix = label;

    // Book method(s)
    for(UInt_t j=0; j<6; ++j) {
      if(Use["BDTCat_BDTG_V0"]) {
        if (j==0) fTMVAReaderICHEP2012[j]->BookMVA( "ElectronIDMVA_V0_BDTG method", "HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat1.weights.xml");                
        if (j==1) fTMVAReaderICHEP2012[j]->BookMVA( "ElectronIDMVA_V0_BDTG method", "HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat2.weights.xml");                
        if (j==2) fTMVAReaderICHEP2012[j]->BookMVA( "ElectronIDMVA_V0_BDTG method", "HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat3.weights.xml");                
        if (j==3) fTMVAReaderICHEP2012[j]->BookMVA( "ElectronIDMVA_V0_BDTG method", "HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat4.weights.xml");                
        if (j==4) fTMVAReaderICHEP2012[j]->BookMVA( "ElectronIDMVA_V0_BDTG method", "HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat5.weights.xml");                
        if (j==5) fTMVAReaderICHEP2012[j]->BookMVA( "ElectronIDMVA_V0_BDTG method", "HiggsAna/HZZ4l/data/ElectronIDMVAWeights/Electrons_BDTG_NonTrigV0_Cat6.weights.xml");                
      }
      if(Use["BDTCat_BDTG_V1"]) {
        //reader->BookMVA( "ElectronIDMVA_BDTG_V0 method", dir + prefix + TString("_BDTCat_BDTG_V0.weights.xml"));              
        if (j==0) fTMVAReader3[j]->BookMVA( "ElectronIDMVA_V1_BDTG method", "/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/weights/ElectronIDMVA_V1_CentralPt5To10_V1_BDTG.weights.xml");
        if (j==1) fTMVAReader3[j]->BookMVA( "ElectronIDMVA_V1_BDTG method", "/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/weights/ElectronIDMVA_V1_TransitionPt5To10_V1_BDTG.weights.xml");
        if (j==2) fTMVAReader3[j]->BookMVA( "ElectronIDMVA_V1_BDTG method", "/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/weights/ElectronIDMVA_V1_EndcapPt5To10_V1_BDTG.weights.xml");
        if (j==3) fTMVAReader3[j]->BookMVA( "ElectronIDMVA_V1_BDTG method", "/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/weights/ElectronIDMVA_V1_CentralPt10ToInf_V1_BDTG.weights.xml");
        if (j==4) fTMVAReader3[j]->BookMVA( "ElectronIDMVA_V1_BDTG method", "/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/weights/ElectronIDMVA_V1_TransitionPt10ToInf_V1_BDTG.weights.xml");
        if (j==5) fTMVAReader3[j]->BookMVA( "ElectronIDMVA_V1_BDTG method", "/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/LeptonSelection/Electrons/IDMVA/weights/ElectronIDMVA_V1_EndcapPt10ToInf_V1_BDTG.weights.xml");        
      }
      if(Use["BDTCat_BDTG_V2"]) {
        if (j==0) fTMVAReader3[j]->BookMVA( "ElectronIDMVA_V2_BDTG method", "weights/ElectronIDMVA_V0_CentralPt5To10_V0_BDTG.weights.xml");                
        if (j==1) fTMVAReader3[j]->BookMVA( "ElectronIDMVA_V2_BDTG method", "weights/ElectronIDMVA_V0_TransitionPt5To10_V0_BDTG.weights.xml");                
        if (j==2) fTMVAReader3[j]->BookMVA( "ElectronIDMVA_V2_BDTG method", "weights/ElectronIDMVA_V0_EndcapPt5To10_V0_BDTG.weights.xml");                
        if (j==3) fTMVAReader3[j]->BookMVA( "ElectronIDMVA_V2_BDTG method", "weights/ElectronIDMVA_V0_CentralPt10ToInf_V0_BDTG.weights.xml");                
        if (j==4) fTMVAReader3[j]->BookMVA( "ElectronIDMVA_V2_BDTG method", "weights/ElectronIDMVA_V0_TransitionPt10ToInf_V0_BDTG.weights.xml");                
        if (j==5) fTMVAReader3[j]->BookMVA( "ElectronIDMVA_V2_BDTG method", "weights/ElectronIDMVA_V0_EndcapPt10ToInf_V0_BDTG.weights.xml");                
      }
    }
 
    //*****************************************************************************************
    // Load input tree
    //*****************************************************************************************
    ElectronTree eleTree;
    eleTree.LoadTree(inputFile.c_str());
    eleTree.InitTree();

    //*****************************************************************************************
    // Output tree
    //*****************************************************************************************
    TFile *EleOutputFile = new TFile(outputFile.c_str(), "RECREATE");
    //EleOutputFile->cd();
    TTree *EleOutputTree = eleTree.tree_->CloneTree(-1, "fast");

    //Add new branches for MVA output
    Float_t fEleIDMVA_BDTG_V0;
    Float_t fEleIDMVA_BDTG_V1;
    Float_t fEleIDMVA_BDTG_V2;
    Float_t fEleIDMVA_BDTG_V3;
    Float_t fEleIDMVA_BDTG_V4;
    Float_t fEleIDMVA_BDTG_V5;
    TBranch* branchEleIDMVA_BDTG_V0 = 0;
    TBranch* branchEleIDMVA_BDTG_V1 = 0;
    TBranch* branchEleIDMVA_BDTG_V2 = 0;
    TBranch* branchEleIDMVA_BDTG_V3 = 0;
    TBranch* branchEleIDMVA_BDTG_V4 = 0;
    TBranch* branchEleIDMVA_BDTG_V5 = 0;
    if(Use["BDTCat_BDTG_V0"])  branchEleIDMVA_BDTG_V0 = EleOutputTree->Branch( "EleIDMVA_BDTG_V0" , &fEleIDMVA_BDTG_V0, "EleIDMVA_BDTG_V0/F" );
    if(Use["BDTCat_BDTG_V1"])  branchEleIDMVA_BDTG_V1 = EleOutputTree->Branch( "EleIDMVA_BDTG_V1" , &fEleIDMVA_BDTG_V1, "EleIDMVA_BDTG_V1/F" );
    if(Use["BDTCat_BDTG_V2"])  branchEleIDMVA_BDTG_V2 = EleOutputTree->Branch( "EleIDMVA_BDTG_V2" , &fEleIDMVA_BDTG_V2, "EleIDMVA_BDTG_V2/F" );
    if(Use["BDTCat_BDTG_V3"])  branchEleIDMVA_BDTG_V3 = EleOutputTree->Branch( "EleIDMVA_BDTG_V3" , &fEleIDMVA_BDTG_V3, "EleIDMVA_BDTG_V3/F" );
    if(Use["BDTCat_BDTG_V4"])  branchEleIDMVA_BDTG_V4 = EleOutputTree->Branch( "EleIDMVA_BDTG_V4" , &fEleIDMVA_BDTG_V4, "EleIDMVA_BDTG_V4/F" );
    if(Use["BDTCat_BDTG_V5"])  branchEleIDMVA_BDTG_V5 = EleOutputTree->Branch( "EleIDMVA_BDTG_V5" , &fEleIDMVA_BDTG_V5, "EleIDMVA_BDTG_V5/F" );

    // --- Event loop
    std::cout << "--- Processing: " << eleTree.tree_->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();

    for (Long64_t ievt=0; ievt<eleTree.tree_->GetEntries();ievt++) {

      if (ievt%10000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      eleTree.tree_->GetEntry(ievt);

      Double_t rho = 0;
      if (!(TMath::IsNaN(eleTree.fRho) || isinf(eleTree.fRho))) rho = eleTree.fRho;

      //--------------------------------------------------------
      // important: here we associate branches to MVA variables
      //--------------------------------------------------------
 
      varfbrem = max(double(eleTree.fEleFBrem),-1.0);
      vardeta = eleTree.fEleDEtaIn;
      vardphi = eleTree.fEleDPhiIn;
      varsee = eleTree.fEleSigmaIEtaIEta;
      varetawidth = eleTree.fEleSCEtaWidth;
      varphiwidth = eleTree.fEleSCPhiWidth;
      varHoE = eleTree.fEleHoverE;
      varEoP = min(double(eleTree.fEleEOverP), 20.0);
      varOneMinuse1x5e5x5 = min(max(double(eleTree.fEleOneMinusSeedE1x5OverE5x5),-1.0),2.0);
      varEoPout = min(double(eleTree.fEleESeedClusterOverPout),20.0);
      vareleEoPout = eleTree.fEleEEleClusterOverPout;
      vardetacalo = eleTree.fEledEtaCalo;
      varkfchi2 = min(double(eleTree.fEleKFTrackChi2OverNDoF),10.0);
      varkfhits = eleTree.fEleKFTrackNHits;
      varkfhitsall = eleTree.fEleKFTrackNHits;
      varspp = eleTree.fEleSigmaIPhiIPhi;
      varIoEmIoP = eleTree.fEleOneOverEMinusOneOverP;
      varR9 = min(double(eleTree.fEleR9), 5.0);
      vargsfchi2 = min(double(eleTree.fEleGsfTrackChi2OverNdof),200.0);
      varPreShowerOverRaw = eleTree.fElePreShowerOverRaw;

      varIoEmIoP_ICHEP2012 = (1./eleTree.fEleEcalEnergy) - (1./(eleTree.fElePt * TMath::CosH(eleTree.fEleEta)));
      vardetaabs = min(fabs(double(eleTree.fEleDEtaIn)),0.06);
      vardphiabs = min(fabs(double(eleTree.fEleDPhiIn)),0.6);
      vardetacaloabs = min(fabs(double(eleTree.fEledEtaCalo)),0.2);
      varPt = eleTree.fElePt;
      varEta = eleTree.fEleEta;
      

      // --- Return the MVA outputs and weights
      TMVA::Reader  *tmpReader = reader;      

      if(Use["BDTCat_BDTG_V0"]) {
        Int_t mvabin = -1;
        if (varPt <= 10 && fabs(varEta) < 0.8) mvabin = 0;
        else if (varPt <= 10 && fabs(varEta) >= 0.8 && fabs(varEta) < 1.479) mvabin = 1;
        else if (varPt <= 10 && fabs(varEta) >= 1.479) mvabin = 2;
        else if (varPt > 10 && fabs(varEta) < 0.8) mvabin = 3;
        else if (varPt > 10 && fabs(varEta) >= 0.8 && fabs(varEta) < 1.479) mvabin = 4;
        else if (varPt > 10 && fabs(varEta) >= 1.479) mvabin = 5;
        assert(mvabin >= 0);

        fEleIDMVA_BDTG_V0            = fTMVAReaderICHEP2012[mvabin]->EvaluateMVA( "ElectronIDMVA_V0_BDTG method" );
        branchEleIDMVA_BDTG_V0->Fill();
      }
      if(Use["BDTCat_BDTG_V1"]) {
        Int_t mvabin = -1;
        if (varPt <= 10 && fabs(varEta) < 0.8) mvabin = 0;
        else if (varPt <= 10 && fabs(varEta) >= 0.8 && fabs(varEta) < 1.479) mvabin = 1;
        else if (varPt <= 10 && fabs(varEta) >= 1.479) mvabin = 2;
        else if (varPt > 10 && fabs(varEta) < 0.8) mvabin = 3;
        else if (varPt > 10 && fabs(varEta) >= 0.8 && fabs(varEta) < 1.479) mvabin = 4;
        else if (varPt > 10 && fabs(varEta) >= 1.479) mvabin = 5;
        assert(mvabin >= 0);

        fEleIDMVA_BDTG_V1            = fTMVAReader3[mvabin]->EvaluateMVA( "ElectronIDMVA_V1_BDTG method" );
        branchEleIDMVA_BDTG_V1->Fill();
      }
      if(Use["BDTCat_BDTG_V2"]) {
        Int_t mvabin = -1;
        if (varPt <= 10 && fabs(varEta) < 0.8) mvabin = 0;
        else if (varPt <= 10 && fabs(varEta) >= 0.8 && fabs(varEta) < 1.479) mvabin = 1;
        else if (varPt <= 10 && fabs(varEta) >= 1.479) mvabin = 2;
        else if (varPt > 10 && fabs(varEta) < 0.8) mvabin = 3;
        else if (varPt > 10 && fabs(varEta) >= 0.8 && fabs(varEta) < 1.479) mvabin = 4;
        else if (varPt > 10 && fabs(varEta) >= 1.479) mvabin = 5;
        assert(mvabin >= 0);

        fEleIDMVA_BDTG_V2            = fTMVAReader3[mvabin]->EvaluateMVA( "ElectronIDMVA_V2_BDTG method" );
        branchEleIDMVA_BDTG_V2->Fill();
      }

    } // End main loop

    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();

    //Write Output file
    EleOutputFile->Write();
    EleOutputFile->Close();

    delete reader;
    
    std::cout << "==> TMVAClassificationApplication is done with sample " << samples.at(i) << endl << std::endl;
  }
}
