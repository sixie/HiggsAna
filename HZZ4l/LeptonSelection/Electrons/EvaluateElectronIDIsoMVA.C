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

void EvaluateElectronIDIsoMVA(
  string inputFile      = "", 
  string outputFile     = "",
  TString label         = "ElectronIDIsoMVA"
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
    for(UInt_t j=0; j<6; ++j) {
      fTMVAReader[j] = new TMVA::Reader( "!Color:!Silent:Error" );  
    }
    TMVA::Reader  *fTMVAReader2[6];
    for(UInt_t j=0; j<6; ++j) {
      fTMVAReader2[j] = new TMVA::Reader( "!Color:!Silent:Error" );  
    }

    Float_t varfbrem;
    Float_t vardeta;
    Float_t vardphi;
    Float_t varsee;
    Float_t varetawidth;
    Float_t varphiwidth;
    Float_t varHoE;
    Float_t varEoP;
    Float_t vare1x5e5x5;
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
    Float_t varPt, varEta, varEventNumberParity;
    Int_t varCiC;
    Float_t                 varChargedIso_DR0p0To0p1;
    Float_t                 varChargedIso_DR0p1To0p2;
    Float_t                 varChargedIso_DR0p2To0p3;
    Float_t                 varChargedIso_DR0p3To0p4;
    Float_t                 varChargedIso_DR0p4To0p5;
    Float_t                 varGammaIso_DR0p0To0p1;
    Float_t                 varGammaIso_DR0p1To0p2;
    Float_t                 varGammaIso_DR0p2To0p3;
    Float_t                 varGammaIso_DR0p3To0p4;
    Float_t                 varGammaIso_DR0p4To0p5;
    Float_t                 varNeutralHadronIso_DR0p0To0p1;
    Float_t                 varNeutralHadronIso_DR0p1To0p2;
    Float_t                 varNeutralHadronIso_DR0p2To0p3;
    Float_t                 varNeutralHadronIso_DR0p3To0p4;
    Float_t                 varNeutralHadronIso_DR0p4To0p5;


    for(UInt_t j=0; j<6; ++j) {
      //Duncan's weights
      fTMVAReader[j]->AddVariable( "fbrem", &varfbrem );
      fTMVAReader[j]->AddVariable( "kfchi2", &varkfchi2 );
      fTMVAReader[j]->AddVariable( "kfhits", &varkfhits );
      //fTMVAReader[j]->AddVariable( "kfhitsall", &varkfhitsall );
      fTMVAReader[j]->AddVariable( "gsfchi2", &vargsfchi2 );
      fTMVAReader[j]->AddVariable( "deta", &vardeta );
      fTMVAReader[j]->AddVariable( "dphi", &vardphi );
      fTMVAReader[j]->AddVariable( "detacalo", &vardetacalo );
      fTMVAReader[j]->AddVariable( "see", &varsee );
      fTMVAReader[j]->AddVariable( "spp", &varspp );
      fTMVAReader[j]->AddVariable( "etawidth", &varetawidth );
      fTMVAReader[j]->AddVariable( "phiwidth", &varphiwidth );
      fTMVAReader[j]->AddVariable( "e1x5e5x5", &vare1x5e5x5 );
      fTMVAReader[j]->AddVariable( "R9", &varR9 );
      fTMVAReader[j]->AddVariable( "HoE", &varHoE );
      fTMVAReader[j]->AddVariable( "EoP", &varEoP );
      fTMVAReader[j]->AddVariable( "IoEmIoP", &varIoEmIoP );
      //fTMVAReader[j]->AddVariable( "eleEoPout", &vareleEoPout );
      fTMVAReader[j]->AddVariable( "EoPout", &varEoPout );
      if (j==2 || j==5) {
        fTMVAReader[j]->AddVariable( "PreShowerOverRaw", &varPreShowerOverRaw );
      }
      fTMVAReader[j]->AddVariable( "ChargedIso_DR0p0To0p1",         &varChargedIso_DR0p0To0p1           );
      fTMVAReader[j]->AddVariable( "ChargedIso_DR0p1To0p2",         &varChargedIso_DR0p1To0p2           );
      fTMVAReader[j]->AddVariable( "ChargedIso_DR0p2To0p3",       &varChargedIso_DR0p2To0p3         );
      fTMVAReader[j]->AddVariable( "ChargedIso_DR0p3To0p4",        &varChargedIso_DR0p3To0p4          );
      fTMVAReader[j]->AddVariable( "ChargedIso_DR0p4To0p5",        &varChargedIso_DR0p4To0p5          );
      fTMVAReader[j]->AddVariable( "GammaIso_DR0p0To0p1",         &varGammaIso_DR0p0To0p1           );
      fTMVAReader[j]->AddVariable( "GammaIso_DR0p1To0p2",         &varGammaIso_DR0p1To0p2           );
      fTMVAReader[j]->AddVariable( "GammaIso_DR0p2To0p3",       &varGammaIso_DR0p2To0p3         );
      fTMVAReader[j]->AddVariable( "GammaIso_DR0p3To0p4",        &varGammaIso_DR0p3To0p4          );
      fTMVAReader[j]->AddVariable( "GammaIso_DR0p4To0p5",        &varGammaIso_DR0p4To0p5          );
      fTMVAReader[j]->AddVariable( "NeutralHadronIso_DR0p0To0p1",         &varNeutralHadronIso_DR0p0To0p1           );
      fTMVAReader[j]->AddVariable( "NeutralHadronIso_DR0p1To0p2",         &varNeutralHadronIso_DR0p1To0p2           );
      fTMVAReader[j]->AddVariable( "NeutralHadronIso_DR0p2To0p3",       &varNeutralHadronIso_DR0p2To0p3         );
      fTMVAReader[j]->AddVariable( "NeutralHadronIso_DR0p3To0p4",        &varNeutralHadronIso_DR0p3To0p4          );
      fTMVAReader[j]->AddVariable( "NeutralHadronIso_DR0p4To0p5",        &varNeutralHadronIso_DR0p4To0p5          );

      fTMVAReader[j]->AddSpectator( "eta",   &varEta );
      fTMVAReader[j]->AddSpectator( "pt" ,   &varPt  );


      //Z+Jets trained
      fTMVAReader2[j]->AddVariable( "fbrem", &varfbrem );
      fTMVAReader2[j]->AddVariable( "kfchi2", &varkfchi2 );
      fTMVAReader2[j]->AddVariable( "kfhits", &varkfhits );
      //fTMVAReader2[j]->AddVariable( "kfhitsall", &varkfhitsall );
      fTMVAReader2[j]->AddVariable( "gsfchi2", &vargsfchi2 );
      fTMVAReader2[j]->AddVariable( "deta", &vardeta );
      fTMVAReader2[j]->AddVariable( "dphi", &vardphi );
      fTMVAReader2[j]->AddVariable( "detacalo", &vardetacalo );
      fTMVAReader2[j]->AddVariable( "see", &varsee );
      fTMVAReader2[j]->AddVariable( "spp", &varspp );
      fTMVAReader2[j]->AddVariable( "etawidth", &varetawidth );
      fTMVAReader2[j]->AddVariable( "phiwidth", &varphiwidth );
      fTMVAReader2[j]->AddVariable( "e1x5e5x5", &vare1x5e5x5 );
      fTMVAReader2[j]->AddVariable( "R9", &varR9 );
      fTMVAReader2[j]->AddVariable( "HoE", &varHoE );
      fTMVAReader2[j]->AddVariable( "EoP", &varEoP );
      fTMVAReader2[j]->AddVariable( "IoEmIoP", &varIoEmIoP );
      //fTMVAReader2[j]->AddVariable( "eleEoPout", &vareleEoPout );
      fTMVAReader2[j]->AddVariable( "EoPout", &varEoPout );
      if (j==2 || j==5) {
        fTMVAReader2[j]->AddVariable( "PreShowerOverRaw", &varPreShowerOverRaw );
      }
      fTMVAReader2[j]->AddVariable( "ChargedIso_DR0p0To0p1",         &varChargedIso_DR0p0To0p1           );
      fTMVAReader2[j]->AddVariable( "ChargedIso_DR0p1To0p2",         &varChargedIso_DR0p1To0p2           );
      fTMVAReader2[j]->AddVariable( "ChargedIso_DR0p2To0p3",       &varChargedIso_DR0p2To0p3         );
      fTMVAReader2[j]->AddVariable( "ChargedIso_DR0p3To0p4",        &varChargedIso_DR0p3To0p4          );
      fTMVAReader2[j]->AddVariable( "ChargedIso_DR0p4To0p5",        &varChargedIso_DR0p4To0p5          );
      fTMVAReader2[j]->AddVariable( "GammaIso_DR0p0To0p1",         &varGammaIso_DR0p0To0p1           );
      fTMVAReader2[j]->AddVariable( "GammaIso_DR0p1To0p2",         &varGammaIso_DR0p1To0p2           );
      fTMVAReader2[j]->AddVariable( "GammaIso_DR0p2To0p3",       &varGammaIso_DR0p2To0p3         );
      fTMVAReader2[j]->AddVariable( "GammaIso_DR0p3To0p4",        &varGammaIso_DR0p3To0p4          );
      fTMVAReader2[j]->AddVariable( "GammaIso_DR0p4To0p5",        &varGammaIso_DR0p4To0p5          );
      fTMVAReader2[j]->AddVariable( "NeutralHadronIso_DR0p0To0p1",         &varNeutralHadronIso_DR0p0To0p1           );
      fTMVAReader2[j]->AddVariable( "NeutralHadronIso_DR0p1To0p2",         &varNeutralHadronIso_DR0p1To0p2           );
      fTMVAReader2[j]->AddVariable( "NeutralHadronIso_DR0p2To0p3",       &varNeutralHadronIso_DR0p2To0p3         );
      fTMVAReader2[j]->AddVariable( "NeutralHadronIso_DR0p3To0p4",        &varNeutralHadronIso_DR0p3To0p4          );
      fTMVAReader2[j]->AddVariable( "NeutralHadronIso_DR0p4To0p5",        &varNeutralHadronIso_DR0p4To0p5          );

      fTMVAReader2[j]->AddSpectator( "eta",   &varEta );
      fTMVAReader2[j]->AddSpectator( "pt" ,   &varPt  );



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
        if (j==0) fTMVAReader[j]->BookMVA( "ElectronIDIsoMVA_V0_BDTG method", "weights/ElectronIDIsoMVA_V0_CentralPt5To10_V0_BDTG.weights.xml");                
        if (j==1) fTMVAReader[j]->BookMVA( "ElectronIDIsoMVA_V0_BDTG method", "weights/ElectronIDIsoMVA_V0_TransitionPt5To10_V0_BDTG.weights.xml");                
        if (j==2) fTMVAReader[j]->BookMVA( "ElectronIDIsoMVA_V0_BDTG method", "weights/ElectronIDIsoMVA_V0_EndcapPt5To10_V0_BDTG.weights.xml");                
        if (j==3) fTMVAReader[j]->BookMVA( "ElectronIDIsoMVA_V0_BDTG method", "weights/ElectronIDIsoMVA_V0_CentralPt10ToInf_V0_BDTG.weights.xml");                
        if (j==4) fTMVAReader[j]->BookMVA( "ElectronIDIsoMVA_V0_BDTG method", "weights/ElectronIDIsoMVA_V0_TransitionPt10ToInf_V0_BDTG.weights.xml");                
        if (j==5) fTMVAReader[j]->BookMVA( "ElectronIDIsoMVA_V0_BDTG method", "weights/ElectronIDIsoMVA_V0_EndcapPt10ToInf_V0_BDTG.weights.xml");                
      }
      if(Use["BDTCat_BDTG_V1"]) {
        if (j==0) fTMVAReader2[j]->BookMVA( "ElectronIDIsoMVA_V1_BDTG method", "weights/ElectronIDIsoMVA_V1_CentralPt5To10_V1_BDTG.weights.xml");                
        if (j==1) fTMVAReader2[j]->BookMVA( "ElectronIDIsoMVA_V1_BDTG method", "weights/ElectronIDIsoMVA_V1_TransitionPt5To10_V1_BDTG.weights.xml");                
        if (j==2) fTMVAReader2[j]->BookMVA( "ElectronIDIsoMVA_V1_BDTG method", "weights/ElectronIDIsoMVA_V1_EndcapPt5To10_V1_BDTG.weights.xml");                
        if (j==3) fTMVAReader2[j]->BookMVA( "ElectronIDIsoMVA_V1_BDTG method", "weights/ElectronIDIsoMVA_V1_CentralPt10ToInf_V1_BDTG.weights.xml");                
        if (j==4) fTMVAReader2[j]->BookMVA( "ElectronIDIsoMVA_V1_BDTG method", "weights/ElectronIDIsoMVA_V1_TransitionPt10ToInf_V1_BDTG.weights.xml");                
        if (j==5) fTMVAReader2[j]->BookMVA( "ElectronIDIsoMVA_V1_BDTG method", "weights/ElectronIDIsoMVA_V1_EndcapPt10ToInf_V1_BDTG.weights.xml");                
      }
    }
 
    //*****************************************************************************************
    //Prepare Tree Variables
    //*****************************************************************************************
    //Variables
    Float_t                 fWeight;
    UInt_t                  fRunNumber;
    UInt_t                  fLumiSectionNumber;
    UInt_t                  fEventNumber;
    Bool_t                  fEleEventNumberParity;
    Float_t                 fElePt; 
    Float_t                 fEleEta; 
    Float_t                 fElePhi; 
    Float_t                 fEleSCEt; 
    Float_t                 fEleSCEta; 
    Float_t                 fEleSCPhi; 
    Bool_t                  fEleIsEcalDriven;
    Float_t                 fElePFIso; 

    //Conversion Variables
    Bool_t                  fEleMatchedConversion;
    Float_t                 fEleConvDist;
    Float_t                 fEleConvDCot;
    UInt_t                  fEleNMissHits;

    //CutBased Variables
    Float_t                 fEleSigmaIEtaIEta; 
    Float_t                 fEleDEtaIn; 
    Float_t                 fEleDPhiIn; 
    Float_t                 fEleHoverE; 
    Float_t                 fEleD0; 
    Float_t                 fEleDZ; 
    Float_t                 fEleFBrem; 
    Float_t                 fEleEOverP; 

    //Additional Vars used in Likelihood
    Float_t                 fEleESeedClusterOverPout; 
    Float_t                 fEleSigmaIPhiIPhi; 
    Float_t                 fEleNBrem; 
    Float_t                 fEleOneOverEMinusOneOverP; 
    Float_t                 fEleESeedClusterOverPIn; 
    Float_t                 fEleIP3d; 
    Float_t                 fEleIP3dSig; 

    Float_t                 fEleHcalDepth1OverEcal;
    Float_t                 fEleHcalDepth2OverEcal;
    Float_t                 fEledEtaCalo;
    Float_t                 fEledPhiCalo;
    Float_t                 fElePreShowerOverRaw;
    Float_t                 fEleSigmaIEtaIPhi;
    Float_t                 fEleSCEtaWidth;
    Float_t                 fEleSCPhiWidth;
    Float_t                 fEleGsfTrackChi2OverNdof;
    UInt_t                  fEleKFTrackNHits;
    Float_t                 fEleKFTrackChi2OverNDoF;
    Float_t                 fEleR9;

    Float_t                 fEleSeedEMaxOverE;
    Float_t                 fEleSeedETopOverE;
    Float_t                 fEleSeedEBottomOverE;
    Float_t                 fEleSeedELeftOverE;
    Float_t                 fEleSeedERightOverE;
    Float_t                 fEleSeedE2ndOverE;
    Float_t                 fEleSeedE2x5RightOverE;
    Float_t                 fEleSeedE2x5LeftOverE;
    Float_t                 fEleSeedE2x5TopOverE;
    Float_t                 fEleSeedE2x5BottomOverE;
    Float_t                 fEleSeedE2x5MaxOverE;
    Float_t                 fEleSeedE1x3OverE;
    Float_t                 fEleSeedE3x1OverE;
    Float_t                 fEleSeedE1x5OverE;
    Float_t                 fEleSeedE2x2OverE;
    Float_t                 fEleSeedE3x2OverE;
    Float_t                 fEleSeedE3x3OverE;
    Float_t                 fEleSeedE4x4OverE;
    Float_t                 fEleSeedE5x5OverE;
    Float_t                 fEleE1x5OverE5x5;
    Float_t                 fEleE3x3OverE5x5;
    Float_t                 fEleE2x5MaxOverE5x5;


    //Isolation Variables
    Float_t                 fEleStandardLikelihood;
    Float_t                 fElePFMVA;
    Float_t                 fEleChargedIso03; 
    Float_t                 fEleNeutralHadronIso03; 
    Float_t                 fEleGammaIso03; 
    Float_t                 fEleChargedIso04; 
    Float_t                 fEleNeutralHadronIso04; 
    Float_t                 fEleGammaIso04; 
    Float_t                 fEleChargedIso04FromOtherVertices; 
    Float_t                 fEleNeutralHadronIso04_10Threshold; 
    Float_t                 fEleGammaIso04_10Threshold; 
    Float_t                 fEleTrkIso03; 
    Float_t                 fEleEMIso03; 
    Float_t                 fEleHadIso03; 
    Float_t                 fEleTrkIso04; 
    Float_t                 fEleEMIso04; 
    Float_t                 fEleHadIso04; 
    Float_t                 fRho; 
    UInt_t                  fNVertices; 

    UInt_t                  fEleTriggerBit;
    Bool_t                  fElePassDenominator;
    Bool_t                  fElePassDenominatorSmurf;

    //More Isolation Variables
    Float_t fChargedIso_DR0p0To0p1;
    Float_t fChargedIso_DR0p1To0p2;
    Float_t fChargedIso_DR0p2To0p3;
    Float_t fChargedIso_DR0p3To0p4;
    Float_t fChargedIso_DR0p4To0p5;
    Float_t fChargedIso_DR0p5To0p7;
    Float_t fChargedIso_DR0p7To1p0;
    Float_t fChargedIso_MeanEta;
    Float_t fChargedIso_MeanPhi;
    Float_t fChargedIso_SigmaEtaEta;
    Float_t fChargedIso_SigmaPhiPhi;
    Float_t fChargedIso_SigmaEtaPhi;
    Float_t fGammaIso_DR0p0To0p1;
    Float_t fGammaIso_DR0p1To0p2;
    Float_t fGammaIso_DR0p2To0p3;
    Float_t fGammaIso_DR0p3To0p4;
    Float_t fGammaIso_DR0p4To0p5;
    Float_t fGammaIso_DR0p5To0p7;
    Float_t fGammaIso_DR0p7To1p0;
    Float_t fGammaIso_MeanEta;
    Float_t fGammaIso_MeanPhi;
    Float_t fGammaIso_SigmaEtaEta;
    Float_t fGammaIso_SigmaPhiPhi;
    Float_t fGammaIso_SigmaEtaPhi;
    Float_t fNeutralHadronIso_DR0p0To0p1;
    Float_t fNeutralHadronIso_DR0p1To0p2;
    Float_t fNeutralHadronIso_DR0p2To0p3;
    Float_t fNeutralHadronIso_DR0p3To0p4;
    Float_t fNeutralHadronIso_DR0p4To0p5;
    Float_t fNeutralHadronIso_DR0p5To0p7;
    Float_t fNeutralHadronIso_DR0p7To1p0;
    Float_t fNeutralHadronIso_MeanEta;
    Float_t fNeutralHadronIso_MeanPhi;
    Float_t fNeutralHadronIso_SigmaEtaEta;
    Float_t fNeutralHadronIso_SigmaPhiPhi;
    Float_t fNeutralHadronIso_SigmaEtaPhi;
    Float_t fDirectionalChargedIso;
    Float_t fDirectionalGammaIso;
    Float_t fDirectionalNeutralHadronIso;

 
    //*****************************************************************************************
    //Prepare Output Tree
    //*****************************************************************************************
    TFile *tmpfile = new TFile(inputFile.c_str(), "READ");
    TTree *tmptree = (TTree*)tmpfile->Get("Electrons");


    TFile *EleOutputFile = new TFile(outputFile.c_str(), "RECREATE");
    TTree *EleOutputTree = tmptree->CloneTree(-1, "fast");

    tmpfile->Close();
    delete tmpfile;

    //Add new branches for MVA output
    Float_t fEleIDIsoMVA_BDTG_V0;
    Float_t fEleIDIsoMVA_BDTG_V1;
    Float_t fEleIDIsoMVA_BDTG_V2;
    Float_t fEleIDIsoMVA_BDTG_V3;
    Float_t fEleIDIsoMVA_BDTG_V4;
    Float_t fEleIDIsoMVA_BDTG_V5;
    TBranch* branchEleIDIsoMVA_BDTG_V0 = 0;
    TBranch* branchEleIDIsoMVA_BDTG_V1 = 0;
    TBranch* branchEleIDIsoMVA_BDTG_V2 = 0;
    TBranch* branchEleIDIsoMVA_BDTG_V3 = 0;
    TBranch* branchEleIDIsoMVA_BDTG_V4 = 0;
    TBranch* branchEleIDIsoMVA_BDTG_V5 = 0;
    if(Use["BDTCat_BDTG_V0"])  branchEleIDIsoMVA_BDTG_V0 = EleOutputTree->Branch( "EleIDIsoMVA_BDTG_V0" , &fEleIDIsoMVA_BDTG_V0, "EleIDIsoMVA_BDTG_V0/F" );
    if(Use["BDTCat_BDTG_V1"])  branchEleIDIsoMVA_BDTG_V1 = EleOutputTree->Branch( "EleIDIsoMVA_BDTG_V1" , &fEleIDIsoMVA_BDTG_V1, "EleIDIsoMVA_BDTG_V1/F" );
    if(Use["BDTCat_BDTG_V2"])  branchEleIDIsoMVA_BDTG_V2 = EleOutputTree->Branch( "EleIDIsoMVA_BDTG_V2" , &fEleIDIsoMVA_BDTG_V2, "EleIDIsoMVA_BDTG_V2/F" );
    if(Use["BDTCat_BDTG_V3"])  branchEleIDIsoMVA_BDTG_V3 = EleOutputTree->Branch( "EleIDIsoMVA_BDTG_V3" , &fEleIDIsoMVA_BDTG_V3, "EleIDIsoMVA_BDTG_V3/F" );
    if(Use["BDTCat_BDTG_V4"])  branchEleIDIsoMVA_BDTG_V4 = EleOutputTree->Branch( "EleIDIsoMVA_BDTG_V4" , &fEleIDIsoMVA_BDTG_V4, "EleIDIsoMVA_BDTG_V4/F" );
    if(Use["BDTCat_BDTG_V5"])  branchEleIDIsoMVA_BDTG_V5 = EleOutputTree->Branch( "EleIDIsoMVA_BDTG_V5" , &fEleIDIsoMVA_BDTG_V5, "EleIDIsoMVA_BDTG_V5/F" );

 

    //*****************************************************************************************
    //Prepare Input Tree
    //*****************************************************************************************
    TFile *EleFile = new TFile(inputFile.c_str(), "READ");
    TTree *eleTree = (TTree*)EleFile->Get("Electrons");

    eleTree->SetBranchAddress("weight",&fWeight);
    eleTree->SetBranchAddress("run",&fRunNumber);
    eleTree->SetBranchAddress("lumi",&fLumiSectionNumber);
    eleTree->SetBranchAddress("event",&fEventNumber);
    eleTree->SetBranchAddress("EventNumberParity",&fEleEventNumberParity);
    eleTree->SetBranchAddress("pt",&fElePt);
    eleTree->SetBranchAddress("eta",&fEleEta);
    eleTree->SetBranchAddress("phi",&fElePhi);
    eleTree->SetBranchAddress("ecaldriven",&fEleIsEcalDriven);
    eleTree->SetBranchAddress("combPFIsoHWW",&fElePFIso);
  
    //Conversion Variables
    eleTree->SetBranchAddress("matchConv",&fEleMatchedConversion);
    eleTree->SetBranchAddress("dist",&fEleConvDist);
    eleTree->SetBranchAddress("dcot",&fEleConvDCot);
    eleTree->SetBranchAddress("missHits",&fEleNMissHits);

    //CutBased Variables
    eleTree->SetBranchAddress("see",&fEleSigmaIEtaIEta);
    eleTree->SetBranchAddress("deta",&fEleDEtaIn);
    eleTree->SetBranchAddress("dphi",&fEleDPhiIn);
    eleTree->SetBranchAddress("HoE",&fEleHoverE);
    eleTree->SetBranchAddress("d0",&fEleD0);
    eleTree->SetBranchAddress("dz",&fEleDZ);
    eleTree->SetBranchAddress("fbrem",&fEleFBrem);
    eleTree->SetBranchAddress("EoP",&fEleEOverP);

    //Additional Vars used in Likelihood
    eleTree->SetBranchAddress("EoPout",&fEleESeedClusterOverPout);
    eleTree->SetBranchAddress("spp",&fEleSigmaIPhiIPhi);
    eleTree->SetBranchAddress("nbrems",&fEleNBrem);
    eleTree->SetBranchAddress("IoEmIoP",&fEleOneOverEMinusOneOverP);
    eleTree->SetBranchAddress("EoPin",&fEleESeedClusterOverPIn);
    eleTree->SetBranchAddress("ip3d",&fEleIP3d);
    eleTree->SetBranchAddress("ip3ds",&fEleIP3dSig);

    eleTree->SetBranchAddress("detacalo",&fEledEtaCalo);


    eleTree->SetBranchAddress("dphicalo",&fEledPhiCalo);
    eleTree->SetBranchAddress("PreShowerOverRaw",&fElePreShowerOverRaw);
    eleTree->SetBranchAddress("sep",&fEleSigmaIEtaIPhi);
    eleTree->SetBranchAddress("etawidth",&fEleSCEtaWidth);
    eleTree->SetBranchAddress("phiwidth",&fEleSCPhiWidth);
    eleTree->SetBranchAddress("gsfchi2",&fEleGsfTrackChi2OverNdof);
    eleTree->SetBranchAddress("kfhits",&fEleKFTrackNHits);
    eleTree->SetBranchAddress("kfchi2",&fEleKFTrackChi2OverNDoF);
    eleTree->SetBranchAddress("R9",&fEleR9);

    eleTree->SetBranchAddress("SeedEMaxOverE",&fEleSeedEMaxOverE);
    eleTree->SetBranchAddress("SeedETopOverE",&fEleSeedETopOverE);
    eleTree->SetBranchAddress("SeedEBottomOverE",&fEleSeedEBottomOverE);
    eleTree->SetBranchAddress("SeedELeftOverE",&fEleSeedELeftOverE);
    eleTree->SetBranchAddress("SeedERightOverE",&fEleSeedERightOverE);
    eleTree->SetBranchAddress("SeedE2ndOverE",&fEleSeedE2ndOverE);
    eleTree->SetBranchAddress("SeedE2x5RightOverE",&fEleSeedE2x5RightOverE);
    eleTree->SetBranchAddress("SeedE2x5LeftOverE",&fEleSeedE2x5LeftOverE);
    eleTree->SetBranchAddress("SeedE2x5TopOverE",&fEleSeedE2x5TopOverE);
    eleTree->SetBranchAddress("SeedE2x5BottompassConversionVetoOverE",&fEleSeedE2x5BottomOverE);
    eleTree->SetBranchAddress("SeedE2x5MaxOverE",&fEleSeedE2x5MaxOverE);
    eleTree->SetBranchAddress("SeedE1x3OverE",&fEleSeedE1x3OverE);
    eleTree->SetBranchAddress("SeedE3x1OverE",&fEleSeedE3x1OverE);
    eleTree->SetBranchAddress("SeedE1x5OverE",&fEleSeedE1x5OverE);
    eleTree->SetBranchAddress("SeedE2x2OverE",&fEleSeedE2x2OverE);
    eleTree->SetBranchAddress("SeedE3x2OverE",&fEleSeedE3x2OverE);
    eleTree->SetBranchAddress("SeedE3x3OverE",&fEleSeedE3x3OverE);
    eleTree->SetBranchAddress("SeedE4x4OverE",&fEleSeedE4x4OverE);
    eleTree->SetBranchAddress("SeedE5x5OverE",&fEleSeedE5x5OverE);
    eleTree->SetBranchAddress("e1x5e5x5",&fEleE1x5OverE5x5);
    eleTree->SetBranchAddress("s9s25",&fEleE3x3OverE5x5);
    eleTree->SetBranchAddress("E2x5MaxOverE5x5",&fEleE2x5MaxOverE5x5);

    //Isolation Variables
    eleTree->SetBranchAddress("PFMVA",&fElePFMVA);
    eleTree->SetBranchAddress("chPFIso03",&fEleChargedIso03);
    eleTree->SetBranchAddress("neuPFIso03",&fEleNeutralHadronIso03);
    eleTree->SetBranchAddress("phoPFIso03",&fEleGammaIso03);
    eleTree->SetBranchAddress("chPFIso04",&fEleChargedIso04);
    eleTree->SetBranchAddress("neuPFIso04",&fEleNeutralHadronIso04);
    eleTree->SetBranchAddress("phoPFIso04",&fEleGammaIso04);

    eleTree->SetBranchAddress("trkIso03",&fEleTrkIso03);
    eleTree->SetBranchAddress("ecalIso03",&fEleEMIso03);
    eleTree->SetBranchAddress("hcalIso03",&fEleHadIso03);
    eleTree->SetBranchAddress("trkIso04",&fEleTrkIso04);
    eleTree->SetBranchAddress("ecalIso04",&fEleEMIso04);
    eleTree->SetBranchAddress("hcalIso04",&fEleHadIso04);
    eleTree->SetBranchAddress("rho",&fRho);
    eleTree->SetBranchAddress("vertices",&fNVertices);

    eleTree->SetBranchAddress("triggerBit",&fEleTriggerBit);
    eleTree->SetBranchAddress("DenomFake",&fElePassDenominator);
    eleTree->SetBranchAddress("DenomFakeSmurf",&fElePassDenominatorSmurf);

    //More Isolation Variables
    eleTree->SetBranchAddress("ChargedIso_DR0p0To0p1",&fChargedIso_DR0p0To0p1);
    eleTree->SetBranchAddress("ChargedIso_DR0p1To0p2",&fChargedIso_DR0p1To0p2);
    eleTree->SetBranchAddress("ChargedIso_DR0p2To0p3",&fChargedIso_DR0p2To0p3);
    eleTree->SetBranchAddress("ChargedIso_DR0p3To0p4",&fChargedIso_DR0p3To0p4);
    eleTree->SetBranchAddress("ChargedIso_DR0p4To0p5",&fChargedIso_DR0p4To0p5);
    eleTree->SetBranchAddress("ChargedIso_DR0p5To0p7",&fChargedIso_DR0p5To0p7);
    eleTree->SetBranchAddress("ChargedIso_DR0p7To1p0",&fChargedIso_DR0p7To1p0);
    eleTree->SetBranchAddress("GammaIso_DR0p0To0p1",&fGammaIso_DR0p0To0p1);
    eleTree->SetBranchAddress("GammaIso_DR0p1To0p2",&fGammaIso_DR0p1To0p2);
    eleTree->SetBranchAddress("GammaIso_DR0p2To0p3",&fGammaIso_DR0p2To0p3);
    eleTree->SetBranchAddress("GammaIso_DR0p3To0p4",&fGammaIso_DR0p3To0p4);
    eleTree->SetBranchAddress("GammaIso_DR0p4To0p5",&fGammaIso_DR0p4To0p5);
    eleTree->SetBranchAddress("GammaIso_DR0p5To0p7",&fGammaIso_DR0p5To0p7);
    eleTree->SetBranchAddress("GammaIso_DR0p7To1p0",&fGammaIso_DR0p7To1p0);
    eleTree->SetBranchAddress("NeutralHadronIso_DR0p0To0p1",&fNeutralHadronIso_DR0p0To0p1);
    eleTree->SetBranchAddress("NeutralHadronIso_DR0p1To0p2",&fNeutralHadronIso_DR0p1To0p2);
    eleTree->SetBranchAddress("NeutralHadronIso_DR0p2To0p3",&fNeutralHadronIso_DR0p2To0p3);
    eleTree->SetBranchAddress("NeutralHadronIso_DR0p3To0p4",&fNeutralHadronIso_DR0p3To0p4);
    eleTree->SetBranchAddress("NeutralHadronIso_DR0p4To0p5",&fNeutralHadronIso_DR0p4To0p5);
    eleTree->SetBranchAddress("NeutralHadronIso_DR0p5To0p7",&fNeutralHadronIso_DR0p5To0p7);
    eleTree->SetBranchAddress("NeutralHadronIso_DR0p7To1p0",&fNeutralHadronIso_DR0p7To1p0);
    eleTree->SetBranchAddress("ChargedIso_MeanEta",&fChargedIso_MeanEta);
    eleTree->SetBranchAddress("ChargedIso_MeanPhi",&fChargedIso_MeanPhi);
    eleTree->SetBranchAddress("ChargedIso_SigmaEtaEta",&fChargedIso_SigmaEtaEta);
    eleTree->SetBranchAddress("ChargedIso_SigmaPhiPhi",&fChargedIso_SigmaPhiPhi);
    eleTree->SetBranchAddress("ChargedIso_SigmaEtaPhi",&fChargedIso_SigmaEtaPhi);
    eleTree->SetBranchAddress("GammaIso_MeanEta",&fGammaIso_MeanEta);
    eleTree->SetBranchAddress("GammaIso_MeanPhi",&fGammaIso_MeanPhi);
    eleTree->SetBranchAddress("GammaIso_SigmaEtaEta",&fGammaIso_SigmaEtaEta);
    eleTree->SetBranchAddress("GammaIso_SigmaPhiPhi",&fGammaIso_SigmaPhiPhi);
    eleTree->SetBranchAddress("GammaIso_SigmaEtaPhi",&fGammaIso_SigmaEtaPhi);
    eleTree->SetBranchAddress("NeutralHadronIso_MeanEta",&fNeutralHadronIso_MeanEta);
    eleTree->SetBranchAddress("NeutralHadronIso_MeanPhi",&fNeutralHadronIso_MeanPhi);
    eleTree->SetBranchAddress("NeutralHadronIso_SigmaEtaEta",&fNeutralHadronIso_SigmaEtaEta);
    eleTree->SetBranchAddress("NeutralHadronIso_SigmaPhiPhi",&fNeutralHadronIso_SigmaPhiPhi);
    eleTree->SetBranchAddress("NeutralHadronIso_SigmaEtaPhi",&fNeutralHadronIso_SigmaEtaPhi);
    eleTree->SetBranchAddress("DirectionalChargedIso",&fDirectionalChargedIso);
    eleTree->SetBranchAddress("DirectionalGammaIso",&fDirectionalGammaIso);
    eleTree->SetBranchAddress("DirectionalNeutralHadronIso",&fDirectionalNeutralHadronIso);




    // --- Event loop
    std::cout << "--- Processing: " << eleTree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();

    for (Long64_t ievt=0; ievt<eleTree->GetEntries();ievt++) {

      if (ievt%10000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      eleTree->GetEntry(ievt);

      Double_t rho = 0;
      if (!(TMath::IsNaN(fRho) || isinf(fRho))) rho = fRho;

      //--------------------------------------------------------
      // important: here we associate branches to MVA variables
      //--------------------------------------------------------
 
      varfbrem = fEleFBrem;
      vardeta = fEleDEtaIn;
      vardphi = fEleDPhiIn;
      varsee = fEleSigmaIEtaIEta;
      varetawidth = fEleSCEtaWidth;
      varphiwidth = fEleSCPhiWidth;
      varHoE = fEleHoverE;
      varEoP = fEleEOverP;
      vare1x5e5x5 = fEleE1x5OverE5x5;
      varEoPout = fEleESeedClusterOverPout;
      vareleEoPout = fEleESeedClusterOverPout;
      vardetacalo = fEledEtaCalo;
      varkfchi2 = fEleKFTrackChi2OverNDoF;
      varkfhits = fEleKFTrackNHits;
      varkfhitsall = fEleKFTrackNHits;
      varspp = fEleSigmaIPhiIPhi;
      varIoEmIoP = fEleOneOverEMinusOneOverP;
      varR9 = fEleR9;
      vargsfchi2 = fEleGsfTrackChi2OverNdof;
      varPreShowerOverRaw = fElePreShowerOverRaw;
      varChargedIso_DR0p0To0p1 	                = fChargedIso_DR0p0To0p1;
      varChargedIso_DR0p1To0p2 	                = fChargedIso_DR0p1To0p2;
      varChargedIso_DR0p2To0p3	        = fChargedIso_DR0p2To0p3;
      varChargedIso_DR0p3To0p4	        = fChargedIso_DR0p3To0p4;
      varChargedIso_DR0p4To0p5	        = fChargedIso_DR0p4To0p5;
      varGammaIso_DR0p0To0p1 	                = fGammaIso_DR0p0To0p1;
      varGammaIso_DR0p1To0p2 	                = fGammaIso_DR0p1To0p2;
      varGammaIso_DR0p2To0p3	                = fGammaIso_DR0p2To0p3;
      varGammaIso_DR0p3To0p4	                = fGammaIso_DR0p3To0p4;
      varGammaIso_DR0p4To0p5	                = fGammaIso_DR0p4To0p5;
      varNeutralHadronIso_DR0p0To0p1 	        = fNeutralHadronIso_DR0p0To0p1;
      varNeutralHadronIso_DR0p1To0p2 	        = fNeutralHadronIso_DR0p1To0p2;
      varNeutralHadronIso_DR0p2To0p3	        = fNeutralHadronIso_DR0p2To0p3;
      varNeutralHadronIso_DR0p3To0p4	        = fNeutralHadronIso_DR0p3To0p4;
      varNeutralHadronIso_DR0p4To0p5	        = fNeutralHadronIso_DR0p4To0p5;
      varPt = fElePt;
      varEta = fEleEta;
      //varEventNumberParity = fEventNumberParity;
      

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

//         cout << varPt << " " << varEta << " : " << mvabin << " : 1\n"; 
        fEleIDIsoMVA_BDTG_V0            = fTMVAReader[mvabin]->EvaluateMVA( "ElectronIDIsoMVA_V0_BDTG method" );
        branchEleIDIsoMVA_BDTG_V0->Fill();
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

//         cout << varPt << " " << varEta << " : " << mvabin << " : 1\n"; 
        fEleIDIsoMVA_BDTG_V1            = fTMVAReader2[mvabin]->EvaluateMVA( "ElectronIDIsoMVA_V1_BDTG method" );
        branchEleIDIsoMVA_BDTG_V1->Fill();
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
