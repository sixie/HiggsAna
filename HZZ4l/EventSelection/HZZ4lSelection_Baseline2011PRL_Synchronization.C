//root -l EWKAna/HZZ4l/Selection/HZZ4lSelection.C+\(\"\"\)

//================================================================================================
//
// HZZ4l selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

// define structures to read in ntuple
#include "HiggsAna/Ntupler/interface/HiggsAnaDefs.hh"
#include "HiggsAna/Ntupler/interface/TEventInfo.hh"
#include "HiggsAna/Ntupler/interface/TElectron.hh"
#include "HiggsAna/Ntupler/interface/TMuon.hh"
#include "HiggsAna/Ntupler/interface/TJet.hh"
#include "HiggsAna/Ntupler/interface/TGenParticle.hh"
#include "HiggsAna/Ntupler/interface/TPFCandidate.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"
#include "HiggsAna/HZZ4l/Utils/LeptonSelection.hh"
#include "HiggsAna/Utils/LeptonIDCuts.hh"

// output data structs
#include "HiggsAna/HZZ4l/interface/HZZKinematics.hh"
#include "HiggsAna/HZZ4l/interface/HZZGenInfo.hh"
#include "HiggsAna/HZZ4l/interface/HZZ4lDefs.hh"

// #include "TMVAGui.C"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================

//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwNtuplerMod");
    if (!dir) {
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}




//=== MAIN MACRO =================================================================================================

void HZZ4lSelection_BaselineSynchronization(const string Label = "", Int_t Option = 0) 
{  
  gBenchmark->Start("HZZTemplate");
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  Bool_t printDebug = kFALSE;
  Double_t lumi = 1092;              // luminosity (pb^-1)


  vector<vector<string> > inputFiles;
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/data/blue/sixie/ntuples/HZZ4l/data/HZZNtuple_r11a-data.1092ipb.FourRecoLeptonSkim.root");

  inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-h115zz4l-gf-v14b-pu_noskim_0000.root");
   inputFiles.back().push_back("/home/sixie/hist/HZZ4lNtuples_NEW/mc/AllNtuple_HZZ4lNtuple_f11-h140zz4l-gf-v14b-pu_noskim_0000.root");
//  inputFiles.back().push_back("/home/sixie/hist/HZZ4lNtuples/mc/AllNtuple_HZZ4lNtuple_f11-h130zz4l-gf-v14b-pu_noskim_0000.root");
  //weight = 16.63*0.000169 / (301479) = 9.32227452e-9 / pb^-1
  //weight(4.6fb^-1) = 4.29e-5

  vector<string> processNames;
  processNames.push_back("HZZ120");
//   processNames.push_back("Data");
//   processNames.push_back("ZZ");
//   processNames.push_back("WZ");
//   processNames.push_back("ttbar");
//   processNames.push_back("WW"); 

  assert(processNames.size() == inputFiles.size());

  //--------------------------------------------------------------------------------------------------------------
  // Yields
  //==============================================================================================================  
  Double_t NGenAccepted_4E = 0;
  Double_t NGenAccepted_4M = 0;
  Double_t NGenAccepted_2E2M = 0;

  Double_t NBaselineRecoID_4E = 0;
  Double_t NBaselineRecoID_4M = 0;
  Double_t NBaselineRecoID_2E2M = 0;

  Double_t NBaselineTrkIsoLoose_4E = 0;
  Double_t NBaselineTrkIsoLoose_4M = 0;
  Double_t NBaselineTrkIsoLoose_2E2M = 0;

  Double_t NBaselineZ1_4E = 0;
  Double_t NBaselineZ1_4M = 0;
  Double_t NBaselineZ1_2E2M = 0;

  Double_t NBaselineBestCand_4E = 0;
  Double_t NBaselineBestCand_4M = 0;
  Double_t NBaselineBestCand_2E2M = 0;

  Double_t NBaselinePairwiseIso_4E = 0;
  Double_t NBaselinePairwiseIso_4M = 0;
  Double_t NBaselinePairwiseIso_2E2M = 0;

  Double_t NBaselineSIP_4E = 0;
  Double_t NBaselineSIP_4M = 0;
  Double_t NBaselineSIP_2E2M = 0;

  Double_t NBaselineSelection_4E = 0;
  Double_t NBaselineSelection_4M = 0;
  Double_t NBaselineSelection_2E2M = 0;

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  vector <TH1F*>  fHZZ4lSelection; 
  vector <TH1F*>  fHZZTo4ESelection; 
  vector <TH1F*>  fHZZTo4MuSelection; 
  vector <TH1F*>  fHZZTo2E2MuSelection; 
  vector <TH1F*>  fZ1Mass; 
  vector <TH1F*>  fZ2Mass; 
  vector <TH1F*>  fZZMass; 
  for (int q=0; q<processNames.size() ; ++q) {
    TH1F *tmpHZZ4lSelection = new TH1F(("hHZZ4lSelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHZZTo4ESelection = new TH1F(("hHZZTo4ESelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHZZTo4MuSelection = new TH1F(("hHZZTo4MuSelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHZZTo2E2MuSelection = new TH1F(("hHZZTo2E2MuSelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpZ1Mass = new TH1F(("hZ1Mass"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 40, 0, 200);
    TH1F *tmpZ2Mass = new TH1F(("hZ2Mass"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 40, 0, 200);
    TH1F *tmpZZMass = new TH1F(("hZZMass"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 40, 0, 600);
    tmpHZZ4lSelection->Sumw2();
    tmpHZZTo4ESelection->Sumw2();
    tmpHZZTo4MuSelection->Sumw2();
    tmpHZZTo2E2MuSelection->Sumw2();   
    tmpZ1Mass->Sumw2();   
    tmpZ2Mass->Sumw2();   
    tmpZZMass->Sumw2();   
    
    fHZZ4lSelection.push_back(tmpHZZ4lSelection);
    fHZZTo4ESelection.push_back(tmpHZZTo4ESelection);
    fHZZTo4MuSelection.push_back(tmpHZZTo4MuSelection);
    fHZZTo2E2MuSelection.push_back(tmpHZZTo2E2MuSelection);
    fZ1Mass.push_back(tmpZ1Mass);
    fZ2Mass.push_back(tmpZ2Mass);
    fZZMass.push_back(tmpZZMass);
  }


  vector<string> CutLabel;
  CutLabel.push_back("Dilepton20/10");

  
  //--------------------------------------------------------------------------------------------------------------
  // output ntuple structure
  //==============================================================================================================  
  HZZKinematics kinematics;
  HZZGenInfo    geninfo;

  TFile *fOutputFile = new TFile(("HZZOutput"+label+".root").c_str(), "RECREATE");
  TTree *fOutputTree = new TTree("hzz4l","hzz4l");
  fOutputTree->Branch("geninfo",&geninfo);
  fOutputTree->Branch("kinematics",&kinematics);



  //--------------------------------------------------------------------------------------------------------------
  // Pileup Reweighting
  //==============================================================================================================  
  TFile *fPUFile = TFile::Open("/data/smurf/sixie/Pileup/weights/PileupReweighting.Fall11_To_Full2011.root");
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;
  Int_t EffectiveAreaVersion = -1;  
  EffectiveAreaVersion = 2;  //Fall11MC



  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  ofstream eventList("BaselineSelection.eventList.txt");
  //
  // Access samples and fill histograms
  TFile *inputFile=0;
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  TClonesArray *genparticleArr = new TClonesArray("mithep::TGenParticle");
  TClonesArray *pfcandidateArr = new TClonesArray("mithep::TPFCandidate");
  
  for (int qq = 0; qq<inputFiles.size() ; ++qq) { 
    for (int f = 0; f < inputFiles[qq].size() ; ++f) {
      //********************************************************
      // Get Tree
      //********************************************************
      cout << "Reading File " << inputFiles[qq][f] << endl;
      eventTree = getTreeFromFile(inputFiles[qq][f].c_str(),"Events"); 
      TBranch *infoBr;
      TBranch *electronBr;
      TBranch *muonBr;
      TBranch *jetBr;
      TBranch *genparticleBr;
      TBranch *pfcandidateBr;


      //*****************************************************************************************
      //Loop over muon Data Tree
      //*****************************************************************************************
      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");
      eventTree->SetBranchAddress("GenParticle", &genparticleArr);         genparticleBr = eventTree->GetBranch("GenParticle");
      eventTree->SetBranchAddress("PFCandidate", &pfcandidateArr);         pfcandidateBr = eventTree->GetBranch("PFCandidate");
  
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
        infoBr->GetEntry(ientry);
        if (ientry % 10000 == 0) cout << "Event " << ientry << endl;
	
        //Use Only Testing events
        //if (info->evtNum % 2 == 0 ) continue;
        //Use Only Training events
        //if (info->evtNum % 2 != 0 ) continue;
       

        //********************************************************
        double mynpu = TMath::Min((double)info->nPUEvents,34.999);
        Int_t npuxbin = fhDPU->GetXaxis()->FindBin(mynpu);
        double npuWeight = fhDPU->GetBinContent(npuxbin);
        //********************************************************
        //********************************************************
        // Pileup Energy Density
        //********************************************************
        Double_t rho = 0;
        if (!(TMath::IsNaN(info->PileupEnergyDensity) || isinf(info->PileupEnergyDensity))) rho = info->PileupEnergyDensity;
        

        //********************************************************
        // Printdebug
        //********************************************************
        printDebug = kFALSE;
        if ((0 == 1) 
            || (info->evtNum == 9316)           
          ) printDebug = kTRUE;
        



        //********************************************************
        // Load the branches
        //********************************************************
        electronArr->Clear(); 
        muonArr->Clear(); 
        jetArr->Clear();  //if (info->evtNum % 2 != 0 ) continue;
       
        genparticleArr->Clear(); 
        pfcandidateArr->Clear(); 
        electronBr->GetEntry(ientry);
        muonBr->GetEntry(ientry);
        jetBr->GetEntry(ientry);
        genparticleBr->GetEntry(ientry);
        pfcandidateBr->GetEntry(ientry);

//         Double_t eventweight = info->eventweight * lumi * npuWeight;
        Double_t eventweight = npuWeight;


        //********************************************************
        // GenInfo
        //********************************************************
        vector<Double_t> pdgid;
        vector<Double_t> pt;
        vector<Double_t> eta;
        vector<Double_t> phi;
        vector<Double_t> pt_ele;
        vector<Double_t> eta_ele;
        vector<Double_t> phi_ele;
        vector<Double_t> pt_mu;
        vector<Double_t> eta_mu;
        vector<Double_t> phi_mu;

        Double_t minGenLeptonDR = 9999;
        for(Int_t k=0; k<genparticleArr->GetEntries(); k++) {
          const mithep::TGenParticle *gen = (mithep::TGenParticle*)((*genparticleArr)[k]);


//           if (gen->status == 1 && (abs(gen->pdgid) == 11 || abs(gen->pdgid) == 13 ) ) {
          if ( (abs(gen->pdgid) == 11 || abs(gen->pdgid) == 13) && gen->motherPdgID == 23 ) {
            if (abs(gen->pdgid) == 11) {
              if (!(fabs(gen->eta) < 2.5 )) continue;
            }
            if (abs(gen->pdgid) == 13) {
              if (!(fabs(gen->eta) < 2.4 )) continue;
            }
            pdgid.push_back(gen->pdgid); 
            pt.push_back(gen->pt); 
            eta.push_back(gen->eta); 
            phi.push_back(gen->phi); 
            if (abs(gen->pdgid) == 11) {
              pt_ele.push_back(gen->pt);
              eta_ele.push_back(gen->eta);
              phi_ele.push_back(gen->phi);
            }
            if (abs(gen->pdgid) == 13) {
              pt_mu.push_back(gen->pt);
              eta_mu.push_back(gen->eta);
              phi_mu.push_back(gen->phi);
            }
          }

          //find pair of gen leptons that are closest in DR
          if (!(gen->status == 1 && (abs(gen->pdgid) == 11 || abs(gen->pdgid) == 13))) continue;
          for(Int_t l=k+1; l<genparticleArr->GetEntries(); l++) {
            const mithep::TGenParticle *p2 = (mithep::TGenParticle*)((*genparticleArr)[l]);
            if (!(p2->status == 1 && (abs(p2->pdgid) == 11 || abs(p2->pdgid) == 13))) continue;
            if (minGenLeptonDR > mithep::MathUtils::DeltaR(gen->phi,gen->eta,p2->phi,p2->eta)) {
              minGenLeptonDR = mithep::MathUtils::DeltaR(gen->phi,gen->eta,p2->phi,p2->eta);
            }
          }
        }

        //******************************************************************
        //Remove events with two leptons that are very close together
        //******************************************************************
        //if (minGenLeptonDR < 0.4) continue;


        geninfo.id_1_a = 0;
        geninfo.id_1_b = 0;
        geninfo.id_2_a = 0;
        geninfo.id_2_b = 0;
        geninfo.pt_1_a = 0;
        geninfo.pt_1_b = 0;
        geninfo.pt_2_a = 0;
        geninfo.pt_2_b = 0;
        geninfo.eta_1_a = 0;
        geninfo.eta_1_b = 0;
        geninfo.eta_2_a = 0;
        geninfo.eta_2_b = 0;

        if (pdgid.size() == 4) {
          if (pt_ele.size() == 4) {
            geninfo.id_1_a = 11;
            geninfo.id_1_b = 11;
            geninfo.id_2_a = 11;
            geninfo.id_2_b = 11;
            geninfo.pt_1_a = pt_ele[0];
            geninfo.pt_1_b = pt_ele[1];
            geninfo.pt_2_a = pt_ele[2];
            geninfo.pt_2_b = pt_ele[3];
            geninfo.eta_1_a = eta_ele[0];
            geninfo.eta_1_b = eta_ele[1];
            geninfo.eta_2_a = eta_ele[2];
            geninfo.eta_2_b = eta_ele[3];
          }
          if (pt_mu.size() == 4) {
            geninfo.id_1_a = 13;
            geninfo.id_1_b = 13;
            geninfo.id_2_a = 13;
            geninfo.id_2_b = 13;
            geninfo.pt_1_a = pt_mu[0];
            geninfo.pt_1_b = pt_mu[1];
            geninfo.pt_2_a = pt_mu[2];
            geninfo.pt_2_b = pt_mu[3];
            geninfo.eta_1_a = eta_mu[0];
            geninfo.eta_1_b = eta_mu[1];
            geninfo.eta_2_a = eta_mu[2];
            geninfo.eta_2_b = eta_mu[3];
          }
          if (pt_ele.size() == 2 && pt_mu.size() == 2) {
            geninfo.id_1_a = 11;
            geninfo.id_1_b = 11;
            geninfo.id_2_a = 13;
            geninfo.id_2_b = 13;
            geninfo.pt_1_a = pt_ele[0];
            geninfo.pt_1_b = pt_ele[1];
            geninfo.pt_2_a = pt_mu[0];
            geninfo.pt_2_b = pt_mu[1];
            geninfo.eta_1_a = eta_ele[0];
            geninfo.eta_1_b = eta_ele[1];
            geninfo.eta_2_a = eta_mu[0];
            geninfo.eta_2_b = eta_mu[1];
          }
        }


        eventList << "Event " << info->runNum << " " << info->lumiSec << " " << info->evtNum << " ";
        if (pt_ele.size() == 4) {
          NGenAccepted_4E += eventweight;
          eventList << "4e ";
        } else if (pt_mu.size() == 4) {
          NGenAccepted_4M += eventweight;
          eventList << "4m ";
        } else if (pt_ele.size() == 2 && pt_mu.size() == 2) {
          NGenAccepted_2E2M += eventweight;
          eventList << "2e2m ";
        } else {
          eventList << "none \n";
        }

        

        //********************************************************
        // Met
        //********************************************************
        TVector3 pfMet;        
        if(info->pfMEx!=0 || info->pfMEy!=0) {       
          pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
        }

        //********************************************************
        // Lepton Selection
        //********************************************************
        Int_t NLeptons = 0;
        vector<Int_t> leptonType;
        vector<Int_t> leptonIndex;
        vector<Double_t> leptonPt;
        vector<Double_t> leptonEta;
        vector<Double_t> leptonPhi;
        vector<Int_t> leptonCharge;

        Int_t NJets = 0;
        const mithep::TJet *leadingJet = 0;
    
        for(Int_t i=0; i<muonArr->GetEntries(); i++) {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
          
          if (info->evtNum == 185882 ) {
            cout << mu->pt << " " << mu->eta << " " << mu->phi << " : " << ((mu->typeBits & kGlobal) == kGlobal) << endl;
          }

          if ( (0==0)
               &&
               mu->pt > 5.0
               &&
               fabs(mu->eta) < 2.4
               && ( ((mu->typeBits & kGlobal) == kGlobal) 
                 )
            ) {
             leptonPt.push_back(mu->pt);
             leptonEta.push_back(mu->eta);
             leptonPhi.push_back(mu->phi);

//             //Test PF muon momentum assignment
//             if (mu->pfPt > 0) {
//               leptonPt.push_back(mu->pfPt);
//               leptonEta.push_back(mu->pfEta);
//               leptonPhi.push_back(mu->pfPhi);
//             } else {
//               leptonPt.push_back(mu->pt);
//               leptonEta.push_back(mu->eta);
//               leptonPhi.push_back(mu->phi);              
//             }

            leptonType.push_back(13);
            leptonIndex.push_back(i);  
            leptonCharge.push_back(mu->q);
          }
        }

        for(Int_t i=0; i<electronArr->GetEntries(); i++) {
          const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);

          //Do we do this?
          Bool_t isMuonOverlap = kFALSE;
          for (int k=0; k<leptonPt.size(); ++k) {
            if ( leptonType[k] == 13 
                 && mithep::MathUtils::DeltaR(ele->phi, ele->eta, leptonPhi[k],leptonEta[k]) < 0.1
              ) {
              isMuonOverlap = kTRUE; 
              break;
            }
          }

          if ( (0==0)
               &&
               ele->pt > 7.0
               && 
               fabs(ele->eta) < 2.5
               &&
               !isMuonOverlap   
            ) {
            leptonPt.push_back(ele->pt);
            leptonEta.push_back(ele->eta);
            leptonPhi.push_back(ele->phi);
            leptonType.push_back(11);
            leptonIndex.push_back(i);
            leptonCharge.push_back(ele->q);
          }
        }

        //sort leptons
        Int_t tempType;
        Int_t tempIndex;
        Double_t tempPt;
        Double_t tempEta;
        Double_t tempPhi;
        Int_t tempCharge;
        for (int l=0; l<leptonIndex.size(); l++) {
          for (int k=0; k < leptonIndex.size() - 1; k++) {
            if (leptonPt[k+1] > leptonPt[k]) {
              tempType = leptonType[k];
              tempIndex = leptonIndex[k];
              tempPt = leptonPt[k];
              tempEta = leptonEta[k];
              tempPhi = leptonPhi[k];
              tempCharge = leptonCharge[k];
          
              leptonType[k] = leptonType[k+1];
              leptonIndex[k] = leptonIndex[k+1];
              leptonPt[k] = leptonPt[k+1];
              leptonEta[k] = leptonEta[k+1];
              leptonPhi[k] = leptonPhi[k+1];
              leptonCharge[k] = leptonCharge[k+1];

              leptonType[k+1] = tempType;
              leptonIndex[k+1] = tempIndex;
              leptonPt[k+1] = tempPt;
              leptonEta[k+1] = tempEta;
              leptonPhi[k+1] = tempPhi;
              leptonCharge[k+1] = tempCharge;
          
            }
          }
        }

        //remove 5lepton events just to be not confused for the pairwise isolation
        if (leptonPt.size() < 4) {
          eventList << " 0 0 0 0 0 0 0\n";
          continue; 
        }
//         if (leptonPt.size() > 4) continue; 

        //******************************************************************************
        //How many events have 4 ID'ed leptons?
        //******************************************************************************
        Int_t NLeptonsPassBaselineID = 0;
        Int_t NLeptonsPassTrkIsoLoose = 0;
        for(int i = 0; i < leptonPt.size(); ++i) {
          if (leptonType[i] == 11) {
            const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[leptonIndex[i]]);

            //Baseline CiC
            if (PassCiCID(ele,1)) {
              NLeptonsPassBaselineID++;
            }
            if(ele->trkIso03 / ele->pt < 0.7) {
              NLeptonsPassTrkIsoLoose++;
            }

          } else if (leptonType[i] == 13) {
            const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[leptonIndex[i]]);  
         
            if (passMuonID_HZZ2011(mu)) {
              NLeptonsPassBaselineID++;
            }
            if(passMuonID_HZZ2011(mu) && mu->trkIso03 / mu->pt < 0.7) {
              NLeptonsPassTrkIsoLoose++;
            }
          } else {
            cout << "weird lepton type : " << leptonType[i] << endl;
          }
        }


        //******************************************************************************
        // BaselineID Only
        //******************************************************************************
        if (NLeptonsPassBaselineID >= 4) {
          eventList << "1 ";
          if (pt_ele.size() >= 4) {
            NBaselineRecoID_4E += eventweight;
          } else if (pt_mu.size() >= 4) {
            NBaselineRecoID_4M += eventweight;
          } else if (pt_ele.size() >= 2 && pt_mu.size() >= 2) {
            NBaselineRecoID_2E2M += eventweight;
          }
        } else {
          eventList << "0 ";
        }
 
        //******************************************************************************
        // TrkIsoLoose
        //******************************************************************************
        if (NLeptonsPassTrkIsoLoose >= 4 && NLeptonsPassBaselineID >= 4) {
          eventList << "1 ";
          if (pt_ele.size() >= 4) {
            NBaselineTrkIsoLoose_4E += eventweight;
          } else if (pt_mu.size() >= 4) {
            NBaselineTrkIsoLoose_4M += eventweight;
          } else if (pt_ele.size() >= 2 && pt_mu.size() >= 2) {
            NBaselineTrkIsoLoose_2E2M += eventweight;
          }
        } else {
          eventList << "0 ";
        }

        //******************************************************************************
        //Event Selection
        //Initialize Event Tree
        //******************************************************************************
        kinematics.l1type = 0;
        kinematics.l2type = 0;
        kinematics.l3type = 0;
        kinematics.l4type = 0;
        kinematics.l1pt = 0;
        kinematics.l2pt = 0;
        kinematics.l3pt = 0;
        kinematics.l4pt = 0;
        kinematics.l1eta = 0;
        kinematics.l2eta = 0;
        kinematics.l3eta = 0;
        kinematics.l4eta = 0;
        kinematics.l1phi = 0;
        kinematics.l2phi = 0;
        kinematics.l3phi = 0;
        kinematics.l4phi = 0;
        kinematics.Z1pt = 0;
        kinematics.Z2pt = 0;
        kinematics.ZZpt = 0;
        kinematics.Z1eta = 0;
        kinematics.Z2eta = 0;
        kinematics.ZZeta = 0;
        kinematics.mZ1 = 0; 
        kinematics.mZ2 = 0; 
        kinematics.m4l = 0;
        kinematics.channel = 0;


        Bool_t passSelection = kTRUE;

        Bool_t passZ1Selection = kFALSE;
        Bool_t passBestCandSelection = kFALSE;
        Bool_t PassPairwiseIso = kTRUE;
        Bool_t passBaselineSelection = kFALSE;
        Bool_t passSIPSelection = kFALSE;


        //trick to break execution if event doesn't satisfy certain requirements
        for (Int_t q=0; q<=0;q++) {

          //******************************************************************************
          //Z1 Selection
          //******************************************************************************
          if (leptonPt.size() < 2) break;

          Int_t Z1LeptonPlusIndex = -1;
          Int_t Z1LeptonMinusIndex = -1;
          Double_t BestZ1Mass = -1;

          for(int i = 0; i < leptonPt.size(); ++i) {
            for(int j = i+1; j < leptonPt.size(); ++j) {

              if (leptonType[i] == 11) {
                const mithep::TElectron *ele1 = (mithep::TElectron*)((*electronArr)[leptonIndex[i]]);  
                if (!(PassCiCID(ele1,1) && ele1->trkIso03/ele1->pt < 0.7)) continue;
              } else if (leptonType[i] == 13) {
                const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr)[leptonIndex[i]]);  
                if (!(passMuonID_HZZ2011(mu1) && mu1->trkIso03 / mu1->pt < 0.7)) continue;
              }
              if (leptonType[j] == 11) {
                const mithep::TElectron *ele2 = (mithep::TElectron*)((*electronArr)[leptonIndex[j]]);  
                if (!(PassCiCID(ele2,1) && ele2->trkIso03/ele2->pt < 0.7)) continue;
              } else if (leptonType[j] == 13) {
                const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr)[leptonIndex[j]]);  
                if (!(passMuonID_HZZ2011(mu2) && mu2->trkIso03 / mu2->pt < 0.7)) continue;
              }
                        
              if (!(leptonPt[i] > 20.0 || leptonPt[j] > 20.0)) continue; //one lepton must be 20 GeV
              if (!(leptonPt[i] > 10.0 && leptonPt[j] > 10.0)) continue; //both leptons must be 10 GeV
              if (leptonCharge[i] == leptonCharge[j]) continue;          //require opp sign
              if (fabs(leptonType[i]) != fabs(leptonType[j])) continue;  //require same flavor

              //Make Z1 hypothesis
              mithep::FourVectorM leptonPlus;
              mithep::FourVectorM leptonMinus;

              Double_t Lepton1Mass = 0.51099892e-3;
              if (leptonType[i] != 11) Lepton1Mass = 105.658369e-3;
              Double_t Lepton2Mass = 0.51099892e-3;
              if (leptonType[j] != 11) Lepton2Mass = 105.658369e-3;

              if (leptonCharge[i] > 0 ) {
                leptonPlus.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], Lepton1Mass );
                leptonMinus.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], Lepton2Mass );
              } else {
                leptonPlus.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], Lepton2Mass );
                leptonMinus.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], Lepton1Mass );
              }
              mithep::FourVectorM dilepton = leptonPlus+leptonMinus;
            
              if (BestZ1Mass < 0 || fabs(dilepton.M() - 91.1876) < fabs(BestZ1Mass - 91.1876)) {
                if (dilepton.M() > 50) {
                  BestZ1Mass = dilepton.M();
                  if (leptonCharge[i] > 0) {
                    Z1LeptonPlusIndex = i;
                    Z1LeptonMinusIndex = j;
                  } else {
                    Z1LeptonPlusIndex = j;
                    Z1LeptonMinusIndex = i;
                  }
                }
              }
            }
          }
        
          if (Z1LeptonPlusIndex == -1) break;

          //******************************************************************************
          //Make Z1 Candidate FourVector
          //******************************************************************************
          mithep::FourVectorM Z1LeptonPlus;
          mithep::FourVectorM Z1LeptonMinus;
          if (leptonType[Z1LeptonPlusIndex] == 11) {
            Z1LeptonPlus.SetCoordinates(leptonPt[Z1LeptonPlusIndex], leptonEta[Z1LeptonPlusIndex], leptonPhi[Z1LeptonPlusIndex], 0.51099892e-3 );
          } else {
            Z1LeptonPlus.SetCoordinates(leptonPt[Z1LeptonPlusIndex], leptonEta[Z1LeptonPlusIndex], leptonPhi[Z1LeptonPlusIndex], 105.658369e-3 );
          }
          if (leptonType[Z1LeptonMinusIndex] == 11) {
            Z1LeptonMinus.SetCoordinates(leptonPt[Z1LeptonMinusIndex], leptonEta[Z1LeptonMinusIndex], leptonPhi[Z1LeptonMinusIndex], 0.51099892e-3 );
          } else {
            Z1LeptonMinus.SetCoordinates(leptonPt[Z1LeptonMinusIndex], leptonEta[Z1LeptonMinusIndex], leptonPhi[Z1LeptonMinusIndex], 105.658369e-3 );
          }
          mithep::FourVectorM Z1Candidate = Z1LeptonPlus+Z1LeptonMinus;
          passZ1Selection = kTRUE;
          

          //******************************************************************************
          //Z2 Selection
          //******************************************************************************
          Int_t Z2LeptonPlusIndex = -1;
          Int_t Z2LeptonMinusIndex = -1;
          Double_t BestZ2Mass = -1;
          for(int i = 0; i < leptonPt.size(); ++i) {
            for(int j = i+1; j < leptonPt.size(); ++j) {

              if (leptonType[i] == 11) {
                const mithep::TElectron *ele1 = (mithep::TElectron*)((*electronArr)[leptonIndex[i]]);  
                if (!(PassCiCID(ele1,1) && ele1->trkIso03/ele1->pt < 0.7)) continue;
              } else if (leptonType[i] == 13) {
                const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr)[leptonIndex[i]]);  
                if (!(passMuonID_HZZ2011(mu1) && mu1->trkIso03 / mu1->pt < 0.7)) continue;
              }
              if (leptonType[j] == 11) {
                const mithep::TElectron *ele2 = (mithep::TElectron*)((*electronArr)[leptonIndex[j]]);  
                if (!(PassCiCID(ele2,1) && ele2->trkIso03/ele2->pt < 0.7)) continue;
              } else if (leptonType[j] == 13) {
                const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr)[leptonIndex[j]]);  
                if (!(passMuonID_HZZ2011(mu2) && mu2->trkIso03 / mu2->pt < 0.7)) continue;
              }

               if (i == Z1LeptonPlusIndex || i == Z1LeptonMinusIndex) continue; //skip Z1 leptons
              if (j == Z1LeptonPlusIndex || j == Z1LeptonMinusIndex) continue; //skip Z1 leptons
              if (leptonCharge[i] == leptonCharge[j]) continue;         //require opp sign
              if (fabs(leptonType[i]) != fabs(leptonType[j])) continue; //require same flavor
            
              //Make Z2 hypothesis
              mithep::FourVectorM leptonPlus;
              mithep::FourVectorM leptonMinus;

              Double_t Lepton1Mass = 0.51099892e-3;
              if (leptonType[i] != 11) Lepton1Mass = 105.658369e-3;
              Double_t Lepton2Mass = 0.51099892e-3;
              if (leptonType[j] != 11) Lepton2Mass = 105.658369e-3;

              if (leptonCharge[i] > 0 ) {
                leptonPlus.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], Lepton1Mass );
                leptonMinus.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], Lepton2Mass );
              } else {
                leptonPlus.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], Lepton2Mass );
                leptonMinus.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], Lepton1Mass );
              }

              mithep::FourVectorM dilepton = leptonPlus+leptonMinus;
              mithep::FourVectorM fourLepton = Z1Candidate + dilepton;


              if (!(dilepton.M() > 12.0)) continue;
              if (!(fourLepton.M() > 100.0)) continue;

              //for 4e and 4mu, require at least 1 of the other opp sign lepton pairs have mass > 12
              if (fabs(leptonType[i]) == fabs(leptonType[Z1LeptonPlusIndex])) {
                mithep::FourVectorM pair1 = Z1LeptonPlus+leptonMinus;
                mithep::FourVectorM pair2 = Z1LeptonMinus+leptonPlus;
                if (!(pair1.M() > 12 || pair2.M() > 12)) continue;
              }
 
              //Disambiguiation is done by choosing the pair with the largest ptMax, and largest ptMin
              if (Z2LeptonPlusIndex < 0) {
                if (leptonCharge[i] > 0) {
                  Z2LeptonPlusIndex = i;
                  Z2LeptonMinusIndex = j;
                } else {
                  Z2LeptonPlusIndex = j;
                  Z2LeptonMinusIndex = i;
                }
              } else {
                Double_t BestPairPtMax = leptonPt[Z2LeptonPlusIndex];               
                Double_t BestPairPtMin = leptonPt[Z2LeptonMinusIndex]; 
                if (leptonPt[Z2LeptonMinusIndex] > BestPairPtMax) {
                  BestPairPtMax = leptonPt[Z2LeptonMinusIndex];
                  BestPairPtMin = leptonPt[Z2LeptonPlusIndex];
                }

                Double_t CurrentPairPtMax = leptonPt[i];               
                Double_t CurrentPairPtMin = leptonPt[j]; 
                if (leptonPt[j] > CurrentPairPtMax) {
                  CurrentPairPtMax = leptonPt[j];
                  CurrentPairPtMin = leptonPt[i];
                }

                if (CurrentPairPtMax > BestPairPtMax) {
                  if (leptonCharge[i] > 0) {
                    Z2LeptonPlusIndex = i;
                    Z2LeptonMinusIndex = j;
                  } else {
                    Z2LeptonPlusIndex = j;
                    Z2LeptonMinusIndex = i;
                  }
                } else if (CurrentPairPtMax  == BestPairPtMax) {
                  if (CurrentPairPtMin > BestPairPtMin) {
                    if (leptonCharge[i] > 0) {
                      Z2LeptonPlusIndex = i;
                      Z2LeptonMinusIndex = j;
                    } else {
                      Z2LeptonPlusIndex = j;
                      Z2LeptonMinusIndex = i;
                    }                  
                  }
                }
              }            
            }
          }
  
          if (Z2LeptonPlusIndex == -1) break;

          //******************************************************************************
          //Make Z2 Candidate FourVector
          //******************************************************************************
          mithep::FourVectorM Z2LeptonPlus;
          mithep::FourVectorM Z2LeptonMinus;
          if (leptonType[Z2LeptonPlusIndex] == 11) {
            Z2LeptonPlus.SetCoordinates(leptonPt[Z2LeptonPlusIndex], leptonEta[Z2LeptonPlusIndex], leptonPhi[Z2LeptonPlusIndex], 0.51099892e-3 );
          } else {
            Z2LeptonPlus.SetCoordinates(leptonPt[Z2LeptonPlusIndex], leptonEta[Z2LeptonPlusIndex], leptonPhi[Z2LeptonPlusIndex], 105.658369e-3 );
          }
          if (leptonType[Z2LeptonMinusIndex] == 11) {
            Z2LeptonMinus.SetCoordinates(leptonPt[Z2LeptonMinusIndex], leptonEta[Z2LeptonMinusIndex], leptonPhi[Z2LeptonMinusIndex], 0.51099892e-3 );
          } else {
            Z2LeptonMinus.SetCoordinates(leptonPt[Z2LeptonMinusIndex], leptonEta[Z2LeptonMinusIndex], leptonPhi[Z2LeptonMinusIndex], 105.658369e-3 );
          }
          mithep::FourVectorM Z2Candidate = Z2LeptonPlus+Z2LeptonMinus;
          mithep::FourVectorM ZZSystem = Z1Candidate + Z2Candidate;
          if (Z2Candidate.M() > 12) {
            passBestCandSelection = kTRUE;
          }

          //***************************************************************
          // remaining kinematic cuts 
          //***************************************************************
          if (Z1Candidate.M() > 120) break;
          if (Z2Candidate.M() > 120) break;
          if (Z2Candidate.M() < 12) break;
          passBaselineSelection = kTRUE;
    
          //***************************************************************
          // Define Channel Type
          //***************************************************************             
          Int_t Channel = -1;
          if (leptonType[Z1LeptonPlusIndex] == 11 && leptonType[Z2LeptonPlusIndex] == 11 ) {
            Channel = kFourEle;
          } else if (leptonType[Z1LeptonPlusIndex] == 13 && leptonType[Z2LeptonPlusIndex] == 13 ) {
            Channel = kFourMu;
          } else if (leptonType[Z1LeptonPlusIndex] == 11 && leptonType[Z2LeptonPlusIndex] == 13 ) {
            Channel = kTwoEleTwoMu;
          } else if (leptonType[Z1LeptonPlusIndex] == 13 && leptonType[Z2LeptonPlusIndex] == 11 ) {
            Channel = kTwoMuTwoEle;
          }
        


          //******************************************************************************
          //Fill Best Candidate Leptons
          //******************************************************************************
          if (info->evtNum == 189969) {

            cout << "Best Candidate: " << endl;
            cout << Z1LeptonPlusIndex << " " << leptonType[Z1LeptonPlusIndex] << " " << leptonIndex[Z1LeptonPlusIndex] << " " << leptonPt[Z1LeptonPlusIndex] 
                 << " " << leptonEta[Z1LeptonPlusIndex] << " " << leptonPhi[Z1LeptonPlusIndex] << " " << leptonCharge[Z1LeptonPlusIndex] << endl;
            cout << Z1LeptonMinusIndex << " " << leptonType[Z1LeptonMinusIndex] << " " << leptonIndex[Z1LeptonMinusIndex] << " " << leptonPt[Z1LeptonMinusIndex] 
                 << " " << leptonEta[Z1LeptonMinusIndex] << " " << leptonPhi[Z1LeptonMinusIndex] << " " << leptonCharge[Z1LeptonMinusIndex] << endl;
            cout << Z2LeptonPlusIndex << " " << leptonType[Z2LeptonPlusIndex] << " " << leptonIndex[Z2LeptonPlusIndex] << " " << leptonPt[Z2LeptonPlusIndex] 
                 << " " << leptonEta[Z2LeptonPlusIndex] << " " << leptonPhi[Z2LeptonPlusIndex] << " " << leptonCharge[Z2LeptonPlusIndex] << endl;
            cout << Z2LeptonMinusIndex << " " << leptonType[Z2LeptonMinusIndex] << " " << leptonIndex[Z2LeptonMinusIndex] << " " << leptonPt[Z2LeptonMinusIndex] 
                 << " " << leptonEta[Z2LeptonMinusIndex] << " " << leptonPhi[Z2LeptonMinusIndex] << " " << leptonCharge[Z2LeptonMinusIndex] << endl;

          }

          vector<Int_t> CandidateLeptonType;
          vector<Int_t> CandidateLeptonIndex;
          vector<Double_t> CandidateLeptonPt;
          vector<Double_t> CandidateLeptonEta;
          vector<Double_t> CandidateLeptonPhi;
          vector<Int_t> CandidateLeptonCharge;
          CandidateLeptonType.push_back(leptonType[Z1LeptonPlusIndex]);
          CandidateLeptonIndex.push_back(leptonIndex[Z1LeptonPlusIndex]);
          CandidateLeptonPt.push_back(leptonPt[Z1LeptonPlusIndex]);
          CandidateLeptonEta.push_back(leptonEta[Z1LeptonPlusIndex]);
          CandidateLeptonPhi.push_back(leptonPhi[Z1LeptonPlusIndex]);
          CandidateLeptonCharge.push_back(leptonCharge[Z1LeptonPlusIndex]);
          
          CandidateLeptonType.push_back(leptonType[Z1LeptonMinusIndex]);
          CandidateLeptonIndex.push_back(leptonIndex[Z1LeptonMinusIndex]);
          CandidateLeptonPt.push_back(leptonPt[Z1LeptonMinusIndex]);
          CandidateLeptonEta.push_back(leptonEta[Z1LeptonMinusIndex]);
          CandidateLeptonPhi.push_back(leptonPhi[Z1LeptonMinusIndex]);
          CandidateLeptonCharge.push_back(leptonCharge[Z1LeptonMinusIndex]);
          
          CandidateLeptonType.push_back(leptonType[Z2LeptonPlusIndex]);
          CandidateLeptonIndex.push_back(leptonIndex[Z2LeptonPlusIndex]);
          CandidateLeptonPt.push_back(leptonPt[Z2LeptonPlusIndex]);
          CandidateLeptonEta.push_back(leptonEta[Z2LeptonPlusIndex]);
          CandidateLeptonPhi.push_back(leptonPhi[Z2LeptonPlusIndex]);
          CandidateLeptonCharge.push_back(leptonCharge[Z2LeptonPlusIndex]);
          
          CandidateLeptonType.push_back(leptonType[Z2LeptonMinusIndex]);
          CandidateLeptonIndex.push_back(leptonIndex[Z2LeptonMinusIndex]);
          CandidateLeptonPt.push_back(leptonPt[Z2LeptonMinusIndex]);
          CandidateLeptonEta.push_back(leptonEta[Z2LeptonMinusIndex]);
          CandidateLeptonPhi.push_back(leptonPhi[Z2LeptonMinusIndex]);
          CandidateLeptonCharge.push_back(leptonCharge[Z2LeptonMinusIndex]);
 
          //******************************************************************************
          //Pairwise Iso Selection
          //******************************************************************************
          for(int i = 0; i < CandidateLeptonPt.size(); ++i) {
            Double_t iso1 = 0;
            if (CandidateLeptonType[i] == 11) {
              const mithep::TElectron *ele1 = (mithep::TElectron*)((*electronArr)[CandidateLeptonIndex[i]]);
              if (!(PassCiCID(ele1,1) && ele1->trkIso03 / ele1->pt < 0.7)) continue;
              Double_t DetIso03PUCorrection = 0;
              
              //****    Effective Area   *************************************************
              //
              //For Data
//             if (fabs(ele1->scEta) < 1.5) DetIso03PUCorrection = rho * (0.078 + 0.026);
//             else DetIso03PUCorrection = rho * (0.046 + 0.072);
              //For MC ( i think this is summer11 though)
              if (fabs(ele1->scEta) < 1.5) DetIso03PUCorrection = rho * (0.101 + 0.021);
              else DetIso03PUCorrection = rho * (0.046 + 0.040);
              //
              //**************************************************************************
              
              iso1 = (ele1->trkIso03 + ele1->emIso03 + ele1->hadIso03 - DetIso03PUCorrection) / ele1->pt;
            } else if (CandidateLeptonType[i] == 13) {
              const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr)[CandidateLeptonIndex[i]]);   
              if (!(passMuonID_HZZ2011(mu1) && mu1->trkIso03 / mu1->pt < 0.7)) continue;
              Double_t DetIso03PUCorrection = 0;

              //****    Effective Area   *************************************************
              //
//             //For Data
//             if (fabs(mu1->eta) < 1.5) DetIso03PUCorrection = rho * (0.087 + 0.042);
//             else DetIso03PUCorrection = rho * (0.049 + 0.059);
              //For MC
              if (fabs(mu1->eta) < 1.5) DetIso03PUCorrection = rho * (0.074 + 0.022);
              else DetIso03PUCorrection = rho * (0.045 + 0.030);
              //
              //**************************************************************************

              iso1 = (mu1->trkIso03 + mu1->emIso03 + mu1->hadIso03 - DetIso03PUCorrection) / mu1->pt;
              if (info->evtNum == 148698) cout << i << " : " << mu1->pt << " " << mu1->eta << " " << mu1->phi << " " << iso1 << endl;
            } else {
              cout << "weird lepton type : " << leptonType[i] << endl;
            }
 
            for(int j = i+1; j < CandidateLeptonPt.size(); ++j) {
              Double_t iso2 = 0;
              if (CandidateLeptonType[j] == 11) {
                const mithep::TElectron *ele2 = (mithep::TElectron*)((*electronArr)[CandidateLeptonIndex[j]]);
                if (!(PassCiCID(ele2,1) && ele2->trkIso03 / ele2->pt < 0.7)) continue;
                Double_t DetIso03PUCorrection = 0;

                //****    Effective Area   *************************************************
                //
                //For Data
//             if (fabs(ele2->scEta) < 1.5) DetIso03PUCorrection = rho * (0.078 + 0.026);
//             else DetIso03PUCorrection = rho * (0.046 + 0.072);
                //For MC ( i think this is summer11 though)
                if (fabs(ele2->scEta) < 1.5) DetIso03PUCorrection = rho * (0.101 + 0.021);
                else DetIso03PUCorrection = rho * (0.046 + 0.040);
                //
                //**************************************************************************

                iso2 = (ele2->trkIso03 + ele2->emIso03 + ele2->hadIso03 - DetIso03PUCorrection) / ele2->pt;
              } else if (CandidateLeptonType[j] == 13) {
                const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr)[CandidateLeptonIndex[j]]);   
                if (!(passMuonID_HZZ2011(mu2) && mu2->trkIso03 / mu2->pt < 0.7)) continue;
                Double_t DetIso03PUCorrection = 0;

                //****    Effective Area   *************************************************
                //
//             //For Data
//             if (fabs(mu2->eta) < 1.5) DetIso03PUCorrection = rho * (0.087 + 0.042);
//             else DetIso03PUCorrection = rho * (0.049 + 0.059);
                //For MC
                if (fabs(mu2->eta) < 1.5) DetIso03PUCorrection = rho * (0.074 + 0.022);
                else DetIso03PUCorrection = rho * (0.045 + 0.030);
                //
                //**************************************************************************

                iso2 = (mu2->trkIso03 + mu2->emIso03 + mu2->hadIso03 - DetIso03PUCorrection) / mu2->pt;
                if (info->evtNum == 148698) cout << j << " : " << mu2->pt << " " << mu2->eta << " " << mu2->phi << " " << iso2 << endl;
              } else {
                cout << "weird lepton type : " << CandidateLeptonType[j] << endl;
              }

              if (info->evtNum == 148698) {
                cout << i << "," << j << iso1 << " " << iso2 << " " << iso1+iso2 << endl;
              } 

              if (iso1 + iso2 > 0.35) PassPairwiseIso = kFALSE;
            }
          }

          if (info->evtNum == 148698) {
            cout << "PairwiseIso: " << PassPairwiseIso << endl;
          }


          //***************************************************************
          // SIP Selection
          //***************************************************************             
          Int_t NLeptonsPassSIP = 0;
          for(int i = 0; i < CandidateLeptonPt.size(); ++i) {
            if (CandidateLeptonType[i] == 11) {
              const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[CandidateLeptonIndex[i]]);
              
              if (info->evtNum == 189969) {
                cout << "Ele " << ele->pt << " " << ele->eta << " " << ele->phi << " : " << ele->ip3dSig << endl;
              }

              if (PassCiCID(ele,1) && ele->trkIso03 / ele->pt < 0.7 && fabs(ele->ip3dSig) < 4) {
                NLeptonsPassSIP++;
              }
              
            } else if (CandidateLeptonType[i] == 13) {
              const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[CandidateLeptonIndex[i]]);  
              
              if (info->evtNum == 189969) {
                cout << "Mu " << mu->pt << " " << mu->eta << " " << mu->phi << " : " << mu->ip3dSig << endl;
              }

              if (passMuonID_HZZ2011(mu) && mu->trkIso03 / mu->pt < 0.7 && passMuonIP_HZZ2011(mu)) {
                NLeptonsPassSIP++;
              }
            } else {
              cout << "weird lepton type : " << CandidateLeptonType[i] << endl;
            }          
          }
          if (NLeptonsPassSIP >= 4) passSIPSelection = kTRUE;

          //***************************************************************
          // Fill Kinematics Ntuple
          //***************************************************************             
          kinematics.l1type = leptonType[Z1LeptonPlusIndex];
          kinematics.l2type = leptonType[Z1LeptonMinusIndex];
          kinematics.l3type = leptonType[Z2LeptonPlusIndex];
          kinematics.l4type = leptonType[Z2LeptonMinusIndex];
          kinematics.l1pt = leptonPt[Z1LeptonPlusIndex];
          kinematics.l2pt = leptonPt[Z1LeptonMinusIndex];
          kinematics.l3pt = leptonPt[Z2LeptonPlusIndex];
          kinematics.l4pt = leptonPt[Z2LeptonMinusIndex];
          kinematics.l1eta = leptonEta[Z1LeptonPlusIndex];
          kinematics.l2eta = leptonEta[Z1LeptonPlusIndex];
          kinematics.l3eta = leptonEta[Z2LeptonMinusIndex];
          kinematics.l4eta = leptonEta[Z2LeptonMinusIndex];
          kinematics.l1phi = leptonEta[Z1LeptonPlusIndex];
          kinematics.l2phi = leptonEta[Z1LeptonMinusIndex];
          kinematics.l3phi = leptonEta[Z2LeptonPlusIndex];
          kinematics.l4phi = leptonEta[Z2LeptonMinusIndex];
          kinematics.Z1pt = Z1Candidate.Pt();
          kinematics.Z2pt = Z2Candidate.Pt();
          kinematics.ZZpt = ZZSystem.Pt();
          kinematics.Z1eta = Z1Candidate.Eta();
          kinematics.Z2eta = Z2Candidate.Eta();
          kinematics.ZZeta = ZZSystem.Eta();
          kinematics.mZ1 = Z1Candidate.M(); 
          kinematics.mZ2 = Z2Candidate.M(); 
          kinematics.m4l = ZZSystem.M();
          kinematics.channel = Channel;
        } // dummy loop over current one event


        //***************************************************************
        // Z1 Selection
        //***************************************************************             

          // BaselineID Only
        if (passZ1Selection 
            && NLeptonsPassTrkIsoLoose >= 4 && NLeptonsPassBaselineID >= 4) {
          eventList << "1 ";
          if (pt_ele.size() >= 4) {
            NBaselineZ1_4E += eventweight;
          } else if (pt_mu.size() >= 4) {
            NBaselineZ1_4M += eventweight;
          } else if (pt_ele.size() >= 2 && pt_mu.size() >= 2) {
            NBaselineZ1_2E2M += eventweight;
          }
        } else {
          eventList << "0 ";
        }

        if (passBestCandSelection 
            && passZ1Selection && NLeptonsPassTrkIsoLoose >= 4 && NLeptonsPassBaselineID >= 4) {
          eventList << "1 ";
          if (pt_ele.size() >= 4) {
            NBaselineBestCand_4E += eventweight;
          } else if (pt_mu.size() >= 4) {
            NBaselineBestCand_4M += eventweight;
          } else if (pt_ele.size() >= 2 && pt_mu.size() >= 2) {
            NBaselineBestCand_2E2M += eventweight;
          }
        } else {
          eventList << "0 ";
        }


        if (PassPairwiseIso 
            && passBestCandSelection && passZ1Selection && NLeptonsPassTrkIsoLoose >= 4 && NLeptonsPassBaselineID >= 4) {
          eventList << "1 ";
          if (pt_ele.size() >= 4) {
            NBaselinePairwiseIso_4E += eventweight;
          } else if (pt_mu.size() >= 4) {
            NBaselinePairwiseIso_4M += eventweight;
          } else if (pt_ele.size() >= 2 && pt_mu.size() >= 2) {
            NBaselinePairwiseIso_2E2M += eventweight;
          }
        } else {
          eventList << "0 ";
        }


        if (passSIPSelection 
            && PassPairwiseIso && passBestCandSelection 
            && passZ1Selection && NLeptonsPassTrkIsoLoose >= 4 && NLeptonsPassBaselineID >= 4) {
          eventList << "1 ";
          if (pt_ele.size() >= 4) {
            NBaselineSIP_4E += eventweight;
          } else if (pt_mu.size() >= 4) {
            NBaselineSIP_4M += eventweight;
          } else if (pt_ele.size() >= 2 && pt_mu.size() >= 2) {
            NBaselineSIP_2E2M += eventweight;
          }
        } else {
          eventList << "0 ";
        }


        if ( passBaselineSelection 
             && passSIPSelection
             && PassPairwiseIso && passBestCandSelection 
             && passZ1Selection && NLeptonsPassTrkIsoLoose >= 4 && NLeptonsPassBaselineID >= 4) {
          eventList << "1 ";
          if (pt_ele.size() >= 4) {
            NBaselineSelection_4E += eventweight;
          } else if (pt_mu.size() >= 4) {
            NBaselineSelection_4M += eventweight;
          } else if (pt_ele.size() >= 2 && pt_mu.size() >= 2) {
            NBaselineSelection_2E2M += eventweight;
          }
        } else {
          eventList << "0 ";
        }

        eventList << endl;


        //Fill Tree
        if (Option == 0) {
          fOutputTree->Fill();
        } 


      } //end loop over data     
    }
  }

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  cout << "Number of GenAccepted Events: " << NGenAccepted_4E << " " << NGenAccepted_4M << " " << NGenAccepted_2E2M << endl;
  cout << "Number of BaselineRecoID Events: " << NBaselineRecoID_4E << " " << NBaselineRecoID_4M << " " << NBaselineRecoID_2E2M << endl;
  cout << "Number of TrkIsoLoose Events: " << NBaselineTrkIsoLoose_4E << " " << NBaselineTrkIsoLoose_4M << " " << NBaselineTrkIsoLoose_2E2M << endl;
  cout << "Number of NBaselineZ1 Events: " << NBaselineZ1_4E << " " << NBaselineZ1_4M << " " << NBaselineZ1_2E2M << endl;
  cout << "Number of NBaselineBestCand Events: " << NBaselineBestCand_4E << " " << NBaselineBestCand_4M << " " << NBaselineBestCand_2E2M << endl;
  cout << "Number of NBaselinePairwiseIso Events: " << NBaselinePairwiseIso_4E << " " << NBaselinePairwiseIso_4M << " " << NBaselinePairwiseIso_2E2M << endl;
  cout << "Number of NBaselineSIP Events: " << NBaselineSIP_4E << " " << NBaselineSIP_4M << " " << NBaselineSIP_2E2M << endl;
  cout << "Number of NBaselineSelection Events: " << NBaselineSelection_4E << " " << NBaselineSelection_4M << " " << NBaselineSelection_2E2M << endl;
  cout << endl;


  cout << "Efficiency of BaselineRecoID: " << NBaselineRecoID_4E/NGenAccepted_4E << " " << NBaselineRecoID_4M/NGenAccepted_4M << " " << NBaselineRecoID_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of NBaselineTrkIsoLoose: " << NBaselineTrkIsoLoose_4E/NGenAccepted_4E << " " << NBaselineTrkIsoLoose_4M/NGenAccepted_4M << " " << NBaselineTrkIsoLoose_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of NBaselineZ1: " << NBaselineZ1_4E/NGenAccepted_4E << " " << NBaselineZ1_4M/NGenAccepted_4M << " " << NBaselineZ1_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of NBaselineBestCand: " << NBaselineBestCand_4E/NGenAccepted_4E << " " << NBaselineBestCand_4M/NGenAccepted_4M << " " << NBaselineBestCand_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of NBaselinePairwiseIso: " << NBaselinePairwiseIso_4E/NGenAccepted_4E << " " << NBaselinePairwiseIso_4M/NGenAccepted_4M << " " << NBaselinePairwiseIso_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of NBaselineSIP: " << NBaselineSIP_4E/NGenAccepted_4E << " " << NBaselineSIP_4M/NGenAccepted_4M << " " << NBaselineSIP_2E2M/NGenAccepted_2E2M << endl;
  cout << "Efficiency of NBaselineSelection: " << NBaselineSelection_4E/NGenAccepted_4E << " " << NBaselineSelection_4M/NGenAccepted_4M << " " << NBaselineSelection_2E2M/NGenAccepted_2E2M << endl;



  fOutputFile->cd();
  fOutputFile->Write();
  fOutputFile->Close();
  

  
  gBenchmark->Show("WWTemplate");       
} 



