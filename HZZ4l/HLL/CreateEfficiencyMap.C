//root -l HiggsAna/HZZ4l/HLL/CreateEfficiencyMap.C+'("HZZEfficiencyMap_HZZ125.root","HZZ125",0)'


//================================================================================================
//
// Create Efficiency Map
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

// data structs
#include "HiggsAna/HZZ4l/interface/HZZEfficiencyMap.hh"

#endif
 

//=== FUNCTION DECLARATIONS ======================================================================================

void initialize3DArray(double ***array, UInt_t NPtBins, UInt_t NEtaBins, UInt_t NPhiBins) {

  array = new double**[NPtBins+2];
  for (uint i=0; i < NPtBins+2; ++i) {
    array[i] = new double*[NEtaBins+2];
    for (uint j=0; j < NEtaBins+2; ++j) {
      array[i][j] = new double[NPhiBins+2];
      for (uint k=0; k < NPhiBins+2; ++k) {
        array[i][j][k] = 0;
      }      
    }
  }
}

void initialize3DArray(vector<vector<vector<double> > > &array, UInt_t NPtBins, UInt_t NEtaBins, UInt_t NPhiBins) {

  array.resize(NPtBins+2);
  for (uint i=0; i < NPtBins+2; ++i) {
    array[i].resize(NEtaBins+2);
    for (uint j=0; j < NEtaBins+2; ++j) {
      array[i][j].resize(NPhiBins+2);
      for (uint k=0; k < NPhiBins+2; ++k) {
        array[i][j][k] = 0;
      }
    }
  }
}

void initialize2DArray(vector<vector<double> > &array, UInt_t NPtBins, UInt_t NEtaBins) {

  array.resize(NPtBins+2);
  for (uint i=0; i < NPtBins+2; ++i) {
    array[i].resize(NEtaBins+2);
    for (uint j=0; j < NEtaBins+2; ++j) {
      array[i][j]= 0;
    }
  }
}
  
UInt_t FindBin( double value, double bins[], UInt_t nbins) {

  UInt_t nbinboundaries = nbins+1;
  UInt_t bin = 0;
  for (uint i=0; i < nbinboundaries; ++i) {
    if (i < nbinboundaries-1) {
      if (value >= bins[i] && value < bins[i+1]) {
        bin = i+1;
        break;
      }
    } else if (i == nbinboundaries-1) {
      if (value >= bins[i]) {
        bin = nbinboundaries;
        break;
      }
    }    
  }
  return bin;
}
  
void computeEfficiencyPtEtaPhi(vector<vector<vector<double> > > &numerator, 
                       vector<vector<vector<double> > > &denominator,
                       vector<vector<vector<double> > > &eff
  ) {

  for (uint i=0; i < numerator.size(); ++i) {
    for (uint j=0; j < numerator[i].size(); ++j) {
      for (uint k=0; k < numerator[i][j].size(); ++k) {
        if (denominator[i][j][k] > 0) {         
          eff[i][j][k] = numerator[i][j][k] / denominator[i][j][k];        
        } else {
          eff[i][j][k] = 0;
        }
      }
    }
  }
}

void computeEfficiencyPtEta(vector<vector<double> > &numerator, 
                       vector<vector<double> > &denominator,
                       vector<vector<double> > &eff
  ) {

  for (uint i=0; i < numerator.size(); ++i) {
    for (uint j=0; j < numerator[i].size(); ++j) {
      if ( denominator[i][j] > 0) {
        eff[i][j] = numerator[i][j] / denominator[i][j];
      } else {
        eff[i][j] = 0;
      }
    }
  }
}


//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = TFile::Open(infname,"read");
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

void CreateEfficiencyMap(const string filename, const string Label = "ZZ", Int_t Option = 0) 
{  
  gBenchmark->Start("HZZTemplate");
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================


//   //--------------------------------------------------------------------------------------------------------------
//   // Pileup Reweighting
//   //==============================================================================================================  
//   TFile *fPUFile = TFile::Open("/data/smurf/sixie/Pileup/weights/PileupReweighting.Fall11_To_Full2011.root");
//   TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
//   assert(fhDPU);
//   fhDPU->SetDirectory(0);
//   delete fPUFile;


  //********************************************************
  // Create Arrays to store the map
  //********************************************************
   const UInt_t NPtBins = 15; 
   const UInt_t NEtaBins = 16;
   const UInt_t NPhiBins = 12;
   double ptBins[NPtBins+1] = { 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 50};
   double etaBins[NEtaBins+1] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.566, 1.8, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6};
   double phiBins[NPhiBins+1] = { -3.2, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5,  2, 2.5, 3.2 };

//   double ***NDenominator_PtEtaPhi;
//   double ***NNumerator_PtEtaPhi;
//   double ***Efficiency_PtEtaPhi;
//   double **NDenominator_PtEta;
//   double **NNumerator_PtEta;
//   double **Efficiency_PtEta;

   vector<vector<vector<double> > > NDenominator_Electrons_PtEtaPhi;
   vector<vector<vector<double> > > NNumerator_Electrons_PtEtaPhi;
   vector<vector<vector<double> > > Efficiency_Electrons_PtEtaPhi;
   vector<vector<double> > NDenominator_Electrons_PtEta;
   vector<vector<double> > NNumerator_Electrons_PtEta;
   vector<vector<double> > Efficiency_Electrons_PtEta;
   vector<vector<vector<double> > > NDenominator_Muons_PtEtaPhi;
   vector<vector<vector<double> > > NNumerator_Muons_PtEtaPhi;
   vector<vector<vector<double> > > Efficiency_Muons_PtEtaPhi;
   vector<vector<double> > NDenominator_Muons_PtEta;
   vector<vector<double> > NNumerator_Muons_PtEta;
   vector<vector<double> > Efficiency_Muons_PtEta;
   
   initialize3DArray(NDenominator_Electrons_PtEtaPhi, NPtBins, NEtaBins, NPhiBins);
   initialize3DArray(NNumerator_Electrons_PtEtaPhi, NPtBins, NEtaBins, NPhiBins);
   initialize3DArray(Efficiency_Electrons_PtEtaPhi, NPtBins, NEtaBins, NPhiBins);
   initialize2DArray(NDenominator_Electrons_PtEta, NPtBins, NEtaBins);
   initialize2DArray(NNumerator_Electrons_PtEta, NPtBins, NEtaBins);
   initialize2DArray(Efficiency_Electrons_PtEta, NPtBins, NEtaBins);
   initialize3DArray(NDenominator_Muons_PtEtaPhi, NPtBins, NEtaBins, NPhiBins);
   initialize3DArray(NNumerator_Muons_PtEtaPhi, NPtBins, NEtaBins, NPhiBins);
   initialize3DArray(Efficiency_Muons_PtEtaPhi, NPtBins, NEtaBins, NPhiBins);
   initialize2DArray(NDenominator_Muons_PtEta, NPtBins, NEtaBins);
   initialize2DArray(NNumerator_Muons_PtEta, NPtBins, NEtaBins);
   initialize2DArray(Efficiency_Muons_PtEta, NPtBins, NEtaBins);
  
   
  //--------------------------------------------------------------------------------------------------------------
  // Read efficiency map ntuple
  //==============================================================================================================  
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  HZZEfficiencyMap *efficiencyMap = new HZZEfficiencyMap();

  //********************************************************
  // Get Tree
  //********************************************************
  cout << "Reading File " << filename << endl;
  eventTree = getTreeFromFile(filename.c_str(),"EfficiencyMap"); 

  TBranch *efficiencyMapBr;

  //*****************************************************************************************
  //Loop over tree
  //*****************************************************************************************
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("efficiencyMap",       &efficiencyMap);      
  efficiencyMapBr = eventTree->GetBranch("efficiencyMap");
  
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    efficiencyMapBr->GetEntry(ientry);
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
    
    //********************************************************
    //double mynpu = TMath::Min((double)info->nPUEvents,34.999);
    //Int_t npuxbin = fhDPU->GetXaxis()->FindBin(mynpu);
    //double npuWeight = fhDPU->GetBinContent(npuxbin);
    //********************************************************
    double weight = 1.0;

    Int_t tmpPtBin = FindBin( efficiencyMap->genpt , ptBins, NPtBins);
    Int_t tmpEtaBin = FindBin( fabs(efficiencyMap->geneta) , etaBins, NEtaBins);
    Int_t tmpPhiBin = FindBin( efficiencyMap->genphi , phiBins, NPhiBins);

//     cout << fabs(efficiencyMap->geneta) << " " << tmpEtaBin << endl;

    if (abs(efficiencyMap->type) == 11) {
      NDenominator_Electrons_PtEtaPhi[tmpPtBin][tmpEtaBin][tmpPhiBin] += weight;
      NDenominator_Electrons_PtEta[tmpPtBin][tmpEtaBin] += weight;
      if (efficiencyMap->pass == kTRUE && efficiencyMap->recopt > 7) {
        NNumerator_Electrons_PtEtaPhi[tmpPtBin][tmpEtaBin][tmpPhiBin] += weight;
        NNumerator_Electrons_PtEta[tmpPtBin][tmpEtaBin] += weight;      
      }
    } 
    if (abs(efficiencyMap->type) == 13) {
      NDenominator_Muons_PtEtaPhi[tmpPtBin][tmpEtaBin][tmpPhiBin] += weight;
      NDenominator_Muons_PtEta[tmpPtBin][tmpEtaBin] += weight;
      if (efficiencyMap->pass == kTRUE && efficiencyMap->recopt > 5) {
        NNumerator_Muons_PtEtaPhi[tmpPtBin][tmpEtaBin][tmpPhiBin] += weight;
        NNumerator_Muons_PtEta[tmpPtBin][tmpEtaBin] += weight;      
      }
    } 


  } //end loop over data     

  computeEfficiencyPtEtaPhi(NNumerator_Electrons_PtEtaPhi, NDenominator_Electrons_PtEtaPhi, Efficiency_Electrons_PtEtaPhi);
  computeEfficiencyPtEta(NNumerator_Electrons_PtEta, NDenominator_Electrons_PtEta, Efficiency_Electrons_PtEta);
  computeEfficiencyPtEtaPhi(NNumerator_Muons_PtEtaPhi, NDenominator_Muons_PtEtaPhi, Efficiency_Muons_PtEtaPhi);
  computeEfficiencyPtEta(NNumerator_Muons_PtEta, NDenominator_Muons_PtEta, Efficiency_Muons_PtEta);

//   for (uint i=0; i < NPtBins+2; ++i) {
//     for (uint j=0; j < NEtaBins+2; ++j) {
//       for (uint k=0; k < NPhiBins+2; ++k) {
//         cout << i << " " << j << " " << k << " : " << Efficiency_Electrons_PtEtaPhi[i][j][k] << endl;
//       }
//     }
//   }


  //********************************************************
  // Produce output lookup table
  //******************************************************** 
  ofstream outf_e("ElectronEfficiencyMap.h");

  outf_e << "UInt_t FindElectronEfficiencyBin( double value, double bins[], UInt_t nbins) {" << endl;
  outf_e << "  UInt_t nbinboundaries = nbins+1;" << endl;
  outf_e << "  UInt_t bin = 0;" << endl;
  outf_e << "  for (uint i=0; i < nbinboundaries; ++i) {" << endl;
  outf_e << "    if (i < nbinboundaries-1) {" << endl;
  outf_e << "      if (value >= bins[i] && value < bins[i+1]) {" << endl;
  outf_e << "        bin = i+1;" << endl;
  outf_e << "        break;" << endl;
  outf_e << "      }" << endl;
  outf_e << "    } else if (i == nbinboundaries-1) {" << endl;
  outf_e << "      if (value >= bins[i]) {" << endl;
  outf_e << "        bin = nbinboundaries;" << endl;
  outf_e << "        break;" << endl;
  outf_e << "      }" << endl;
  outf_e << "    }    " << endl;
  outf_e << "  }" << endl;
  outf_e << "  return bin;" << endl;
  outf_e << "}" << endl;

  outf_e << endl;
  outf_e << endl;

  outf_e << "Double_t GetElectronEfficiencyPtEtaPhi(Double_t Pt, Double_t Eta, Double_t Phi) {" << endl;

  outf_e << endl;
  outf_e << "  Double_t ptBins[" << NPtBins+1 << "] = {";
  for (uint i=0; i < NPtBins+1; ++i) {
    outf_e << ptBins[i];
    if (i < NPtBins) {
      outf_e << ",";
    }
  }
  outf_e << "};\n";

  outf_e << "  Double_t etaBins[" << NEtaBins+1 << "] = {";
  for (uint i=0; i < NEtaBins+1; ++i) {
    outf_e << etaBins[i];
    if (i < NEtaBins) {
      outf_e << ",";
    }
  }
  outf_e << "};\n";

  outf_e << "  Double_t phiBins[" << NPhiBins+1 << "] = {";
  for (uint i=0; i < NPhiBins+1; ++i) {
    outf_e << phiBins[i];
    if (i < NPhiBins) {
      outf_e << ",";
    }
  }
  outf_e << "};\n";

  outf_e << endl;
  outf_e << endl;

  outf_e << "  Double_t Efficiency[" << NPtBins+2 << "][" << NEtaBins+2 << "][" << NPhiBins+2 << "]  = {";
  outf_e << endl;

  for (uint i=0; i < NPtBins+2; ++i) {
    outf_e << "    {" << endl;
    for (uint j=0; j < NEtaBins+2; ++j) {
      outf_e << "      {";
      for (uint k=0; k < NPhiBins+2; ++k) {
        outf_e << Efficiency_Electrons_PtEtaPhi[i][j][k];
        if (k< NPhiBins+1) {
          outf_e << ",";
        }        
      }
      if (j< NEtaBins+1) {
        outf_e << "},";
      } else {
        outf_e << "}";
      }
      outf_e << endl;
    }
    if (i< NPtBins+1) {
      outf_e << "    },";
    } else {
      outf_e << "}";
    }
    outf_e << endl;
  }

  outf_e << "};" << endl;

  outf_e << endl;
  outf_e << endl;

  outf_e << "  Int_t tmpPtBin = FindElectronEfficiencyBin( Pt , ptBins, " << NPtBins << ");" << endl;
  outf_e << "  Int_t tmpEtaBin = FindElectronEfficiencyBin( Eta , etaBins, " << NEtaBins << ");" << endl;
  outf_e << "  Int_t tmpPhiBin = FindElectronEfficiencyBin( Phi , phiBins, " << NPhiBins << ");" << endl;
  outf_e << "  return Efficiency[tmpPtBin][tmpEtaBin][tmpPhiBin];" << endl;
  outf_e << "}" << endl;


  outf_e << endl;
  outf_e << endl;
  outf_e << endl;
  outf_e << endl;


  outf_e << "Double_t GetElectronEfficiencyPtEta(Double_t Pt, Double_t Eta) {" << endl;

  outf_e << endl;
  outf_e << "  Double_t ptBins[" << NPtBins+1 << "] = {";
  for (uint i=0; i < NPtBins+1; ++i) {
    outf_e << ptBins[i];
    if (i < NPtBins) {
      outf_e << ",";
    }
  }
  outf_e << "};\n";

  outf_e << "  Double_t etaBins[" << NEtaBins+1 << "] = {";
  for (uint i=0; i < NEtaBins+1; ++i) {
    outf_e << etaBins[i];
    if (i < NEtaBins) {
      outf_e << ",";
    }
  }
  outf_e << "};\n";


  outf_e << endl;
  outf_e << endl;

  outf_e << "  Double_t Efficiency[" << NPtBins+2 << "][" << NEtaBins+2 << "] = {";
  outf_e << endl;

  for (uint i=0; i < NPtBins+2; ++i) {
    outf_e << "    {";
    for (uint j=0; j < NEtaBins+2; ++j) {
      outf_e << Efficiency_Electrons_PtEta[i][j];
      if (j< NEtaBins+1) {
        outf_e << ",";
      }
    }
    if (i< NPtBins+1) {
      outf_e << "    },";
    } else {
      outf_e << "}";
    }
    outf_e << endl;
  }
  
  outf_e << "  };" << endl;

  outf_e << endl;
  outf_e << endl;

  outf_e << "  Int_t tmpPtBin = FindElectronEfficiencyBin( Pt , ptBins, " << NPtBins << ");" << endl;
  outf_e << "  Int_t tmpEtaBin = FindElectronEfficiencyBin( Eta , etaBins, " << NEtaBins << ");" << endl;
  outf_e << "  return Efficiency[tmpPtBin][tmpEtaBin];" << endl;
  outf_e << "}" << endl;


  outf_e.close();




  ofstream outf_m("MuonEfficiencyMap.h");

  outf_m << "UInt_t FindMuonEfficiencyBin( double value, double bins[], UInt_t nbins) {" << endl;
  outf_m << "  UInt_t nbinboundaries = nbins+1;" << endl;
  outf_m << "  UInt_t bin = 0;" << endl;
  outf_m << "  for (uint i=0; i < nbinboundaries; ++i) {" << endl;
  outf_m << "    if (i < nbinboundaries-1) {" << endl;
  outf_m << "      if (value >= bins[i] && value < bins[i+1]) {" << endl;
  outf_m << "        bin = i+1;" << endl;
  outf_m << "        break;" << endl;
  outf_m << "      }" << endl;
  outf_m << "    } else if (i == nbinboundaries-1) {" << endl;
  outf_m << "      if (value >= bins[i]) {" << endl;
  outf_m << "        bin = nbinboundaries;" << endl;
  outf_m << "        break;" << endl;
  outf_m << "      }" << endl;
  outf_m << "    }    " << endl;
  outf_m << "  }" << endl;
  outf_m << "  return bin;" << endl;
  outf_m << "}" << endl;

  outf_m << endl;
  outf_m << endl;

  outf_m << "Double_t GetMuonEfficiencyPtEtaPhi(Double_t Pt, Double_t Eta, Double_t Phi) {" << endl;

  outf_m << endl;
  outf_m << "  Double_t ptBins[" << NPtBins+1 << "] = {";
  for (uint i=0; i < NPtBins+1; ++i) {
    outf_m << ptBins[i];
    if (i < NPtBins) {
      outf_m << ",";
    }
  }
  outf_m << "};\n";

  outf_m << "  Double_t etaBins[" << NEtaBins+1 << "] = {";
  for (uint i=0; i < NEtaBins+1; ++i) {
    outf_m << etaBins[i];
    if (i < NEtaBins) {
      outf_m << ",";
    }
  }
  outf_m << "};\n";

  outf_m << "  Double_t phiBins[" << NPhiBins+1 << "] = {";
  for (uint i=0; i < NPhiBins+1; ++i) {
    outf_m << phiBins[i];
    if (i < NPhiBins) {
      outf_m << ",";
    }
  }
  outf_m << "};\n";

  outf_m << endl;
  outf_m << endl;

  outf_m << "  Double_t Efficiency[" << NPtBins+2 << "][" << NEtaBins+2 << "][" << NPhiBins+2 << "]  = {";
  outf_m << endl;

  for (uint i=0; i < NPtBins+2; ++i) {
    outf_m << "    {" << endl;
    for (uint j=0; j < NEtaBins+2; ++j) {
      outf_m << "      {";
      for (uint k=0; k < NPhiBins+2; ++k) {
        outf_m << Efficiency_Muons_PtEtaPhi[i][j][k];
        if (k< NPhiBins+1) {
          outf_m << ",";
        }        
      }
      if (j< NEtaBins+1) {
        outf_m << "},";
      } else {
        outf_m << "}";
      }
      outf_m << endl;
    }
    if (i< NPtBins+1) {
      outf_m << "    },";
    } else {
      outf_m << "}";
    }
    outf_m << endl;
  }

  outf_m << "};" << endl;

  outf_m << endl;
  outf_m << endl;

  outf_m << "  Int_t tmpPtBin = FindMuonEfficiencyBin( Pt , ptBins, " << NPtBins << ");" << endl;
  outf_m << "  Int_t tmpEtaBin = FindMuonEfficiencyBin( Eta , etaBins, " << NEtaBins << ");" << endl;
  outf_m << "  Int_t tmpPhiBin = FindMuonEfficiencyBin( Phi , phiBins, " << NPhiBins << ");" << endl;
  outf_m << "  return Efficiency[tmpPtBin][tmpEtaBin][tmpPhiBin];" << endl;
  outf_m << "}" << endl;


  outf_m << endl;
  outf_m << endl;
  outf_m << endl;
  outf_m << endl;


  outf_m << "Double_t GetMuonEfficiencyPtEta(Double_t Pt, Double_t Eta) {" << endl;

  outf_m << endl;
  outf_m << "  Double_t ptBins[" << NPtBins+1 << "] = {";
  for (uint i=0; i < NPtBins+1; ++i) {
    outf_m << ptBins[i];
    if (i < NPtBins) {
      outf_m << ",";
    }
  }
  outf_m << "};\n";

  outf_m << "  Double_t etaBins[" << NEtaBins+1 << "] = {";
  for (uint i=0; i < NEtaBins+1; ++i) {
    outf_m << etaBins[i];
    if (i < NEtaBins) {
      outf_m << ",";
    }
  }
  outf_m << "};\n";


  outf_m << endl;
  outf_m << endl;

  outf_m << "  Double_t Efficiency[" << NPtBins+2 << "][" << NEtaBins+2 << "] = {";
  outf_m << endl;

  for (uint i=0; i < NPtBins+2; ++i) {
    outf_m << "    {";
    for (uint j=0; j < NEtaBins+2; ++j) {
      outf_m << Efficiency_Muons_PtEta[i][j];
      if (j< NEtaBins+1) {
        outf_m << ",";
      }
    }
    if (i< NPtBins+1) {
      outf_m << "    },";
    } else {
      outf_m << "}";
    }
    outf_m << endl;
  }
  
  outf_m << "  };" << endl;

  outf_m << endl;
  outf_m << endl;

  outf_m << "  Int_t tmpPtBin = FindMuonEfficiencyBin( Pt , ptBins, " << NPtBins << ");" << endl;
  outf_m << "  Int_t tmpEtaBin = FindMuonEfficiencyBin( Eta , etaBins, " << NEtaBins << ");" << endl;
  outf_m << "  return Efficiency[tmpPtBin][tmpEtaBin];" << endl;
  outf_m << "}" << endl;


  outf_m.close();




  gBenchmark->Show("WWTemplate");       
} 


