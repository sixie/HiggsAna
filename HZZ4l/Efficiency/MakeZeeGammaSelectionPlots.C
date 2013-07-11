//root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt",0)'
//================================================================================================
//
// HWW selection macro
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
#include <TH2F.h>                   // 2D histograms
#include <TGraphAsymmErrors.h>     
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <MyStyle.h>
#include "TLegend.h"
#include "TEfficiency.h"

#include "CITCommon/CommonData/interface/EffDataZGamma.hh" 
#include "HiggsAna/Utils/EfficiencyUtils.hh"


#endif


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
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("eleIDdir");
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

//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}

//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void MakeZeeGammaSelectionPlots(string inputfile = "/afs/cern.ch/work/s/sixie/public/HZZ4l/Efficiency/Electrons/Summer12_EleHZZICHEP2012WPWithZeeGammaStudy/probes.root")
{  


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *Mass_NoSelection_All = new TH1F("Mass_NoSelection_All", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  75, 50 , 200);
  TH1F *Mass_NoSelection_Bkg = new TH1F("Mass_NoSelection_Bkg", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  75, 50 , 200);
  TH1F *Mass_NoSelection_Sig = new TH1F("Mass_NoSelection_Sig", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  75, 50 , 200);

  TH1F *Mass_AfterPixelVeto_All = new TH1F("Mass_AfterPixelVeto_All", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  75, 50 , 200);
  TH1F *Mass_AfterPixelVeto_Bkg = new TH1F("Mass_AfterPixelVeto_Bkg", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  75, 50 , 200);
  TH1F *Mass_AfterPixelVeto_Sig = new TH1F("Mass_AfterPixelVeto_Sig", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  75, 50 , 200);

  TH1F *MassPlusDileptonmass_AfterPixelVeto_All = new TH1F("MassPlusDileptonmass_AfterPixelVeto_All", "; m_{ee #gamma} + m_{ee} [GeV/c^{2}] ; Number of Events ",  75, 50 , 300);
  TH1F *MassPlusDileptonmass_AfterPixelVeto_Bkg = new TH1F("MassPlusDileptonmass_AfterPixelVeto_Bkg", "; m_{ee #gamma} + m_{ee} [GeV/c^{2}] ; Number of Events ",  75, 50 , 300);
  TH1F *MassPlusDileptonmass_AfterPixelVeto_Sig = new TH1F("MassPlusDileptonmass_AfterPixelVeto_Sig", "; m_{ee #gamma} + m_{ee} [GeV/c^{2}] ; Number of Events ",  125, 50 , 300);


  TH1F *Mass_AfterISRRemoval_All = new TH1F("Mass_AfterISRRemoval_All", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  75, 50 , 200);
  TH1F *Mass_AfterISRRemoval_Bkg = new TH1F("Mass_AfterISRRemoval_Bkg", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  75, 50 , 200);
  TH1F *Mass_AfterISRRemoval_Sig = new TH1F("Mass_AfterISRRemoval_Sig", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  75, 50 , 200);

  TH1F *Dileptonmass_AfterISRRemoval_All = new TH1F("Dileptonmass_AfterISRRemoval_All", "; m_{ee} [GeV/c^{2}] ; Number of Events ",  100, 0 , 200);
  TH1F *Dileptonmass_AfterISRRemoval_Bkg = new TH1F("Dileptonmass_AfterISRRemoval_Bkg", "; m_{ee} [GeV/c^{2}] ; Number of Events ",  100, 0 , 200);
  TH1F *Dileptonmass_AfterISRRemoval_Sig = new TH1F("Dileptonmass_AfterISRRemoval_Sig", "; m_{ee} [GeV/c^{2}] ; Number of Events ",  100, 0 , 200);

  TH1F *Mass_AfterLowPtSelection_All = new TH1F("Mass_AfterLowPtSelection_All", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  75, 50 , 200);
  TH1F *Mass_AfterLowPtSelection_Bkg = new TH1F("Mass_AfterLowPtSelection_Bkg", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  75, 50 , 200);
  TH1F *Mass_AfterLowPtSelection_Sig = new TH1F("Mass_AfterLowPtSelection_Sig", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  75, 50 , 200);

  TH1F *R9_AfterLowPtSelection_All = new TH1F("R9_AfterLowPtSelection_All", "; Photon R9 ; Number of Events ",  50, 0 , 1);
  TH1F *R9_AfterLowPtSelection_Bkg = new TH1F("R9_AfterLowPtSelection_Bkg", "; Photon R9 ; Number of Events ",  50, 0 , 1);
  TH1F *R9_AfterLowPtSelection_Sig = new TH1F("R9_AfterLowPtSelection_Sig", "; Photon R9 ; Number of Events ",  50, 0 , 1);

  TH1F *SigmaIEtaIEtaBarrel_AfterLowPtSelection_All = new TH1F("SigmaIEtaIEtaBarrel_AfterLowPtSelection_All", "; Photon #sigma_{i#eta i#eta} ; Number of Events ",  50, 0 , 0.015);
  TH1F *SigmaIEtaIEtaBarrel_AfterLowPtSelection_Bkg = new TH1F("SigmaIEtaIEtaBarrel_AfterLowPtSelection_Bkg", "; Photon #sigma_{i#eta i#eta} ; Number of Events ",  50, 0 , 0.015);
  TH1F *SigmaIEtaIEtaBarrel_AfterLowPtSelection_Sig = new TH1F("SigmaIEtaIEtaBarrel_AfterLowPtSelection_Sig", "; Photon #sigma_{i#eta i#eta} ; Number of Events ",  50, 0 , 0.015);
  TH1F *SigmaIEtaIEtaEndcap_AfterLowPtSelection_All = new TH1F("SigmaIEtaIEtaEndcap_AfterLowPtSelection_All", "; Photon #sigma_{i#eta i#eta} ; Number of Events ",  50, 0.01 , 0.03);
  TH1F *SigmaIEtaIEtaEndcap_AfterLowPtSelection_Bkg = new TH1F("SigmaIEtaIEtaEndcap_AfterLowPtSelection_Bkg", "; Photon #sigma_{i#eta i#eta} ; Number of Events ",  50, 0.01 , 0.03);
  TH1F *SigmaIEtaIEtaEndcap_AfterLowPtSelection_Sig = new TH1F("SigmaIEtaIEtaEndcap_AfterLowPtSelection_Sig", "; Photon #sigma_{i#eta i#eta} ; Number of Events ",  50, 0.01 , 0.03);


  TH1F *DeltaRProbePhoton_AfterPhotonSelection_All = new TH1F("DeltaRProbePhoton_AfterPhotonSelection_All", "; #Delta R(Probe,#gamma) ; Number of Events ",  50, 0 , 5);
  TH1F *DeltaRProbePhoton_AfterPhotonSelection_Bkg = new TH1F("DeltaRProbePhoton_AfterPhotonSelection_Bkg", "; #Delta R(Probe,#gamma) ; Number of Events ",  50, 0 , 5);
  TH1F *DeltaRProbePhoton_AfterPhotonSelection_Sig = new TH1F("DeltaRProbePhoton_AfterPhotonSelection_Sig", "; #Delta R(Probe,#gamma) ; Number of Events ",  50, 0 , 5);

  TH1F *Mass_AfterAllSelection_All = new TH1F("Mass_AfterAllSelection_All", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  50, 50 , 200);
  TH1F *Mass_AfterAllSelection_Bkg = new TH1F("Mass_AfterAllSelection_Bkg", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  50, 50 , 200);
  TH1F *Mass_AfterAllSelection_Sig = new TH1F("Mass_AfterAllSelection_Sig", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  50, 50 , 200);

  TH1F *Mass_HighDPhiRegion_Bkg = new TH1F("Mass_HighDPhiRegion_Bkg", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  50, 50 , 200);

  TH1F *Mass_AfterAllSelectionProbeFail_All = new TH1F("Mass_AfterAllSelectionProbeFail_All", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  50, 50 , 200);
  TH1F *Mass_AfterAllSelectionProbeFail_Bkg = new TH1F("Mass_AfterAllSelectionProbeFail_Bkg", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  50, 50 , 200);
  TH1F *Mass_AfterAllSelectionProbeFail_Sig = new TH1F("Mass_AfterAllSelectionProbeFail_Sig", "; m_{ee #gamma} [GeV/c^{2}] ; Number of Events ",  50, 50 , 200);

  TH1F *MassTagPhoton_AfterAllSelection_All = new TH1F("MassTagPhoton_AfterAllSelection_All", "; m_{tag #gamma} [GeV/c^{2}] ; Number of Events ",  50, 0 , 200);
  TH1F *MassTagPhoton_AfterAllSelection_Bkg = new TH1F("MassTagPhoton_AfterAllSelection_Bkg", "; m_{tag #gamma} [GeV/c^{2}] ; Number of Events ",  50, 0 , 200);
  TH1F *MassTagPhoton_AfterAllSelection_Sig = new TH1F("MassTagPhoton_AfterAllSelection_Sig", "; m_{tag #gamma} [GeV/c^{2}] ; Number of Events ",  50, 0 , 200);


  TH1F *DeltaRProbePhoton_AfterPhotonSelection_Sig_Numerator = new TH1F("DeltaRProbePhoton_AfterPhotonSelection_Sig_Numerator", "; #Delta R(Probe,#gamma) ; Number of Events ",  10, 0 , 2.0);
  TH1F *DeltaRProbePhoton_AfterPhotonSelection_Sig_Denominator = new TH1F("DeltaRProbePhoton_AfterPhotonSelection_Sig_Denominator", "; #Delta R(Probe,#gamma) ; Number of Events ",  10, 0 , 2.0);




  //*****************************************************************************************
  //Load Branches
  //*****************************************************************************************
  TTree *tree = getTreeFromFile(inputfile.c_str(),"Events");

  EffDataZGamma data;
  tree->SetBranchAddress("Events",&data);

  for(UInt_t ientry=0; ientry < tree->GetEntries(); ientry++) {       	
    tree->GetEntry(ientry);
    
    Mass_NoSelection_All->Fill(data.mass);
    if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
      Mass_NoSelection_Bkg->Fill(data.mass);      
    } else {
      Mass_NoSelection_Sig->Fill(data.mass);      
    }

    //********************************************************************************
    //apply pixel veto
    //********************************************************************************
    if (!data.phopasspixelveto) continue;

    Mass_AfterPixelVeto_All->Fill(data.mass);
    if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
      Mass_AfterPixelVeto_Bkg->Fill(data.mass);      
    } else {
      Mass_AfterPixelVeto_Sig->Fill(data.mass);      
    }

    MassPlusDileptonmass_AfterPixelVeto_All->Fill(data.mass + data.massll);
    if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
      MassPlusDileptonmass_AfterPixelVeto_Bkg->Fill(data.mass + data.massll);      
    } else {
      MassPlusDileptonmass_AfterPixelVeto_Sig->Fill(data.mass + data.massll);      
    }

    //********************************************************************************
    //apply mass+massll cut
    //********************************************************************************
    if (!(data.mass + data.massll < 180)) continue;

    Mass_AfterISRRemoval_All->Fill(data.mass);
    if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
      Mass_AfterISRRemoval_Bkg->Fill(data.mass);      
    } else {
      Mass_AfterISRRemoval_Sig->Fill(data.mass);      
    }

    Dileptonmass_AfterISRRemoval_All->Fill(data.massll);
    if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
      Dileptonmass_AfterISRRemoval_Bkg->Fill(data.massll);      
    } else {
      Dileptonmass_AfterISRRemoval_Sig->Fill(data.massll);      
    }

    //********************************************************************************
    //apply massll cut
    //********************************************************************************
    if (!(data.massll > 40)) continue;

    //********************************************************************************
    //low pT leptons
    //********************************************************************************
    if (!(data.pt > 7 && data.pt < 10)) continue;

    Mass_AfterLowPtSelection_All->Fill(data.mass);
    if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
      Mass_AfterLowPtSelection_Bkg->Fill(data.mass);      
    } else {
      Mass_AfterLowPtSelection_Sig->Fill(data.mass);      
    }

    R9_AfterLowPtSelection_All->Fill(data.phoR9);
    if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
      R9_AfterLowPtSelection_Bkg->Fill(data.phoR9);      
    } else {
      R9_AfterLowPtSelection_Sig->Fill(data.phoR9);      
    }

    //********************************************************************************
    //R9 selection
    //********************************************************************************
    if (!(data.phoR9 > 0.9)) continue;

    if (fabs(data.phoeta) < 1.479) {
      SigmaIEtaIEtaBarrel_AfterLowPtSelection_All->Fill(data.phosigiEtaiEta);
      if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
        SigmaIEtaIEtaBarrel_AfterLowPtSelection_Bkg->Fill(data.phosigiEtaiEta);      
      } else {
        SigmaIEtaIEtaBarrel_AfterLowPtSelection_Sig->Fill(data.phosigiEtaiEta);      
      }
    } else {
      SigmaIEtaIEtaEndcap_AfterLowPtSelection_All->Fill(data.phosigiEtaiEta);
      if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
        SigmaIEtaIEtaEndcap_AfterLowPtSelection_Bkg->Fill(data.phosigiEtaiEta);      
      } else {
        SigmaIEtaIEtaEndcap_AfterLowPtSelection_Sig->Fill(data.phosigiEtaiEta);      
      }
    }

    //********************************************************************************
    //sigmaIEtaIEta cut
    //********************************************************************************
    if (fabs(data.phoeta) < 1.479) {
      if (!(data.phosigiEtaiEta < 0.0105)) continue;
    } else {
      if (!(data.phosigiEtaiEta < 0.028)) continue;
    }

    DeltaRProbePhoton_AfterPhotonSelection_All->Fill(data.drprobepho);
    if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
      DeltaRProbePhoton_AfterPhotonSelection_Bkg->Fill(data.drprobepho);      
    } else {
      DeltaRProbePhoton_AfterPhotonSelection_Sig->Fill(data.drprobepho);      
    }


    //********************************************************************************
    //Efficiency Vs DeltaR
    //********************************************************************************
    if ((data.phoisreal && data.probeisreal && data.tagisreal)) {
      DeltaRProbePhoton_AfterPhotonSelection_Sig_Denominator->Fill(data.drprobepho);
      if (data.pass) DeltaRProbePhoton_AfterPhotonSelection_Sig_Numerator->Fill(data.drprobepho);
    }


    //********************************************************************************
    //massll cut
    //********************************************************************************
    if (!(data.drprobepho < 1.5)) continue;


    //Modelling background without dphi cut
    if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
      Mass_HighDPhiRegion_Bkg->Fill(data.mass);  
    }



    //********************************************************************************
    //delta phi cut
    //********************************************************************************
    if (!(data.drprobepho < 0.7)) continue;

    Mass_AfterAllSelection_All->Fill(data.mass);
    if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
      Mass_AfterAllSelection_Bkg->Fill(data.mass);      
    } else {
      Mass_AfterAllSelection_Sig->Fill(data.mass);      
    }

    //Mass TagPhoton
    MassTagPhoton_AfterAllSelection_All->Fill(data.masstagpho);
    if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
      MassTagPhoton_AfterAllSelection_Bkg->Fill(data.masstagpho);      
    } else {
      MassTagPhoton_AfterAllSelection_Sig->Fill(data.masstagpho);      
    }


    if (data.pass) continue;

    Mass_AfterAllSelectionProbeFail_All->Fill(data.mass);
    if (!(data.phoisreal && data.probeisreal && data.tagisreal)) {
      Mass_AfterAllSelectionProbeFail_Bkg->Fill(data.mass);      
    } else {
      Mass_AfterAllSelectionProbeFail_Sig->Fill(data.mass);      
    }

   


  } 
  



  TCanvas *cv = 0;
  TLegend *legend = 0;

  //**********************************************************************************
  //No Selection
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.6,0.7,0.94,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(Mass_NoSelection_All,"Total", "L");
  legend->AddEntry(Mass_NoSelection_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(Mass_NoSelection_Sig,"Z#rightarrow ee#gamma", "L");

  Mass_NoSelection_All->SetLineColor(kBlue);
  Mass_NoSelection_Bkg->SetLineColor(kRed);
  Mass_NoSelection_Sig->SetLineColor(kGreen+2);

  Mass_NoSelection_All->GetXaxis()->SetTitleOffset(1.05);
  Mass_NoSelection_All->Draw("hist");
  Mass_NoSelection_Bkg->Draw("same,hist");
  Mass_NoSelection_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("MassEEGamma_NoSelection.gif");



  //**********************************************************************************
  //After PixelVeto Selection
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.6,0.7,0.94,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(Mass_AfterPixelVeto_All,"Total", "L");
  legend->AddEntry(Mass_AfterPixelVeto_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(Mass_AfterPixelVeto_Sig,"Z#rightarrow ee#gamma", "L");

  Mass_AfterPixelVeto_All->SetLineColor(kBlue);
  Mass_AfterPixelVeto_Bkg->SetLineColor(kRed);
  Mass_AfterPixelVeto_Sig->SetLineColor(kGreen+2);

  Mass_AfterPixelVeto_All->GetXaxis()->SetTitleOffset(1.05);
  Mass_AfterPixelVeto_All->Draw("hist");
  Mass_AfterPixelVeto_Bkg->Draw("same,hist");
  Mass_AfterPixelVeto_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("MassEEGamma_AfterPixelVeto.gif");


  //**********************************************************************************
  //mass + massll
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.2,0.7,0.6,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(MassPlusDileptonmass_AfterPixelVeto_All,"Total", "L");
  legend->AddEntry(MassPlusDileptonmass_AfterPixelVeto_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(MassPlusDileptonmass_AfterPixelVeto_Sig,"Z#rightarrow ee#gamma", "L");

  MassPlusDileptonmass_AfterPixelVeto_All->SetLineColor(kBlue);
  MassPlusDileptonmass_AfterPixelVeto_Bkg->SetLineColor(kRed);
  MassPlusDileptonmass_AfterPixelVeto_Sig->SetLineColor(kGreen+2);

  MassPlusDileptonmass_AfterPixelVeto_All->GetXaxis()->SetTitleOffset(1.05);
  MassPlusDileptonmass_AfterPixelVeto_All->Draw("hist");
  MassPlusDileptonmass_AfterPixelVeto_Bkg->Draw("same,hist");
  MassPlusDileptonmass_AfterPixelVeto_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("MassPlusDileptonmass_AfterPixelVeto.gif");


  //**********************************************************************************
  //After ISRRemoval
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.6,0.7,0.94,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(Mass_AfterISRRemoval_All,"Total", "L");
  legend->AddEntry(Mass_AfterISRRemoval_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(Mass_AfterISRRemoval_Sig,"Z#rightarrow ee#gamma", "L");

  Mass_AfterISRRemoval_All->SetLineColor(kBlue);
  Mass_AfterISRRemoval_Bkg->SetLineColor(kRed);
  Mass_AfterISRRemoval_Sig->SetLineColor(kGreen+2);

  Mass_AfterISRRemoval_All->GetXaxis()->SetTitleOffset(1.05);
  Mass_AfterISRRemoval_All->Draw("hist");
  Mass_AfterISRRemoval_Bkg->Draw("same,hist");
  Mass_AfterISRRemoval_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("MassEEGamma_AfterISRRemoval.gif");

  //**********************************************************************************
  //Dileptonmass After ISR Removal
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.6,0.7,0.94,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(Dileptonmass_AfterISRRemoval_All,"Total", "L");
  legend->AddEntry(Dileptonmass_AfterISRRemoval_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(Dileptonmass_AfterISRRemoval_Sig,"Z#rightarrow ee#gamma", "L");

  Dileptonmass_AfterISRRemoval_All->SetLineColor(kBlue);
  Dileptonmass_AfterISRRemoval_Bkg->SetLineColor(kRed);
  Dileptonmass_AfterISRRemoval_Sig->SetLineColor(kGreen+2);

  Dileptonmass_AfterISRRemoval_All->GetXaxis()->SetTitleOffset(1.05);
  Dileptonmass_AfterISRRemoval_All->Draw("hist");
  Dileptonmass_AfterISRRemoval_Bkg->Draw("same,hist");
  Dileptonmass_AfterISRRemoval_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("Dileptonmass_AfterISRRemoval.gif");

  //**********************************************************************************
  //After Low Pt Selection
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.6,0.7,0.94,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(Mass_AfterLowPtSelection_All,"Total", "L");
  legend->AddEntry(Mass_AfterLowPtSelection_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(Mass_AfterLowPtSelection_Sig,"Z#rightarrow ee#gamma", "L");

  Mass_AfterLowPtSelection_All->SetLineColor(kBlue);
  Mass_AfterLowPtSelection_Bkg->SetLineColor(kRed);
  Mass_AfterLowPtSelection_Sig->SetLineColor(kGreen+2);

  Mass_AfterLowPtSelection_All->GetXaxis()->SetTitleOffset(1.05);
  Mass_AfterLowPtSelection_All->Draw("hist");
  Mass_AfterLowPtSelection_Bkg->Draw("same,hist");
  Mass_AfterLowPtSelection_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("MassEEGamma_AfterLowPtSelection.gif");


  //**********************************************************************************
  //Photon R9
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.6,0.7,0.94,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(R9_AfterLowPtSelection_All,"Total", "L");
  legend->AddEntry(R9_AfterLowPtSelection_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(R9_AfterLowPtSelection_Sig,"Z#rightarrow ee#gamma", "L");

  R9_AfterLowPtSelection_All->SetLineColor(kBlue);
  R9_AfterLowPtSelection_Bkg->SetLineColor(kRed);
  R9_AfterLowPtSelection_Sig->SetLineColor(kGreen+2);

  R9_AfterLowPtSelection_All->GetXaxis()->SetTitleOffset(1.05);
  R9_AfterLowPtSelection_All->Draw("hist");
  R9_AfterLowPtSelection_Bkg->Draw("same,hist");
  R9_AfterLowPtSelection_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("R9_AfterLowPtSelection.gif");

  //**********************************************************************************
  //Photon SigmaIEtaIEta
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.2,0.7,0.6,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(SigmaIEtaIEtaBarrel_AfterLowPtSelection_All,"Total", "L");
  legend->AddEntry(SigmaIEtaIEtaBarrel_AfterLowPtSelection_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(SigmaIEtaIEtaBarrel_AfterLowPtSelection_Sig,"Z#rightarrow ee#gamma", "L");

  SigmaIEtaIEtaBarrel_AfterLowPtSelection_All->SetLineColor(kBlue);
  SigmaIEtaIEtaBarrel_AfterLowPtSelection_Bkg->SetLineColor(kRed);
  SigmaIEtaIEtaBarrel_AfterLowPtSelection_Sig->SetLineColor(kGreen+2);

  SigmaIEtaIEtaBarrel_AfterLowPtSelection_All->GetXaxis()->SetTitleOffset(1.05);
  SigmaIEtaIEtaBarrel_AfterLowPtSelection_All->Draw("hist");
  SigmaIEtaIEtaBarrel_AfterLowPtSelection_Bkg->Draw("same,hist");
  SigmaIEtaIEtaBarrel_AfterLowPtSelection_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("SigmaIEtaIEtaBarrel_AfterLowPtSelection.gif");

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.2,0.7,0.6,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(SigmaIEtaIEtaEndcap_AfterLowPtSelection_All,"Total", "L");
  legend->AddEntry(SigmaIEtaIEtaEndcap_AfterLowPtSelection_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(SigmaIEtaIEtaEndcap_AfterLowPtSelection_Sig,"Z#rightarrow ee#gamma", "L");

  SigmaIEtaIEtaEndcap_AfterLowPtSelection_All->SetLineColor(kBlue);
  SigmaIEtaIEtaEndcap_AfterLowPtSelection_Bkg->SetLineColor(kRed);
  SigmaIEtaIEtaEndcap_AfterLowPtSelection_Sig->SetLineColor(kGreen+2);

  SigmaIEtaIEtaEndcap_AfterLowPtSelection_All->GetXaxis()->SetTitleOffset(1.05);
  SigmaIEtaIEtaEndcap_AfterLowPtSelection_All->Draw("hist");
  SigmaIEtaIEtaEndcap_AfterLowPtSelection_Bkg->Draw("same,hist");
  SigmaIEtaIEtaEndcap_AfterLowPtSelection_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("SigmaIEtaIEtaEndcap_AfterLowPtSelection.gif");


  //**********************************************************************************
  //DeltaRProbePhoton After Photon Selection
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.6,0.7,0.94,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(DeltaRProbePhoton_AfterPhotonSelection_All,"Total", "L");
  legend->AddEntry(DeltaRProbePhoton_AfterPhotonSelection_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(DeltaRProbePhoton_AfterPhotonSelection_Sig,"Z#rightarrow ee#gamma", "L");

  DeltaRProbePhoton_AfterPhotonSelection_All->SetLineColor(kBlue);
  DeltaRProbePhoton_AfterPhotonSelection_Bkg->SetLineColor(kRed);
  DeltaRProbePhoton_AfterPhotonSelection_Sig->SetLineColor(kGreen+2);

  DeltaRProbePhoton_AfterPhotonSelection_All->GetXaxis()->SetTitleOffset(1.05);
  DeltaRProbePhoton_AfterPhotonSelection_All->Draw("hist");
  DeltaRProbePhoton_AfterPhotonSelection_Bkg->Draw("same,hist");
  DeltaRProbePhoton_AfterPhotonSelection_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("DeltaRProbePhoton_AfterPhotonSelection.gif");


  //**********************************************************************************
  //After All Selection
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.6,0.7,0.94,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(Mass_AfterAllSelection_All,"Total", "L");
  legend->AddEntry(Mass_AfterAllSelection_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(Mass_AfterAllSelection_Sig,"Z#rightarrow ee#gamma", "L");

  Mass_AfterAllSelection_All->SetLineColor(kBlue);
  Mass_AfterAllSelection_Bkg->SetLineColor(kRed);
  Mass_AfterAllSelection_Sig->SetLineColor(kGreen+2);

  Mass_AfterAllSelection_All->GetXaxis()->SetTitleOffset(1.05);
  Mass_AfterAllSelection_All->Draw("hist");
  Mass_AfterAllSelection_Bkg->Draw("same,hist");
  Mass_AfterAllSelection_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("MassEEGamma_AfterAllSelection.gif");


  //**********************************************************************************
  //After All Selection
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.6,0.7,0.94,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(MassTagPhoton_AfterAllSelection_All,"Total", "L");
  legend->AddEntry(MassTagPhoton_AfterAllSelection_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(MassTagPhoton_AfterAllSelection_Sig,"Z#rightarrow ee#gamma", "L");

  MassTagPhoton_AfterAllSelection_All->SetLineColor(kBlue);
  MassTagPhoton_AfterAllSelection_Bkg->SetLineColor(kRed);
  MassTagPhoton_AfterAllSelection_Sig->SetLineColor(kGreen+2);

  MassTagPhoton_AfterAllSelection_All->GetXaxis()->SetTitleOffset(1.05);
  MassTagPhoton_AfterAllSelection_All->Draw("hist");
  MassTagPhoton_AfterAllSelection_Bkg->Draw("same,hist");
  MassTagPhoton_AfterAllSelection_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("MassTagPhoton_AfterAllSelection.gif");



  //**********************************************************************************
  //After All Selection - Fail Sample
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.6,0.7,0.94,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(Mass_AfterAllSelectionProbeFail_All,"Total", "L");
  legend->AddEntry(Mass_AfterAllSelectionProbeFail_Bkg,"Bkg (ele fakes photon)", "L");
  legend->AddEntry(Mass_AfterAllSelectionProbeFail_Sig,"Z#rightarrow ee#gamma", "L");

  Mass_AfterAllSelectionProbeFail_All->SetLineColor(kBlue);
  Mass_AfterAllSelectionProbeFail_Bkg->SetLineColor(kRed);
  Mass_AfterAllSelectionProbeFail_Sig->SetLineColor(kGreen+2);

  Mass_AfterAllSelectionProbeFail_All->GetXaxis()->SetTitleOffset(1.05);
  Mass_AfterAllSelectionProbeFail_All->Draw("hist");
  Mass_AfterAllSelectionProbeFail_Bkg->Draw("same,hist");
  Mass_AfterAllSelectionProbeFail_Sig->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("MassEEGamma_AfterAllSelectionProbeFail.gif");



  //**********************************************************************************
  //Modeling of the background
  //**********************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.6,0.7,0.94,0.9);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  legend->AddEntry(Mass_AfterAllSelection_Bkg,"Nominal Selection", "L");
  legend->AddEntry(Mass_HighDPhiRegion_Bkg,"0.7 < #Delta R(probe,#gamma) < 1.5", "L");

//   Mass_AfterAllSelection_Bkg->Rebin(2);
//   Mass_HighDPhiRegion_Bkg->Rebin(2);

  NormalizeHist(Mass_AfterAllSelection_Bkg);
  NormalizeHist(Mass_HighDPhiRegion_Bkg);

  Mass_AfterAllSelection_Bkg->SetLineColor(kBlack);
  Mass_HighDPhiRegion_Bkg->SetLineColor(kRed);

  Mass_AfterAllSelection_Bkg->GetXaxis()->SetTitleOffset(1.05);
  Mass_AfterAllSelection_Bkg->Draw("e1");
  Mass_HighDPhiRegion_Bkg->Draw("same,hist");
  
  legend->Draw();
  cv->SaveAs("MassEEGamma_BackgroundShapeModeling.gif");


  //**********************************************************************************
  //Efficiency vs DeltaR
  //**********************************************************************************

  TGraphAsymmErrors *efficiency_deltaR = higgsana::createEfficiencyGraph(DeltaRProbePhoton_AfterPhotonSelection_Sig_Numerator, DeltaRProbePhoton_AfterPhotonSelection_Sig_Denominator, "", vector<double>(0),  -99, -99, 0, 1);
  cv = new TCanvas("cv", "cv", 800, 600);
  efficiency_deltaR->Draw("AP");
  efficiency_deltaR->GetXaxis()->SetTitleOffset(1.05);
  efficiency_deltaR->GetYaxis()->SetTitleOffset(1.2);
  cv->SaveAs("Efficiency_DeltaRProbePhoton.gif");



} 

