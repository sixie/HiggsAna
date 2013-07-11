
#include <string.h>
#include <iostream>
#include "TMath.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"




void MakeEffPullPlots(string filename, string Label = "") {

  string label = "";
  if (Label != "") label = "_" + Label;  

  gStyle->SetOptStat(0);

  //Histograms
  TH1F *histEff = new TH1F("histEff"," ;Eff; Number of Toys", 50, 0, 1);
  TH1F *histEffPull = new TH1F("histEffPull"," ;Eff / #Delta Eff; Number of Toys", 50, -5, 5);

  TFile *file = new TFile(filename.c_str(),"read");  
  
  Float_t                 varEff;
  Float_t                 varEffErrL;
  Float_t                 varEffErrH;

  TTree *tree = (TTree*)file->Get("eff");
  tree->SetBranchAddress( "eff", &varEff);
  tree->SetBranchAddress( "efferrl", &varEffErrL);
  tree->SetBranchAddress( "efferrh", &varEffErrH);

  double inputEff = 0.5502;

  for (UInt_t i=0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    //**********************************************************
    //fill fit result and pull plots
    //**********************************************************
    histEff->Fill(varEff);
    if (varEff > inputEff) {
      histEffPull->Fill((varEff - inputEff) / varEffErrL);
    } else {
      histEffPull->Fill((varEff - inputEff) / varEffErrH);
    }

  }
  
  TCanvas *cv = 0;
    TPaveText *pt;

  cv = new TCanvas ("cv", "cv", 800, 600);
  histEff->Draw("hist");
  histEff->GetYaxis()->SetTitleOffset(1.2);
  histEff->GetXaxis()->SetTitleOffset(1.05);
  pt = new TPaveText(0.1, histEff->GetMaximum()*0.4, 0.4, histEff->GetMaximum()*0.8);

  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->AddText(Form("Mean : %.3f #pm %.3f", histEff->GetMean(), histEff->GetMeanError()));
  pt->AddText(Form("RMS : %.3f #pm %.3f", histEff->GetRMS(), histEff->GetRMSError()));
  pt->SetTextSize(0.04);
  pt->Draw("same");
  cv->SaveAs(("Eff" + label + ".gif").c_str());


  cv = new TCanvas ("cv", "cv", 800, 600);
  histEffPull->Draw("hist");
  pt = new TPaveText(1.5, histEffPull->GetMaximum()*0.4, 2.5, histEffPull->GetMaximum()*0.8);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->AddText(Form("Mean : %.3f #pm %.3f", histEffPull->GetMean(), histEffPull->GetMeanError()));
  pt->AddText(Form("RMS : %.3f #pm %.3f", histEffPull->GetRMS(), histEffPull->GetRMSError()));
  pt->SetTextSize(0.04);
  pt->Draw("same");
  cv->SaveAs(("EffPull" + label + ".gif").c_str());
  


}

void ZeeGammaMassFitSystematicPlotPulls() {

  MakeEffPullPlots("/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/ZeeGammaMassFitSystematicStudyToys/EffToyResults_Option0.root", "Option0");
    MakeEffPullPlots("/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/ZeeGammaMassFitSystematicStudyToys/EffToyResults_Option1.root", "Option1");
     MakeEffPullPlots("/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/ZeeGammaMassFitSystematicStudyToys/EffToyResults_Option2.root", "Option2");
//   MakeEffPullPlots("/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/ZeeGammaMassFitSystematicStudyToys/EffToyResults_Option3.root", "Option3");
//   MakeEffPullPlots("/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/ZeeGammaMassFitSystematicStudyToys/EffToyResults_Option4.root", "Option4");
//    MakeEffPullPlots("/afs/cern.ch/work/s/sixie/public/releases/CMSSW_5_2_3_patch1/src/HiggsAna/HZZ4l/Efficiency/results/ZeeGammaMassFitSystematicStudyToys/EffToyResults_Option10.root", "Option10");

}
