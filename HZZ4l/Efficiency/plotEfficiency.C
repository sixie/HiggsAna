//root -l plotEfficiencyComparisons.C+'("results/Data2012_EleHZZICHEP2012WPWithTightTag_190456-196531/basic2_76_106_FineBins/eff.root","results/Data2012_EleHZZICHEP2012WPWithTightTag/basic2_76_106_FineBins/eff.root","Run2012A/B","Run2012C","ElectronHZZICHEP2012Selection")'

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TGraphAsymmErrors.h>      // graphs
#include <TH2F.h>                   // 2D histograms
#include <TLegend.h>                   // 2D histograms
#include <TCanvas.h>                   // 2D histograms
#include <TMath.h>                  // ROOT math library
#include <TLatex.h>                
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#endif


void plotEfficiencyComparisonsVsSeparateObservables(string filename1, 
                              string filename2, 
                              string Label1,
                              string Label2,                               
                              string Label)
{
  gBenchmark->Start("printMuonWPEff");
    
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TFile efffile1(filename1.c_str());
  TFile efffile2(filename2.c_str());
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   
  
  TGraphAsymmErrors* EffVsNVtx1   = ( TGraphAsymmErrors*)efffile1.Get("grEffNPV");
  TGraphAsymmErrors* EffVsPt1   = ( TGraphAsymmErrors*)efffile1.Get("grEffPt");
  TGraphAsymmErrors* EffVsEta1   = ( TGraphAsymmErrors*)efffile1.Get("grEffEta");
  TGraphAsymmErrors* EffVsNVtx2   = ( TGraphAsymmErrors*)efffile2.Get("grEffNPV");
  TGraphAsymmErrors* EffVsPt2   = ( TGraphAsymmErrors*)efffile2.Get("grEffPt");
  TGraphAsymmErrors* EffVsEta2   = ( TGraphAsymmErrors*)efffile2.Get("grEffEta");

//   assert(EffVsNVtx1);
//   assert(EffVsPt1);
//   assert(EffVsEta1);
//   assert(EffVsNVtx2);
//   assert(EffVsPt2);
//   assert(EffVsEta2);
 
  TCanvas *cv = 0;
  TLegend *legend = 0;

//   //***************************************************************
//   //Plot Eff Vs NVtx
//   //***************************************************************
//   cv = new TCanvas("cv", "Canvas", 800, 600);
//   legend = new TLegend( 0.6, 0.7, 0.9, 0.9);
//   legend->SetTextSize(0.03);
//   legend->SetBorderSize(0);
//   legend->SetFillStyle(0);
//   legend->AddEntry(EffVsNVtx1, Label1.c_str(), "LP");
//   legend->AddEntry(EffVsNVtx2, Label2.c_str(), "LP");
//   EffVsNVtx1->SetTitle("");
//   EffVsNVtx1->GetXaxis()->SetTitle("Number of Reconstructed Primary Vertices");
//   EffVsNVtx1->GetXaxis()->SetTitleOffset(1.05);
//   EffVsNVtx1->GetYaxis()->SetTitle("Efficiency");
//   EffVsNVtx1->GetYaxis()->SetTitleOffset(1.2);
//   EffVsNVtx1->GetYaxis()->SetRangeUser(0.75,1.0);
//   EffVsNVtx1->SetLineColor(kRed);
//   EffVsNVtx1->SetMarkerColor(kRed);
//   EffVsNVtx2->SetLineColor(kBlue);
//   EffVsNVtx2->SetMarkerColor(kBlue);
//   EffVsNVtx1->Draw("AP");
//   EffVsNVtx2->Draw("P");
//   legend->Draw();
//   cv->SaveAs( ("EfficiencyVsNVtx"+label+".gif").c_str());

//   //***************************************************************
//   //Plot Eff Vs Pt
//   //***************************************************************
//   cv = new TCanvas("cv", "Canvas", 800, 600);
//   legend = new TLegend( 0.6, 0.4, 0.9, 0.7);
//   legend->SetTextSize(0.03);
//   legend->SetBorderSize(0);
//   legend->SetFillStyle(0);
//   legend->AddEntry(EffVsPt1, Label1.c_str(), "LP");
//   legend->AddEntry(EffVsPt2, Label2.c_str(), "LP");
//   EffVsPt1->SetTitle("");
//   EffVsPt1->GetXaxis()->SetTitle("p_{T} [GeV/c]");
//   EffVsPt1->GetXaxis()->SetTitleOffset(1.05);
//   EffVsPt1->GetYaxis()->SetTitle("Efficiency");
//   EffVsPt1->GetYaxis()->SetTitleOffset(1.2);
//   EffVsPt1->GetYaxis()->SetRangeUser(0.0,1.0);
//   EffVsPt1->SetLineColor(kRed);
//   EffVsPt1->SetMarkerColor(kRed);
//   EffVsPt2->SetLineColor(kBlue);
//   EffVsPt2->SetMarkerColor(kBlue);
//   EffVsPt1->Draw("AP");
//   EffVsPt2->Draw("P");
//   legend->Draw();
//   cv->SaveAs( ("EfficiencyVsPt"+label+".gif").c_str());

  //***************************************************************
  //Plot Eff Vs Eta
  //***************************************************************
  cv = new TCanvas("cv", "Canvas", 800, 600);
  legend = new TLegend( 0.2, 0.2, 0.6, 0.5);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(EffVsEta1, Label1.c_str(), "LP");
  legend->AddEntry(EffVsEta2, Label2.c_str(), "LP");
  EffVsEta1->SetTitle("");
  EffVsEta1->GetXaxis()->SetTitle("#eta");
  EffVsEta1->GetXaxis()->SetTitleOffset(1.05);
  EffVsEta1->GetYaxis()->SetTitle("Efficiency");
  EffVsEta1->GetYaxis()->SetTitleOffset(1.2);
  EffVsEta1->GetYaxis()->SetRangeUser(0.0,1.0);
  EffVsEta1->SetLineColor(kRed);
  EffVsEta1->SetMarkerColor(kRed);
  EffVsEta2->SetLineColor(kBlue);
  EffVsEta2->SetMarkerColor(kBlue);
  EffVsEta1->Draw("AP");
  EffVsEta2->Draw("P");
  legend->Draw();
  cv->SaveAs( ("EfficiencyVsEta"+label+".gif").c_str());


}



void plotEfficiencyComparison(string filenameData, 
                              string filenameMC, 
                              string legendLabelData,
                              string legendLabelMC,
                              Int_t dataset,
                              string Label)
{
  gBenchmark->Start("printMuonWPEff");
    
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TFile efffileData(filenameData.c_str());
  TFile efffileMC(filenameMC.c_str());
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   
  
  TH2D* EffEtaPtData   = ( TH2D*)efffileData.Get("hEffEtaPt");
  TH2D* EffErrLowEtaPtData   = ( TH2D*)efffileData.Get("hErrlEtaPt");
  TH2D* EffErrHighEtaPtData   = ( TH2D*)efffileData.Get("hErrhEtaPt");
  TH2D* EffEtaPtMC   = ( TH2D*)efffileMC.Get("hEffEtaPt");
  TH2D* EffErrLowEtaPtMC   = ( TH2D*)efffileMC.Get("hErrlEtaPt");
  TH2D* EffErrHighEtaPtMC   = ( TH2D*)efffileMC.Get("hErrhEtaPt");

  assert(EffEtaPtData);
  assert(EffErrLowEtaPtData);
  assert(EffErrHighEtaPtData);
  assert(EffEtaPtMC);
  assert(EffErrLowEtaPtMC);
  assert(EffErrHighEtaPtMC);
 
  TCanvas *cv = 0;
  TLegend *legend = 0;

  const UInt_t n = 7;
  Double_t xval[n], xerr[n];
  xval[0] = 8.5; xerr[0] = 1.5;
  xval[1] = 12.5; xerr[1] = 2.5;
  xval[2] = 17.5; xerr[2] = 2.5;
  xval[3] = 25; xerr[3] = 5;
  xval[4] = 35; xerr[4] = 5;
  xval[5] = 45; xerr[5] = 5;
  xval[6] = 75; xerr[6] = 25;

  Double_t yvalData[n], yerrlData[n], yerrhData[n];
  Double_t yvalMC[n], yerrlMC[n], yerrhMC[n];

  for (UInt_t b=1; b < 6; ++b) {

    string binLabel;
    if (b==1) binLabel = "Eta0p0To0p8";
    if (b==2) binLabel = "Eta0p8To1p4442";
    if (b==3) binLabel = "Eta1p4442To1p566";
    if (b==4) binLabel = "Eta1p566To2p0";
    if (b==5) binLabel = "Eta2p0To2p5";

    string binTexLabel;
    if (b==1) binTexLabel = "Probe in Barrel: 0 < |#eta_{probe}| < 0.8";
    if (b==2) binTexLabel = "Probe in Barrel: 0.8 < |#eta_{probe}| < 1.4442";
    if (b==3) binTexLabel = "Probe in Crack: 1.4442 < |#eta_{probe}| < 1.566";
    if (b==4) binTexLabel = "Probe in Endcap: 1.566 < |#eta_{probe}| < 2.0";
    if (b==5) binTexLabel = "Probe in Endcap: 2.0 < |#eta_{probe}| < 2.5";

    yvalData[0] = EffEtaPtData->GetBinContent(1,b); 
    yvalData[1] = EffEtaPtData->GetBinContent(2,b);
    yvalData[2] = EffEtaPtData->GetBinContent(3,b);
    yvalData[3] = EffEtaPtData->GetBinContent(4,b);
    yvalData[4] = EffEtaPtData->GetBinContent(5,b);
    yvalData[5] = EffEtaPtData->GetBinContent(6,b);
    yvalData[6] = EffEtaPtData->GetBinContent(7,b);
    yerrlData[0] = EffErrLowEtaPtData->GetBinContent(1,b); 
    yerrlData[1] = EffErrLowEtaPtData->GetBinContent(2,b);
    yerrlData[2] = EffErrLowEtaPtData->GetBinContent(3,b);
    yerrlData[3] = EffErrLowEtaPtData->GetBinContent(4,b);
    yerrlData[4] = EffErrLowEtaPtData->GetBinContent(5,b);
    yerrlData[5] = EffErrLowEtaPtData->GetBinContent(6,b);
    yerrlData[6] = EffErrLowEtaPtData->GetBinContent(7,b);
    yerrhData[0] = EffErrHighEtaPtData->GetBinContent(1,b); 
    yerrhData[1] = EffErrHighEtaPtData->GetBinContent(2,b);
    yerrhData[2] = EffErrHighEtaPtData->GetBinContent(3,b);
    yerrhData[3] = EffErrHighEtaPtData->GetBinContent(4,b);
    yerrhData[4] = EffErrHighEtaPtData->GetBinContent(5,b);
    yerrhData[5] = EffErrHighEtaPtData->GetBinContent(6,b);
    yerrhData[6] = EffErrHighEtaPtData->GetBinContent(7,b);

    yvalMC[0] = EffEtaPtMC->GetBinContent(1,b); 
    yvalMC[1] = EffEtaPtMC->GetBinContent(2,b);
    yvalMC[2] = EffEtaPtMC->GetBinContent(3,b);
    yvalMC[3] = EffEtaPtMC->GetBinContent(4,b);
    yvalMC[4] = EffEtaPtMC->GetBinContent(5,b);
    yvalMC[5] = EffEtaPtMC->GetBinContent(6,b);
    yvalMC[6] = EffEtaPtMC->GetBinContent(7,b);
    yerrlMC[0] = EffErrLowEtaPtMC->GetBinContent(1,b); 
    yerrlMC[1] = EffErrLowEtaPtMC->GetBinContent(2,b);
    yerrlMC[2] = EffErrLowEtaPtMC->GetBinContent(3,b);
    yerrlMC[3] = EffErrLowEtaPtMC->GetBinContent(4,b);
    yerrlMC[4] = EffErrLowEtaPtMC->GetBinContent(5,b);
    yerrlMC[5] = EffErrLowEtaPtMC->GetBinContent(6,b);
    yerrlMC[6] = EffErrLowEtaPtMC->GetBinContent(7,b);
    yerrhMC[0] = EffErrHighEtaPtMC->GetBinContent(1,b); 
    yerrhMC[1] = EffErrHighEtaPtMC->GetBinContent(2,b);
    yerrhMC[2] = EffErrHighEtaPtMC->GetBinContent(3,b);
    yerrhMC[3] = EffErrHighEtaPtMC->GetBinContent(4,b);
    yerrhMC[4] = EffErrHighEtaPtMC->GetBinContent(5,b);
    yerrhMC[5] = EffErrHighEtaPtMC->GetBinContent(6,b);
    yerrhMC[6] = EffErrHighEtaPtMC->GetBinContent(7,b);

    //make tgraphs
    TGraphAsymmErrors *EfficiencyGraphMC = new TGraphAsymmErrors(n,xval,yvalMC,xerr,xerr,yerrlMC,yerrhMC);
    TGraphAsymmErrors *EfficiencyGraphData = new TGraphAsymmErrors(n,xval,yvalData,xerr,xerr,yerrlData,yerrhData);

    cv = new TCanvas("cv", "Canvas", 800, 600);
    legend = new TLegend( 0.5, 0.20, 0.8, 0.40);
    legend->SetTextSize(0.03);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(EfficiencyGraphData, legendLabelData.c_str(), "LP");
    legend->AddEntry(EfficiencyGraphMC, legendLabelMC.c_str(), "F");
    EfficiencyGraphData->SetTitle("");
    EfficiencyGraphData->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    EfficiencyGraphData->GetXaxis()->SetTitleOffset(1.05);
    EfficiencyGraphData->GetYaxis()->SetTitle("RECO+ID+ISO+SIP Efficiency");
    EfficiencyGraphData->GetYaxis()->SetTitleOffset(1.2);
    EfficiencyGraphData->GetXaxis()->SetRangeUser(0.0,100);
    EfficiencyGraphData->GetYaxis()->SetRangeUser(0.0,1.0);
    EfficiencyGraphData->SetLineColor(kBlack);
    EfficiencyGraphData->SetLineWidth(2);
    EfficiencyGraphData->SetMarkerStyle(20);
    EfficiencyGraphData->SetMarkerColor(kBlack);
    EfficiencyGraphData->SetMarkerSize(1.25);
    EfficiencyGraphMC->SetLineColor(kBlue);
    EfficiencyGraphMC->SetMarkerColor(kBlue);
    EfficiencyGraphMC->SetFillColor(kAzure+1);
    EfficiencyGraphData->Draw("AP");
    EfficiencyGraphMC->Draw("2");
    EfficiencyGraphData->Draw("P");
    legend->Draw();


    // Print Fit Values
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(.04);
    tex->SetTextFont(2);
    tex->Draw();
    tex->SetTextSize(0.030);
    tex->SetTextColor(kBlack);
    if (dataset == 2011) tex->DrawLatex(0.400, 0.92, "CMS Preliminary 2011 #sqrt{s} = 7 TeV, L=5.1 fb^{-1}");
    if(dataset == 2012) tex->DrawLatex(0.400, 0.92, "CMS Preliminary 2012 #sqrt{s} = 8 TeV, L=19.6 fb^{-1}");
    tex->DrawLatex(0.500, 0.42, binTexLabel.c_str() );
    


    cv->SaveAs( ("ElectronEfficiencyVsPt"+binLabel+label+".gif").c_str());

  }


}


void plotEfficiencyDataVsMC(string filenameData, 
                            string filenameMC, 
                            Int_t dataset,
                            string Label)
{
  plotEfficiencyComparison(filenameData,filenameMC,"Data","MC Simulation",dataset,Label);
}




void plotEfficiencyScaleFactor(string filename, 
                              string filename2, 
                              string legendLabel,
                              string legendLabel2,
                               Int_t dataset,
                              string Label)
{
  gBenchmark->Start("printMuonWPEff");
    
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TFile file(filename.c_str(), "READ");
  TFile *file2 = 0;
  if (filename2 != "") file2 = new TFile(filename2.c_str(), "READ");
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   
  
  TH2D* ScaleFactors   = ( TH2D*)file.Get("heff_electron_selection");
  TH2D* ScaleFactors2  = 0;
  if (file2) ScaleFactors2 = ( TH2D*)file2->Get("heff_electron_selection");

  assert(ScaleFactors);
  if (file2) assert(ScaleFactors2);

  TCanvas *cv = 0;
  TLegend *legend = 0;

  const UInt_t n = 7;
  Double_t xval[n], xerr[n];
  xval[0] = 8.5; xerr[0] = 1.5;
  xval[1] = 12.5; xerr[1] = 2.5;
  xval[2] = 17.5; xerr[2] = 2.5;
  xval[3] = 25; xerr[3] = 5;
  xval[4] = 35; xerr[4] = 5;
  xval[5] = 45; xerr[5] = 5;
  xval[6] = 75; xerr[6] = 25;

  Double_t yval[n], yerrl[n], yerrh[n];
  Double_t yval2[n], yerrl2[n], yerrh2[n];

  for (UInt_t b=1; b < 6; ++b) {

    string binLabel;
    if (b==1) binLabel = "Eta0p0To0p8";
    if (b==2) binLabel = "Eta0p8To1p4442";
    if (b==3) binLabel = "Eta1p4442To1p566";
    if (b==4) binLabel = "Eta1p566To2p0";
    if (b==5) binLabel = "Eta2p0To2p5";

    string binTexLabel;
    if (b==1) binTexLabel = "Probe in Barrel: 0 < |#eta_{probe}| < 0.8";
    if (b==2) binTexLabel = "Probe in Barrel: 0.8 < |#eta_{probe}| < 1.4442";
    if (b==3) binTexLabel = "Probe in Crack: 1.4442 < |#eta_{probe}| < 1.566";
    if (b==4) binTexLabel = "Probe in Endcap: 1.566 < |#eta_{probe}| < 2.0";
    if (b==5) binTexLabel = "Probe in Endcap: 2.0 < |#eta_{probe}| < 2.5";

    yval[0] = ScaleFactors->GetBinContent(1,b); 
    yval[1] = ScaleFactors->GetBinContent(2,b);
    yval[2] = ScaleFactors->GetBinContent(3,b);
    yval[3] = ScaleFactors->GetBinContent(4,b);
    yval[4] = ScaleFactors->GetBinContent(5,b);
    yval[5] = ScaleFactors->GetBinContent(6,b);
    yval[6] = ScaleFactors->GetBinContent(7,b);
    yerrl[0] = ScaleFactors->GetBinError(1,b) + 0.00; 
    yerrl[1] = ScaleFactors->GetBinError(2,b) + 0.00;
    yerrl[2] = ScaleFactors->GetBinError(3,b) + 0.00;
    yerrl[3] = ScaleFactors->GetBinError(4,b) + 0.00;
    yerrl[4] = ScaleFactors->GetBinError(5,b) + 0.00;
    yerrl[5] = ScaleFactors->GetBinError(6,b) + 0.00;
    yerrl[6] = ScaleFactors->GetBinError(7,b) + 0.00;
    yerrh[0] = ScaleFactors->GetBinError(1,b) + 0.00; 
    yerrh[1] = ScaleFactors->GetBinError(2,b) + 0.00;
    yerrh[2] = ScaleFactors->GetBinError(3,b) + 0.00;
    yerrh[3] = ScaleFactors->GetBinError(4,b) + 0.00;
    yerrh[4] = ScaleFactors->GetBinError(5,b) + 0.00;
    yerrh[5] = ScaleFactors->GetBinError(6,b) + 0.00;
    yerrh[6] = ScaleFactors->GetBinError(7,b) + 0.00;

    if (ScaleFactors2) {
      yval2[0] = ScaleFactors2->GetBinContent(1,b); 
      yval2[1] = ScaleFactors2->GetBinContent(2,b);
      yval2[2] = ScaleFactors2->GetBinContent(3,b);
      yval2[3] = ScaleFactors2->GetBinContent(4,b);
      yval2[4] = ScaleFactors2->GetBinContent(5,b);
      yval2[5] = ScaleFactors2->GetBinContent(6,b);
      yval2[6] = ScaleFactors2->GetBinContent(7,b);
      yerrl2[0] = ScaleFactors2->GetBinError(1,b) + 0.00; 
      yerrl2[1] = ScaleFactors2->GetBinError(2,b) + 0.00;
      yerrl2[2] = ScaleFactors2->GetBinError(3,b) + 0.00;
      yerrl2[3] = ScaleFactors2->GetBinError(4,b) + 0.00;
      yerrl2[4] = ScaleFactors2->GetBinError(5,b) + 0.00;
      yerrl2[5] = ScaleFactors2->GetBinError(6,b) + 0.00;
      yerrl2[6] = ScaleFactors2->GetBinError(7,b) + 0.00;
      yerrh2[0] = ScaleFactors2->GetBinError(1,b) + 0.00; 
      yerrh2[1] = ScaleFactors2->GetBinError(2,b) + 0.00;
      yerrh2[2] = ScaleFactors2->GetBinError(3,b) + 0.00;
      yerrh2[3] = ScaleFactors2->GetBinError(4,b) + 0.00;
      yerrh2[4] = ScaleFactors2->GetBinError(5,b) + 0.00;
      yerrh2[5] = ScaleFactors2->GetBinError(6,b) + 0.00;
      yerrh2[6] = ScaleFactors2->GetBinError(7,b) + 0.00;

    }

    //make tgraphs
    TGraphAsymmErrors *ScaleFactorGraph = new TGraphAsymmErrors(n,xval,yval,xerr,xerr,yerrl,yerrh);
    TGraphAsymmErrors *ScaleFactorGraph2 = 0;
    if (ScaleFactors2) {
      ScaleFactorGraph2 = new TGraphAsymmErrors(n,xval,yval2,xerr,xerr,yerrl2,yerrh2);
    }


    cv = new TCanvas("cv", "Canvas", 800, 600);
    legend = new TLegend( 0.5, 0.20, 0.8, 0.40);
    legend->SetTextSize(0.03);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(ScaleFactorGraph, legendLabel.c_str(), "F");
    if (ScaleFactors2) legend->AddEntry(ScaleFactorGraph2, legendLabel2.c_str(), "LP");
    ScaleFactorGraph->SetTitle("");
    ScaleFactorGraph->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    ScaleFactorGraph->GetXaxis()->SetTitleOffset(1.05);
    ScaleFactorGraph->GetYaxis()->SetTitle("RECO+ID+ISO+SIP Efficiency Scale Factor");
    ScaleFactorGraph->GetYaxis()->SetTitleOffset(1.2);
    ScaleFactorGraph->GetXaxis()->SetRangeUser(0.0,100);
    ScaleFactorGraph->GetYaxis()->SetRangeUser(0.7,1.1);
    ScaleFactorGraph->SetLineColor(kBlue);
    ScaleFactorGraph->SetMarkerColor(kBlue);
    ScaleFactorGraph->SetFillColor(kAzure+1);
    ScaleFactorGraph->Draw("A2");

    if (ScaleFactors2)ScaleFactorGraph2->SetLineColor(kBlack);
    if (ScaleFactors2)ScaleFactorGraph2->SetLineWidth(2);
    if (ScaleFactors2)ScaleFactorGraph2->SetMarkerStyle(20);
    if (ScaleFactors2)ScaleFactorGraph2->SetMarkerColor(kBlack);
    if (ScaleFactors2)ScaleFactorGraph2->SetMarkerSize(1.25);
    if (ScaleFactors2) ScaleFactorGraph2->Draw("P");
    legend->Draw();


    // Print Fit Values
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(.04);
    tex->SetTextFont(2);
    tex->Draw();
    tex->SetTextSize(0.030);
    tex->SetTextColor(kBlack);
    if (dataset == 2011) tex->DrawLatex(0.400, 0.92, "CMS Preliminary 2011 #sqrt{s} = 7 TeV, L=5.1 fb^{-1}");
    if (dataset == 2012) tex->DrawLatex(0.400, 0.92, "CMS Preliminary 2012 #sqrt{s} = 8 TeV, L=19.6 fb^{-1}");

    tex->DrawLatex(0.500, 0.42, binTexLabel.c_str() );
    


    cv->SaveAs( ("ElectronEfficiencyScaleFactorVsPt"+binLabel+label+".gif").c_str());

  }


}





void plotEfficiency() {

//   plotEfficiencyComparisons("results/Data2012_EleHZZMoriond2013IDGivenIsoWithTightTag/basic2_76_106_PtAbove20_2012DGood/eff.root",
//                             "results/Data2012_EleHZZMoriond2013IDGivenIsoWithTightTag/basic2_76_106_PtAbove20_2012DBad/eff.root",
//                             "Run2012D Good",
//                             "Run2012D Bad",
//                             "ElectronHZZHCP2012IDGivenIso");


//   plotEfficiencyDataVsMC( "results/Data2012_EleHZZMoriond2013WPWithCombinedMethods/effData.root",
//                           "results/Data2012_EleHZZMoriond2013WPWithCombinedMethods/effMC.root",
//                           "Moriond2013" );
  
//   plotEfficiencyScaleFactor( "efficiency_results_EleHZZMoriond2013WPMixed_Moriond2013.root", 
//                              "/afs/cern.ch/user/s/sixie/work/public/HZZ4l/auxiliar/2012/efficiency_results_EleHZZHCP2012WPMixed_HCP2012.UseFactorizedSFFor7To10Bin.root",
//                              "Moriond2013");

  //******************************************************************************************
  // Run1 Legacy Paper vs HCP WP
  //******************************************************************************************

  //2012 MC
  plotEfficiencyComparison( "results/Data2012_EleHZZRun1LegacyPaperWPWithCombinedMethods/effMC.root",
                          "results/Data2012_EleHZZMoriond2013WPWithCombinedMethods/effMC.root",
                            "Run1 Legacy Paper WP",
                            "HCP 2012 WP",
                            2012,
                            "Run1LegacyPaperVsHCP2012WP_2012_MC" );
  //2012 Data
  plotEfficiencyComparison( "results/Data2012_EleHZZRun1LegacyPaperWPWithCombinedMethods/effData.root",
                          "results/Data2012_EleHZZMoriond2013WPWithCombinedMethods/effData.root",
                            "Run1 Legacy Paper WP (Rereco)",
                            "HCP 2012 WP (PromptReco)",
                            2012,
                            "Run1LegacyPaperVsHCP2012WP_2012_Data" );


  //2011 MC
  plotEfficiencyComparison( "results/Data2011_EleHZZRun1LegacyPaperWPWithCombinedMethods/effMC.root",
                            "results/Data2011_EleHZZHCP2012WPWithCombinedMethods/effMC.root",
                            "Run1 Legacy Paper WP",
                            "HCP 2012 WP",
                            2011,
                            "Run1LegacyPaperVsHCP2012WP_2011_MC" );
  //2011 Data
  plotEfficiencyComparison( "results/Data2011_EleHZZRun1LegacyPaperWPWithCombinedMethods/effData.root",
                          "results/Data2011_EleHZZHCP2012WPWithCombinedMethods/effData.root",
                            "Run1 Legacy Paper WP",
                            "HCP 2012 WP",
                            2011,
                            "Run1LegacyPaperVsHCP2012WP_2011_Data" );


  //******************************************************************************************
  // Run1 Legacy Paper Data vs MC
  //******************************************************************************************
  //2012
  plotEfficiencyComparison( "results/Data2012_EleHZZRun1LegacyPaperWPWithCombinedMethods/effData.root",
                            "results/Data2012_EleHZZRun1LegacyPaperWPWithCombinedMethods/effMC.root",
                            "Data",
                            "MC Simulation",
                            2012,
                            "Run1LegacyPaperWPDataVsMC_2012" );

  plotEfficiencyComparison( "results/Data2011_EleHZZRun1LegacyPaperWPWithCombinedMethods/effData.root",
                            "results/Data2011_EleHZZRun1LegacyPaperWPWithCombinedMethods/effMC.root",
                            "Data",
                            "MC Simulation",
                            2011,
                            "Run1LegacyPaperWPDataVsMC_2011" );

  //******************************************************************************************
  // Run1 Legacy Paper Scale Factors
  //******************************************************************************************
  plotEfficiencyScaleFactor( "efficiency_results_EleHZZRun1LegacyPaperWPFromZeeWithLowestPtBinFactorized_Run1LegacyPaper_2011Data.root", 
                             "efficiency_results_EleHZZHCP2012WPFromZeeWithLowestPtBinFactorized_HCP2012_2011Data.root",
                             "Run1 Legacy Paper WP",
                             "HCP 2012 WP",
                             2011,
                             "Run1LegacyPaperWPVsHCP2012WP_2011");
  
  plotEfficiencyScaleFactor( "efficiency_results_EleHZZRun1LegacyPaperWPFromZeeWithLowestPtBinFactorized_Run1LegacyPaper_2012Data.root", 
                             "efficiency_results_EleHZZMoriond2013WPFromZeeWithLowestPtBinFactorized_Moriond2013_2012Data.root",
                             "Run1 Legacy Paper WP",
                             "HCP 2012 WP",
                             2012,
                             "Run1LegacyPaperWPVsHCP2012WP_2012");


  plotEfficiencyScaleFactor( "efficiency_results_EleHZZRun1LegacyPaperWPFromZeeWithLowestPtBinFactorized_Run1LegacyPaper_2011Data.root", 
                             "",
                             "Run1 Legacy Paper WP",
                             "HCP 2012 WP",
                             2011,
                             "Run1LegacyPaperWP_2011");
  
  plotEfficiencyScaleFactor( "efficiency_results_EleHZZRun1LegacyPaperWPFromZeeWithLowestPtBinFactorized_Run1LegacyPaper_2012Data.root", 
                             "",
                             "Run1 Legacy Paper WP",
                             "HCP 2012 WP",
                             2012,
                             "Run1LegacyPaperWP_2012");


  
}

