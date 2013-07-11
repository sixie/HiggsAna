
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TGraphAsymmErrors.h>      // graphs
#include <TH2F.h>                   // 2D histograms
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#endif

//=== Functions =================================================================================================
void doEfficiencyScaleFactors(string StandardEffBinsFilename,
                              string DataFilename, 
                              string DataLowPtBinsFilename, 
                              string DataLowPtBinsZeeGammaFilename, 
                              string MCFilename, 
                              string MCLowPtBinsFilename, 
                              string MCLowPtBinsZeeGammaFilenam, 
                              string outputDir, 
                              string histName,
                              string Label);

void doEfficiencyScaleFactors_LowestPtBinUsesZeeFactorized(string StandardEffBinsFilename,
                                                           string DataFilename, 
                                                           string DataLowPtBinsFilename, 
                                                           string MCFilename, 
                                                           string MCLowPtBinsFilename,
                                                           Int_t dataset,
                                                           string outputDir, 
                                                           string histName,
                                                           string Label);

void doEfficiencyScaleFactorsOneSample(string DataFilename = "Data_EleWPEffTP/basic2_76_106/eff.root", 
                                       string MCFilename = "Summer11_Zee_EleWPEffTP/basic_76_106/eff.root" , 
                                       string outputDir = "Data_EleWPEffTP/basic2_76_106/", 
                                       string histName = "h2_results_electron_selection",
                                       string Label = "SmurfV6");

//=== MAIN MACRO =================================================================================================
void makeEfficiencyScaleFactors(string DataFilename, 
                                string MCFilename , 
                                string outputDir, 
                                string histName,
                                string Label)
{
  doEfficiencyScaleFactorsOneSample(DataFilename, MCFilename, outputDir, histName, Label);
}

void makeEfficiencyScaleFactors(string StandardEffBinsFilename,
                                string DataFilename, 
                                string DataLowPtBinsFilename, 
                                string DataLowPtBinsZeeGammaFilename,
                                string MCFilename, 
                                string MCLowPtBinsFilename, 
                                string MCLowPtBinsZeeGammaFilename,
                                string outputDir, 
                                string histName,
                                string Label)
{
  doEfficiencyScaleFactors(StandardEffBinsFilename, DataFilename, DataLowPtBinsFilename, DataLowPtBinsZeeGammaFilename,
                           MCFilename, MCLowPtBinsFilename, MCLowPtBinsZeeGammaFilename,
                           outputDir, histName, Label);
}


void makeEfficiencyScaleFactors(string StandardEffBinsFilename,
                                string DataFilename, 
                                string DataLowPtBinsFilename, 
                                string MCFilename, 
                                string MCLowPtBinsFilename, 
                                Int_t dataset,
                                string outputDir, 
                                string histName,
                                string Label)
{
  doEfficiencyScaleFactors_LowestPtBinUsesZeeFactorized(StandardEffBinsFilename, 
                                                        DataFilename, DataLowPtBinsFilename,
                                                        MCFilename, MCLowPtBinsFilename,
                                                        dataset, outputDir, histName, Label);
}


void doEfficiencyScaleFactors( string StandardEffBinsFilename,
                               string DataFilename, 
                               string DataLowPtBinsFilename, 
                               string DataLowPtBinsZeeGammaFilename,
                               string MCFilename, 
                               string MCLowPtBinsFilename, 
                               string MCLowPtBinsZeeGammaFilename,
                               string outputDir, 
                               string histName,
                               string Label)
{
  gBenchmark->Start("printMuonWPEff");
    
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TString outfname = (outputDir + "/eff_table.txt").c_str();
  TFile standardeffbinsfile(StandardEffBinsFilename.c_str());
  TFile mcfile(MCFilename.c_str());
  TFile datafile(DataFilename.c_str());
  TFile lowpt_mcfile(MCLowPtBinsFilename.c_str());
  TFile lowpt_datafile(DataLowPtBinsFilename.c_str());
  TFile lowptZeeGamma_mcfile(MCLowPtBinsZeeGammaFilename.c_str());
  TFile lowptZeeGamma_datafile(DataLowPtBinsZeeGammaFilename.c_str());
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   
  
  //Load Standard Efficiency Bins histogram
  TH2F *hStandardEff=0;
  hStandardEff  = (TH2F*)standardeffbinsfile.Get("hEffEtaPt");

  //Load Standard Tag and Probe efficiencies
  TH2F *hMCEff=0,    *hMCErrl=0,    *hMCErrh=0;
  TH2F *hDataEff=0, *hDataErrl=0, *hDataErrh=0;
  hMCEff  = (TH2F*)mcfile.Get("hEffEtaPt");
  hMCErrl = (TH2F*)mcfile.Get("hErrlEtaPt");
  hMCErrh = (TH2F*)mcfile.Get("hErrhEtaPt");
  hDataEff  = (TH2F*)datafile.Get("hEffEtaPt");
  hDataErrl = (TH2F*)datafile.Get("hErrlEtaPt");
  hDataErrh = (TH2F*)datafile.Get("hErrhEtaPt");

  //Load Low Pt Tag and Probe efficiencies
  TH2F *hMCEff_LowPt=0,    *hMCErrl_LowPt=0,    *hMCErrh_LowPt=0;
  TH2F *hDataEff_LowPt=0, *hDataErrl_LowPt=0, *hDataErrh_LowPt=0;
  hMCEff_LowPt  = (TH2F*)lowpt_mcfile.Get("hEffEtaPt");
  hMCErrl_LowPt = (TH2F*)lowpt_mcfile.Get("hErrlEtaPt");
  hMCErrh_LowPt = (TH2F*)lowpt_mcfile.Get("hErrhEtaPt");
  hDataEff_LowPt  = (TH2F*)lowpt_datafile.Get("hEffEtaPt");
  hDataErrl_LowPt = (TH2F*)lowpt_datafile.Get("hErrlEtaPt");
  hDataErrh_LowPt = (TH2F*)lowpt_datafile.Get("hErrhEtaPt");

  //Load Low Pt ZeeGamma Tag and Probe efficiencies
  TH2F *hMCEff_LowPtZeeGamma=0,    *hMCErrl_LowPtZeeGamma=0,    *hMCErrh_LowPtZeeGamma=0;
  TH2F *hDataEff_LowPtZeeGamma=0, *hDataErrl_LowPtZeeGamma=0, *hDataErrh_LowPtZeeGamma=0;
  hMCEff_LowPtZeeGamma  = (TH2F*)lowptZeeGamma_mcfile.Get("hEffEtaPt");
  hMCErrl_LowPtZeeGamma = (TH2F*)lowptZeeGamma_mcfile.Get("hErrlEtaPt");
  hMCErrh_LowPtZeeGamma = (TH2F*)lowptZeeGamma_mcfile.Get("hErrhEtaPt");
  hDataEff_LowPtZeeGamma  = (TH2F*)lowptZeeGamma_datafile.Get("hEffEtaPt");
  hDataErrl_LowPtZeeGamma = (TH2F*)lowptZeeGamma_datafile.Get("hErrlEtaPt");
  hDataErrh_LowPtZeeGamma = (TH2F*)lowptZeeGamma_datafile.Get("hErrhEtaPt");
//   if (!hMCEff_LowPtZeeGamma) {
//     hMCEff_LowPtZeeGamma  = (TH2F*)((TH2D*)lowptZeeGamma_mcfile.Get("hEffEtaPt"));
//     hMCErrl_LowPtZeeGamma = (TH2F*)((TH2D*)lowptZeeGamma_mcfile.Get("hErrlEtaPt"));
//     hMCErrh_LowPtZeeGamma = (TH2F*)((TH2D*)lowptZeeGamma_mcfile.Get("hErrhEtaPt"));
//   }
//   if (!hDataEff_LowPtZeeGamma) {
//     hDataEff_LowPtZeeGamma  = (TH2F*)((TH2D*)lowptZeeGamma_datafile.Get("hEffEtaPt"));
//     hDataErrl_LowPtZeeGamma = (TH2F*)((TH2D*)lowptZeeGamma_datafile.Get("hErrlEtaPt"));
//     hDataErrh_LowPtZeeGamma = (TH2F*)((TH2D*)lowptZeeGamma_datafile.Get("hErrhEtaPt"));    
//   }

  TGraphAsymmErrors* MCEffVsRho   = ( TGraphAsymmErrors*)mcfile.Get("grEffRho");
  TGraphAsymmErrors* DataEffVsRho = ( TGraphAsymmErrors*)datafile.Get("grEffRho");
 
  TGraphAsymmErrors* MCEffVsNVtx   = ( TGraphAsymmErrors*)mcfile.Get("grEffNPV");
  TGraphAsymmErrors* DataEffVsNVtx = ( TGraphAsymmErrors*)datafile.Get("grEffNPV");


  //--------------------------------------------------------------------------------------------------------------
  // Update root file histograms
  //==============================================================================================================   
  const Int_t nx = hStandardEff->GetNbinsX();
  const Int_t ny = hStandardEff->GetNbinsY();

  TFile *outputFile = new TFile(("efficiency_results"+label+".root").c_str(), "UPDATE");

  //Do Binning
  Double_t *ptbins = new Double_t[ny+1];
  Double_t *etabins = new Double_t[nx+1];
  for(Int_t iy=1; iy<=ny; iy++) {
      ptbins[iy-1] = hStandardEff->GetYaxis()->GetBinLowEdge(iy);
  }
  for(Int_t ix=1; ix<=nx; ix++) {
    etabins[ix-1]= hStandardEff->GetXaxis()->GetBinLowEdge(ix);
  }
  ptbins[ny] = 60;
  etabins[nx] = 2.5;
    
  TH2F *h2_results_selection = new TH2F(histName.c_str(),"",ny,ptbins,nx,etabins);
  for(Int_t iy=0; iy<=ny+2; iy++) {
    for(Int_t ix=0; ix<=nx+2; ix++) {
      h2_results_selection->SetCellContent(iy,ix, 1.0);
      h2_results_selection->SetCellError(iy,ix, 0.0);
    }
  }

  //save output efficiencies into root file
  TH2F *hMCEffCombined=new TH2F("hEffEtaPt","",ny,ptbins,nx,etabins); 
  TH2F *hMCErrlCombined=new TH2F("hErrlEtaPt","",ny,ptbins,nx,etabins);
  TH2F *hMCErrhCombined=new TH2F("hErrhEtaPt","",ny,ptbins,nx,etabins);
  TH2F *hDataEffCombined=new TH2F("hEffEtaPt","",ny,ptbins,nx,etabins);
  TH2F *hDataErrlCombined=new TH2F("hErrlEtaPt","",ny,ptbins,nx,etabins);
  TH2F *hDataErrhCombined=new TH2F("hErrhEtaPt","",ny,ptbins,nx,etabins);

  cout << "ptbins : ";
  for(Int_t i=0; i<ny+1;++i) {
    cout << ptbins[i] << " ";
  }
  cout << endl;
  cout << "etabins : ";
  for(Int_t i=0; i<nx+1;++i) {
    cout << etabins[i] << " ";
  }
  cout << endl;

  //--------------------------------------------------------------------------------------------------------------
  // Produce Text file table
  //==============================================================================================================   

  ofstream txtfile;
  txtfile.open(outfname.Data());
  assert(txtfile.is_open());
    
  
  txtfile << " pT        ";
  txtfile << " eta           ";
  txtfile << "    MC T&P           ";
  txtfile << "    Data T&P         ";  
  txtfile << "    Scale factor     ";
  txtfile << endl;
  txtfile << "----------------------------------------------------------------------------------------------------------------------------------------------------";
  txtfile << endl;
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {
      txtfile << "[" << setw(4) << hStandardEff->GetYaxis()->GetBinLowEdge(iy) << "," << setw(4) << hStandardEff->GetYaxis()->GetBinLowEdge(iy+1) << "]";
      txtfile << "[" << setw(6) << hStandardEff->GetXaxis()->GetBinLowEdge(ix) << "," << setw(6) << hStandardEff->GetXaxis()->GetBinLowEdge(ix+1) << "]";      
      ios_base::fmtflags flags = txtfile.flags();
      txtfile.precision(4);
      
      Double_t binXValue = hStandardEff->GetXaxis()->GetBinCenter(ix);
      Double_t binYValue = hStandardEff->GetYaxis()->GetBinCenter(iy);
      //cout << "Bin : " << binXValue << " " << binYValue << endl;

      TH2F *tmphistMCEff=0,    *tmphistMCErrl=0,    *tmphistMCErrh=0;
      TH2F *tmphistDataEff=0, *tmphistDataErrl=0, *tmphistDataErrh=0;
      
      //use LowPtZeeGamma
      if (hStandardEff->GetYaxis()->GetBinCenter(iy) < 10) {
        if (!(hDataEff_LowPtZeeGamma && hMCEff_LowPtZeeGamma)) { cout << "ZeeGamma eff hist not found\n"; assert(false);}
        tmphistMCEff = hMCEff_LowPtZeeGamma;
        tmphistMCErrl = hMCErrl_LowPtZeeGamma;
        tmphistMCErrh = hMCErrh_LowPtZeeGamma;
        tmphistDataEff= hDataEff_LowPtZeeGamma;
        tmphistDataErrl= hDataErrl_LowPtZeeGamma;
        tmphistDataErrh= hDataErrh_LowPtZeeGamma;
      } else if(hDataEff_LowPt && hMCEff_LowPt && hStandardEff->GetYaxis()->GetBinCenter(iy) < 20)  {
        if (!(hDataEff_LowPt && hMCEff_LowPt)) { cout << "Zee low pt eff hist not found\n"; assert(false);}
        tmphistMCEff = hMCEff_LowPt;
        tmphistMCErrl = hMCErrl_LowPt;
        tmphistMCErrh = hMCErrh_LowPt;
        tmphistDataEff= hDataEff_LowPt;
        tmphistDataErrl= hDataErrl_LowPt;
        tmphistDataErrh= hDataErrh_LowPt;
      } else {
        tmphistMCEff = hMCEff;
        tmphistMCErrl = hMCErrl;
        tmphistMCErrh = hMCErrh;
        tmphistDataEff= hDataEff;
        tmphistDataErrl= hDataErrl;
        tmphistDataErrh= hDataErrh;
      }

      Double_t mceff  = tmphistMCEff->GetCellContent(tmphistMCEff->GetXaxis()->FindFixBin(binXValue),tmphistMCEff->GetYaxis()->FindFixBin(binYValue));
      Double_t mcerrl = tmphistMCErrl->GetCellContent(tmphistMCErrl->GetXaxis()->FindFixBin(binXValue),tmphistMCErrl->GetYaxis()->FindFixBin(binYValue));
      Double_t mcerrh = tmphistMCErrh->GetCellContent(tmphistMCErrh->GetXaxis()->FindFixBin(binXValue),tmphistMCErrh->GetYaxis()->FindFixBin(binYValue));
      txtfile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);

      Double_t dataeff  = tmphistDataEff->GetCellContent(tmphistDataEff->GetXaxis()->FindFixBin(binXValue),tmphistDataEff->GetYaxis()->FindFixBin(binYValue));
      Double_t dataerrl = tmphistDataErrl->GetCellContent(tmphistDataErrl->GetXaxis()->FindFixBin(binXValue),tmphistDataErrl->GetYaxis()->FindFixBin(binYValue));
      Double_t dataerrh = tmphistDataErrh->GetCellContent(tmphistDataErrh->GetXaxis()->FindFixBin(binXValue),tmphistDataErrh->GetYaxis()->FindFixBin(binYValue));
      txtfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      
      Double_t scale     = (tmphistDataEff->GetCellContent(tmphistDataEff->GetXaxis()->FindFixBin(binXValue),tmphistDataEff->GetYaxis()->FindFixBin(binYValue)))/(tmphistMCEff->GetCellContent(tmphistMCEff->GetXaxis()->FindFixBin(binXValue),tmphistMCEff->GetYaxis()->FindFixBin(binYValue)));
      Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
      txtfile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);

      txtfile << endl;
      txtfile.flags(flags);

      cout << "Set " << iy << " " << ix << " " << scale << endl;
      h2_results_selection->SetCellContent(iy,ix, scale);
      h2_results_selection->SetCellError(iy,ix, (scaleErrl + scaleErrh)/2);
      hDataEffCombined->SetCellContent(iy,ix, dataeff);
      hDataErrlCombined->SetCellContent(iy,ix, dataerrl);
      hDataErrhCombined->SetCellContent(iy,ix, dataerrh);
      hMCEffCombined->SetCellContent(iy,ix, mceff);
      hMCErrlCombined->SetCellContent(iy,ix, mcerrl);
      hMCErrhCombined->SetCellContent(iy,ix, mcerrh);

      //fill overflow bins with the same values as last bin
      if (ix == nx) {
      cout << "Set " << iy << " " << ix+1 << " " << scale << endl;
        h2_results_selection->SetCellContent(iy,nx+1, scale);
        h2_results_selection->SetCellError(iy,nx+1, (scaleErrl + scaleErrh)/2);
        hDataEffCombined->SetCellContent(iy,nx+1, dataeff);
        hDataErrlCombined->SetCellContent(iy,nx+1, dataerrl);
        hDataErrhCombined->SetCellContent(iy,nx+1, dataerrh);
        hMCEffCombined->SetCellContent(iy,nx+1, mceff);
        hMCErrlCombined->SetCellContent(iy,nx+1, mcerrl);
        hMCErrhCombined->SetCellContent(iy,nx+1, mcerrh);
      }
    }
    if (iy == ny) {
      for(Int_t ix=1; ix<=nx; ix++) {
        Double_t binXValue = hStandardEff->GetXaxis()->GetBinCenter(ix);
        Double_t binYValue = hStandardEff->GetYaxis()->GetBinCenter(iy);

        TH2F *tmphistMCEff=0,    *tmphistMCErrl=0,    *tmphistMCErrh=0;
        TH2F *tmphistDataEff=0, *tmphistDataErrl=0, *tmphistDataErrh=0;
      
        //use LowPtZeeGamma
        if (hDataEff_LowPtZeeGamma && hMCEff_LowPtZeeGamma && hStandardEff->GetYaxis()->GetBinCenter(iy) < 10) {
          tmphistMCEff = hMCEff_LowPtZeeGamma;
          tmphistMCErrl = hMCErrl_LowPtZeeGamma;
          tmphistMCErrh = hMCErrh_LowPtZeeGamma;
          tmphistDataEff= hDataEff_LowPtZeeGamma;
          tmphistDataErrl= hDataErrl_LowPtZeeGamma;
          tmphistDataErrh= hDataErrh_LowPtZeeGamma;
        } else if(hDataEff_LowPt && hMCEff_LowPt && hStandardEff->GetYaxis()->GetBinCenter(iy) < 20)  {
          tmphistMCEff = hMCEff_LowPt;
          tmphistMCErrl = hMCErrl_LowPt;
          tmphistMCErrh = hMCErrh_LowPt;
          tmphistDataEff= hDataEff_LowPt;
          tmphistDataErrl= hDataErrl_LowPt;
          tmphistDataErrh= hDataErrh_LowPt;
        } else {
          tmphistMCEff = hMCEff;
          tmphistMCErrl = hMCErrl;
          tmphistMCErrh = hMCErrh;
          tmphistDataEff= hDataEff;
          tmphistDataErrl= hDataErrl;
          tmphistDataErrh= hDataErrh;
        }

        Double_t mceff  = tmphistMCEff->GetCellContent(tmphistMCEff->GetXaxis()->FindFixBin(binXValue),tmphistMCEff->GetYaxis()->FindFixBin(binYValue));
        Double_t mcerrl = tmphistMCErrl->GetCellContent(tmphistMCErrl->GetXaxis()->FindFixBin(binXValue),tmphistMCErrl->GetYaxis()->FindFixBin(binYValue));
        Double_t mcerrh = tmphistMCErrh->GetCellContent(tmphistMCErrh->GetXaxis()->FindFixBin(binXValue),tmphistMCErrh->GetYaxis()->FindFixBin(binYValue));
        Double_t dataeff  = tmphistDataEff->GetCellContent(tmphistDataEff->GetXaxis()->FindFixBin(binXValue),tmphistDataEff->GetYaxis()->FindFixBin(binYValue));
        Double_t dataerrl = tmphistDataErrl->GetCellContent(tmphistDataErrl->GetXaxis()->FindFixBin(binXValue),tmphistDataErrl->GetYaxis()->FindFixBin(binYValue));
        Double_t dataerrh = tmphistDataErrh->GetCellContent(tmphistDataErrh->GetXaxis()->FindFixBin(binXValue),tmphistDataErrh->GetYaxis()->FindFixBin(binYValue));
        Double_t scale     = (tmphistDataEff->GetCellContent(tmphistDataEff->GetXaxis()->FindFixBin(binXValue),tmphistDataEff->GetYaxis()->FindFixBin(binYValue)))/(tmphistMCEff->GetCellContent(tmphistMCEff->GetXaxis()->FindFixBin(binXValue),tmphistMCEff->GetYaxis()->FindFixBin(binYValue)));
        Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
        Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
        
        cout << "Set " << iy+1 << " " << ix << " " << scale << endl;
        h2_results_selection->SetCellContent(iy+1,ix, scale);
        h2_results_selection->SetCellError(iy+1,ix, (scaleErrl + scaleErrh)/2);
        hDataEffCombined->SetCellContent(iy+1,ix, dataeff);
        hDataErrlCombined->SetCellContent(iy+1,ix, dataerrl);
        hDataErrhCombined->SetCellContent(iy+1,ix, dataerrh);
        hMCEffCombined->SetCellContent(iy+1,ix, mceff);
        hMCErrlCombined->SetCellContent(iy+1,ix, mcerrl);
        hMCErrhCombined->SetCellContent(iy+1,ix, mcerrh);
        
        //fill overflow bins with the same values as last bin
        if (ix == nx) {
          cout << "Set " << iy+1 << " " << ix+1 << " " << scale << endl;
          h2_results_selection->SetCellContent(iy+1,nx+1, scale);
          h2_results_selection->SetCellError(iy+1,nx+1, (scaleErrl + scaleErrh)/2);
          hDataEffCombined->SetCellContent(iy+1,nx+1, dataeff);
          hDataErrlCombined->SetCellContent(iy+1,nx+1, dataerrl);
          hDataErrhCombined->SetCellContent(iy+1,nx+1, dataerrh);
          hMCEffCombined->SetCellContent(iy+1,nx+1, mceff);
          hMCErrlCombined->SetCellContent(iy+1,nx+1, mcerrl);
          hMCErrhCombined->SetCellContent(iy+1,nx+1, mcerrh);
        }
      }
    }

    txtfile << endl;
  }
  txtfile.close();
  
  cout << outfname << " created!" << endl;
  



  //--------------------------------------------------------------------------------------------------------------
  // Create TEX table
  //==============================================================================================================   

  ofstream texfile;
  texfile.open((outputDir + "/eff_table.tex").c_str());
  assert(texfile.is_open());

  texfile << " \\begin{table}[!ht]" << endl;
  texfile << " \\begin{center} " << endl;
  texfile << " \\begin{tabular}{|c|c|c|c|}" << endl;
  texfile << " \\hline\n";


  texfile << " $p_{T}$ / $\\eta$ bin    &  Monte Carlo Efficiency    &  Data Efficiency   &  MC to Data Scale Factor \\\\  ";
  texfile << " \\hline           ";
  texfile << endl;
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {

      Double_t binXValue = hStandardEff->GetXaxis()->GetBinCenter(ix);
      Double_t binYValue = hStandardEff->GetYaxis()->GetBinCenter(iy);

      string binLabel = Form("$%5.1f < p_{T} \\le %5.1f$ , $%5.1f  \\le |\\eta| < %5.1f$", 
                             hStandardEff->GetYaxis()->GetBinLowEdge(iy), hStandardEff->GetYaxis()->GetBinLowEdge(iy+1),
                             hStandardEff->GetXaxis()->GetBinLowEdge(ix), hStandardEff->GetXaxis()->GetBinLowEdge(ix+1));
      if (iy == ny) {
        binLabel = Form("$%5.1f < p_{T} $ , $%5.1f  \\le |\\eta| < %5.1f$", 
                        hStandardEff->GetYaxis()->GetBinLowEdge(iy), 
                        hStandardEff->GetXaxis()->GetBinLowEdge(ix), hStandardEff->GetXaxis()->GetBinLowEdge(ix+1));
      }
      
      TH2F *tmphistMCEff=0,    *tmphistMCErrl=0,    *tmphistMCErrh=0;
      TH2F *tmphistDataEff=0, *tmphistDataErrl=0, *tmphistDataErrh=0;
      
      //use LowPtZeeGamma
      if (hDataEff_LowPtZeeGamma && hMCEff_LowPtZeeGamma && hStandardEff->GetYaxis()->GetBinCenter(iy) < 10) {
        tmphistMCEff = hMCEff_LowPtZeeGamma;
        tmphistMCErrl = hMCErrl_LowPtZeeGamma;
        tmphistMCErrh = hMCErrh_LowPtZeeGamma;
        tmphistDataEff= hDataEff_LowPtZeeGamma;
        tmphistDataErrl= hDataErrl_LowPtZeeGamma;
        tmphistDataErrh= hDataErrh_LowPtZeeGamma;
      } else if(hDataEff_LowPt && hMCEff_LowPt && hStandardEff->GetYaxis()->GetBinCenter(iy) < 20)  {
        tmphistMCEff = hMCEff_LowPt;
        tmphistMCErrl = hMCErrl_LowPt;
        tmphistMCErrh = hMCErrh_LowPt;
        tmphistDataEff= hDataEff_LowPt;
        tmphistDataErrl= hDataErrl_LowPt;
        tmphistDataErrh= hDataErrh_LowPt;
      } else {
        tmphistMCEff = hMCEff;
        tmphistMCErrl = hMCErrl;
        tmphistMCErrh = hMCErrh;
        tmphistDataEff= hDataEff;
        tmphistDataErrl= hDataErrl;
        tmphistDataErrh= hDataErrh;
      }



      texfile << binLabel;
      texfile << "   &   ";

      ios_base::fmtflags flags = texfile.flags();
      texfile.precision(4);
      
      Double_t mceff  = tmphistMCEff->GetCellContent(tmphistMCEff->GetXaxis()->FindFixBin(binXValue),tmphistMCEff->GetYaxis()->FindFixBin(binYValue));
      Double_t mcerrl = tmphistMCErrl->GetCellContent(tmphistMCErrl->GetXaxis()->FindFixBin(binXValue),tmphistMCErrl->GetYaxis()->FindFixBin(binYValue));
      Double_t mcerrh = tmphistMCErrh->GetCellContent(tmphistMCErrh->GetXaxis()->FindFixBin(binXValue),tmphistMCErrh->GetYaxis()->FindFixBin(binYValue));
      texfile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);
      texfile << "   &   ";

      Double_t dataeff  = tmphistDataEff->GetCellContent(tmphistDataEff->GetXaxis()->FindFixBin(binXValue),tmphistDataEff->GetYaxis()->FindFixBin(binYValue));
      Double_t dataerrl = tmphistDataErrl->GetCellContent(tmphistDataErrl->GetXaxis()->FindFixBin(binXValue),tmphistDataErrl->GetYaxis()->FindFixBin(binYValue));
      Double_t dataerrh = tmphistDataErrh->GetCellContent(tmphistDataErrh->GetXaxis()->FindFixBin(binXValue),tmphistDataErrh->GetYaxis()->FindFixBin(binYValue));
      texfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      texfile << "   &   ";
     
      Double_t scale     = (tmphistDataEff->GetCellContent(tmphistDataEff->GetXaxis()->FindFixBin(binXValue),tmphistDataEff->GetYaxis()->FindFixBin(binYValue)))/(tmphistMCEff->GetCellContent(tmphistMCEff->GetXaxis()->FindFixBin(binXValue),tmphistMCEff->GetYaxis()->FindFixBin(binYValue)));
      Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
      texfile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);
      texfile << "   \\\\   ";
      texfile << endl;
      texfile << "\\hline";
      texfile << endl;

    }
  }
  texfile << "\\end{tabular}" << endl;
  texfile << "\\caption{CAPTION.}" << endl;
  texfile << "\\label{tab:eff_ele_offline}" << endl;
  texfile << "\\end{center}" << endl;
  texfile << "\\end{table}" << endl;
  texfile.close();
  cout << outputDir + "/eff_table.tex" << " created!" << endl;




  //--------------------------------------------------------------------------------------------------------------
  // Create TWIKI table
  //==============================================================================================================   

  ofstream twikifile;
  twikifile.open((outputDir + "/eff_table.twiki").c_str());
  assert(twikifile.is_open());
  

  twikifile << "| *pT bin* | *eta bin* | *MC Efficiency* | *Data Efficiency* | *ScaleFactor* | \n";
  
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {

      Double_t binXValue = hStandardEff->GetXaxis()->GetBinCenter(ix);
      Double_t binYValue = hStandardEff->GetYaxis()->GetBinCenter(iy);

      string binLabel = Form("| *%4.1f < pT \\le %4.1f* | *%3.1f  <= eta < %3.1f* | ", 
                             hStandardEff->GetYaxis()->GetBinLowEdge(iy), hStandardEff->GetYaxis()->GetBinLowEdge(iy+1),
                             hStandardEff->GetXaxis()->GetBinLowEdge(ix), hStandardEff->GetXaxis()->GetBinLowEdge(ix+1));
      if (iy == ny) {
        binLabel = Form("| *%4.1f < pT* | *%3.1f  <= eta < %3.1f* | ", 
                        hStandardEff->GetYaxis()->GetBinLowEdge(iy), 
                        hStandardEff->GetXaxis()->GetBinLowEdge(ix), hStandardEff->GetXaxis()->GetBinLowEdge(ix+1));
      }
      
      TH2F *tmphistMCEff=0,    *tmphistMCErrl=0,    *tmphistMCErrh=0;
      TH2F *tmphistDataEff=0, *tmphistDataErrl=0, *tmphistDataErrh=0;
      
      //use LowPtZeeGamma
      if (hDataEff_LowPtZeeGamma && hMCEff_LowPtZeeGamma && hStandardEff->GetYaxis()->GetBinCenter(iy) < 10) {
        tmphistMCEff = hMCEff_LowPtZeeGamma;
        tmphistMCErrl = hMCErrl_LowPtZeeGamma;
        tmphistMCErrh = hMCErrh_LowPtZeeGamma;
        tmphistDataEff= hDataEff_LowPtZeeGamma;
        tmphistDataErrl= hDataErrl_LowPtZeeGamma;
        tmphistDataErrh= hDataErrh_LowPtZeeGamma;
      } else if(hDataEff_LowPt && hMCEff_LowPt && hStandardEff->GetYaxis()->GetBinCenter(iy) < 20)  {
        tmphistMCEff = hMCEff_LowPt;
        tmphistMCErrl = hMCErrl_LowPt;
        tmphistMCErrh = hMCErrh_LowPt;
        tmphistDataEff= hDataEff_LowPt;
        tmphistDataErrl= hDataErrl_LowPt;
        tmphistDataErrh= hDataErrh_LowPt;
      } else {
        tmphistMCEff = hMCEff;
        tmphistMCErrl = hMCErrl;
        tmphistMCErrh = hMCErrh;
        tmphistDataEff= hDataEff;
        tmphistDataErrl= hDataErrl;
        tmphistDataErrh= hDataErrh;
      }

      twikifile << binLabel ;

      ios_base::fmtflags flags = twikifile.flags();
      twikifile.precision(4);
      
      Double_t mceff  = tmphistMCEff->GetCellContent(tmphistMCEff->GetXaxis()->FindFixBin(binXValue),tmphistMCEff->GetYaxis()->FindFixBin(binYValue));
      Double_t mcerrl = tmphistMCErrl->GetCellContent(tmphistMCErrl->GetXaxis()->FindFixBin(binXValue),tmphistMCErrl->GetYaxis()->FindFixBin(binYValue));
      Double_t mcerrh = tmphistMCErrh->GetCellContent(tmphistMCErrh->GetXaxis()->FindFixBin(binXValue),tmphistMCErrh->GetYaxis()->FindFixBin(binYValue));
      twikifile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);
      twikifile << " | ";

      Double_t dataeff  = tmphistDataEff->GetCellContent(tmphistDataEff->GetXaxis()->FindFixBin(binXValue),tmphistDataEff->GetYaxis()->FindFixBin(binYValue));
      Double_t dataerrl = tmphistDataErrl->GetCellContent(tmphistDataErrl->GetXaxis()->FindFixBin(binXValue),tmphistDataErrl->GetYaxis()->FindFixBin(binYValue));
      Double_t dataerrh = tmphistDataErrh->GetCellContent(tmphistDataErrh->GetXaxis()->FindFixBin(binXValue),tmphistDataErrh->GetYaxis()->FindFixBin(binYValue));
      twikifile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      twikifile << " | ";
     
      Double_t scale     = (tmphistDataEff->GetCellContent(tmphistDataEff->GetXaxis()->FindFixBin(binXValue),tmphistDataEff->GetYaxis()->FindFixBin(binYValue)))/(tmphistMCEff->GetCellContent(tmphistMCEff->GetXaxis()->FindFixBin(binXValue),tmphistMCEff->GetYaxis()->FindFixBin(binYValue)));
      Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
      twikifile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);
      twikifile << " |\n";

    }
  }
  twikifile.close();
  cout << outputDir + "/eff_table.twiki" << " created!" << endl;


  cout << "ptbins : ";
  for(Int_t i=0; i<ny+1;++i) {
    cout << ptbins[i] << " ";
  }
  cout << endl;
  cout << "etabins : ";
  for(Int_t i=0; i<nx+1;++i) {
    cout << etabins[i] << " ";
  }
  cout << endl;

  //Turn around histogram axes
  TH2F *h2_results_selection_rotated = new TH2F( "h2_results_selection", h2_results_selection->GetTitle(), ny, ptbins, nx,  etabins);
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {
      cout << ix << " " << iy << " : " << h2_results_selection->GetCellContent(iy,ix) << endl;
      h2_results_selection_rotated->SetCellContent(iy,ix,h2_results_selection->GetCellContent(iy,ix));
      h2_results_selection_rotated->SetCellError(iy,ix,h2_results_selection->GetCellError(iy,ix));
    }
  }

//   outputFile->WriteTObject(h2_results_selection, h2_results_selection->GetName(), "WriteDelete");
  outputFile->WriteTObject(h2_results_selection_rotated, h2_results_selection->GetName(), "WriteDelete");

  TFile *DataEfficiencyOutputFile = new TFile((outputDir + "/effData.root").c_str(), "UPDATE");
  TFile *MCEfficiencyOutputFile = new TFile((outputDir + "/effMC.root").c_str(), "UPDATE");

  DataEfficiencyOutputFile->WriteTObject(hDataEffCombined, hDataEffCombined->GetName(), "WriteDelete");
  DataEfficiencyOutputFile->WriteTObject(hDataErrlCombined, hDataErrlCombined->GetName(), "WriteDelete");
  DataEfficiencyOutputFile->WriteTObject(hDataErrhCombined, hDataErrhCombined->GetName(), "WriteDelete");
  MCEfficiencyOutputFile->WriteTObject(hMCEffCombined, hMCEffCombined->GetName(), "WriteDelete");
  MCEfficiencyOutputFile->WriteTObject(hMCErrlCombined, hMCErrlCombined->GetName(), "WriteDelete");
  MCEfficiencyOutputFile->WriteTObject(hMCErrhCombined, hMCErrhCombined->GetName(), "WriteDelete");

  outputFile->Close();
  DataEfficiencyOutputFile->Close();
  MCEfficiencyOutputFile->Close();



  gBenchmark->Show("printMuonWPEff"); 
}




void doEfficiencyScaleFactors_LowestPtBinUsesZeeFactorized( string StandardEffBinsFilename,
                                                            string DataFilename, 
                                                            string DataLowPtBinsFilename, 
                                                            string MCFilename, 
                                                            string MCLowPtBinsFilename, 
                                                            Int_t dataset,
                                                            string outputDir, 
                                                            string histName,
                                                            string Label)
{
  gBenchmark->Start("printMuonWPEff");
    
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TString outfname = (outputDir + "/eff_table.txt").c_str();
  TFile standardeffbinsfile(StandardEffBinsFilename.c_str());
  TFile mcfile(MCFilename.c_str());
  TFile datafile(DataFilename.c_str());
  TFile lowpt_mcfile(MCLowPtBinsFilename.c_str());
  TFile lowpt_datafile(DataLowPtBinsFilename.c_str());
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   
  
  //Load Standard Efficiency Bins histogram
  TH2F *hStandardEff=0;
  hStandardEff  = (TH2F*)standardeffbinsfile.Get("hEffEtaPt");

  //Load Standard Tag and Probe efficiencies
  TH2F *hMCEff=0,    *hMCErrl=0,    *hMCErrh=0;
  TH2F *hDataEff=0, *hDataErrl=0, *hDataErrh=0;
  hMCEff  = (TH2F*)mcfile.Get("hEffEtaPt");
  hMCErrl = (TH2F*)mcfile.Get("hErrlEtaPt");
  hMCErrh = (TH2F*)mcfile.Get("hErrhEtaPt");
  hDataEff  = (TH2F*)datafile.Get("hEffEtaPt");
  hDataErrl = (TH2F*)datafile.Get("hErrlEtaPt");
  hDataErrh = (TH2F*)datafile.Get("hErrhEtaPt");

  //Load Low Pt Tag and Probe efficiencies
  TH2F *hMCEff_LowPt=0,    *hMCErrl_LowPt=0,    *hMCErrh_LowPt=0;
  TH2F *hDataEff_LowPt=0, *hDataErrl_LowPt=0, *hDataErrh_LowPt=0;
  hMCEff_LowPt  = (TH2F*)lowpt_mcfile.Get("hEffEtaPt");
  hMCErrl_LowPt = (TH2F*)lowpt_mcfile.Get("hErrlEtaPt");
  hMCErrh_LowPt = (TH2F*)lowpt_mcfile.Get("hErrhEtaPt");
  hDataEff_LowPt  = (TH2F*)lowpt_datafile.Get("hEffEtaPt");
  hDataErrl_LowPt = (TH2F*)lowpt_datafile.Get("hErrlEtaPt");
  hDataErrh_LowPt = (TH2F*)lowpt_datafile.Get("hErrhEtaPt");


  TGraphAsymmErrors* MCEffVsRho   = ( TGraphAsymmErrors*)mcfile.Get("grEffRho");
  TGraphAsymmErrors* DataEffVsRho = ( TGraphAsymmErrors*)datafile.Get("grEffRho");
 
  TGraphAsymmErrors* MCEffVsNVtx   = ( TGraphAsymmErrors*)mcfile.Get("grEffNPV");
  TGraphAsymmErrors* DataEffVsNVtx = ( TGraphAsymmErrors*)datafile.Get("grEffNPV");


  //--------------------------------------------------------------------------------------------------------------
  // Update root file histograms
  //==============================================================================================================   
  const Int_t nx = hStandardEff->GetNbinsX();
  const Int_t ny = hStandardEff->GetNbinsY();

  TFile *outputFile = new TFile(("efficiency_results"+label+".root").c_str(), "UPDATE");

  //Do Binning
  Double_t *ptbins = new Double_t[ny+1];
  Double_t *etabins = new Double_t[nx+1];
  for(Int_t iy=1; iy<=ny; iy++) {
      ptbins[iy-1] = hStandardEff->GetYaxis()->GetBinLowEdge(iy);
  }
  for(Int_t ix=1; ix<=nx; ix++) {
    etabins[ix-1]= hStandardEff->GetXaxis()->GetBinLowEdge(ix);
  }
  ptbins[ny] = 60;
  etabins[nx] = 2.5;
    
  TH2F *h2_results_selection = new TH2F(histName.c_str(),"",ny,ptbins,nx,etabins);
  for(Int_t iy=0; iy<=ny+2; iy++) {
    for(Int_t ix=0; ix<=nx+2; ix++) {
      h2_results_selection->SetCellContent(iy,ix, 1.0);
      h2_results_selection->SetCellError(iy,ix, 0.0);
    }
  }

  //save output efficiencies into root file
  TH2F *hMCEffCombined=new TH2F("hEffEtaPt","",ny,ptbins,nx,etabins); 
  TH2F *hMCErrlCombined=new TH2F("hErrlEtaPt","",ny,ptbins,nx,etabins);
  TH2F *hMCErrhCombined=new TH2F("hErrhEtaPt","",ny,ptbins,nx,etabins);
  TH2F *hDataEffCombined=new TH2F("hEffEtaPt","",ny,ptbins,nx,etabins);
  TH2F *hDataErrlCombined=new TH2F("hErrlEtaPt","",ny,ptbins,nx,etabins);
  TH2F *hDataErrhCombined=new TH2F("hErrhEtaPt","",ny,ptbins,nx,etabins);

  cout << "ptbins : ";
  for(Int_t i=0; i<ny+1;++i) {
    cout << ptbins[i] << " ";
  }
  cout << endl;
  cout << "etabins : ";
  for(Int_t i=0; i<nx+1;++i) {
    cout << etabins[i] << " ";
  }
  cout << endl;

  //--------------------------------------------------------------------------------------------------------------
  // Produce Text file table
  //==============================================================================================================   

  ofstream txtfile;
  txtfile.open(outfname.Data());
  assert(txtfile.is_open());
    
  
  txtfile << " pT        ";
  txtfile << " eta           ";
  txtfile << "    MC T&P           ";
  txtfile << "    Data T&P         ";  
  txtfile << "    Scale factor     ";
  txtfile << endl;
  txtfile << "----------------------------------------------------------------------------------------------------------------------------------------------------";
  txtfile << endl;
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {
      txtfile << "[" << setw(4) << hStandardEff->GetYaxis()->GetBinLowEdge(iy) << "," << setw(4) << hStandardEff->GetYaxis()->GetBinLowEdge(iy+1) << "]";
      txtfile << "[" << setw(6) << hStandardEff->GetXaxis()->GetBinLowEdge(ix) << "," << setw(6) << hStandardEff->GetXaxis()->GetBinLowEdge(ix+1) << "]";      
      ios_base::fmtflags flags = txtfile.flags();
      txtfile.precision(4);
      
      Double_t binXValue = hStandardEff->GetXaxis()->GetBinCenter(ix);
      Double_t binYValue = hStandardEff->GetYaxis()->GetBinCenter(iy);
      //cout << "Bin : " << binXValue << " " << binYValue << endl;

      TH2F *tmphistMCEff=0,    *tmphistMCErrl=0,    *tmphistMCErrh=0;
      TH2F *tmphistDataEff=0, *tmphistDataErrl=0, *tmphistDataErrh=0;
      
      if(hDataEff_LowPt && hMCEff_LowPt && hStandardEff->GetYaxis()->GetBinCenter(iy) < 20)  {
        if (!(hDataEff_LowPt && hMCEff_LowPt)) { cout << "Zee low pt eff hist not found\n"; assert(false);}
        tmphistMCEff = hMCEff_LowPt;
        tmphistMCErrl = hMCErrl_LowPt;
        tmphistMCErrh = hMCErrh_LowPt;
        tmphistDataEff= hDataEff_LowPt;
        tmphistDataErrl= hDataErrl_LowPt;
        tmphistDataErrh= hDataErrh_LowPt;
      } else {
        tmphistMCEff = hMCEff;
        tmphistMCErrl = hMCErrl;
        tmphistMCErrh = hMCErrh;
        tmphistDataEff= hDataEff;
        tmphistDataErrl= hDataErrl;
        tmphistDataErrh= hDataErrh;
      }

      Double_t mceff  = tmphistMCEff->GetCellContent(tmphistMCEff->GetXaxis()->FindFixBin(binXValue),tmphistMCEff->GetYaxis()->FindFixBin(binYValue));
      Double_t mcerrl = tmphistMCErrl->GetCellContent(tmphistMCErrl->GetXaxis()->FindFixBin(binXValue),tmphistMCErrl->GetYaxis()->FindFixBin(binYValue));
      Double_t mcerrh = tmphistMCErrh->GetCellContent(tmphistMCErrh->GetXaxis()->FindFixBin(binXValue),tmphistMCErrh->GetYaxis()->FindFixBin(binYValue));

      Double_t dataeff  = tmphistDataEff->GetCellContent(tmphistDataEff->GetXaxis()->FindFixBin(binXValue),tmphistDataEff->GetYaxis()->FindFixBin(binYValue));
      Double_t dataerrl = tmphistDataErrl->GetCellContent(tmphistDataErrl->GetXaxis()->FindFixBin(binXValue),tmphistDataErrl->GetYaxis()->FindFixBin(binYValue));
      Double_t dataerrh = tmphistDataErrh->GetCellContent(tmphistDataErrh->GetXaxis()->FindFixBin(binXValue),tmphistDataErrh->GetYaxis()->FindFixBin(binYValue));
      
      Double_t scale     = (tmphistDataEff->GetCellContent(tmphistDataEff->GetXaxis()->FindFixBin(binXValue),tmphistDataEff->GetYaxis()->FindFixBin(binYValue)))/(tmphistMCEff->GetCellContent(tmphistMCEff->GetXaxis()->FindFixBin(binXValue),tmphistMCEff->GetYaxis()->FindFixBin(binYValue)));
      Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);


      //*************************************************************************************
      //For 7-10 bin we take the factorized IDGivenIso x IsoGivenID scale factors
      //*************************************************************************************
      if (iy==1 && (ix == 1 || ix == 2)) {
        if (dataset == 2011) {
          //For Run1 Legacy Paper 2011 Data
          dataeff = 0.8917*0.8479; 
          dataerrl = 0.8917*0.8479 * sqrt( pow(0.0276/0.8917,2) + pow(0.0224/0.8479,2) ); 
          dataerrh = 0.8917*0.8479 * sqrt( pow(0.0276/0.8917,2) + pow(0.0224/0.8479,2) );
          mceff = 0.9054*0.8705; 
          mcerrl = 0.9054*0.8705 * sqrt( pow(0.0035/0.9054,2) + pow(0.0040/0.8705,2) ); 
          mcerrh = 0.9054*0.8705 * sqrt( pow(0.0035/0.9054,2) + pow(0.0040/0.8705,2) );
          scale = 0.9593;
          scaleErrl = 0.0395;       
          scaleErrh = 0.0395;               
        } else if (dataset == 2012) {
          //For Run1 Legacy Paper 2012 Data
          dataeff = 0.8452*0.7846; 
          dataerrl = 0.8452*0.7846 * sqrt( pow(0.0342/0.8452,2) + pow(0.0029/0.7846,2) ); 
          dataerrh = 0.8452*0.7846 * sqrt( pow(0.0342/0.8452,2) + pow(0.0029/0.7846,2) );
          mceff = 0.8832*0.8712; 
          mcerrl = 0.8832*0.8712 * sqrt( pow(0.0063/0.8832,2) + pow(0.0063/0.8712,2) ); 
          mcerrh = 0.8832*0.8712 * sqrt( pow(0.0063/0.8832,2) + pow(0.0063/0.8712,2) );
          scale = 0.8618;
          scaleErrl = 0.0361;       
          scaleErrh = 0.0361;               
        } else {
          cout << "Error: Dataset " << dataset << " is not supported.\n";
        }
      }
      if (iy==1 && (ix == 3 || ix == 4 || ix == 5)) {
        if (dataset == 2011) {
          //For Run1 Legacy Paper 2011 Data
          dataeff = 0.7024*0.8519; 
          dataerrl = 0.7024*0.8519 * sqrt( pow(0.0280/0.7024,2) + pow(0.0272/0.8519,2) ); 
          dataerrh = 0.7024*0.8519 * sqrt( pow(0.0280/0.7024,2) + pow(0.0272/0.8519,2) );
          mceff = 0.7776*0.8834; 
          mcerrl = 0.7776*0.8834 * sqrt( pow(0.0044/0.7776,2) + pow(0.0036/0.8834,2) ); 
          mcerrh = 0.7776*0.8834 * sqrt( pow(0.0044/0.7776,2) + pow(0.0036/0.8834,2) );
          scale = 0.8710;
          scaleErrl = 0.0449;  
          scaleErrh = 0.0449;  
        } else if (dataset == 2012) {
          //For Run1 Legacy Paper 2012 Data
          //I artificially made the uncertainty on the data fit for IDGivenIso 3%, because
          //the fit error was abnormally small
          dataeff = 0.5307*0.8809; 
          dataerrl = 0.5307*0.8809 * sqrt( pow(0.03/0.5307,2) + pow(0.0161/0.8809,2) ); 
          dataerrh = 0.5307*0.8809 * sqrt( pow(0.03/0.5307,2) + pow(0.0161/0.8809,2) );
          mceff = 0.6695*0.8913; 
          mcerrl = 0.6695*0.8913 * sqrt( pow(0.0086/0.6695,2) + pow(0.0068/0.8913,2) ); 
          mcerrh = 0.6695*0.8913 * sqrt( pow(0.0086/0.6695,2) + pow(0.0068/0.8913,2) );
          scale = 0.7835;
          scaleErrl = 0.0480;  
          scaleErrh = 0.0480;  
        } else {
          cout << "Error: Dataset " << dataset << " is not supported.\n";
        }
      }
      //*************************************************************************************


      txtfile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);
      txtfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      txtfile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);

      txtfile << endl;
      txtfile.flags(flags);

      cout << "Set " << iy << " " << ix+1 << " " << dataeff << " / " << mceff << " = " << scale << endl;
      h2_results_selection->SetCellContent(iy,ix, scale);
      h2_results_selection->SetCellError(iy,ix, (scaleErrl + scaleErrh)/2);
      hDataEffCombined->SetCellContent(iy,ix, dataeff);
      hDataErrlCombined->SetCellContent(iy,ix, dataerrl);
      hDataErrhCombined->SetCellContent(iy,ix, dataerrh);
      hMCEffCombined->SetCellContent(iy,ix, mceff);
      hMCErrlCombined->SetCellContent(iy,ix, mcerrl);
      hMCErrhCombined->SetCellContent(iy,ix, mcerrh);

      //fill overflow bins with the same values as last bin
      if (ix == nx) {
        cout << "Set " << iy << " " << ix+1 << " " << dataeff << " / " << mceff << " = " << scale << endl;
        h2_results_selection->SetCellContent(iy,nx+1, scale);
        h2_results_selection->SetCellError(iy,nx+1, (scaleErrl + scaleErrh)/2);
        hDataEffCombined->SetCellContent(iy,nx+1, dataeff);
        hDataErrlCombined->SetCellContent(iy,nx+1, dataerrl);
        hDataErrhCombined->SetCellContent(iy,nx+1, dataerrh);
        hMCEffCombined->SetCellContent(iy,nx+1, mceff);
        hMCErrlCombined->SetCellContent(iy,nx+1, mcerrl);
        hMCErrhCombined->SetCellContent(iy,nx+1, mcerrh);
      }
    }
    if (iy == ny) {
      for(Int_t ix=1; ix<=nx; ix++) {
        Double_t binXValue = hStandardEff->GetXaxis()->GetBinCenter(ix);
        Double_t binYValue = hStandardEff->GetYaxis()->GetBinCenter(iy);

        TH2F *tmphistMCEff=0,    *tmphistMCErrl=0,    *tmphistMCErrh=0;
        TH2F *tmphistDataEff=0, *tmphistDataErrl=0, *tmphistDataErrh=0;
      
        if(hDataEff_LowPt && hMCEff_LowPt && hStandardEff->GetYaxis()->GetBinCenter(iy) < 20)  {
          tmphistMCEff = hMCEff_LowPt;
          tmphistMCErrl = hMCErrl_LowPt;
          tmphistMCErrh = hMCErrh_LowPt;
          tmphistDataEff= hDataEff_LowPt;
          tmphistDataErrl= hDataErrl_LowPt;
          tmphistDataErrh= hDataErrh_LowPt;
        } else {
          tmphistMCEff = hMCEff;
          tmphistMCErrl = hMCErrl;
          tmphistMCErrh = hMCErrh;
          tmphistDataEff= hDataEff;
          tmphistDataErrl= hDataErrl;
          tmphistDataErrh= hDataErrh;
        }

        Double_t mceff  = tmphistMCEff->GetCellContent(tmphistMCEff->GetXaxis()->FindFixBin(binXValue),tmphistMCEff->GetYaxis()->FindFixBin(binYValue));
        Double_t mcerrl = tmphistMCErrl->GetCellContent(tmphistMCErrl->GetXaxis()->FindFixBin(binXValue),tmphistMCErrl->GetYaxis()->FindFixBin(binYValue));
        Double_t mcerrh = tmphistMCErrh->GetCellContent(tmphistMCErrh->GetXaxis()->FindFixBin(binXValue),tmphistMCErrh->GetYaxis()->FindFixBin(binYValue));
        Double_t dataeff  = tmphistDataEff->GetCellContent(tmphistDataEff->GetXaxis()->FindFixBin(binXValue),tmphistDataEff->GetYaxis()->FindFixBin(binYValue));
        Double_t dataerrl = tmphistDataErrl->GetCellContent(tmphistDataErrl->GetXaxis()->FindFixBin(binXValue),tmphistDataErrl->GetYaxis()->FindFixBin(binYValue));
        Double_t dataerrh = tmphistDataErrh->GetCellContent(tmphistDataErrh->GetXaxis()->FindFixBin(binXValue),tmphistDataErrh->GetYaxis()->FindFixBin(binYValue));
        Double_t scale     = (tmphistDataEff->GetCellContent(tmphistDataEff->GetXaxis()->FindFixBin(binXValue),tmphistDataEff->GetYaxis()->FindFixBin(binYValue)))/(tmphistMCEff->GetCellContent(tmphistMCEff->GetXaxis()->FindFixBin(binXValue),tmphistMCEff->GetYaxis()->FindFixBin(binYValue)));
        Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
        Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
        
        cout << "Set " << iy+1 << " " << ix << " " << scale << endl;
        h2_results_selection->SetCellContent(iy+1,ix, scale);
        h2_results_selection->SetCellError(iy+1,ix, (scaleErrl + scaleErrh)/2);
        hDataEffCombined->SetCellContent(iy+1,ix, dataeff);
        hDataErrlCombined->SetCellContent(iy+1,ix, dataerrl);
        hDataErrhCombined->SetCellContent(iy+1,ix, dataerrh);
        hMCEffCombined->SetCellContent(iy+1,ix, mceff);
        hMCErrlCombined->SetCellContent(iy+1,ix, mcerrl);
        hMCErrhCombined->SetCellContent(iy+1,ix, mcerrh);
        
        //fill overflow bins with the same values as last bin
        if (ix == nx) {
          cout << "Set " << iy+1 << " " << ix+1 << " " << scale << endl;
          h2_results_selection->SetCellContent(iy+1,nx+1, scale);
          h2_results_selection->SetCellError(iy+1,nx+1, (scaleErrl + scaleErrh)/2);
          hDataEffCombined->SetCellContent(iy+1,nx+1, dataeff);
          hDataErrlCombined->SetCellContent(iy+1,nx+1, dataerrl);
          hDataErrhCombined->SetCellContent(iy+1,nx+1, dataerrh);
          hMCEffCombined->SetCellContent(iy+1,nx+1, mceff);
          hMCErrlCombined->SetCellContent(iy+1,nx+1, mcerrl);
          hMCErrhCombined->SetCellContent(iy+1,nx+1, mcerrh);
        }
      }
    }

    txtfile << endl;
  }
  txtfile.close();
  
  cout << outfname << " created!" << endl;
  



  //--------------------------------------------------------------------------------------------------------------
  // Create TEX table
  //==============================================================================================================   

  ofstream texfile;
  texfile.open((outputDir + "/eff_table.tex").c_str());
  assert(texfile.is_open());

  texfile << " \\begin{table}[!ht]" << endl;
  texfile << " \\begin{center} " << endl;
  texfile << " \\begin{tabular}{|c|c|c|c|}" << endl;
  texfile << " \\hline\n";


  texfile << " $p_{T}$ / $\\eta$ bin    &  Monte Carlo Efficiency    &  Data Efficiency   &  MC to Data Scale Factor \\\\  ";
  texfile << " \\hline           ";
  texfile << endl;
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {


      Double_t binXValue = hStandardEff->GetYaxis()->GetBinCenter(iy);
      Double_t binYValue = hStandardEff->GetXaxis()->GetBinCenter(ix);

      string binLabel = Form("$%5.1f < p_{T} \\le %5.1f$ , $%5.1f  \\le |\\eta| < %5.1f$", 
                             hStandardEff->GetYaxis()->GetBinLowEdge(iy), hStandardEff->GetYaxis()->GetBinLowEdge(iy+1),
                             hStandardEff->GetXaxis()->GetBinLowEdge(ix), hStandardEff->GetXaxis()->GetBinLowEdge(ix+1));
      if (iy == ny) {
        binLabel = Form("$%5.1f < p_{T} $ , $%5.1f  \\le |\\eta| < %5.1f$", 
                        hStandardEff->GetYaxis()->GetBinLowEdge(iy), 
                        hStandardEff->GetXaxis()->GetBinLowEdge(ix), hStandardEff->GetXaxis()->GetBinLowEdge(ix+1));
      }
      
      texfile << binLabel;
      texfile << "   &   ";

      ios_base::fmtflags flags = texfile.flags();
      texfile.precision(4);
      
      Double_t mceff  = hMCEffCombined->GetCellContent(hMCEffCombined->GetXaxis()->FindFixBin(binXValue),hMCEffCombined->GetYaxis()->FindFixBin(binYValue));
      Double_t mcerrl = hMCErrlCombined->GetCellContent(hMCErrlCombined->GetXaxis()->FindFixBin(binXValue),hMCErrlCombined->GetYaxis()->FindFixBin(binYValue));
      Double_t mcerrh = hMCErrhCombined->GetCellContent(hMCErrhCombined->GetXaxis()->FindFixBin(binXValue),hMCErrhCombined->GetYaxis()->FindFixBin(binYValue));
      texfile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);
      texfile << "   &   ";

      Double_t dataeff  = hDataEffCombined->GetCellContent(hDataEffCombined->GetXaxis()->FindFixBin(binXValue),hDataEffCombined->GetYaxis()->FindFixBin(binYValue));
      Double_t dataerrl = hDataErrlCombined->GetCellContent(hDataErrlCombined->GetXaxis()->FindFixBin(binXValue),hDataErrlCombined->GetYaxis()->FindFixBin(binYValue));
      Double_t dataerrh = hDataErrhCombined->GetCellContent(hDataErrhCombined->GetXaxis()->FindFixBin(binXValue),hDataErrhCombined->GetYaxis()->FindFixBin(binYValue));
      texfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      texfile << "   &   ";
     
      Double_t scale     = h2_results_selection->GetCellContent(h2_results_selection->GetXaxis()->FindFixBin(binXValue),h2_results_selection->GetYaxis()->FindFixBin(binYValue));
      Double_t scaleErrl = h2_results_selection->GetCellError(h2_results_selection->GetXaxis()->FindFixBin(binXValue),h2_results_selection->GetYaxis()->FindFixBin(binYValue));
      Double_t scaleErrh = h2_results_selection->GetCellError(h2_results_selection->GetXaxis()->FindFixBin(binXValue),h2_results_selection->GetYaxis()->FindFixBin(binYValue));      
      texfile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);
      texfile << "   \\\\   ";
      texfile << endl;
      texfile << "\\hline";
      texfile << endl;

    }
  }
  texfile << "\\end{tabular}" << endl;
  texfile << "\\caption{CAPTION.}" << endl;
  texfile << "\\label{tab:eff_ele_offline}" << endl;
  texfile << "\\end{center}" << endl;
  texfile << "\\end{table}" << endl;
  texfile.close();
  cout << outputDir + "/eff_table.tex" << " created!" << endl;




  //--------------------------------------------------------------------------------------------------------------
  // Create TWIKI table
  //==============================================================================================================   

  ofstream twikifile;
  twikifile.open((outputDir + "/eff_table.twiki").c_str());
  assert(twikifile.is_open());
  

  twikifile << "| *pT bin* | *eta bin* | *MC Efficiency* | *Data Efficiency* | *ScaleFactor* | \n";
  
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {

      Double_t binXValue = hStandardEff->GetYaxis()->GetBinCenter(iy);
      Double_t binYValue = hStandardEff->GetXaxis()->GetBinCenter(ix);

      string binLabel = Form("| *%4.1f < pT \\le %4.1f* | *%3.1f  <= eta < %3.1f* | ", 
                             hStandardEff->GetYaxis()->GetBinLowEdge(iy), hStandardEff->GetYaxis()->GetBinLowEdge(iy+1),
                             hStandardEff->GetXaxis()->GetBinLowEdge(ix), hStandardEff->GetXaxis()->GetBinLowEdge(ix+1));
      if (iy == ny) {
        binLabel = Form("| *%4.1f < pT* | *%3.1f  <= eta < %3.1f* | ", 
                        hStandardEff->GetYaxis()->GetBinLowEdge(iy), 
                        hStandardEff->GetXaxis()->GetBinLowEdge(ix), hStandardEff->GetXaxis()->GetBinLowEdge(ix+1));
      }
      
      TH2F *tmphistMCEff=0,    *tmphistMCErrl=0,    *tmphistMCErrh=0;
      TH2F *tmphistDataEff=0, *tmphistDataErrl=0, *tmphistDataErrh=0;
      
      if(hDataEff_LowPt && hMCEff_LowPt && hStandardEff->GetYaxis()->GetBinCenter(iy) < 20)  {
        tmphistMCEff = hMCEff_LowPt;
        tmphistMCErrl = hMCErrl_LowPt;
        tmphistMCErrh = hMCErrh_LowPt;
        tmphistDataEff= hDataEff_LowPt;
        tmphistDataErrl= hDataErrl_LowPt;
        tmphistDataErrh= hDataErrh_LowPt;
      } else {
        tmphistMCEff = hMCEff;
        tmphistMCErrl = hMCErrl;
        tmphistMCErrh = hMCErrh;
        tmphistDataEff= hDataEff;
        tmphistDataErrl= hDataErrl;
        tmphistDataErrh= hDataErrh;
      }

      twikifile << binLabel ;

      ios_base::fmtflags flags = twikifile.flags();
      twikifile.precision(4);
      
      Double_t mceff  = hMCEffCombined->GetCellContent(hMCEffCombined->GetXaxis()->FindFixBin(binXValue),hMCEffCombined->GetYaxis()->FindFixBin(binYValue));
      Double_t mcerrl = hMCErrlCombined->GetCellContent(hMCErrlCombined->GetXaxis()->FindFixBin(binXValue),hMCErrlCombined->GetYaxis()->FindFixBin(binYValue));
      Double_t mcerrh = hMCErrhCombined->GetCellContent(hMCErrhCombined->GetXaxis()->FindFixBin(binXValue),hMCErrhCombined->GetYaxis()->FindFixBin(binYValue));
      twikifile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);
      twikifile << " | ";

      Double_t dataeff  = hDataEffCombined->GetCellContent(hDataEffCombined->GetXaxis()->FindFixBin(binXValue),hDataEffCombined->GetYaxis()->FindFixBin(binYValue));
      Double_t dataerrl = hDataErrlCombined->GetCellContent(hDataErrlCombined->GetXaxis()->FindFixBin(binXValue),hDataErrlCombined->GetYaxis()->FindFixBin(binYValue));
      Double_t dataerrh = hDataErrhCombined->GetCellContent(hDataErrhCombined->GetXaxis()->FindFixBin(binXValue),hDataErrhCombined->GetYaxis()->FindFixBin(binYValue));
      twikifile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      twikifile << " | ";
     
      Double_t scale     = h2_results_selection->GetCellContent(h2_results_selection->GetXaxis()->FindFixBin(binXValue),h2_results_selection->GetYaxis()->FindFixBin(binYValue));
      Double_t scaleErrl = h2_results_selection->GetCellError(h2_results_selection->GetXaxis()->FindFixBin(binXValue),h2_results_selection->GetYaxis()->FindFixBin(binYValue));
      Double_t scaleErrh = h2_results_selection->GetCellError(h2_results_selection->GetXaxis()->FindFixBin(binXValue),h2_results_selection->GetYaxis()->FindFixBin(binYValue));      
      twikifile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);
      twikifile << " |\n";

    }
  }
  twikifile.close();
  cout << outputDir + "/eff_table.twiki" << " created!" << endl;


  cout << "ptbins : ";
  for(Int_t i=0; i<ny+1;++i) {
    cout << ptbins[i] << " ";
  }
  cout << endl;
  cout << "etabins : ";
  for(Int_t i=0; i<nx+1;++i) {
    cout << etabins[i] << " ";
  }
  cout << endl;

  //Turn around histogram axes
  TH2F *h2_results_selection_rotated = new TH2F( "h2_results_selection", h2_results_selection->GetTitle(), ny, ptbins, nx,  etabins);
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {
      cout << ix << " " << iy << " : " << h2_results_selection->GetCellContent(iy,ix) << endl;
      h2_results_selection_rotated->SetCellContent(iy,ix,h2_results_selection->GetCellContent(iy,ix));
      h2_results_selection_rotated->SetCellError(iy,ix,h2_results_selection->GetCellError(iy,ix));
    }
  }

//   outputFile->WriteTObject(h2_results_selection, h2_results_selection->GetName(), "WriteDelete");
  outputFile->WriteTObject(h2_results_selection_rotated, h2_results_selection->GetName(), "WriteDelete");

  TFile *DataEfficiencyOutputFile = new TFile((outputDir + "/effData.root").c_str(), "UPDATE");
  TFile *MCEfficiencyOutputFile = new TFile((outputDir + "/effMC.root").c_str(), "UPDATE");

  DataEfficiencyOutputFile->WriteTObject(hDataEffCombined, hDataEffCombined->GetName(), "WriteDelete");
  DataEfficiencyOutputFile->WriteTObject(hDataErrlCombined, hDataErrlCombined->GetName(), "WriteDelete");
  DataEfficiencyOutputFile->WriteTObject(hDataErrhCombined, hDataErrhCombined->GetName(), "WriteDelete");
  MCEfficiencyOutputFile->WriteTObject(hMCEffCombined, hMCEffCombined->GetName(), "WriteDelete");
  MCEfficiencyOutputFile->WriteTObject(hMCErrlCombined, hMCErrlCombined->GetName(), "WriteDelete");
  MCEfficiencyOutputFile->WriteTObject(hMCErrhCombined, hMCErrhCombined->GetName(), "WriteDelete");

  outputFile->Close();
  DataEfficiencyOutputFile->Close();
  MCEfficiencyOutputFile->Close();



  gBenchmark->Show("printMuonWPEff"); 
}





void doEfficiencyScaleFactorsOneSample(string DataFilename,
                   string MCFilename, 
                   string outputDir, 
                   string histName,
                   string Label)
{
    
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
//   TString outfname = "Data_EleWPEffTP/basic2_76_106/eff_table.txt";
//   TFile mcfile("Summer11_Zee_EleWPEffTP/basic_76_106/eff.root");
// //  TFile datafile1("Data_PR_MuonWPEffTP/basic1_60_120/eff.root");
// //   TFile datafile2("Data_PR_MuonIPEffTP/test2_76_106/eff.root");
//   TFile datafile2("Data_EleWPEffTP/basic2_76_106/eff.root");
  
//   TString outfname = "Data_EleLHEffTP/basic2_76_106/eff_table.txt";
//   TFile mcfile("Summer11_Zee_EleLHEffTP/basic_76_106/eff.root");
//   TFile datafile2("Data_EleLHEffTP/basic2_76_106/eff.root");
  
  TString outfname = (outputDir + "/eff_table.txt").c_str();
  TFile mcfile(MCFilename.c_str());
  TFile datafile2(DataFilename.c_str());
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   
  
  TH2F *hMCEff=0,    *hMCErrl=0,    *hMCErrh=0;
  TH2F *hDataEff1=0, *hDataErrl1=0, *hDataErrh1=0;
  TH2F *hDataEff2=0, *hDataErrl2=0, *hDataErrh2=0;
  
  hMCEff  = (TH2F*)mcfile.Get("hEffEtaPt");
  hMCErrl = (TH2F*)mcfile.Get("hErrlEtaPt");
  hMCErrh = (TH2F*)mcfile.Get("hErrhEtaPt");
  
//  hDataEff1  = (TH2F*)datafile1.Get("hEffEtaPt");
//  hDataErrl1 = (TH2F*)datafile1.Get("hErrlEtaPt");
//  hDataErrh1 = (TH2F*)datafile1.Get("hErrhEtaPt");
  
  hDataEff2  = (TH2F*)datafile2.Get("hEffEtaPt");
  hDataErrl2 = (TH2F*)datafile2.Get("hErrlEtaPt");
  hDataErrh2 = (TH2F*)datafile2.Get("hErrhEtaPt");
  
  TGraphAsymmErrors* MCEffVsRho   = ( TGraphAsymmErrors*)mcfile.Get("grEffRho");
  TGraphAsymmErrors* DataEffVsRho = ( TGraphAsymmErrors*)datafile2.Get("grEffRho");
 
  TGraphAsymmErrors* MCEffVsNVtx   = ( TGraphAsymmErrors*)mcfile.Get("grEffNPV");
  TGraphAsymmErrors* DataEffVsNVtx = ( TGraphAsymmErrors*)datafile2.Get("grEffNPV");


  //--------------------------------------------------------------------------------------------------------------
  // Update root file histograms
  //==============================================================================================================   
  const Int_t nx = hMCEff->GetNbinsX();
  const Int_t ny = hMCEff->GetNbinsY();

  TFile *outputFile = new TFile(("efficiency_results"+label+".root").c_str(), "UPDATE");

  //Do Binning
  Double_t *ptbins = new Double_t[ny+1];
  Double_t *etabins = new Double_t[nx+1];
  for(Int_t iy=1; iy<=ny; iy++) {
      ptbins[iy-1] = hMCEff->GetYaxis()->GetBinLowEdge(iy);
  }
  for(Int_t ix=1; ix<=nx; ix++) {
    etabins[ix-1]= hMCEff->GetXaxis()->GetBinLowEdge(ix);
  }
  ptbins[ny] = 60;
  etabins[nx] = 2.5;
    
  TH2F *h2_results_selection = new TH2F(histName.c_str(),"",ny,ptbins,nx,etabins);
  for(Int_t iy=0; iy<=ny+2; iy++) {
    for(Int_t ix=0; ix<=nx+2; ix++) {
      h2_results_selection->SetCellContent(iy,ix, 1.0);
      h2_results_selection->SetCellError(iy,ix, 0.0);
    }
  }

  cout << "ptbins : ";
  for(Int_t i=0; i<ny+1;++i) {
    cout << ptbins[i] << " ";
  }
  cout << endl;
  cout << "etabins : ";
  for(Int_t i=0; i<nx+1;++i) {
    cout << etabins[i] << " ";
  }
  cout << endl;

  //--------------------------------------------------------------------------------------------------------------
  // Produce Text file table
  //==============================================================================================================   

  ofstream txtfile;
  txtfile.open(outfname.Data());
  assert(txtfile.is_open());
    
  
  txtfile << " pT        ";
  txtfile << " eta           ";
  txtfile << "    MC T&P           ";
  txtfile << "    Data T&P         ";  
  txtfile << "    Scale factor     ";
  txtfile << endl;
  txtfile << "----------------------------------------------------------------------------------------------------------------------------------------------------";
  txtfile << endl;
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {
      txtfile << "[" << setw(4) << hMCEff->GetYaxis()->GetBinLowEdge(iy) << "," << setw(4) << hMCEff->GetYaxis()->GetBinLowEdge(iy+1) << "]";
      txtfile << "[" << setw(6) << hMCEff->GetXaxis()->GetBinLowEdge(ix) << "," << setw(6) << hMCEff->GetXaxis()->GetBinLowEdge(ix+1) << "]";      
      ios_base::fmtflags flags = txtfile.flags();
      txtfile.precision(4);
      
      Double_t mceff  = hMCEff->GetCellContent(ix,iy);
      Double_t mcerrl = hMCErrl->GetCellContent(ix,iy);
      Double_t mcerrh = hMCErrh->GetCellContent(ix,iy);
      txtfile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);
/*      
      Double_t dataeff  = hDataEff1->GetCellContent(ix,iy);
      Double_t dataerrl = hDataErrl1->GetCellContent(ix,iy);
      Double_t dataerrh = hDataErrh1->GetCellContent(ix,iy);
      txtfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      
      Double_t scale     = (hDataEff1->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
      Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
      txtfile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);
      
      dataeff  = hDataEff2->GetCellContent(ix,iy);
      dataerrl = hDataErrl2->GetCellContent(ix,iy);
      dataerrh = hDataErrh2->GetCellContent(ix,iy);
      txtfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
            
      scale	= (hDataEff2->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
      scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
      txtfile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);
*/
      Double_t dataeff  = hDataEff2->GetCellContent(ix,iy);
      Double_t dataerrl = hDataErrl2->GetCellContent(ix,iy);
      Double_t dataerrh = hDataErrh2->GetCellContent(ix,iy);
      txtfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      
      Double_t scale     = (hDataEff2->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
      Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);

//       //For RECO muon efficiency: adhoc run2011B averaging
//       Double_t scale     = (2.1/4.7)*1 + (2.6/4.7)*(hDataEff2->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
//       Double_t scaleErrl = (2.6/4.7)*scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
//       Double_t scaleErrh = (2.6/4.7)*scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
 
      txtfile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);

      txtfile << endl;
      txtfile.flags(flags);

      cout << "Set " << iy << " " << ix << " " << scale << endl;
      h2_results_selection->SetCellContent(iy,ix, scale);
      h2_results_selection->SetCellError(iy,ix, (scaleErrl + scaleErrh)/2);
      //fill overflow bins with the same values as last bin
      if (ix == nx) {
      cout << "Set " << iy << " " << ix+1 << " " << scale << endl;
        h2_results_selection->SetCellContent(iy,nx+1, scale);
        h2_results_selection->SetCellError(iy,nx+1, (scaleErrl + scaleErrh)/2);
      }
    }
    if (iy == ny) {
      for(Int_t ix=1; ix<=nx; ix++) {
        Double_t mceff  = hMCEff->GetCellContent(ix,iy);
        Double_t mcerrl = hMCErrl->GetCellContent(ix,iy);
        Double_t mcerrh = hMCErrh->GetCellContent(ix,iy);
        Double_t dataeff  = hDataEff2->GetCellContent(ix,iy);
        Double_t dataerrl = hDataErrl2->GetCellContent(ix,iy);
        Double_t dataerrh = hDataErrh2->GetCellContent(ix,iy);
        Double_t scale     = (hDataEff2->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
        Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
        Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
        
        cout << "Set " << iy+1 << " " << ix << " " << scale << endl;
        h2_results_selection->SetCellContent(iy+1,ix, scale);
        h2_results_selection->SetCellError(iy+1,ix, (scaleErrl + scaleErrh)/2);
        //fill overflow bins with the same values as last bin
        if (ix == nx) {
          cout << "Set " << iy+1 << " " << ix+1 << " " << scale << endl;
          h2_results_selection->SetCellContent(iy+1,nx+1, scale);
          h2_results_selection->SetCellError(iy+1,nx+1, (scaleErrl + scaleErrh)/2);
        }
      }
    }

    txtfile << endl;
  }
  txtfile.close();
  
  cout << outfname << " created!" << endl;
  



  //--------------------------------------------------------------------------------------------------------------
  // Create TEX table
  //==============================================================================================================   

  ofstream texfile;
  texfile.open((outputDir + "/eff_table.tex").c_str());
  assert(texfile.is_open());

  texfile << " \\begin{table}[!ht]" << endl;
  texfile << " \\begin{center} " << endl;
  texfile << " \\begin{tabular}{|c|c|c|c|}" << endl;
  texfile << " \\hline\n";


  texfile << " $p_{T}$ / $\\eta$ bin    &  Monte Carlo Efficiency    &  Data Efficiency   &  MC to Data Scale Factor \\\\  ";
  texfile << " \\hline           ";
  texfile << endl;
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {

      string binLabel = Form("$%5.1f < p_{T} \\le %5.1f$ , $%5.1f  \\le |\\eta| < %5.1f$", 
                             hMCEff->GetYaxis()->GetBinLowEdge(iy), hMCEff->GetYaxis()->GetBinLowEdge(iy+1),
                             hMCEff->GetXaxis()->GetBinLowEdge(ix), hMCEff->GetXaxis()->GetBinLowEdge(ix+1));
      if (iy == ny) {
        binLabel = Form("$%5.1f < p_{T} $ , $%5.1f  \\le |\\eta| < %5.1f$", 
                        hMCEff->GetYaxis()->GetBinLowEdge(iy), 
                        hMCEff->GetXaxis()->GetBinLowEdge(ix), hMCEff->GetXaxis()->GetBinLowEdge(ix+1));
      }
      
      texfile << binLabel;
      texfile << "   &   ";

      ios_base::fmtflags flags = texfile.flags();
      texfile.precision(4);
      
      Double_t mceff  = hMCEff->GetCellContent(ix,iy);
      Double_t mcerrl = hMCErrl->GetCellContent(ix,iy);
      Double_t mcerrh = hMCErrh->GetCellContent(ix,iy);
      texfile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);
      texfile << "   &   ";

      Double_t dataeff  = hDataEff2->GetCellContent(ix,iy);
      Double_t dataerrl = hDataErrl2->GetCellContent(ix,iy);
      Double_t dataerrh = hDataErrh2->GetCellContent(ix,iy);
      texfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      texfile << "   &   ";
     
      Double_t scale     = (hDataEff2->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
      Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
      texfile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);
      texfile << "   \\\\   ";
      texfile << endl;
      texfile << "\\hline";
      texfile << endl;

    }
  }
  texfile << "\\end{tabular}" << endl;
  texfile << "\\caption{CAPTION.}" << endl;
  texfile << "\\label{tab:eff_ele_offline}" << endl;
  texfile << "\\end{center}" << endl;
  texfile << "\\end{table}" << endl;
  texfile.close();
  cout << outputDir + "/eff_table.tex" << " created!" << endl;




  //--------------------------------------------------------------------------------------------------------------
  // Create TWIKI table
  //==============================================================================================================   

  ofstream twikifile;
  twikifile.open((outputDir + "/eff_table.twiki").c_str());
  assert(twikifile.is_open());
  

  twikifile << "| *pT bin* | *eta bin* | *MC Efficiency* | *Data Efficiency* | *ScaleFactor* | \n";
  
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {

      string binLabel = Form("| *%4.1f < pT \\le %4.1f* | *%3.1f  <= eta < %3.1f* | ", 
                             hMCEff->GetYaxis()->GetBinLowEdge(iy), hMCEff->GetYaxis()->GetBinLowEdge(iy+1),
                             hMCEff->GetXaxis()->GetBinLowEdge(ix), hMCEff->GetXaxis()->GetBinLowEdge(ix+1));
      if (iy == ny) {
        binLabel = Form("| *%4.1f < pT* | *%3.1f  <= eta < %3.1f* | ", 
                        hMCEff->GetYaxis()->GetBinLowEdge(iy), 
                        hMCEff->GetXaxis()->GetBinLowEdge(ix), hMCEff->GetXaxis()->GetBinLowEdge(ix+1));
      }
      
      twikifile << binLabel ;

      ios_base::fmtflags flags = twikifile.flags();
      twikifile.precision(4);
      
      Double_t mceff  = hMCEff->GetCellContent(ix,iy);
      Double_t mcerrl = hMCErrl->GetCellContent(ix,iy);
      Double_t mcerrh = hMCErrh->GetCellContent(ix,iy);
      twikifile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);
      twikifile << " | ";

      Double_t dataeff  = hDataEff2->GetCellContent(ix,iy);
      Double_t dataerrl = hDataErrl2->GetCellContent(ix,iy);
      Double_t dataerrh = hDataErrh2->GetCellContent(ix,iy);
      twikifile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      twikifile << " | ";
     
      Double_t scale     = (hDataEff2->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
      Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
      twikifile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);
      twikifile << " |\n";

    }
  }
  twikifile.close();
  cout << outputDir + "/eff_table.twiki" << " created!" << endl;

  outputFile->WriteTObject(h2_results_selection, h2_results_selection->GetName(), "WriteDelete");
  outputFile->Close();

}
