#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TStyle.h>
#include "PlotStyle.hh"

void PlotStyle() {
  SetStyle();
}

TCanvas* MakeCanvas(const char* name, const char *title, int dX, int dY)
{
  // Start with a canvas
  TCanvas *canvas = new TCanvas(name,title,0,0,dX,dY);
  // General overall stuff
  canvas->SetFillColor      (0);
  canvas->SetBorderMode     (0);
  canvas->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canvas->SetLeftMargin     (0.18);
  canvas->SetRightMargin    (0.05);
  canvas->SetTopMargin      (0.08);
  canvas->SetBottomMargin   (0.15);
  // Setup a frame which makes sense
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);

  return canvas;
}

void InitSubPad(TPad* pad, int i)
{
  //printf("Pad: %p, index: %d\n",pad,i);

  pad->cd(i);
  TPad *tmpPad = (TPad*) pad->GetPad(i);
  tmpPad->SetLeftMargin  (0.18);
  tmpPad->SetTopMargin   (0.05);
  tmpPad->SetRightMargin (0.07);
  tmpPad->SetBottomMargin(0.15);
  return;
}

void InitHist(TH1 *hist, const char *xtit, const char *ytit, EColor color)
{
  hist->SetXTitle(xtit);
  hist->SetYTitle(ytit);
  hist->SetLineColor(color);
  hist->SetTitleSize  (0.055,"Y");
  hist->SetTitleOffset(1.600,"Y");
  hist->SetLabelOffset(0.014,"Y");
  hist->SetLabelSize  (0.050,"Y");
  hist->SetLabelFont  (42   ,"Y");
  hist->SetTitleSize  (0.055,"X");
  hist->SetTitleOffset(1.300,"X");
  hist->SetLabelOffset(0.014,"X");
  hist->SetLabelSize  (0.050,"X");
  hist->SetLabelFont  (42   ,"X");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize (0.6);
  // Strangely enough this cannot be set anywhere else??
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->SetTitle("");  
  return;
}

void SetStyle()
{
  TStyle *PlotStyle = new TStyle("PlotStyle","PlotStyle");
  gStyle = PlotStyle;

  // Canvas
  PlotStyle->SetCanvasColor     (0);
  PlotStyle->SetCanvasBorderSize(10);
  PlotStyle->SetCanvasBorderMode(0);
  PlotStyle->SetCanvasDefH      (700);
  PlotStyle->SetCanvasDefW      (700);
  PlotStyle->SetCanvasDefX      (100);
  PlotStyle->SetCanvasDefY      (100);

  // Pads
  PlotStyle->SetPadColor       (0);
  PlotStyle->SetPadBorderSize  (10);
  PlotStyle->SetPadBorderMode  (0);
  PlotStyle->SetPadBottomMargin(0.13);
  PlotStyle->SetPadTopMargin   (0.08);
  PlotStyle->SetPadLeftMargin  (0.15);
  PlotStyle->SetPadRightMargin (0.05);
  PlotStyle->SetPadGridX       (0);
  PlotStyle->SetPadGridY       (0);
  PlotStyle->SetPadTickX       (0);
  PlotStyle->SetPadTickY       (0);

  // Frames
  PlotStyle->SetFrameFillStyle ( 0);
  PlotStyle->SetFrameFillColor ( 0);
  PlotStyle->SetFrameLineColor ( 1);
  PlotStyle->SetFrameLineStyle ( 0);
  PlotStyle->SetFrameLineWidth ( 1);
  PlotStyle->SetFrameBorderSize(10);
  PlotStyle->SetFrameBorderMode( 0);

  // Histograms
  PlotStyle->SetHistFillColor(2);
  PlotStyle->SetHistFillStyle(0);
  PlotStyle->SetHistLineColor(1);
  PlotStyle->SetHistLineStyle(0);
  PlotStyle->SetHistLineWidth(2);
  PlotStyle->SetNdivisions(505);

  // Functions
  PlotStyle->SetFuncColor(1);
  PlotStyle->SetFuncStyle(0);
  PlotStyle->SetFuncWidth(2);

  // Various
  PlotStyle->SetMarkerStyle(20);
  PlotStyle->SetMarkerColor(kBlack);
  PlotStyle->SetMarkerSize (1.2);

  PlotStyle->SetTitleBorderSize(0);
  PlotStyle->SetTitleFillColor (0);
  PlotStyle->SetTitleX         (0.2);

  PlotStyle->SetTitleSize  (0.055,"X");
  PlotStyle->SetTitleOffset(1.200,"X");
  PlotStyle->SetLabelOffset(0.005,"X");
  PlotStyle->SetLabelSize  (0.050,"X");
  PlotStyle->SetLabelFont  (42   ,"X");

  PlotStyle->SetStripDecimals(kFALSE);

  PlotStyle->SetTitleSize  (0.055,"Y");
  PlotStyle->SetTitleOffset(1.600,"Y");
  PlotStyle->SetLabelOffset(0.010,"Y");
  PlotStyle->SetLabelSize  (0.050,"Y");
  PlotStyle->SetLabelFont  (42   ,"Y");

  PlotStyle->SetTextSize   (0.055);
  PlotStyle->SetTextFont   (42);

  PlotStyle->SetStatFont   (42);
  PlotStyle->SetTitleFont  (42);
  PlotStyle->SetTitleFont  (42,"X");
  PlotStyle->SetTitleFont  (42,"Y");

  PlotStyle->SetOptStat    (0);
  return;
}
