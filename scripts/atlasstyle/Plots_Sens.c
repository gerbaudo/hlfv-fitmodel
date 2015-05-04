
#include "AtlasUtils.h"
#include "AtlasLabels.h"
#include "AtlasStyle.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TGraphErrors.h"
#include "TPad.h"

void Plots_Sens()
{
  SetAtlasStyle();
#ifdef __CINT__
  gROOT->LoadMacro("AtlasUtils.C");
#endif

	TFile* f = new TFile("sens.root");
	TH1D* hsens = (TH1D*)f->Get("sens");
	TH1D* hgr1 = (TH1D*)f->Get("greenquantile");
	TH1D* hgr2 = (TH1D*)f->Get("greenquantile2");
	TH1D* hyl1 = (TH1D*)f->Get("yellowquantile");
	TH1D* hyl2 = (TH1D*)f->Get("yellowquantile2");

	TCanvas* c = new TCanvas("sensitivity","sensitivity",600,600);
	hsens->GetXaxis()->SetTitle("#mu (%BR)"); hsens->Draw();
	hgr1->Draw("sames");hgr2->Draw("sames");hyl1->Draw("sames");hyl2->Draw("sames");
	myText(0.7,0.85,1,"#sqrt{s} = 8 TeV");
	myText(0.5,0.6,1,"median 2#sigma sensitivity (1.37#pm 0.34)%");
	#ifdef __CINT__
	  gROOT->LoadMacro("AtlasLabels.C");
	#endif
  	ATLASLabel(0.6,0.5,"Internal");



  return;
}

#ifndef __CINT__
int main() { 
  BasicExample();
//  gPad->Print("basic.png");
  return 0;
}
#endif
