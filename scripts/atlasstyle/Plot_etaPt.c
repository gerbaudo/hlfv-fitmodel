
#include "AtlasUtils.h"
#include "AtlasLabels.h"
#include "AtlasStyle.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TString.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TGraphErrors.h"
#include "TPad.h"


  SetAtlasStyle();
#ifdef __CINT__
  gROOT->LoadMacro("AtlasUtils.C");
#endif


void Plot_etaPt(TString file)
{

	SetAtlasStyle();
	#ifdef __CINT__
	  gROOT->LoadMacro("AtlasUtils.C");
	#endif

  	double xmin = 15;
  	double xmax = 80;


  	TFile* f = new TFile(file);

  	TH2D* h_mu =(TH2D*)f->Get("nom/mu_pt_vs_eta");
  	TH2D* h_el =(TH2D*)f->Get("nom/el_pt_vs_eta");

  	h_mu->GetXaxis()->SetRangeUser(xmin,xmax);
  	h_mu->GetXaxis()->SetTitle("p_{T}^{#mu}");
	h_mu->GetYaxis()->SetTitle("#eta^{#mu}");
  	h_mu->GetYaxis()->SetRangeUser(-2.8,2.8);
  	h_el->GetXaxis()->SetRangeUser(xmin,xmax);
  	h_el->GetXaxis()->SetTitle("p_{T}^{e}");
  	h_el->GetYaxis()->SetTitle("#eta^{e}");
  	h_el->GetYaxis()->SetRangeUser(-2.8,2.8);

//  	h_mu->RebinX(4); h_mu->RebinY(4);
//  	h_el->RebinX(4); h_el->RebinY(4);

  	gStyle->SetPalette(1);
  	gStyle->SetPadRightMargin(0.2);
  	TString op("colz");

  	TCanvas* c_el = new TCanvas("el pt vs. eta","el pt vs. eta",600,600); c_el=c_el;
  	h_el->Draw(op);
	#ifdef __CINT__
  		gROOT->LoadMacro("AtlasLabels.C");
	#endif

  	myText(0.3,0.65,1,"#int L dt = 20.3 fb^{-1}  #sqrt{s} = 8 TeV");
  	ATLASLabel(0.55,0.7,"Internal");


  	TCanvas* c_mu = new TCanvas("mu pt vs. eta","mu pt vs. eta",600,600); c_mu=c_mu;
  	h_mu->Draw(op);

	#ifdef __CINT__
  		gROOT->LoadMacro("AtlasLabels.C");
	#endif

  	myText(0.3,0.65,1,"#int L dt = 20.3 fb^{-1}  #sqrt{s} = 8 TeV");
  	ATLASLabel(0.55,0.7,"Internal");

  	c_el->Update(); c_mu->Update();


  return;
}

#ifndef __CINT__
int main() { 
  BasicExample();
//  gPad->Print("basic.png");
  return 0;
}
#endif
