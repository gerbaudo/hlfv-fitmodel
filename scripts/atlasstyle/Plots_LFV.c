
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

void Plots_LFV()
{
  SetAtlasStyle();
#ifdef __CINT__
  gROOT->LoadMacro("AtlasUtils.C");
#endif

	TFile* f = new TFile("histos.root");
	TH1D* hEM = (TH1D*)f->Get("EM_McollHiggs_Unblind");
        TH1D* hME = (TH1D*)f->Get("ME_McollHiggs_Unblind");
	TH1D* hEM_blind = (TH1D*)f->Get("EM_McollHiggs_Unblind_blind");
	TH1D* hME_blind = (TH1D*)f->Get("ME_McollHiggs_Unblind_blind");
	TH1D* hsignal = (TH1D*)f->Get("h_signal");
	TH1D* h_b = (TH1D*)f->Get("hb");

	//hsignal->Add(h_b);
	hsignal->SetLineStyle(2);
	hsignal->SetMarkerStyle(0);
	double nbins = hEM_blind->GetXaxis()->GetNbins();
	Bins = new double[nbins+1]; 
	const TArrayD* array = hEM->GetXaxis()->GetXbins();
	const Double_t *bins = array->GetArray();
	memcpy(Bins,bins,(nbins+1)*sizeof(double));

	double xmax = 400;
	double xmin = 0;
	gStyle->SetErrorX(0.1); gStyle->SetOptStat(0);
	TCanvas* c1 = new TCanvas("BG Estimation","BG Estimation",600,600); c1=c1;
	TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
	TPad *pad2 = new TPad("pad2", "ratio",0.0,0.0,1.0,0.2,22);
	pad1->SetFillColor(0);pad2->SetFillColor(0);
	pad1->Draw(); pad2->Draw();
	pad1->cd(); //gPad-> SetLogy();
	TLegend* leg = new TLegend(0.5,0.7,0.7,0.9);
	leg->SetFillColor(kWhite); leg->SetBorderSize(1); leg->SetLineColor(0); leg->SetTextFont(42);

	h_b->GetXaxis()->SetTitle("M_{Collinear} (GeV)"); h_b->GetXaxis()->SetRangeUser(xmin,xmax);
	h_b->GetYaxis()->SetTitle("Events/10 GeV");
	h_b->SetLineColor(kBlack); h_b->SetLineWidth(2);
	h_b->GetXaxis()->SetRangeUser(xmin,xmax); h_b->SetLineColor(kGreen+2); h_b->Draw("e1");
	h_b->SetMarkerSize(0.2);
	hME_blind->GetXaxis()->SetTitle("M_{Collinear} (GeV)");hME_blind->SetLineColor(kRed); 
	hME_blind->SetLineWidth(2); hME_blind->GetXaxis()->SetRangeUser(xmin,xmax);
	hEM_blind->SetMarkerStyle(2);hEM_blind->SetMarkerColor(kBlue);
	hME_blind->SetMarkerStyle(2);hME_blind->SetMarkerColor(kRed);
	h_b->SetMarkerStyle(6); 
	hME_blind->Draw("e1 sames");
	hEM_blind->SetLineColor(kBlue); hEM_blind->SetLineWidth(2); 
	hEM_blind->GetXaxis()->SetRangeUser(xmin,xmax);
	hEM_blind->Draw("e1 sames");
	hsignal->Draw("hist sames");
	h_b->Draw("hist sames");


	TLine *vline1 = new TLine(100,0,100,1000);
	TLine *vline2 = new TLine(170,0,170,1000);
	vline1->SetLineStyle(2); vline2->SetLineStyle(2); vline1->Draw(); vline2->Draw();

	myText(0.7,0.85,1,"#sqrt{s} = 8 TeV");
	#ifdef __CINT__
	  gROOT->LoadMacro("AtlasLabels.C");
	#endif
  	ATLASLabel(0.6,0.5,"Internal");

	//draw uncertainty
	TH1D* Buncertainty = new TH1D("Buncertainty","Buncertainty",nbins, Bins);
	for (int i=1; i<=nbins; i++){
		double sigb = h_b->GetBinError(i);
		double b = h_b->GetBinContent(i);
		Buncertainty->SetBinError(i,sigb); Buncertainty->SetBinContent(i,b);
	}
	Buncertainty->SetFillStyle(3004); Buncertainty->SetFillColor(kBlack);
	Buncertainty->SetLineColor(kBlack); Buncertainty->SetLineStyle(2);
	Buncertainty->SetMarkerStyle(6);
	Buncertainty->Draw("E3 sames");


	leg->AddEntry(hME_blind,"#mue","l");
	leg->AddEntry(hEM_blind,"e#mu","l");
	leg->AddEntry(h_b,"BG estimation","l");
	leg->AddEntry(Buncertainty,"BG uncertainty","f");
	leg->AddEntry(hsignal,"BR(h#rightarrow#tau#mu) = 1%");
	leg->Draw();

	pad2->cd();
	TH1D* ratio1 = (TH1D*)hEM_blind->Clone("ratio1");
	ratio1->Divide(h_b);
	TH1D* ratio2 = (TH1D*)hME_blind->Clone("ratio2");
	ratio2->Divide(h_b);
	ratio1->GetYaxis()->SetLabelSize(0.1); ratio1->GetXaxis()->SetLabelSize(0.1);
	ratio1->SetLineColor(kBlack); ratio1->SetMarkerStyle(kFullCircle); ratio1->SetMarkerSize(0.8);
	ratio1->SetMarkerColor(kBlue); ratio2->SetMarkerColor(kRed);ratio2->GetYaxis()->SetLabelSize(0.1); ratio2->GetXaxis()->SetLabelSize(0.1);
	ratio2->SetLineColor(kBlack); ratio2->SetMarkerStyle(kFullCircle); ratio2->SetMarkerSize(0.8);
	TH1D* uncertainty = new TH1D("uncertainty","uncertainty",nbins, Bins);
	for (int i=1; i<=nbins; i++){
		double sigb = h_b->GetBinError(i);
		double b = h_b->GetBinContent(i);
		double n = hEM->GetBinContent(i);
		double m = hME->GetBinContent(i);
		ratio1->SetBinError(i,(1/b)*sqrt(n));
		ratio2->SetBinError(i,(1/b)*sqrt(m));
		uncertainty->SetBinError(i,sigb/b); uncertainty->SetBinContent(i,1);
		if (ratio1->GetBinContent(i)==0) {ratio1->SetBinContent(i,1); ratio1->SetBinError(i,0);}
		if (ratio2->GetBinContent(i)==0) {ratio2->SetBinContent(i,1); ratio2->SetBinError(i,0);}
	}
	uncertainty->SetFillStyle(3004); uncertainty->SetFillColor(kBlack);
	uncertainty->SetMarkerStyle(6); uncertainty->SetLineColor(kBlack); uncertainty->SetLineStyle(2);
	ratio1->SetLineWidth(1); ratio2->SetLineWidth(1);
	ratio1->SetMinimum(0); ratio1->SetMaximum(2);
	ratio1->SetLineColor(kBlue); ratio2->SetLineColor(kRed);
	ratio1->Draw("ep"); ratio2->Draw("ep sames");
	uncertainty->Draw("E3 sames");
	TLine *line = new TLine(0,1,700,1);
	TLine *vline12 = new TLine(100,0,100,2);
	TLine *vline22 = new TLine(150,0,150,2);
	vline12->SetLineStyle(2); vline22->SetLineStyle(2); vline12->Draw(); vline22->Draw();
	line->Draw();  




  return;
}

#ifndef __CINT__
int main() { 
  BasicExample();
//  gPad->Print("basic.png");
  return 0;
}
#endif
