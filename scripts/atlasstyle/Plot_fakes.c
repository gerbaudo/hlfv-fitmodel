
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
#include "TRandom.h"
#include "TGraphErrors.h"
#include "TPad.h"

void Plot_fakes(TString filename)
{
  SetAtlasStyle();
#ifdef __CINT__
  gROOT->LoadMacro("AtlasUtils.C");
#endif

	TFile* f = new TFile(filename);
	TH1D* hEM = (TH1D*)f->Get("nom/EM_McollHiggs_Unblind");
        TH1D* hME = (TH1D*)f->Get("nom/ME_McollHiggs_Unblind");

//	double nbins = hEM->GetXaxis()->GetNbins();
//	Bins = new double[nbins+1]; 
//	const TArrayD* array = hEM->GetXaxis()->GetXbins();
//	const Double_t *bins = array->GetArray();
//	memcpy(Bins,bins,(nbins+1)*sizeof(double));

	double xmax = 400;
	double xmin = 0;
	gStyle->SetErrorX(0.1); gStyle->SetOptStat(0);
	TCanvas* c1 = new TCanvas("Fake Estimation","Fake Estimation",600,600); c1=c1;
	TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0,21);
	TPad *pad2 = new TPad("pad2", "ratio",0.0,0.0,1.0,0.2,22);
	pad1->SetFillColor(0);pad2->SetFillColor(0);
	pad1->Draw(); pad2->Draw();
	pad1->cd(); //gPad-> SetLogy();
	TLegend* leg = new TLegend(0.5,0.7,0.7,0.9);
	leg->SetFillColor(kWhite); leg->SetBorderSize(1); leg->SetLineColor(0); leg->SetTextFont(42);
	leg->SetTextSize(.05);

	hEM->GetXaxis()->SetTitle("M_{Collinear} (GeV)"); hEM->GetXaxis()->SetRangeUser(xmin,xmax);
	hEM->GetYaxis()->SetTitle("Events/10 GeV");
	hEM->SetLineColor(kGreen+2); hEM->SetMarkerStyle(2); hEM->SetMarkerColor(kGreen+2);hEM->Draw("e1");
	hME->GetXaxis()->SetTitle("M_{Collinear} (GeV)");hME->SetLineColor(kBlue); 
	hME->SetLineWidth(2); hME->GetXaxis()->SetRangeUser(xmin,xmax);
	hME->SetMarkerStyle(2); hME->SetMarkerColor(kBlue);hME->Draw("e1 sames");

	TLine *vline1 = new TLine(100,0,100,1000);
	TLine *vline2 = new TLine(170,0,170,1000);
//	vline1->SetLineStyle(2); vline2->SetLineStyle(2); vline1->Draw(); vline2->Draw();

	TLine *hline = new TLine(0,0,400,0); hline->SetLineStyle(2);hline->Draw();

	myText(0.7,0.85,1,"#sqrt{s} = 8 TeV");
	#ifdef __CINT__
	  gROOT->LoadMacro("AtlasLabels.C");
	#endif
  	ATLASLabel(0.6,0.5,"Internal");

	//draw uncertainty

	leg->AddEntry(hME,"#mue","l");
	leg->AddEntry(hEM,"e#mu","l");
	//leg->AddEntry(hsignal,"BR(h#rightarrow#tau#mu) = 1%");
	leg->Draw();

	pad2->cd();
	TH1D* ratio1 = (TH1D*)hEM->Clone("ratio1");
	ratio1->Divide(hME);
//	TH1D* ratio2 = (TH1D*)hME_blind->Clone("ratio2");
//	ratio2->Divide(h_b);
	ratio1->GetYaxis()->SetLabelSize(0.1); ratio1->GetXaxis()->SetLabelSize(0.1);
	ratio1->SetLineColor(kBlack); ratio1->SetMarkerStyle(kFullCircle); ratio1->SetMarkerSize(0.8);
	ratio1->GetYaxis()->SetTitle("EM/ME");
	//ratio1-GetYaxis()->SetTitleSize(0.04);
	ratio1->SetMarkerColor(kBlack);// ratio2->SetMarkerColor(kRed);ratio2->GetYaxis()->SetLabelSize(0.1); ratio2->GetXaxis()->SetLabelSize(0.1);
//	ratio2->SetLineColor(kBlack); ratio2->SetMarkerStyle(kFullCircle); ratio2->SetMarkerSize(0.8);
//	TH1D* uncertainty = new TH1D("uncertainty","uncertainty",nbins, Bins);
//	for (int i=1; i<=nbins; i++){
//		double sigb = h_b->GetBinError(i);
//		double b = h_b->GetBinContent(i);
//		double n = hEM->GetBinContent(i);
//		double m = hME->GetBinContent(i);
//		ratio1->SetBinError(i,(1/b)*sqrt(n));
//		ratio2->SetBinError(i,(1/b)*sqrt(m));
//		uncertainty->SetBinError(i,sigb/b); uncertainty->SetBinContent(i,1);
//		if (ratio1->GetBinContent(i)==0) {ratio1->SetBinContent(i,1); ratio1->SetBinError(i,0);}
//		if (ratio2->GetBinContent(i)==0) {ratio2->SetBinContent(i,1); ratio2->SetBinError(i,0);}
//	}
//	uncertainty->SetFillStyle(3004); uncertainty->SetFillColor(kBlack);
//	uncertainty->SetMarkerStyle(6); uncertainty->SetLineColor(kBlack); uncertainty->SetLineStyle(2);
	ratio1->SetLineWidth(1); //ratio2->SetLineWidth(1);
//	ratio1->SetMinimum(0); ratio1->SetMaximum(2);
//	ratio1->SetLineColor(kBlue); ratio2->SetLineColor(kRed);
	ratio1->Draw("ep"); //ratio2->Draw("ep sames");
//	uncertainty->Draw("E3 sames");
	TLine *line = new TLine(0,1,700,1);
//	TLine *vline12 = new TLine(100,0,100,2);
//	TLine *vline22 = new TLine(150,0,150,2);
//	vline12->SetLineStyle(2); vline22->SetLineStyle(2); vline12->Draw(); vline22->Draw();
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
