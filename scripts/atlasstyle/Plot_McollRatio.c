
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


  SetAtlasStyle();
#ifdef __CINT__
  gROOT->LoadMacro("AtlasUtils.C");
#endif
//
// TH1D* GetEM(TString filename,TString title,double xmin, double xmax)
//{
//  	double xbins[]= {0,100,150,200,300,1000};
//  	double xbins2[]= {0,100,150,200,300,400};
//  	TFile* f = new TFile(filename);
//  	TH1D* hEM =(TH1D*)f->Get("nom/EM_McollHiggs_Unblind");
////  	hEM->SetLineColor(kGreen+2); hEM->SetLineWidth(2);
//  	TH1D* htemp = (TH1D*)hEM->Rebin(5,"hEMRebinned",xbins);
//  	TH1D* hEMrebinned = new TH1D("hEMrebinned","hEMrebinned",5,xbins2);
//  	for (int bin=1; bin<6; bin++){
//  		hEMrebinned->SetBinContent(bin,htemp->GetBinContent(bin));
//  		hEMrebinned->SetBinError(bin,htemp->GetBinError(bin));
//  	}
//  	hEMrebinned->GetYaxis()->SetTitleOffset(0.2);
//  	hEMrebinned->GetXaxis()->SetRangeUser(xmin,xmax);
//  	hEMrebinned->SetMarkerStyle(8); hEMrebinned->SetMarkerSize(0.7);
//  	hEMrebinned->SetLineColor(kPink + 8);hEMrebinned->SetOption("e1");
//  	hEMrebinned->SetMarkerColor(kPink + 8);
//  	hEMrebinned->SetTitle(title); //hEMrebinned->GetXaxis()->SetTitle("M_{coll} (GeV)");
//  	hEMrebinned->GetXaxis()->SetLabelOffset(0);
//  	hEMrebinned->GetXaxis()->SetLabelSize(0);
//    delete htemp;
//  	return hEMrebinned;
// }
//
//  TH1D* GetME(TString filename,TString title,double xmin, double xmax)
//  {
//  	double xbins1[]= {0,100,150,200,300,1000};
//  	double xbins2[]= {0,100,150,200,300,400};
//  	TFile* f = new TFile(filename);
//  	TH1D* hME =(TH1D*)f->Get("nom/ME_McollHiggs_Unblind");
//  	hME->SetLineColor(kBlue); hME->SetLineWidth(2);
// 	TH1D* htemp = (TH1D*)hME->Rebin(5,"hMERebinned",xbins1);
//  	TH1D* hMErebinned = new TH1D("hMERebinned","hMERebinned",5,xbins2);
//  	for (int bin=1; bin<6; bin++){
//  		hMERebinned->SetBinContent(bin,htemp->GetBinContent(bin));
//  		hMERebinned->SetBinError(bin,htemp->GetBinError(bin));
//  	}
//  	hMErebinned->GetYaxis()->SetTitleOffset(0.2);
//  	hMErebinned->GetXaxis()->SetRangeUser(xmin,xmax);
//  	hMErebinned->SetMarkerStyle(8); hMErebinned->SetMarkerSize(0.7);
//  	hMErebinned->SetLineColor(kTeal - 6);hMErebinned->SetOption("e1");
//  	hMErebinned->SetMarkerColor(kTeal - 6);
//  	hMErebinned->SetTitle(title);
//  	//hMErebinned->GetXaxis()->SetTitle("M_{coll} (GeV)");
//  	hMErebinned->GetXaxis()->SetLabelOffset(0);
//  	hMErebinned->GetXaxis()->SetLabelSize(0);
//  	delete htemp;
//  	return hMErebinned;
//  }
//
//  TH1D* GetEMblind(TString filename,TString title,double xmin, double xmax)
//  {
//  	double xbins[]= {0,100,150,200,300,1000};
//  	double xbins2[]= {0,100,150,200,300,400};
//  	TFile* f = new TFile(filename);
//  	TH1D* hEM =(TH1D*)f->Get("nom/EM_McollHiggs");
//  	hEM->SetLineColor(kGreen+2); hEM->SetLineWidth(2);
//  	TH1D* htemp = (TH1D*)hEM->Rebin(5,"hEMRebinned",xbins);
//  	TH1D* hEMrebinned = new TH1D("hEMrebinned","hEMrebinned",5,xbins2);
//  	for (int bin=1; bin<6; bin++){
//  		hEMrebinned->SetBinContent(bin,htemp->GetBinContent(bin));
//  		hEMrebinned->SetBinError(bin,htemp->GetBinError(bin));
//  	}
//  	hEMrebinned->GetYaxis()->SetTitleOffset(0.2);
//  	hEMrebinned->GetXaxis()->SetRangeUser(xmin,xmax);
//  	hEMrebinned->SetMarkerStyle(8); hEMrebinned->SetMarkerSize(0.7);
//  	hEMrebinned->SetLineColor(kPink + 8);hEMrebinned->SetOption("e1");
//  	hEMrebinned->SetMarkerColor(kPink + 8);
//  	hEMrebinned->SetTitle(title); hEMrebinned->SetTitleSize(0.2); hEMrebinned->SetTitleOffset(0.4);
//  	hEMrebinned->GetXaxis()->SetLabelOffset(0);
//  	hEMrebinned->GetXaxis()->SetLabelSize(0);
//  	delete htemp;
//  	return hEMrebinned;
//  }
//
//  TH1D* GetMEblind(TString filename,TString title,double xmin, double xmax)
//  {
//  	double xbins[]= {0,100,150,200,300,1000};
//  	double xbins2[]= {0,100,150,200,300,400};
//  	TFile* f = new TFile(filename);
//  	TH1D* hME =(TH1D*)f->Get("nom/ME_McollHiggs");
//  	TH1D* htemp = (TH1D*)hME->Rebin(5,"hMERebinned",xbins);
//  	TH1D* hMErebinned = new TH1D("hMERebinned","hMERebinned",5,xbins2);
//  	for (int bin=1; bin<6; bin++){
//  		hMERebinned->SetBinContent(bin,htemp->GetBinContent(bin));
//  		hMERebinned->SetBinError(bin,htemp->GetBinError(bin));
//  	}
//  	hMErebinned->GetYaxis()->SetTitleOffset(0.2);
//  	hMErebinned->GetXaxis()->SetRangeUser(xmin,xmax);
//  	hMErebinned->SetMarkerStyle(8); hMErebinned->SetMarkerSize(0.7);
//  	hMErebinned->SetLineColor(kTeal - 6);hMErebinned->SetOption("e1");
//  	hMErebinned->SetMarkerColor(kTeal - 6);
//  	hMErebinned->SetTitle(title); hMErebinned->SetTitleSize(0.02); hMErebinned->SetTitleOffset(0.4);
//  	hMErebinned->GetXaxis()->SetLabelOffset(0);
//  	hMErebinned->GetXaxis()->SetLabelSize(0);
//  	delete htemp;
//  	return hMErebinned;
//  }

  TH1D* GetRatio(TH1D* hEM,TH1D* hME,double xmin, double xmax)
  {
  	TH1D* ratio = (TH1D*)hEM->Clone("ratio");
  	ratio->Divide(hME);
//  	ratio->GetYaxis()->SetTitle("Ratio"); //ratio->SetTitleFont(64);
  	ratio->GetYaxis()->SetTitle("Ratio"); ratio->GetYaxis()->SetTitleSize(0.1);
  	ratio->GetYaxis()->SetTitleOffset(0.3);
  	ratio->GetYaxis()->CenterTitle();
//  	ratio->SetTitleSize(0.1);
  	ratio->GetXaxis()->SetRangeUser(xmin,xmax); ratio->SetLineColor(kBlack);
  	ratio->GetYaxis()->SetRangeUser(0,2);
  	ratio->SetLineWidth(2); ratio->SetMarkerStyle(8); ratio->SetMarkerSize(0.7);
  	ratio->SetMarkerColor(kBlack);
  	ratio->GetXaxis()->SetTitle("M_{coll} (GeV)"); ratio->GetXaxis()->SetTitleSize(0.15);
  	ratio->GetXaxis()->SetTitleOffset(0.8);
	ratio->GetXaxis()->SetLabelOffset();
  	ratio->GetYaxis()->SetLabelSize(0.1);
  	ratio->GetYaxis()->SetNdivisions(5);
  	ratio->GetXaxis()->SetLabelSize(0.1);

  	ratio->GetXaxis()->SetRangeUser(xmin,xmax);
  	for (int i=1; i<=hEM->GetXaxis()->FindBin(99); i++){
  		double n = hEM->GetBinContent(i);
  		double m = hME->GetBinContent(i);
  		double deltaN = hEM->GetBinError(i);
  		double deltaM = hME->GetBinError(i);
  		double err = (1./m)*TMath::Sqrt(TMath::Power(deltaN,2)+TMath::Power(n*deltaM/m,2));
  		ratio->SetBinError(i,err);
  	}
	for (int i=hEM->GetXaxis()->FindBin(151); i<=hEM->GetXaxis()->GetNbins(); i++){
  		double n = hEM->GetBinContent(i);
  		double m = hME->GetBinContent(i);
  		double deltaN = hEM->GetBinError(i);
  		double deltaM = hME->GetBinError(i);
  		double err = (1./m)*TMath::Sqrt(TMath::Power(deltaN,2)+TMath::Power(n*deltaM/m,2));
  		ratio->SetBinError(i,err);
  	}
  	return ratio;
  }

  TH1D* GetDiff(TH1D* hEM,TH1D* hME,double xmin, double xmax)
  {
  	TH1D* diff = (TH1D*)hEM->Clone("diff");
  	diff->Add(hME,-1);
  	diff->GetYaxis()->SetTitle("Diff"); diff->GetYaxis()->SetTitleSize(0.1);
  	diff->GetYaxis()->SetTitleOffset(0.3);
  	diff->GetYaxis()->CenterTitle();
  	diff->GetXaxis()->SetRangeUser(xmin,xmax); diff->SetLineColor(kBlack);
  	diff->GetYaxis()->SetLabelSize(0.08);
  	diff->SetLineWidth(2); diff->SetMarkerStyle(8); diff->SetMarkerSize(0.7);
  	diff->SetMarkerColor(kBlack);
  	diff->GetXaxis()->SetTitle("M_{coll} (GeV)"); diff->GetXaxis()->SetTitleSize(0.15);
  	diff->GetXaxis()->SetTitleOffset(0.8);
  	diff->GetXaxis()->SetLabelOffset();
  	diff->GetXaxis()->SetLabelSize(0.1);
  	return diff;
  }


void Plot_McollRatio(TString file, TString fakefile, TString sample, TString SR, bool subtractFakes)
{

	SetAtlasStyle();
	#ifdef __CINT__
	  gROOT->LoadMacro("AtlasUtils.C");
	#endif

  	double xmin = 50;
  	double xmax = 400;


  	TFile* f = new TFile(file);
  	TFile* fakef = new TFile(fakefile);

  	TH1D* hME;
  	TH1D* hEM;

  	//if (isScaled){
//  		hME =(TH1D*)f->Get("ME_scaled");//nom/ME_McollHiggs_Unblind");
//  		hEM =(TH1D*)f->Get("EM");//nom/EM_McollHiggs_Unblind");
//  		hME_orig = (TH1D*)f->Get("ME_original");
//  		hEM_orig = (TH1D*)f->Get("EM_original");
  	//}
  	//else{
  		hME =(TH1D*)f->Get("nom/ME_McollHiggs_Unblind");
  		hEM =(TH1D*)f->Get("nom/EM_McollHiggs_Unblind");
  	//}

  	if (subtractFakes){
  		TH1D* hME_fake = (TH1D*)fakef->Get("nom/ME_McollHiggs_Unblind");
  		TH1D* hEM_fake = (TH1D*)fakef->Get("nom/EM_McollHiggs_Unblind");

  		hME->Add(hME_fake,-1);
  		hEM->Add(hEM_fake,-1);
 	}
  	hEM->SetMarkerStyle(8); hEM->SetMarkerSize(0.7);
  	hEM->SetLineColor(kPink + 8);hEM->SetOption("e1");
  	hEM->SetMarkerColor(kPink + 8);
  	hME->SetMarkerStyle(8); hME->SetMarkerSize(0.7);
  	hME->SetLineColor(kTeal - 6);hME->SetOption("e1");
  	hME->SetMarkerColor(kTeal - 6);
  	hME->GetXaxis()->SetLabelOffset(0); hME->GetXaxis()->SetLabelSize(0);
  	hEM->GetXaxis()->SetLabelOffset(0); hEM->GetXaxis()->SetLabelSize(0);
  	hEM->GetXaxis()->SetRangeUser(xmin,xmax);   hME->GetXaxis()->SetRangeUser(xmin,xmax);
  	hEM->GetYaxis()->SetTitleOffset(1.0);
  	hME->GetYaxis()->SetTitleOffset(1.0);

//  	hEM_orig->SetMarkerStyle(8); hEM_orig->SetMarkerSize(0.7);
//  	hEM_orig->SetLineColor(kPink + 8);hEM_orig->SetOption("e1");
//  	hEM_orig->SetMarkerColor(kPink + 8);
//  	hME_orig->SetMarkerStyle(8); hME_orig->SetMarkerSize(0.7);
//  	hME_orig->SetLineColor(kTeal - 6);hME_orig->SetOption("e1");
//  	hME_orig->SetMarkerColor(kTeal - 6);
//  	hME_orig->GetXaxis()->SetLabelOffset(0); hME_orig->GetXaxis()->SetLabelSize(0);
//  	hEM_orig->GetXaxis()->SetLabelOffset(0); hEM_orig->GetXaxis()->SetLabelSize(0);
//  	hEM_orig->GetXaxis()->SetRangeUser(xmin,xmax);   hME_orig->GetXaxis()->SetRangeUser(xmin,xmax);
//  	hEM_orig->GetYaxis()->SetTitleOffset(1.0);
//  	hME_orig->GetYaxis()->SetTitleOffset(1.0);

  	double numEM = hEM->Integral();
  	double numME = hME->Integral();

  	cout<<"# entries EM = " << numEM << endl;
  	cout<<"# entries ME = " << numME << endl;

//  	TH1D* ratio2 = GetRatio(hEM_orig,hME_orig,xmin,xmax);
//  	TH1D* diff2 = GetDiff(hEM_orig,hME_orig,xmin,xmax);

  	TH1D* ratio = GetRatio(hEM,hME,xmin,xmax);
  	TH1D* diff = GetDiff(hEM,hME,xmin,xmax);

  	TLegend* leg = new TLegend(0.6,0.45,0.75,0.55);
  	leg->SetFillColor(kWhite); leg->SetBorderSize(1); leg->SetLineColor(0); leg->SetTextFont(42);
  	leg->SetTextSize(.05);

  	TLine* line1 = new TLine(xmin,1,xmax,1); line1->SetLineColor(kRed); line1->SetLineStyle(2);
  	TLine* line14 = new TLine(xmin,1.5,xmax,1.5); line14->SetLineStyle(9);
  	line14->SetLineColor(kBlack);
  	TLine* line13 = new TLine(xmin,2./3,xmax,2./3); line13->SetLineStyle(2);
  	line13->SetLineColor(kRed);
  	TLine* line12 = new TLine(xmin,2,xmax,2); line12->SetLineStyle(9);
  	line12->SetLineColor(kBlack);
  	TLine* line15 = new TLine(xmin,0.5,xmax,0.5); line15->SetLineStyle(9);
  	line15->SetLineColor(kBlack);
  	TLine* line2 = new TLine(xmin,0,xmax,0); line2->SetLineColor(kBlack); line2->SetLineStyle(9);

  	TLine* vline1 = new TLine(100,0,100,3000); vline1->SetLineStyle(2);
  	TLine* vline2 = new TLine(150,0,150,3000); vline2->SetLineStyle(2);

  	TCanvas* c0 = new TCanvas("mcoll ","mcoll ",600,600); c0=c0;
  	TPad *pad1 =  new TPad("pad1", "12<L1<15",0.0,0.2,1,1.0,21); pad1->SetMargin(0.1,0.1,0.02,0.2);
  	TPad *pad2 =  new TPad("pad2", "ratio",   0.0,0,1,0.2,21); pad2->SetMargin(0.1,0.1,0.3,0.02);
//  	TPad *pad3 =  new TPad("pad3", "diff",    0.0,0,1,0.2,21); pad3->SetMargin(0.1,0.1,0.3,0.02);

  	pad1->SetFillColor(0);pad2->SetFillColor(0);//pad3->SetFillColor(0);
  	pad1->Draw(); pad2->Draw();//pad3->Draw();


  	pad1->cd();

	#ifdef __CINT__
  		gROOT->LoadMacro("AtlasLabels.C");
	#endif

  	leg->AddEntry(hME,"#mue","le");
  	leg->AddEntry(hEM,"e#mu","le");


  	hEM->Draw("e1"); hME->Draw("e1 sames");
  	vline1->Draw(); vline2->Draw();


  	pad2->cd();
  	ratio->Draw();
  	line14->Draw();line15->Draw();line1->Draw();

//  	pad3->cd();
//  	diff->Draw();
//  	line2->Draw();

  	pad1->cd();
  	myText(0.5,0.65,1,"#int L dt = 20.3 fb^{-1}  #sqrt{s} = 8 TeV");
  	myText(0.6,0.35,1,sample);
  	myText(0.6,0.25,1,SR);

  	ATLASLabel(0.62,0.75,"Internal");

  	leg->Draw();


  	c0->Update();

//  	TCanvas* c1 = new TCanvas("mcoll original","mcoll original",600,600); c0=c0;
//  	TPad *pad1_2 =  new TPad("pad1", "12<L1<15",0.0,0.4,1,1.0,21); pad1_2->SetMargin(0.1,0.1,0.02,0.2);
//  	TPad *pad2_2 =  new TPad("pad2", "ratio",   0.0,0.2,1,0.4,21); pad2_2->SetMargin(0.1,0.1,0.02,0.02);
//  	TPad *pad3_2 =  new TPad("pad3", "diff",    0.0,0,1,0.2,21); pad3_2->SetMargin(0.1,0.1,0.3,0.02);
//
//  	pad1_2->SetFillColor(0);pad2_2->SetFillColor(0);pad3_2->SetFillColor(0);
//  	pad1_2->Draw(); pad2_2->Draw();pad3_2->Draw();
//
//
//  	pad1_2->cd();

//  	hEM_orig->Draw("e1"); hME_orig->Draw("e1 sames");
//  	vline1->Draw(); vline2->Draw();


//  	pad2_2->cd();
//  	ratio2->Draw();
//  	line14->Draw();line15->Draw();line1->Draw();

//  	pad3_2->cd();
//  	diff2->Draw();
//  	line2->Draw();

//  	pad1_2->cd();
//  	myText(0.5,0.65,1,"#int L dt = 20.3 fb^{-1}  #sqrt{s} = 8 TeV");
//  	myText(0.6,0.35,1,sample);
//
//  	ATLASLabel(0.62,0.75,"Internal");
//
//  	leg->Draw();
//
//
//  	c1->Update();



  return;
}

#ifndef __CINT__
int main() { 
  BasicExample();
//  gPad->Print("basic.png");
  return 0;
}
#endif
