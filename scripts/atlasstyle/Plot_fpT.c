
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

 TH1D* GetEM(TString filename,TString title,double xmin, double xmax)
{
  	double xbins[]= {0,100,150,200,300,1000};
  	double xbins2[]= {0,100,150,200,300,400};
  	TFile* f = new TFile(filename);
  	TH1D* hEM =(TH1D*)f->Get("nom/EM_McollHiggs_Unblind");
//  	hEM->SetLineColor(kGreen+2); hEM->SetLineWidth(2);
  	TH1D* htemp = (TH1D*)hEM->Rebin(5,"hEMRebinned",xbins);
  	TH1D* hEMrebinned = new TH1D("hEMrebinned","hEMrebinned",5,xbins2);
  	for (int bin=1; bin<6; bin++){
  		hEMrebinned->SetBinContent(bin,htemp->GetBinContent(bin));
  		hEMrebinned->SetBinError(bin,htemp->GetBinError(bin));
  	}
  	hEMrebinned->GetXaxis()->SetRangeUser(xmin,xmax);
  	hEMrebinned->SetMarkerStyle(8); hEMrebinned->SetMarkerSize(0.7);
  	hEMrebinned->SetLineColor(kPink + 8);hEMrebinned->SetOption("e1");
  	hEMrebinned->SetMarkerColor(kPink + 8);
  	hEMrebinned->SetTitle(title); //hEMrebinned->GetXaxis()->SetTitle("M_{coll} (GeV)");
  	hEMrebinned->GetXaxis()->SetLabelOffset(0);
  	hEMrebinned->GetXaxis()->SetLabelSize(0);
    delete htemp;
  	return hEMrebinned;
 }

  TH1D* GetME(TString filename,TString title,double xmin, double xmax)
  {
  	double xbins1[]= {0,100,150,200,300,1000};
  	double xbins2[]= {0,100,150,200,300,400};
  	TFile* f = new TFile(filename);
  	TH1D* hME =(TH1D*)f->Get("nom/ME_McollHiggs_Unblind");
  	hME->SetLineColor(kBlue); hME->SetLineWidth(2);
 	TH1D* htemp = (TH1D*)hME->Rebin(5,"hMERebinned",xbins1);
  	TH1D* hMErebinned = new TH1D("hMERebinned","hMERebinned",5,xbins2);
  	for (int bin=1; bin<6; bin++){
  		hMERebinned->SetBinContent(bin,htemp->GetBinContent(bin));
  		hMERebinned->SetBinError(bin,htemp->GetBinError(bin));
  	}
  	hMErebinned->GetXaxis()->SetRangeUser(xmin,xmax);
  	hMErebinned->SetMarkerStyle(8); hMErebinned->SetMarkerSize(0.7);
  	hMErebinned->SetLineColor(kTeal - 6);hMErebinned->SetOption("e1");
  	hMErebinned->SetMarkerColor(kTeal - 6);
  	hMErebinned->SetTitle(title);
  	//hMErebinned->GetXaxis()->SetTitle("M_{coll} (GeV)");
  	hMErebinned->GetXaxis()->SetLabelOffset(0);
  	hMErebinned->GetXaxis()->SetLabelSize(0);
  	delete htemp;
  	return hMErebinned;
  }

  TH1D* GetEMblind(TString filename,TString title,double xmin, double xmax)
  {
  	double xbins[]= {0,100,150,200,300,1000};
  	double xbins2[]= {0,100,150,200,300,400};
  	TFile* f = new TFile(filename);
  	TH1D* hEM =(TH1D*)f->Get("nom/EM_McollHiggs");
  	hEM->SetLineColor(kGreen+2); hEM->SetLineWidth(2);
  	TH1D* htemp = (TH1D*)hEM->Rebin(5,"hEMRebinned",xbins);
  	TH1D* hEMrebinned = new TH1D("hEMrebinned","hEMrebinned",5,xbins2);
  	for (int bin=1; bin<6; bin++){
  		hEMrebinned->SetBinContent(bin,htemp->GetBinContent(bin));
  		hEMrebinned->SetBinError(bin,htemp->GetBinError(bin));
  	}
  	hEMrebinned->GetXaxis()->SetRangeUser(xmin,xmax);
  	hEMrebinned->SetMarkerStyle(8); hEMrebinned->SetMarkerSize(0.7);
  	hEMrebinned->SetLineColor(kPink + 8);hEMrebinned->SetOption("e1");
  	hEMrebinned->SetMarkerColor(kPink + 8);
  	hEMrebinned->SetTitle(title); hEMrebinned->SetTitleSize(0.2); hEMrebinned->SetTitleOffset(0.4);
  	hEMrebinned->GetXaxis()->SetLabelOffset(0);
  	hEMrebinned->GetXaxis()->SetLabelSize(0);
  	delete htemp;
  	return hEMrebinned;
  }

  TH1D* GetMEblind(TString filename,TString title,double xmin, double xmax)
  {
  	double xbins[]= {0,100,150,200,300,1000};
  	double xbins2[]= {0,100,150,200,300,400};
  	TFile* f = new TFile(filename);
  	TH1D* hME =(TH1D*)f->Get("nom/ME_McollHiggs");
  	TH1D* htemp = (TH1D*)hME->Rebin(5,"hMERebinned",xbins);
  	TH1D* hMErebinned = new TH1D("hMERebinned","hMERebinned",5,xbins2);
  	for (int bin=1; bin<6; bin++){
  		hMERebinned->SetBinContent(bin,htemp->GetBinContent(bin));
  		hMERebinned->SetBinError(bin,htemp->GetBinError(bin));
  	}
  	hMErebinned->GetXaxis()->SetRangeUser(xmin,xmax);
  	hMErebinned->SetMarkerStyle(8); hMErebinned->SetMarkerSize(0.7);
  	hMErebinned->SetLineColor(kTeal - 6);hMErebinned->SetOption("e1");
  	hMErebinned->SetMarkerColor(kTeal - 6);
  	hMErebinned->SetTitle(title); hMErebinned->SetTitleSize(0.02); hMErebinned->SetTitleOffset(0.4);
  	hMErebinned->GetXaxis()->SetLabelOffset(0);
  	hMErebinned->GetXaxis()->SetLabelSize(0);
  	delete htemp;
  	return hMErebinned;
  }

  TH1D* GetRatio(TH1D* hEM,TH1D* hME,double xmin, double xmax)
  {
  	TH1D* ratio = (TH1D*)hEM->Clone("ratio");
  	ratio->Divide(hME);
//  	ratio->GetYaxis()->SetTitle("Ratio"); //ratio->SetTitleFont(64);
  	ratio->GetYaxis()->SetTitle("Ratio"); ratio->GetYaxis()->SetTitleSize(0.1);
  	ratio->GetYaxis()->SetTitleOffset(0.5);
  	ratio->GetYaxis()->CenterTitle();
//  	ratio->SetTitleSize(0.1);
  	ratio->GetXaxis()->SetRangeUser(xmin,xmax); ratio->SetLineColor(kBlack);
  	ratio->GetYaxis()->SetRangeUser(0,2);
  	ratio->SetLineWidth(2); ratio->SetMarkerStyle(8); ratio->SetMarkerSize(0.7);
  	ratio->SetMarkerColor(kBlack);
  	ratio->GetYaxis()->SetLabelSize(0.1);
  	ratio->GetYaxis()->SetNdivisions(5);
  	return ratio;
  }

  TH1D* GetDiff(TH1D* hEM,TH1D* hME,double xmin, double xmax)
  {
  	TH1D* diff = (TH1D*)hEM->Clone("diff");
  	diff->Add(hME,-1);
  	diff->GetYaxis()->SetTitle("Diff"); diff->GetYaxis()->SetTitleSize(0.1);
  	diff->GetYaxis()->SetTitleOffset(0.5);
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


void Plot_fpT(TString dir, TString l0pt))
{

	SetAtlasStyle();
	#ifdef __CINT__
	  gROOT->LoadMacro("AtlasUtils.C");
	#endif

  	double xmin = 0;
  	double xmax = 400;

//  	TString l0pt("45");

  	TH1D* h_EM1 = GetEM(dir+l0pt+"_12_5.0_0.0_0.0_0.0.root", " 12<L1<15",xmin,xmax);
  	TH1D* h_ME1 = GetME(dir+l0pt+"_12_5.0_0.0_0.0_0.0.root"," 12<L1<15",xmin,xmax);
  	TH1D* ratio1 = GetRatio(h_EM1,h_ME1,xmin,xmax);
  	TH1D* diff1 = GetDiff(h_EM1,h_ME1,xmin,xmax);

  	TH1D* h_EM2 = GetEM(dir+l0pt+"_15_5.0_0.0_0.0_0.0.root"," 15<L1<18",xmin,xmax);
  	TH1D* h_ME2 = GetME(dir+l0pt+"_15_5.0_0.0_0.0_0.0.root"," 15<L1<18",xmin,xmax);
  	TH1D* ratio2 = GetRatio(h_EM2,h_ME2,xmin,xmax);
  	TH1D* diff2 = GetDiff(h_EM2,h_ME2,xmin,xmax);

  	TH1D* h_EM3 = GetEM(dir+l0pt+"_18_5.0_0.0_0.0_0.0.root"," 18<L1<20",xmin,xmax);
  	TH1D* h_ME3 = GetME(dir+l0pt+"_18_5.0_0.0_0.0_0.0.root"," 18<L1<20",xmin,xmax);
  	TH1D* ratio3 = GetRatio(h_EM3,h_ME3,xmin,xmax);
  	TH1D* diff3 = GetDiff(h_EM3,h_ME3,xmin,xmax);

  	TH1D* h_EM4 = GetEM(dir+l0pt+"_20_5.0_0.0_0.0_0.0.root"," 20<L1<25",xmin,xmax);
  	TH1D* h_ME4 = GetME(dir+l0pt+"_20_5.0_0.0_0.0_0.0.root"," 20<L1<25",xmin,xmax);
  	TH1D* ratio4 = GetRatio(h_EM4,h_ME4,xmin,xmax);
  	TH1D* diff4 = GetDiff(h_EM4,h_ME4,xmin,xmax);

  	TH1D* h_EM5 = GetEM(dir+l0pt+"_25_5.0_0.0_0.0_0.0.root"," 25<L1",xmin,xmax);
  	TH1D* h_ME5 = GetME(dir+l0pt+"_25_5.0_0.0_0.0_0.0.root"," 25<L1",xmin,xmax);
  	TH1D* ratio5 = GetRatio(h_EM5,h_ME5,xmin,xmax);
  	TH1D* diff5 = GetDiff(h_EM5,h_ME5,xmin,xmax);

	TH1D* h_EM6 = GetEM(dir+l0pt+"_30_5.0_0.0_0.0_0.0.root"," 25<L1",xmin,xmax);
  	TH1D* h_ME6 = GetME(dir+l0pt+"_30_5.0_0.0_0.0_0.0.root"," 25<L1",xmin,xmax);
  	TH1D* ratio6 = GetRatio(h_EM6,h_ME6,xmin,xmax);
  	TH1D* diff6 = GetDiff(h_EM6,h_ME6,xmin,xmax);

//  	TH1D* h_EM12 = GetEM(dir2+l0pt+"_12_5.0_0.0_0.0_0.0.root",name2+" 12<L1<15",xmin,xmax);
//  	TH1D* h_ME12 = GetME(dir2+l0pt+"_12_5.0_0.0_0.0_0.0.root",name2+" 12<L1<15",xmin,xmax);
//  	TH1D* ratio12 = GetRatio(h_EM12,h_ME12,xmin,xmax);
//  	TH1D* diff12 = GetDiff(h_EM12,h_ME12,xmin,xmax);
//
//  	TH1D* h_EM22 = GetEM(dir2+l0pt+"_15_5.0_0.0_0.0_0.0.root",name2+" 15<L1<18",xmin,xmax);
//  	TH1D* h_ME22 = GetME(dir2+l0pt+"_15_5.0_0.0_0.0_0.0.root",name2+" 15<L1<18",xmin,xmax);
//  	TH1D* ratio22 = GetRatio(h_EM22,h_ME22,xmin,xmax);
//  	TH1D* diff22 = GetDiff(h_EM22,h_ME22,xmin,xmax);
//
//  	TH1D* h_EM32 = GetEM(dir2+l0pt+"_18_5.0_0.0_0.0_0.0.root",name2+" 18<L1<20",xmin,xmax);
//  	TH1D* h_ME32 = GetME(dir2+l0pt+"_18_5.0_0.0_0.0_0.0.root",name2+" 18<L1<20",xmin,xmax);
//  	TH1D* ratio32 = GetRatio(h_EM32,h_ME32,xmin,xmax);
//  	TH1D* diff32= GetDiff(h_EM32,h_ME32,xmin,xmax);
//
//  	TH1D* h_EM42 = GetEM(dir2+l0pt+"_20_5.0_0.0_0.0_0.0.root",name2+" 20<L1<25",xmin,xmax);
//  	TH1D* h_ME42 = GetME(dir2+l0pt+"_20_5.0_0.0_0.0_0.0.root",name2+" 20<L1<25",xmin,xmax);
//  	TH1D* ratio42 = GetRatio(h_EM42,h_ME42,xmin,xmax);
//  	TH1D* diff42 = GetDiff(h_EM42,h_ME42,xmin,xmax);
//
//  	TH1D* h_EM52 = GetEM(dir2+l0pt+"_25_5.0_0.0_0.0_0.0.root",name2+" 25<L1",xmin,xmax);
//  	TH1D* h_ME52 = GetME(dir2+l0pt+"_25_5.0_0.0_0.0_0.0.root",name2+" 25<L1",xmin,xmax);
//  	TH1D* ratio52 = GetRatio(h_EM52,h_ME52,xmin,xmax);
//  	TH1D* diff52 = GetDiff(h_EM52,h_ME52,xmin,xmax);

  	TLegend* leg = new TLegend(0.7,0.55,0.85,0.65);
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


  	TCanvas* c0 = new TCanvas("mcoll","mcoll",1800,600); c0=c0;
  	TPad *pad1 =  new TPad("pad1", "12<L1<15",0.0,0.4,0.2,1.0,21); pad1->SetMargin(0.1,0.1,0.02,0.2);
  	TPad *pad12 = new TPad("pad12","15<L1<18",1./6,0.4,2./6,1.0,21); pad12->SetMargin(0.1,0.1,0.02,0.2);
  	TPad *pad13 = new TPad("pad13","18<L1<20",2./6,0.4,3./6,1.0,21); pad13->SetMargin(0.1,0.1,0.02,0.2);
  	TPad *pad14 = new TPad("pad14","20<L1<25",3./6,0.4,4./6,1.0,21); pad14->SetMargin(0.1,0.1,0.02,0.2);
  	TPad *pad15 = new TPad("pad15","25<L1<30", 4./6,0.4,5./6,1.0,21); pad15->SetMargin(0.1,0.1,0.02,0.2);
  	TPad *pad16 = new TPad("pad16","L1>30",   5./6,0.4,1.0,1.0,21); pad16->SetMargin(0.1,0.1,0.02,0.2);
  	TPad *pad2 =  new TPad("pad2", "ratio",   0.0,0.2,1./6,0.4,21); pad2->SetMargin(0.1,0.1,0.02,0.02);
  	TPad *pad22 = new TPad("pad22", "ratio",  1./6,0.2,2./6,0.4,21); pad22->SetMargin(0.1,0.1,0.02,0.02);
  	TPad *pad23 = new TPad("pad23", "ratio",  2./6,0.2,3./6,0.4,21); pad23->SetMargin(0.1,0.1,0.02,0.02);
  	TPad *pad24 = new TPad("pad24","ratio",   3./6,0.2,4./6,0.4,21); pad24->SetMargin(0.1,0.1,0.02,0.02);
  	TPad *pad25 = new TPad("pad25","ratio",   4./6,0.2,5./6,0.4,21); pad25->SetMargin(0.1,0.1,0.02,0.02);
  	TPad *pad26 = new TPad("pad26","ratio",   5./6,0.2,1,0.4,21); pad26->SetMargin(0.1,0.1,0.02,0.02);
  	TPad *pad3 =  new TPad("pad3", "diff",    0.0,0,1./6,0.2,21); pad3->SetMargin(0.1,0.1,0.3,0.02);
  	TPad *pad32 = new TPad("pad32", "diff",   1./6,0,2./6,0.2,21); pad32->SetMargin(0.1,0.1,0.3,0.02);
  	TPad *pad33 = new TPad("pad33", "diff",   2./6,0,3./6,0.2,21); pad33->SetMargin(0.1,0.1,0.3,0.02);
  	TPad *pad34 = new TPad("pad34", "diff",   3./6,0,4./6,0.2,21); pad34->SetMargin(0.1,0.1,0.3,0.02);
  	TPad *pad35 = new TPad("pad35", "diff",   4./6,0,5./6,0.2,21); pad35->SetMargin(0.1,0.1,0.3,0.02);
  	TPad *pad36 = new TPad("pad36", "diff",   5./6,0,1.0,0.2,21); pad36->SetMargin(0.1,0.1,0.3,0.02);
//  	TPad *pad4 =  new TPad("pad4", "12<L1<15",0.0,0.2,0.2,0.5,21); pad4->SetMargin(0.1,0.1,0.02,0.2);
//  	TPad *pad42 = new TPad("pad42","15<L1<18",0.2,0.2,0.4,0.5,21); pad42->SetMargin(0.1,0.1,0.02,0.2);
//  	TPad *pad43 = new TPad("pad43","18<L1<20",0.4,0.2,0.6,0.5,21); pad43->SetMargin(0.1,0.1,0.02,0.2);
//  	TPad *pad44 = new TPad("pad44","20<L1<25",0.6,0.2,0.8,0.5,21); pad44->SetMargin(0.1,0.1,0.02,0.2);
//  	TPad *pad45 = new TPad("pad45","25<L1",   0.8,0.2,1.0,0.5,21); pad45->SetMargin(0.1,0.1,0.02,0.2);
//  	TPad *pad5 =  new TPad("pad5", "ratio",   0.0,0.1,0.2,0.2,21); pad5->SetMargin(0.1,0.1,0.02,0.02);
//  	TPad *pad52 = new TPad("pad52", "ratio",  0.2,0.1,0.4,0.2,21); pad52->SetMargin(0.1,0.1,0.02,0.02);
//  	TPad *pad53 = new TPad("pad53", "ratio",  0.4,0.1,0.6,0.2,21); pad53->SetMargin(0.1,0.1,0.02,0.02);
//  	TPad *pad54 = new TPad("pad54","ratio",   0.6,0.1,0.8,0.2,21); pad54->SetMargin(0.1,0.1,0.02,0.02);
//  	TPad *pad55 = new TPad("pad55","ratio",   0.8,0.1,1.0,0.2,21); pad55->SetMargin(0.1,0.1,0.02,0.02);
//  	TPad *pad6 =  new TPad("pad6", "diff",    0.0,0.0,0.2,0.1,21); pad6->SetMargin(0.1,0.1,0.3,0.02);
//  	TPad *pad62 = new TPad("pad62", "diff",   0.2,0.0,0.4,0.1,21); pad62->SetMargin(0.1,0.1,0.3,0.02);
//  	TPad *pad63 = new TPad("pad63", "diff",   0.4,0.0,0.6,0.1,21); pad63->SetMargin(0.1,0.1,0.3,0.02);
//  	TPad *pad64 = new TPad("pad64", "diff",   0.6,0.0,0.8,0.1,21); pad64->SetMargin(0.1,0.1,0.3,0.02);
//  	TPad *pad65 = new TPad("pad65", "diff",   0.8,0.0,1.0,0.1,21); pad65->SetMargin(0.1,0.1,0.3,0.02);
  	pad1->SetFillColor(0);pad2->SetFillColor(0);pad3->SetFillColor(0);
  	pad1->Draw(); pad2->Draw();pad3->Draw();
  	pad12->SetFillColor(0);pad22->SetFillColor(0);pad32->SetFillColor(0);
  	pad13->SetFillColor(0);pad23->SetFillColor(0);pad33->SetFillColor(0);
  	pad14->SetFillColor(0);pad24->SetFillColor(0);pad34->SetFillColor(0);
  	pad15->SetFillColor(0);pad25->SetFillColor(0);pad35->SetFillColor(0);
  	pad16->SetFillColor(0);pad26->SetFillColor(0);pad36->SetFillColor(0);
  	pad12->Draw(); pad13->Draw(); pad22->Draw(); pad23->Draw();pad32->Draw();pad33->Draw();
  	pad14->Draw(); pad15->Draw(); pad24->Draw(); pad25->Draw(); pad36->Draw(); pad34->Draw(); pad35->Draw(); pad36->Draw();

//  	pad4->SetFillColor(0);pad5->SetFillColor(0);pad6->SetFillColor(0);
//  	pad4->Draw(); pad5->Draw();pad6->Draw();
//  	pad42->SetFillColor(0);pad52->SetFillColor(0);pad62->SetFillColor(0);
//  	pad43->SetFillColor(0);pad53->SetFillColor(0);pad63->SetFillColor(0);
//  	pad44->SetFillColor(0);pad54->SetFillColor(0);pad64->SetFillColor(0);
//  	pad45->SetFillColor(0);pad55->SetFillColor(0);pad65->SetFillColor(0);
//  	pad42->Draw(); pad43->Draw(); pad52->Draw(); pad53->Draw();pad62->Draw();pad63->Draw();
//  	pad44->Draw(); pad45->Draw(); pad54->Draw(); pad55->Draw(); pad64->Draw(); pad65->Draw();


  	pad1->cd();

	#ifdef __CINT__
  		gROOT->LoadMacro("AtlasLabels.C");
	#endif

  	leg->AddEntry(h_ME1,"#mue","le");
  	leg->AddEntry(h_EM1,"e#mu","le");




  	h_EM1->Draw(); h_ME1->Draw("sames");
  	pad12->cd(); h_EM2->Draw(); h_ME2->Draw("sames");
  	pad13->cd(); h_EM3->Draw(); h_ME3->Draw("sames");
  	pad14->cd(); h_EM4->Draw(); h_ME4->Draw("sames");
  	pad15->cd(); h_EM5->Draw(); h_ME5->Draw("sames");
  	pad15->cd(); h_EM6->Draw(); h_ME6->Draw("sames");



  	pad2->cd();
  	ratio1->Draw();
//  	line12->Draw();line13->Draw();
  	line14->Draw();line15->Draw();line1->Draw();

  	pad22->cd();
  	ratio2->Draw();
//  	line12->Draw();line13->Draw();
  	line14->Draw();line15->Draw();line1->Draw();

  	pad23->cd();
  	ratio3->Draw();
//  	line12->Draw();line13->Draw();
  	line14->Draw();line15->Draw();line1->Draw();

  	pad24->cd();
  	ratio4->Draw();
//  	line12->Draw();line13->Draw();
  	line14->Draw();line15->Draw();line1->Draw();

  	pad25->cd();
  	ratio5->Draw("e1");
//  	line12->Draw();line13->Draw();
  	line14->Draw();line15->Draw();line1->Draw();

  	pad26->cd();
  	ratio6->Draw("e1");
  	line14->Draw();line15->Draw();line1->Draw();



  	pad3->cd();
  	diff1->Draw();
//  	line2->Draw();

  	pad32->cd();
  	diff2->Draw();
//  	line2->Draw();

  	pad33->cd();
  	diff3->Draw();
//  	line2->Draw();

  	pad34->cd();
  	diff4->Draw();
//  	line2->Draw();

  	pad35->cd();
  	diff5->Draw();
//  	line2->Draw();

  	pad36->cd();
  	diff6->Draw();

//  	pad4->cd();
//
//  	h_EM12->Draw(); h_ME12->Draw("sames");
//  	pad42->cd(); h_EM22->Draw(); h_ME22->Draw("sames");
//  	pad43->cd(); h_EM32->Draw(); h_ME32->Draw("sames");
//  	pad44->cd(); h_EM42->Draw(); h_ME42->Draw("sames");
//  	pad45->cd(); h_EM52->Draw(); h_ME52->Draw("sames");
//
//  	pad5->cd();
//  	ratio12->Draw();
////  	line12->Draw();line13->Draw();
//  	line14->Draw();line15->Draw();line1->Draw();
//
//  	pad52->cd();
//  	ratio22->Draw();
////  	line12->Draw();line13->Draw();
//  	line14->Draw();line15->Draw();line1->Draw();
//
//  	pad53->cd();
//  	ratio32->Draw();
////  	line12->Draw();line13->Draw();
//  	line14->Draw();line15->Draw();line1->Draw();
//
//  	pad54->cd();
//  	ratio42->Draw();
////  	line12->Draw();line13->Draw();
//  	line14->Draw();line15->Draw();line1->Draw();
//
//  	pad55->cd();
//  	ratio52->Draw();
////  	line12->Draw();line13->Draw();
//  	line14->Draw();line15->Draw();line1->Draw();
//
//  	pad6->cd();
//  	diff12->Draw();
////  	line2->Draw();
//
//  	pad62->cd();
//  	diff22->Draw();
////  	line2->Draw();
//
//  	pad63->cd();
//  	diff32->Draw();
////  	line2->Draw();
//
//  	pad64->cd();
//  	diff42->Draw();
////  	line2->Draw();
//
//  	pad65->cd();
//  	diff52->Draw();
////  	line2->Draw();

  	pad1->cd();
  	myText(0.3,0.65,1,"#int L dt = 20.3 fb^{-1}  #sqrt{s} = 8 TeV");
  	myText(0.4,0.9,1,"12<p^{T}_{l1}<15");
  	ATLASLabel(0.55,0.7,"Internal");
  	leg->Draw();

  	pad12->cd(); leg->Draw(); myText(0.4,0.9,1,"15<p^{T}_{l1}<18");
  	pad13->cd(); leg->Draw(); myText(0.4,0.9,1,"18<p^{T}_{l1}<20");
  	pad14->cd(); leg->Draw(); myText(0.4,0.9,1,"20<p^{T}_{l1}<25");
  	pad15->cd(); leg->Draw(); myText(0.4,0.9,1,"25<p^{T}_{l1},30");
  	pad16->cd(); leg->Draw(); myText(0.4,0.9,1,"30<p^{T}_{l1}");

//  	pad4->cd(); leg->Draw();  myText(0.4,0.9,1,"12<p^{T}_{l1}<15");
//  	pad42->cd(); leg->Draw();  myText(0.4,0.9,1,"15<p^{T}_{l1}<15");
//  	pad43->cd(); leg->Draw();  myText(0.4,0.9,1,"18<p^{T}_{l1}<20");
//  	pad44->cd(); leg->Draw();  myText(0.4,0.9,1,"20<p^{T}_{l1}<25");
//  	pad45->cd(); leg->Draw(); myText(0.4,0.9,1,"25<p^{T}_{l1}");




  	c0->Update();


  return;
}

#ifndef __CINT__
int main() { 
  BasicExample();
//  gPad->Print("basic.png");
  return 0;
}
#endif
