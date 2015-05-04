
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

#include <cmath>



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
  	ratio->GetYaxis()->SetTitle("Trigger Efficiency"); //ratio->SetTitleFont(64);
//  	ratio->GetYaxis()->SetTitle("Ratio"); ratio->GetYaxis()->SetTitleSize(0.1);
//  	ratio->GetYaxis()->SetTitleOffset(0.3);
//  	ratio->GetYaxis()->CenterTitle();
//  	ratio->SetTitleSize(0.1);
  	ratio->GetXaxis()->SetTitle("L_{1} p_{T}");
  	ratio->GetXaxis()->SetRangeUser(xmin,xmax); ratio->SetLineColor(kBlack);
  	ratio->GetYaxis()->SetRangeUser(0,1);
  	ratio->SetLineWidth(2); ratio->SetMarkerStyle(8); ratio->SetMarkerSize(0.7);
  	ratio->SetMarkerColor(kBlack);
//  	ratio->GetYaxis()->SetLabelSize(0.1);
//  	ratio->GetYaxis()->SetNdivisions(5);
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


void Plot_turnOn(TString fileAll, TString fileTrig,TString fileAll2, TString fileTrig2,
		TString fileAll3, TString fileTrig3)
{

	SetAtlasStyle();
	#ifdef __CINT__
	  gROOT->LoadMacro("AtlasUtils.C");
	#endif

  	double xmin = 0;
  	double xmax = 100;


  	TFile* fAll = new TFile(fileAll);
  	TFile* fTrig = new TFile(fileTrig);
  	TH1D* hMEAll1 =(TH1D*)fAll->Get("nom/ME_l1_pt");
  	TH1D* hEMAll1 =(TH1D*)fAll->Get("nom/EM_l1_pt");
  	TH1D* hMETrig1 =(TH1D*)fTrig->Get("nom/ME_l1_pt");
  	TH1D* hEMTrig1 =(TH1D*)fTrig->Get("nom/EM_l1_pt");

  	TFile* fAll2 = new TFile(fileAll2);
  	TFile* fTrig2 = new TFile(fileTrig2);
  	TH1D* hMEAll2 =(TH1D*)fAll2->Get("nom/ME_l1_pt");
  	TH1D* hEMAll2 =(TH1D*)fAll2->Get("nom/EM_l1_pt");
  	TH1D* hMETrig2 =(TH1D*)fTrig2->Get("nom/ME_l1_pt");
  	TH1D* hEMTrig2 =(TH1D*)fTrig2->Get("nom/EM_l1_pt");

  	TFile* fAll3 = new TFile(fileAll3);
  	TFile* fTrig3 = new TFile(fileTrig3);
  	TH1D* hMEAll3 =(TH1D*)fAll3->Get("nom/ME_l1_pt");
  	TH1D* hEMAll3 =(TH1D*)fAll3->Get("nom/EM_l1_pt");
  	TH1D* hMETrig3 =(TH1D*)fTrig3->Get("nom/ME_l1_pt");
  	TH1D* hEMTrig3 =(TH1D*)fTrig3->Get("nom/EM_l1_pt");

  	double xbins[]= {0,5,10,15,20,25,30,35,40,50,60,70,100};

  	TH1D* hMEAll = (TH1D*)hMEAll1->Rebin(12,"hMEAllrebin",xbins);
  	TH1D* hEMAll = (TH1D*)hEMAll1->Rebin(12,"hEMAllrebin",xbins);
  	TH1D* hMETrig = (TH1D*)hMETrig1->Rebin(12,"hMETrigrebin",xbins);
  	TH1D* hEMTrig = (TH1D*)hEMTrig1->Rebin(12,"hEMTrigrebin",xbins);
  	TH1D* hMEAll_2 = (TH1D*)hMEAll2->Rebin(12,"hMEAllrebin_2",xbins);
  	TH1D* hEMAll_2 = (TH1D*)hEMAll2->Rebin(12,"hEMAllrebin_2",xbins);
  	TH1D* hMETrig_2 = (TH1D*)hMETrig2->Rebin(12,"hMETrigrebin_2",xbins);
  	TH1D* hEMTrig_2 = (TH1D*)hEMTrig2->Rebin(12,"hEMTrigrebin_2",xbins);
  	TH1D* hMEAll_3 = (TH1D*)hMEAll3->Rebin(12,"hMEAllrebin_3",xbins);
  	TH1D* hEMAll_3 = (TH1D*)hEMAll3->Rebin(12,"hEMAllrebin_3",xbins);
  	TH1D* hMETrig_3 = (TH1D*)hMETrig3->Rebin(12,"hMETrigrebin_3",xbins);
  	TH1D* hEMTrig_3 = (TH1D*)hEMTrig3->Rebin(12,"hEMTrigrebin_3",xbins);

  	TH1D* ratioEM = GetRatio(hEMTrig,hEMAll,xmin,xmax);
  	TH1D* ratioME = GetRatio(hMETrig,hMEAll,xmin,xmax);
  	TH1D* ratioEM2 = GetRatio(hEMTrig_2,hEMAll_2,xmin,xmax);
  	TH1D* ratioME2 = GetRatio(hMETrig_2,hMEAll_2,xmin,xmax);
  	TH1D* ratioEM3 = GetRatio(hEMTrig_3,hEMAll_3,xmin,xmax);
  	TH1D* ratioME3 = GetRatio(hMETrig_3,hMEAll_3,xmin,xmax);

  	for (int i=1; i<=12; i++)
  	{
  		double EM1err = ratioEM->GetBinError(i);
  		double EM1Trig = hEMTrig->GetBinContent(i);
  		double EM1All = hEMAll->GetBinContent(i);
  		ratioEM->SetBinError(i,EM1err*std::sqrt((EM1All-EM1Trig)/(EM1All+EM1Trig)));
  		double EM2err = ratioEM2->GetBinError(i);
  		double EM2Trig = hEMTrig_2->GetBinContent(i);
  		double EM2All = hEMAll_2->GetBinContent(i);
  		ratioEM2->SetBinError(i,EM2err*std::sqrt((EM2All-EM2Trig)/(EM2All+EM2Trig)));
  		double EM3err = ratioEM3->GetBinError(i);
  		double EM3Trig = hEMTrig_3->GetBinContent(i);
  		double EM3All = hEMAll_3->GetBinContent(i);
  		ratioEM3->SetBinError(i,EM3err*std::sqrt((EM3All-EM3Trig)/(EM3All+EM3Trig)));
  		double ME1err = ratioME->GetBinError(i);
  		double ME1Trig = hMETrig->GetBinContent(i);
  		double ME1All = hMEAll->GetBinContent(i);
  		ratioME->SetBinError(i,ME1err*std::sqrt((ME1All-ME1Trig)/(ME1All+ME1Trig)));
  		double ME2err = ratioME2->GetBinError(i);
  		double ME2Trig = hMETrig_2->GetBinContent(i);
  		double ME2All = hMEAll_2->GetBinContent(i);
  		ratioME2->SetBinError(i,ME2err*std::sqrt((ME2All-ME2Trig)/(ME2All+ME2Trig)));
  		double ME3err = ratioME3->GetBinError(i);
  		double ME3Trig = hMETrig_3->GetBinContent(i);
  		double ME3All = hMEAll_3->GetBinContent(i);
  		ratioME3->SetBinError(i,ME3err*std::sqrt((ME3All-ME3Trig)/(ME3All+ME3Trig)));
  	}

  	ratioME->SetMarkerStyle(8); ratioME->SetMarkerSize(0.7);
  	ratioME->SetLineColor(kBlue);ratioME->SetOption("e1");
  	ratioME->SetMarkerColor(kBlue);
  	ratioEM->SetMarkerStyle(8); ratioEM->SetMarkerSize(0.7);
  	ratioEM->SetLineColor(kBlue);ratioEM->SetOption("e1");
  	ratioEM->SetMarkerColor(kBlue);
  	ratioME2->SetMarkerStyle(8); ratioME2->SetMarkerSize(0.7);
  	ratioME2->SetLineColor(kRed);ratioME2->SetOption("e1");
  	ratioME->SetLineStyle(2); ratioME2->SetLineStyle(2); ratioME3->SetLineStyle(2);
  	ratioME2->SetMarkerColor(kRed);
  	ratioME3->SetMarkerStyle(8); ratioME3->SetMarkerSize(0.7);
  	ratioME3->SetLineColor(kGreen);ratioME3->SetOption("e1");
  	ratioME3->SetMarkerColor(kGreen);
  	ratioEM2->SetMarkerStyle(8); ratioEM2->SetMarkerSize(0.7);
  	ratioEM2->SetLineColor(kRed);ratioEM2->SetOption("e1");
  	ratioEM2->SetMarkerColor(kRed);
//  	ratioEM2->SetLineStyle(2);
  	ratioEM3->SetMarkerStyle(8); ratioEM3->SetMarkerSize(0.7);
  	ratioEM3->SetLineColor(kGreen);ratioEM3->SetOption("e1");
  	ratioEM3->SetMarkerColor(kGreen);



  	TLegend* leg = new TLegend(0.6,0.40,0.85,0.65);
  	leg->SetFillColor(kWhite); leg->SetBorderSize(1); leg->SetLineColor(0); leg->SetTextFont(42);
  	leg->SetTextSize(.05);

  	TCanvas* c0 = new TCanvas("turnOn","turnOn",600,600); c0=c0;


	#ifdef __CINT__
  		gROOT->LoadMacro("AtlasLabels.C");
	#endif

  	leg->AddEntry(ratioME,"WW nJets=0 #mue","le");
  	leg->AddEntry(ratioEM,"WW nJets=0 e#mu","le");
  	leg->AddEntry(ratioME2,"WW nJets=1 #mue","le");
  	leg->AddEntry(ratioEM2,"WW nJets=1 e#mu","le");
  	leg->AddEntry(ratioME3,"WW nJets=2 #mue","le");
  	leg->AddEntry(ratioEM3,"WW nJets=2 e#mu","le");

  	ratioEM->Draw(); ratioME->Draw("sames");
  	ratioEM2->Draw("sames"); ratioME2->Draw("sames");
  	ratioEM3->Draw("sames"); ratioME3->Draw("sames");

  	myText(0.3,0.65,1,"#int L dt = 20.3 fb^{-1}  #sqrt{s} = 8 TeV");
  	ATLASLabel(0.55,0.7,"Internal");
  	leg->Draw();

  	c0->Update();

  	TCanvas* c1 = new TCanvas("ratio","ratio",600,600); c1=c1;
  	TH1D* ratioRatio = GetRatio(ratioEM,ratioME,xmin,xmax);
  	TH1D* ratioRatio2 = GetRatio(ratioEM2,ratioME2,xmin,xmax);
  	TH1D* ratioRatio3 = GetRatio(ratioEM3,ratioME3,xmin,xmax);

  	ratioRatio->SetLineColor(kBlue); ratioRatio2->SetLineColor(kRed);
  	ratioRatio3->SetLineColor(kGreen);
  	ratioRatio->GetYaxis()->SetRangeUser(0.8,1.5);
  	ratioRatio->GetYaxis()->SetTitle("e#mu eff. / #mue eff.");

  	TLegend* leg2 = new TLegend(0.6,0.40,0.85,0.65);
  	leg2->SetFillColor(kWhite); leg2->SetBorderSize(1); leg2->SetLineColor(0); leg2->SetTextFont(42);
  	leg2->SetTextSize(.05);

  	leg2->AddEntry(ratioRatio,"nJets = 0","le");
  	leg2->AddEntry(ratioRatio2,"nJets = 1","le");
  	leg2->AddEntry(ratioRatio3,"nJets = 2","le");

  	ratioRatio->Draw("E1");
  	ratioRatio2->Draw("e1 sames");
  	ratioRatio3->Draw("e1 sames");
  	leg2->Draw();
  	myText(0.3,0.65,1,"#int L dt = 20.3 fb^{-1}  #sqrt{s} = 8 TeV");
  	ATLASLabel(0.55,0.7,"Internal");


  return;
}

#ifndef __CINT__
int main() { 
  BasicExample();
//  gPad->Print("basic.png");
  return 0;
}
#endif
