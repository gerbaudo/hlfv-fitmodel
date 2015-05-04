
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

TH2D* Rebin(TH2D* old){

  int nx=11;
  int ny=36;
  double xbins[]={12,15,20,25,30,35,40,45,50,55,60,1000};
  double ybins[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,400};
  TH2D *h = new TH2D("oldrebin",old->GetTitle(),nx,xbins,ny,ybins);
  TAxis *xaxis = old->GetXaxis();
  TAxis *yaxis = old->GetYaxis();
  for (int j=1;j<=yaxis->GetNbins();j++) {
     for (int i=1;i<=xaxis->GetNbins();i++) {
        h->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),old->GetBinContent(i,j));
     }
  }
return h;
}

TH1D* GetRatio(TH1D* hEM,TH1D* hME,double xmin, double xmax)
{
  TH1D* ratio = (TH1D*)hEM->Clone("ratio");
  ratio->Divide(hME);
  ratio->GetYaxis()->SetTitle("Ratio"); ratio->GetYaxis()->SetTitleSize(0.1);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->GetXaxis()->SetRangeUser(xmin,xmax); ratio->SetLineColor(kBlack);
  ratio->GetYaxis()->SetRangeUser(0,2);
  ratio->SetLineWidth(2); ratio->SetMarkerStyle(8); ratio->SetMarkerSize(0.7);
  ratio->SetMarkerColor(kBlack);
  ratio->GetYaxis()->SetLabelSize(0.1);
  ratio->GetYaxis()->SetNdivisions(5);
  return ratio;
}

void fpT(TString file1, TString file2, TString file3, TString fakefile, TString fakefile2, bool fakes)
{

  SetAtlasStyle();
#ifdef __CINT__
  gROOT->LoadMacro("AtlasUtils.C");
#endif

#ifdef __CINT__
  gROOT->LoadMacro("AtlasLabels.C");
#endif
  gStyle->SetOptStat(0);

  	double xmin = 0;
  	double xmax = 400;

 	TFile* f_0 = new TFile(file1);//"~/workspace/CutFlowOptimization/DataBase/SR_noJets/AllData/SR_noJets.40_12_5.0_0.0_0.0_0.0.root");
  	
  	TFile* f_1 = new TFile(file2);//"~/workspace/CutFlowOptimization/DataBase/SR_with1Jets/AllData/SR_with1Jets.40_12_5.0_0.0_0.0_0.0.root");
  	
	TFile* f_2 = new TFile(file3);//"~/workspace/CutFlowOptimization/DataBase/SR_with2Jets/AllData/SR_with2Jets.40_12_5.0_0.0_0.0_0.0.root");

//	TFile* f_3 = new TFile("~/workspace/CutFlowOptimization/DataBase/SR_All/AllData/SR_All.40_12_5.0_0.0_0.0_0.0.root");

  	
//  	TFile* f_3 = new TFile("~/workspace/CutFlowOptimization/40_40_40/SR_All/AllData/SR_All.40_12_9.9_0.0_0.0_0.0.root");
//  	TFile* f_2 = new TFile("~/workspace/CutFlowOptimization/45_45_45/SR_All/AllData/SR_All.45_12_9.9_0.0_0.0_0.0.root");
//  	TFile* f_1 = new TFile("~/workspace/CutFlowOptimization/50_50_50/SR_All/AllData/SR_All.50_12_9.9_0.0_0.0_0.0.root");
//  	TFile* f_0 = new TFile("~/workspace/CutFlowOptimization/55_55_55/SR_All/AllData/SR_All.55_12_9.9_0.0_0.0_0.0.root");
  	
  	TH1D* hMEl1_0 =(TH1D*)f_0->Get("nom/ME_l1_pt");
  	TH1D* hEMl1_0 =(TH1D*)f_0->Get("nom/EM_l1_pt");
  	TH1D* hMEl1_1 =(TH1D*)f_1->Get("nom/ME_l1_pt");
  	TH1D* hEMl1_1 =(TH1D*)f_1->Get("nom/EM_l1_pt");
  	TH1D* hMEl1_2 =(TH1D*)f_2->Get("nom/ME_l1_pt");
  	TH1D* hEMl1_2 =(TH1D*)f_2->Get("nom/EM_l1_pt");
//  	TH1D* hMEl1_3 =(TH1D*)f_3->Get("nom/ME_l1_pt");
//  	TH1D* hEMl1_3 =(TH1D*)f_3->Get("nom/EM_l1_pt");

        Double_t xbins[25] = {12,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
	double Min=12.; double Max=60.;


		if (fakes){
			TFile* ff = new TFile(fakefile);
			TH1D* h_l1pt_MEfake = (TH1D*)ff->Get("nom/ME_l1_pt");
			TH1D* h_l1pt_EMfake = (TH1D*)ff->Get("nom/EM_l1_pt");

			TFile* ff2 = new TFile(fakefile2);
			TH1D* h_l1pt_MEfake2 = (TH1D*)ff2->Get("nom/ME_l1_pt");
			TH1D* h_l1pt_EMfake2 = (TH1D*)ff2->Get("nom/EM_l1_pt");

			hMEl1_0->Add(h_l1pt_MEfake,-1);
			hEMl1_0->Add(h_l1pt_EMfake,-1);

			hMEl1_1->Add(h_l1pt_MEfake2,-1);
			hEMl1_1->Add(h_l1pt_EMfake2,-1);
		}


  	TH1D* NewhEMl1_0 = (TH1D*) hEMl1_0->Rebin(24,"NewhEMl1_0",xbins);
	TH1D* NewhMEl1_0 = (TH1D*) hMEl1_0->Rebin(24,"NewhMEl1_0",xbins);

 	TH1D* r_0 = GetRatio(NewhEMl1_0,NewhMEl1_0,Min,Max);
 	r_0->SetLineColor(kRed);         r_0->SetMarkerColor(kRed); 

	TH1D* NewhEMl1_1 = (TH1D*) hEMl1_1->Rebin(24,"NewhEMl1_1",xbins);
	TH1D* NewhMEl1_1 = (TH1D*) hMEl1_1->Rebin(24,"NewhMEl1_1",xbins);

 	TH1D* r_1 = GetRatio(NewhEMl1_1,NewhMEl1_1,Min,Max);
 	r_1->SetLineColor(kGreen);         r_1->SetMarkerColor(kGreen);

	TH1D* NewhEMl1_2 = (TH1D*) hEMl1_2->Rebin(24,"NewhEMl1_2",xbins);
	TH1D* NewhMEl1_2 = (TH1D*) hMEl1_2->Rebin(24,"NewhMEl1_2",xbins);

 	TH1D* r_2 = GetRatio(NewhEMl1_2,NewhMEl1_2,Min,Max);
 	r_2->SetLineColor(6);         r_2->SetMarkerColor(6);

//	TH1D* NewhEMl1_3 = (TH1D*) hEMl1_3->Rebin(24,"NewhEMl1_3",xbins);
//	TH1D* NewhMEl1_3 = (TH1D*) hMEl1_3->Rebin(24,"NewhMEl1_3",xbins);

// 	TH1D* r_3 = GetRatio(NewhEMl1_3,NewhMEl1_3,Min,Max);
// 	r_3->SetLineColor(kBlue);         r_3->SetMarkerColor(kBlue);

  	TCanvas* cr = new TCanvas("c","c",50,50,600,600);
 	TPad* thePad = (TPad*)cr->cd();

 	TH1F *h = thePad->DrawFrame(Min,0.0,Max,2.5);
 	h->SetYTitle("e#mu/#mue");
 	h->SetXTitle("l_{1} p_{T} [GeV]");
 	h->GetYaxis()->SetTitleOffset(1.4);
 	h->GetXaxis()->SetTitleOffset(1.4);
 	h->Draw();         
        r_0->Draw("e1 sames");     
	r_1->Draw("e1 sames");     
	r_2->Draw("e1 sames");
//	r_3->Draw("e1 sames");
 	TLine* line = new TLine(Min,1,Max,1); line->SetLineColor(kBlack); line->SetLineStyle(2);
 	line->Draw();
 	myText(0.35,0.75,1,"#int L dt = 20.3 fb^{-1}  #sqrt{s} = 8 TeV"); 
	myText(0.2,  0.2, 1, "Data");
  	ATLASLabel(0.6,0.85,"Internal");
 	TLegend* leg = new TLegend(0.7,0.2,0.9,0.4);
	leg->SetFillColor(kWhite); leg->SetBorderSize(1); leg->SetLineColor(0); leg->SetTextFont(42);
	leg->SetTextSize(.05);
//	leg->AddEntry(r_0,"40-45","le");
//	leg->AddEntry(r_1,"45-50","le");
//	leg->AddEntry(r_2,"50-55","le");
//	leg->AddEntry(r_3,"#geq55","le");
	leg->AddEntry(r_0,"WW","le");
	leg->AddEntry(r_1,"ttbar","le");
	leg->AddEntry(r_2,"Zj","le");
//	leg->AddEntry(r_3,"All jets","le");
        leg->Draw();

  return;
}

#ifndef __CINT__
int main() { 
  BasicExample();
//  gPad->Print("basic.png");
  return 0;
}
#endif
