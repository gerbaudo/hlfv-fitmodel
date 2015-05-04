
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

TString SR = "N_{j}=0";
TString L = "#int L dt = 20.3 fb^{-1}  #sqrt{s} = 8 TeV";
TString SAMPLE = "Data VR BTB";
TString SAMPLEVR = "Data VR notBTB";
bool BLIND = 1;
Double_t l1bins[] = {12,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
int l1nbins = 24;

TH2D* Rebin2DHisto(TH2D* old, int ny, Double_t* ybins){
  
  TAxis *xaxis = old->GetXaxis();
  TAxis *yaxis = old->GetYaxis();

//  const Double_t* ybins = yaxis->GetXbins()->GetArray();
//  int ny = yaxis->GetNbins();

  TH2D *h = new TH2D("oldrebin",old->GetTitle(),l1nbins,l1bins,ny,ybins);

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

TH1D* SetEMHists(TH1D* hem)
{
  hem->SetMarkerStyle(8); hem->SetMarkerSize(0.7); hem->SetLineStyle(1);
  hem->SetLineColor(kPink + 8);hem->SetOption("e1");
  hem->SetMarkerColor(kPink + 8);
  hem->GetXaxis()->SetLabelOffset(0); hem->GetXaxis()->SetLabelSize(0);
  hem->GetXaxis()->SetRangeUser(0,400);  hem->GetYaxis()->SetRangeUser(0,2500);
  hem->GetYaxis()->SetTitleOffset(1.0);
  
 return hem;
}

TH1D* SetMEHists(TH1D* hme)
{
  hme->SetMarkerStyle(8); hme->SetMarkerSize(0.7); hme->SetLineStyle(1);
  hme->SetLineColor(kTeal - 6);hme->SetOption("e1");
  hme->SetMarkerColor(kTeal - 6);
  hme->GetXaxis()->SetLabelOffset(0); hme->GetXaxis()->SetLabelSize(0);
  hme->GetXaxis()->SetRangeUser(0,400);  hme->GetYaxis()->SetRangeUser(0,2500);
  hme->GetYaxis()->SetTitleOffset(1.0);

  return hme;
}

TH1D* GetBlindHistos(TH1D* h)
{
  TH1D* h_blind=(TH1D*)h->Clone("h_blind");
  TAxis *xaxis = h->GetXaxis();
  for (int i=1;i<=xaxis->GetNbins();i++) {
    if (100<=h->GetBinCenter(i)&&h->GetBinCenter(i)<=150&&BLIND){h_blind->SetBinContent(i,0);h_blind->SetBinError(i,0);}
  }
  return h_blind;
}

void DrawRatio(TH1D* hem, TH1D* hme, double Min, double Max, TString Xtitle, TString c, bool ValReg)
{

  TH1D* r = GetRatio(hem,hme,Min,Max);
  r->SetLineColor(kBlue);         r->SetMarkerColor(kBlue); 

  TCanvas* cr = new TCanvas(c,c,50,50,600,600);
  TPad* thePad = (TPad*)cr->cd();

  TH1F *h = thePad->DrawFrame(Min,0.0,Max,2.5);
  h->SetYTitle("e#mu/#mue");
  h->SetXTitle(Xtitle);
  h->GetYaxis()->SetTitleOffset(1.4);
  h->GetXaxis()->SetTitleOffset(1.4);
  h->Draw();         r->Draw("sames");     
  TLine* line = new TLine(Min,1,Max,1); line->SetLineColor(kBlack); line->SetLineStyle(2);
  line->Draw();
  myText(0.35,0.75,1,L); 
  myText(0.2,  0.2, 1, SR);
  if(ValReg){myText(0.2,  0.3, 1, SAMPLEVR);} else{myText(0.2,  0.3, 1, SAMPLE);}
  ATLASLabel(0.6,0.85,"Internal");
  return;
}

void DrawHistos(TH1D* hem, TH1D* hme, TH1D* hmeScaled, double Min, double Max, TString Xtitle, TString c)
{

  hem->SetMarkerColor(kBlue); hem->SetLineColor(kBlue); hem->SetMarkerSize(0.5); 
  hme->SetMarkerColor(kRed);  hme->SetLineColor(kRed);  hme->SetMarkerSize(0.5); 
  hmeScaled->SetMarkerColor(kGreen); hmeScaled->SetLineColor(kGreen); hmeScaled->SetMarkerSize(0.5); 

  TCanvas* cr = new TCanvas(c,c,50,50,600,600);
  TPad* thePad = (TPad*)cr->cd();
  TH1F *h = thePad->DrawFrame(Min,0.0,Max,2000.);
  h->SetYTitle("");
  h->SetXTitle(Xtitle);
  h->GetYaxis()->SetTitleOffset(1.4);
  h->GetXaxis()->SetTitleOffset(1.4);
  TLegend* leg = new TLegend(0.6,0.45,0.75,0.55);
  leg->SetFillColor(kWhite); leg->SetBorderSize(1); leg->SetLineColor(0); leg->SetTextFont(42);
  leg->SetTextSize(.05);
  leg->AddEntry(hem,"e#mu","le");
  leg->AddEntry(hme,"#mue","le");
  leg->AddEntry(hmeScaled,"#mue scaled","le");

  h->Draw();         hem->Draw("e1 sames");     hme->Draw("e1 sames");      hmeScaled->Draw("e1 sames");     leg->Draw("sames");
  myText(0.35,0.75,1,L); 
  myText(0.2,  0.2, 1, SR);
  myText(0.2,  0.3, 1, SAMPLE);
  ATLASLabel(0.6,0.85,"Internal");
  ATLASLabel(0.6,0.85,"Internal");
  return;
}

void DrawMcoll(TH1D* hem, TH1D* hme,TString c)
{
  hem=SetEMHists(hem);
  hme=SetMEHists(hme);
//  hem=GetBlindHistos(hem);
//  hme=GetBlindHistos(hme);
  Double_t ymin=0.0;  Double_t ymax=2500.;  Double_t xmin=0.;  Double_t xmax=400.; 
  TH1D* ratio = GetRatio(hem,hme,xmin,xmax);
  TH1D* diff = GetDiff(hem,hme,xmin,xmax);
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
  TLine* vline1 = new TLine(100,0,100,2500); vline1->SetLineStyle(2);
  TLine* vline2 = new TLine(150,0,150,2500); vline2->SetLineStyle(2);

  TCanvas* c0 = new TCanvas(c,c,600,600); c0=c0;
  TPad *pad1 =  new TPad("pad1", "Events",0.0,0.4,1,1.0,21); pad1->SetMargin(0.1,0.1,0.02,0.2);
  TPad *pad2 =  new TPad("pad2", "ratio",   0.0,0.2,1,0.4,21); pad2->SetMargin(0.1,0.1,0.02,0.02);
  TPad *pad3 =  new TPad("pad3", "diff",    0.0,0,1,0.2,21); pad3->SetMargin(0.1,0.1,0.3,0.02);
  pad1->SetFillColor(0);pad2->SetFillColor(0);pad3->SetFillColor(0);
  pad1->Draw(); pad2->Draw();pad3->Draw();

  pad1->cd();

  leg->AddEntry(hme,"#mue","le");
  leg->AddEntry(hem,"e#mu","le");

  hem->Draw("e1"); hme->Draw("e1 sames");
  vline1->Draw(); vline2->Draw();

  pad2->cd();
  ratio->SetLineStyle(1); ratio->Draw(); 
  line1->Draw();

  pad3->cd();
  diff->SetLineStyle(1);  diff->Draw("");
  line2->Draw();

  pad1->cd();
  myText(0.5,0.65,1,L); 
  myText(0.56,  0.2, 1, SR);
  myText(0.56,  0.3, 1, SAMPLE);
  ATLASLabel(0.55,0.7,"Internal");

  leg->Draw();

  c0->Update();
}

void Scale_VR(TString filename, TString VRfilename)
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

        Double_t l0bins[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,400};
        int l0nbins = 36;
        Double_t mcollbins[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,400};
        int mcollnbins = 36;
        Double_t dphibins[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4};
        int dphinbins = 17;
        Double_t METbins[] = {0, 10, 20, 30, 40, 50, 60, 70,  80, 90, 100, 120, 150, 200};
        int METnbins = 13;

// get opt and VR files

  	TFile* f = new TFile(filename);
  	TFile* fVR = new TFile(VRfilename);

// get l1 pT histograms frop opt and VR and rebin them

  	TH1D* hMEl1_temp =(TH1D*)f->Get("nom/ME_l1_pt");
  	TH1D* hEMl1_temp =(TH1D*)f->Get("nom/EM_l1_pt");
  	TH1D* hMEl1VR_temp =(TH1D*)fVR->Get("nom/ME_l1_pt");
  	TH1D* hEMl1VR_temp =(TH1D*)fVR->Get("nom/EM_l1_pt");

        TH1D* hMEl1 = hMEl1_temp->Rebin(l1nbins,"hMEl1",l1bins);
        TH1D* hEMl1 = hEMl1_temp->Rebin(l1nbins,"hEMl1",l1bins);
        TH1D* hMEl1VR = hMEl1VR_temp->Rebin(l1nbins,"hMEl1VR",l1bins);
        TH1D* hEMl1VR = hEMl1VR_temp->Rebin(l1nbins,"hEMl1VR",l1bins);

 // calculate fpT

	TH1D* rfine = GetRatio(hEMl1VR,hMEl1VR,0,400);
	for (int bin = rfine->FindBin(50); bin< rfine->GetXaxis()->GetNbins(); bin++){
		rfine->SetBinContent(bin,1);
		rfine->SetBinError(bin,0);
	}

// get 2D histograms and rebin them

	TH2D* h2EMl1Mcoll_temp =(TH2D*)f->Get("nom/EM_l1pt_vs_Mcoll_Unblind");	
        TH2D* h2MEl1Mcoll_temp =(TH2D*)f->Get("nom/ME_l1pt_vs_Mcoll_Unblind");
	TH2D* h2EMl1l0_temp =(TH2D*)f->Get("nom/EM_l1pt_vs_l0pt");	 	
        TH2D* h2MEl1l0_temp =(TH2D*)f->Get("nom/ME_l1pt_vs_l0pt");		
	TH2D* h2EMl1MET_temp =(TH2D*)f->Get("nom/EM_l1pt_vs_MET");		
        TH2D* h2MEl1MET_temp =(TH2D*)f->Get("nom/ME_l1pt_vs_MET");		
//	TH2D* h2EMl1dPhi_temp =(TH2D*)f->Get("nom/EM_l1pt_vs_dPhil0l1");
//        TH2D* h2MEl1dPhi_temp =(TH2D*)f->Get("nom/ME_l1pt_vs_dPhil0l1");

	TH2D* h2EMl1Mcoll=Rebin2DHisto(h2EMl1Mcoll_temp,mcollnbins,mcollbins);	
	TH2D* h2MEl1Mcoll=Rebin2DHisto(h2MEl1Mcoll_temp,mcollnbins,mcollbins);
	TH2D* h2EMl1l0=Rebin2DHisto(h2EMl1l0_temp,l0nbins,l0bins);
	TH2D* h2MEl1l0=Rebin2DHisto(h2MEl1l0_temp,l0nbins,l0bins);
	TH2D* h2EMl1MET=Rebin2DHisto(h2EMl1MET_temp,METnbins,METbins);
	TH2D* h2MEl1MET=Rebin2DHisto(h2MEl1MET_temp,METnbins,METbins);
//	TH2D* h2EMl1dPhi=Rebin2DHisto(h2EMl1dPhi_temp,dphinbins,dphibins);
//	TH2D* h2MEl1dPhi=Rebin2DHisto(h2MEl1dPhi_temp,dphinbins,dphibins);

// scaled histogram (scaling taken from VR)

	TH2D *h2MEl1McollScaled = (TH2D*)h2MEl1Mcoll->Clone("h2MEScaled");    
      	TH2D *h2MEl1l0Scaled = (TH2D*)h2MEl1l0->Clone("h2MEl1l0Scaled");
      	TH2D *h2MEl1METScaled = (TH2D*)h2MEl1MET->Clone("h2MEl1METScaled");
//      	TH2D *h2MEl1dPhiScaled = (TH2D*)h2MEl1dPhi->Clone("h2MEl1dPhiScaled");

        for (int i=1;i<=h2MEl1Mcoll->GetXaxis()->GetNbins();i++) {
         for (int j=1;j<=h2MEl1Mcoll->GetYaxis()->GetNbins();j++) {
            h2MEl1McollScaled->SetBinContent(i,j,h2MEl1Mcoll->GetBinContent(i,j)*rfine->GetBinContent(i));
          }
         for (int j=1;j<=h2MEl1l0->GetYaxis()->GetNbins();j++) {
            h2MEl1l0Scaled->SetBinContent(i,j,h2MEl1l0->GetBinContent(i,j)*rfine->GetBinContent(i));
          }
         for (int j=1;j<=h2MEl1MET->GetYaxis()->GetNbins();j++) {
            h2MEl1METScaled->SetBinContent(i,j,h2MEl1MET->GetBinContent(i,j)*rfine->GetBinContent(i));
          }
//         for (int j=1;j<=h2MEl1dPhi->GetYaxis()->GetNbins();j++) {
//            h2MEl1dPhiScaled->SetBinContent(i,j,h2MEl1dPhi->GetBinContent(i,j)*rfine->GetBinContent(i));
//          }
        }

// project 2D histos into 1D

 	TH1D* hMEl0 = h2MEl1l0->ProjectionY("ME l0 pT", 1,  h2MEl1l0->GetXaxis()->GetNbins(),"e");
  	TH1D* hEMl0 = h2EMl1l0->ProjectionY("EM l0 pT", 1,  h2MEl1l0->GetXaxis()->GetNbins(),"e");
 	TH1D* hMEMET = h2MEl1MET->ProjectionY("ME MET", 1,  h2MEl1MET->GetXaxis()->GetNbins(),"e");
  	TH1D* hEMMET = h2EMl1MET->ProjectionY("EM MET", 1,  h2MEl1MET->GetXaxis()->GetNbins(),"e");
  	TH1D* hME = h2MEl1Mcoll->ProjectionY("ME Mcoll", 1, h2MEl1Mcoll->GetXaxis()->GetNbins(),"e");
  	TH1D* hEM = h2EMl1Mcoll->ProjectionY("EM Mcoll", 1, h2MEl1Mcoll->GetXaxis()->GetNbins(),"e");
//        TH1D *hMEdPhi = h2MEl1dPhi->ProjectionY("ME dPhi", 1,  h2EMl1dPhi->GetXaxis()->GetNbins(),"e");
//        TH1D *hEMdPhi = h2EMl1dPhi->ProjectionY("EM dPhi", 1,  h2EMl1dPhi->GetXaxis()->GetNbins(),"e");

        TH1D* hMEl1Scaled = h2MEl1McollScaled->ProjectionX("ME l1 pT scaled", 1,  h2MEl1McollScaled->GetYaxis()->GetNbins(),"e");
        TH1D *hMEScaled = h2MEl1McollScaled->ProjectionY("ME Mcoll scaled", 1, h2MEl1McollScaled->GetXaxis()->GetNbins(),"e");
        TH1D *hMEl0Scaled = h2MEl1l0Scaled->ProjectionY("ME l0 pT Scaled", 1,  h2MEl1l0Scaled->GetXaxis()->GetNbins(),"e");
        TH1D *hMEMETScaled = h2MEl1METScaled->ProjectionY("ME MET Scaled", 1,  h2MEl1METScaled->GetXaxis()->GetNbins(),"e");
//        TH1D *hMEdPhiScaled = h2MEl1dPhiScaled->ProjectionY("ME dPhi Scaled", 1,  h2MEl1dPhiScaled->GetXaxis()->GetNbins(),"e");

// plots

//        DrawRatio(hEMl1VR,hMEl1VR,12.,100.,"l_{1} p_{T} [GeV]","l1 VR",1);
        TCanvas* crfine = new TCanvas("l1 VR","l1 VR",600,600);
        rfine->Draw("e1");
        DrawRatio(hEMl1,hMEl1,12.,100.,"l_{1} p_{T} [GeV]","l1",0);
        DrawRatio(hEMl1,hMEl1Scaled,12.,100.,"l_{1} p_{T} [GeV]","l1 scaled",0);

        DrawMcoll(hEM,hME,"Mcoll");
        DrawMcoll(hEM,hMEScaled,"McollScaled");

        DrawRatio(hEMl0,hMEl0,45.,200.,"l_{0} p_{T} [GeV]","l0",0);
        DrawRatio(hEMl0,hMEl0Scaled,45.,200.,"l_{0} p_{T} [GeV]","l0 scaled",0);        

        DrawRatio(hEMMET,hMEMET,0.,200.,"MET [GeV]","MET",0);
        DrawRatio(hEMMET,hMEMETScaled,0.,200.,"MET [GeV]","MET scaled",0);        
 
//        DrawRatio(hEMdPhi,hMEdPhi,0.,3.3,"#Delta#phi","dPhi",0);
//        DrawRatio(hEMdPhi,hMEdPhiScaled,0.,3.3,"#Delta#phi","dPhi scaled",0);
  
  return;
}

#ifndef __CINT__
int main() { 
  BasicExample();
//  gPad->Print("basic.png");
  return 0;
}
#endif
