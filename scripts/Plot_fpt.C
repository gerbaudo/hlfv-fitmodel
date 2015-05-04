#include <fstream>

#include "TMath.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "TFile.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TString.h"

#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooNLLVar.h"
#include "RooStats/ModelConfig.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "Math/MinimizerOptions.h"
#include "TStopwatch.h"
#include "RooMinimizerFcn.h"
#include "RooMinimizer.h"
#include "RooCategory.h"
#include "RooRealSumPdf.h"
#include "RooSimultaneous.h"
#include "RooProduct.h"
#include "TParameter.h"


#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

using namespace RooStats;
using namespace RooFit;


template <typename T>
  string NumberToString ( T Number )
{
     ostringstream ss;
     ss << Number;
     return ss.str();
}

template <typename T>
  char NumberToChar ( T Number)
{
     ostringstream ss;
     ss << Number;
     string st = ss.str();
     return st.c_str();

}

double FindError(double a,double da, double b, double db){
            double e=10000000;
            if(b>0.01&&a>0.01){
                        e=sqrt((a*a*db*db+b*b*da*da)/(b*b*b*b));
            }
            return e;
}


TH1D* GetRatio(TH1D* hEM,TH1D* hME)
{
  TH1D* ratio = (TH1D*)hEM->Clone("ratio");
  ratio->Divide(hME);
  ratio->GetYaxis()->SetTitle("Ratio"); ratio->GetYaxis()->SetTitleSize(0.05);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->CenterTitle();
  ratio->SetLineColor(kBlack);
  ratio->GetYaxis()->SetRangeUser(0,2);
  ratio->SetLineWidth(2); ratio->SetMarkerStyle(8); ratio->SetMarkerSize(0.7);
  ratio->SetMarkerColor(kBlack);
  ratio->GetYaxis()->SetLabelSize(0.05);
  ratio->GetYaxis()->SetNdivisions(5);
 
  for(int j=1; j<=ratio->GetNbinsX(); j++){
            ratio->SetBinError(j,FindError(hEM->GetBinContent(j),hEM->GetBinError(j),hME->GetBinContent(j),hME->GetBinError(j)));
  }
  return ratio;
}
 

void Plot_fpt()
{
	TFile *fz1 = new TFile("../../SRHistos/zjets_NOM_sr_emu_os.root");
	TFile *fz2 = new TFile("../../SRHistos/zjets_NOM_sr_mue_os.root");
	TH1D* hz_EM = (TH1D*)fz1->Get("h_pt1_zjets_NOM_sr_emu_os");
	TH1D* hz_ME = (TH1D*)fz2->Get("h_pt1_zjets_NOM_sr_mue_os");	
	
	TFile *ft1 = new TFile("../../SRHistos/top_NOM_sr_emu_os.root");
        TFile *ft2 = new TFile("../../SRHistos/top_NOM_sr_mue_os.root");
        TH1D* ht_EM = (TH1D*)ft1->Get("h_pt1_top_NOM_sr_emu_os");
        TH1D* ht_ME = (TH1D*)ft2->Get("h_pt1_top_NOM_sr_mue_os");
	
	TFile *fh1 = new TFile("../../SRHistos/higgs_NOM_sr_emu_os.root");
        TFile *fh2 = new TFile("../../SRHistos/higgs_NOM_sr_mue_os.root");
        TH1D* hh_EM = (TH1D*)fh1->Get("h_pt1_higgs_NOM_sr_emu_os");
        TH1D* hh_ME = (TH1D*)fh2->Get("h_pt1_higgs_NOM_sr_mue_os");

	TFile *fd1 = new TFile("../../SRHistos/diboson_NOM_sr_emu_os.root");
        TFile *fd2 = new TFile("../../SRHistos/diboson_NOM_sr_mue_os.root");
        TH1D* hd_EM = (TH1D*)fd1->Get("h_pt1_diboson_NOM_sr_emu_os");
        TH1D* hd_ME = (TH1D*)fd2->Get("h_pt1_diboson_NOM_sr_mue_os");

	TH1D* hEM = (TH1D*)hz_EM->Clone("hEM_sum");
	hEM->Add(ht_EM);
	hEM->Add(hh_EM);
	hEM->Add(hd_EM);
	TH1D* hME = (TH1D*)hz_ME->Clone("hME_sum");
        hME->Add(ht_ME);
        hME->Add(hh_ME);
        hME->Add(hd_ME);

	TH1D* fpt = GetRatio(hEM,hME);

	// Draw
	TCanvas* c1 = new TCanvas("l1pT em/me","l1pT em/me",600,600);
	fpt->SetLineWidth(2);
	fpt->Draw();
	
	TFile *fs1 = new TFile("../../SRHistos/signaltaumu_NOM_sr_mue_os.root");
        TFile *fs2 = new TFile("../../SRHistos/signaltaumu_NOM_sr_emu_os.root");
        TH1D* hs_EM = (TH1D*)fs2->Get("h_pt1_signaltaumu_NOM_sr_emu_os");
        TH1D* hs_ME = (TH1D*)fs1->Get("h_pt1_signaltaumu_NOM_sr_mue_os");

	TH1D* hMEplusSig = (TH1D*)hME->Clone("hME_sum_plus_signal");
	hMEplusSig->Add(hs_ME);

	TH1D* hEMplusSig = (TH1D*)hEM->Clone("hEM_sum_plus_wrong_signal");
        hEMplusSig->Add(hs_EM);

	TH1D* fpt_Sig1 = GetRatio(hEMplusSig,hMEplusSig);

	fpt_Sig1->SetLineColor(kGreen+2); fpt_Sig1->SetLineWidth(2);
	fpt_Sig1->Draw("sames");

	TH1D* hMEplusSig2 = (TH1D*)hME->Clone("hME_sum_plus_2signal");
	hMEplusSig2->Add(hs_ME,2);

	TH1D* hEMplusSig2 = (TH1D*)hEM->Clone("hEM_sum_plus_2signal");
        hEMplusSig2->Add(hs_EM,2);

        TH1D* fpt_Sig2 = GetRatio(hEMplusSig2,hMEplusSig2);

        fpt_Sig2->SetLineColor(kMagenta);
        fpt_Sig2->Draw("sames");

	TLine* line1 = new TLine(0,1,100,1); line1->SetLineColor(kRed); line1->SetLineStyle(2);

	line1->Draw();	

	TCanvas* c2 = new TCanvas("c2","c2",600,600);
        hs_ME->Draw();

	TCanvas* c3 = new TCanvas("c3","c3",600,600);	
	hEM->Draw();
	hME->Draw("sames");
	hMEplusSig->Draw("sames");

}
