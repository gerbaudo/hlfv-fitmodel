#include <fstream>

#include "TMath.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "TFile.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TString.h"

#include "RooWorkspace.h"
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

template <typename T>
  string NumberToString ( T Number )
{
     ostringstream ss;
     ss << Number;
     return ss.str();
}




void CreateToyData(TString outfile)
{
	double Fpt[6] = {1.00891, 0.995867, 1.07463, 1.17698, 1.07954,1.0};

	TFile* inputF = TFile::Open("atlasstyle-00-03-05/Data_BGEstimations.root");
	TH1F* hBG0 = (TH1F*)inputF->Get("BG0");
	TH1F* hBG1 = (TH1F*)inputF->Get("BG1");
	TH1F* hBG2 = (TH1F*)inputF->Get("BG2");
	TH1F* hBG3 = (TH1F*)inputF->Get("BG3");
	TH1F* hBG4 = (TH1F*)inputF->Get("BG4");
	TH1F* hBG5 = (TH1F*)inputF->Get("BG5");
	TH1F* hFakeEM0 = (TH1F*)inputF->Get("Mcoll_Fakes_EM_l1pt0_rebin");
	TH1F* hFakeME0 = (TH1F*)inputF->Get("Mcoll_Fakes_ME_l1pt0_rebin");
	TH1F* hFakeEM1 = (TH1F*)inputF->Get("Mcoll_Fakes_EM_l1pt1_rebin");
        TH1F* hFakeME1 = (TH1F*)inputF->Get("Mcoll_Fakes_ME_l1pt1_rebin");
	TH1F* hFakeEM2 = (TH1F*)inputF->Get("Mcoll_Fakes_EM_l1pt2_rebin");
        TH1F* hFakeME2 = (TH1F*)inputF->Get("Mcoll_Fakes_ME_l1pt2_rebin");
	TH1F* hFakeEM3 = (TH1F*)inputF->Get("Mcoll_Fakes_EM_l1pt3_rebin");
        TH1F* hFakeME3 = (TH1F*)inputF->Get("Mcoll_Fakes_ME_l1pt3_rebin");
	TH1F* hFakeEM4 = (TH1F*)inputF->Get("Mcoll_Fakes_EM_l1pt4_rebin");
        TH1F* hFakeME4 = (TH1F*)inputF->Get("Mcoll_Fakes_ME_l1pt4_rebin");
	TH1F* hFakeEM5 = (TH1F*)inputF->Get("Mcoll_Fakes_EM_l1pt5_rebin");
        TH1F* hFakeME5 = (TH1F*)inputF->Get("Mcoll_Fakes_ME_l1pt5_rebin");


	TRandom rand(12);

	TH1F* hEMtoy0 = (TH1F*)hBG0->Clone("hEMtoy0");
	TH1F* hMEtoy0 = (TH1F*)hBG0->Clone("hMEtoy0");
	TH1F* hEMtoy1 = (TH1F*)hBG1->Clone("hEMtoy1");
        TH1F* hMEtoy1 = (TH1F*)hBG1->Clone("hMEtoy1");
	TH1F* hEMtoy2 = (TH1F*)hBG2->Clone("hEMtoy2");
        TH1F* hMEtoy2 = (TH1F*)hBG2->Clone("hMEtoy2");
	TH1F* hEMtoy3 = (TH1F*)hBG3->Clone("hEMtoy3");
        TH1F* hMEtoy3 = (TH1F*)hBG3->Clone("hMEtoy3");
	TH1F* hEMtoy4 = (TH1F*)hBG4->Clone("hEMtoy4");
        TH1F* hMEtoy4 = (TH1F*)hBG4->Clone("hMEtoy4");
	TH1F* hEMtoy5 = (TH1F*)hBG5->Clone("hEMtoy5");
        TH1F* hMEtoy5 = (TH1F*)hBG5->Clone("hMEtoy5");


	//testing
	//double a = 5;
	//double b = 550.3;
	//cout << "random around 5 = " <<rand.Poisson(a) << ", random around 550.3 = " << rand.Poisson(b) << endl;
	//cout << "TAKE 2:" << endl;
	//cout << "random around 5 = " <<rand.Poisson(a) << ", random around 550.3 = " << rand.Poisson(b) << endl;

	
	for (int i=1; i<=hBG0->GetNbinsX(); i++)
	{
		double EMbase_0 = Fpt[0]*hBG0->GetBinContent(i)+hFakeEM0->GetBinContent(i);
		hEMtoy0->SetBinContent(i,rand.Poisson(EMbase_0));
		double MEbase_0 = hBG0->GetBinContent(i)+hFakeME0->GetBinContent(i);
		hMEtoy0->SetBinContent(i,rand.Poisson(MEbase_0));
		double EMbase_1 = Fpt[1]*hBG1->GetBinContent(i)+hFakeEM1->GetBinContent(i);
                hEMtoy1->SetBinContent(i,rand.Poisson(EMbase_1));
                double MEbase_1 = hBG1->GetBinContent(i)+hFakeME1->GetBinContent(i);
                hMEtoy1->SetBinContent(i,rand.Poisson(MEbase_1));
		double EMbase_2 = Fpt[2]*hBG2->GetBinContent(i)+hFakeEM2->GetBinContent(i);
                hEMtoy2->SetBinContent(i,rand.Poisson(EMbase_2));
                double MEbase_2 = hBG2->GetBinContent(i)+hFakeME2->GetBinContent(i);
                hMEtoy2->SetBinContent(i,rand.Poisson(MEbase_2));
		double EMbase_3 = Fpt[3]*hBG3->GetBinContent(i)+hFakeEM3->GetBinContent(i);
                hEMtoy3->SetBinContent(i,rand.Poisson(EMbase_3));
                double MEbase_3 = hBG3->GetBinContent(i)+hFakeME3->GetBinContent(i);
                hMEtoy3->SetBinContent(i,rand.Poisson(MEbase_3));
		double EMbase_4 = Fpt[4]*hBG4->GetBinContent(i)+hFakeEM4->GetBinContent(i);
                hEMtoy4->SetBinContent(i,rand.Poisson(EMbase_4));
                double MEbase_4 = hBG4->GetBinContent(i)+hFakeME4->GetBinContent(i);
                hMEtoy4->SetBinContent(i,rand.Poisson(MEbase_4));
		double EMbase_5 = Fpt[5]*hBG5->GetBinContent(i)+hFakeEM5->GetBinContent(i);
                hEMtoy5->SetBinContent(i,rand.Poisson(EMbase_5));
                double MEbase_5 = hBG5->GetBinContent(i)+hFakeME5->GetBinContent(i);
                hMEtoy5->SetBinContent(i,rand.Poisson(MEbase_5));
		
	}

	TFile* outputF = new TFile(outfile,"RECREATE");
	
	hEMtoy0->Write();
	hMEtoy0->Write();
	hEMtoy1->Write();
        hMEtoy1->Write();
	hEMtoy2->Write();
        hMEtoy2->Write();
	hEMtoy3->Write();
        hMEtoy3->Write();
	hEMtoy4->Write();
        hMEtoy4->Write();
	hEMtoy5->Write();
        hMEtoy5->Write();


}
