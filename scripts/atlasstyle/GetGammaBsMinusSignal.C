#include "AtlasUtils.h"
#include "AtlasLabels.h"
#include "AtlasStyle.h"

#include <fstream>

#include "TMath.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "TFile.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TROOT.h"
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


//  SetAtlasStyle();
//#ifdef __CINT__
//  gROOT->LoadMacro("AtlasUtils.C");
//#endif

//char* binLabels[19] = {"60","70","80","90","100","110","120","130","140","150","160","170","180","190","200","250","300","400","1000"};

char* binLabels[33] = {"0","60","65","70","75", "80", "85", "90", "95", "100", "105", "110", "115", "120", "125", "130", "135", "140", "145", "150", "155", "160", "165", "170", "180", "200","220", "240","260", "300", "350", "400", "450"};

//char* binLabels[57] = {"55","60","65", "70","75", "80", "85", "90", "95", "100", "105", "110", "115", "120", "125", "130", "135", "140", "145", "150", "155", "160", "165", "170", "175", "180","185", "190","195", "200","205", "210","215", "220","225", "230","235", "240","245", "250","255","260","265","270","275","280","285","290","295", "300","310","320","330","340", "350", "400", "450"};


void makeBlind(TH1* inhist, TH1* outhist)
{
        for (int i=1; i<10; i++){
                outhist->SetBinContent(i,inhist->GetBinContent(i));
                outhist->SetBinError(i,inhist->GetBinError(i));
        }
        for (int i=20; i<=outhist->GetNbinsX(); i++){
                outhist->SetBinContent(i,inhist->GetBinContent(i));
                outhist->SetBinError(i,inhist->GetBinError(i));
        }
}

void SetLabels(TH1* hist)
{
 	for (int k=1; k<=hist->GetNbinsX(); k++){
        	hist->GetXaxis()->SetBinLabel(k,binLabels[k-1]);
        	hist->GetXaxis()->SetLabelSize(0.03);
        	hist->GetXaxis()->SetTitle("M_{coll} (GeV)"); hist->GetXaxis()->SetTitleSize(0.04);
        	hist->GetYaxis()->SetTitle("Events"); hist->GetYaxis()->SetTitleSize(0.04);
 	}
}

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


int minimize(RooAbsReal* fcn)
{
  static int nrItr = 0;
int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
  int save_strat = strat;
  RooMinimizer minim(*fcn);
  minim.setStrategy(strat);
  minim.setPrintLevel(printLevel);


  int status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());


//up the strategy
if (status != 0 && status != 1 && strat < 2)
  {
    strat++;
    cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
    minim.setStrategy(strat);
    status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  }

  if (status != 0 && status != 1 && strat < 2)
  {
    strat++;
    cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
    minim.setStrategy(strat);
    status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  }
 if (status != 0 && status != 1)
  {
    string minType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
    string newMinType;
    if (minType == "Minuit2") newMinType = "Minuit";
    else newMinType = "Minuit2";

    cout << "Switching minuit type from " << minType << " to " << newMinType << endl;

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(newMinType.c_str());
    strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
    minim.setStrategy(strat);

    status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());


    if (status != 0 && status != 1 && strat < 2)
    {
      strat++;
      cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
      minim.setStrategy(strat);
      status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }

    if (status != 0 && status != 1 && strat < 2)
    {
      strat++;
      cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
      minim.setStrategy(strat);
      status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minType.c_str());
  }

  if (status != 0 && status != 1)
  {
    nrItr++;
    if (nrItr > 3)
    {
      nrItr = 0;
      cout << "WARNING::Fit failure unresolved with status " << status << endl;
      return status;
    }
}

  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(msglevel);
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(save_strat);


  if (nrItr != 0) cout << "Successful fit" << endl;
  nrItr=0;
  return status;
}







void GetGammaBsMinusSignal(TString wsname,TString inputfile,TString outputfile)
{

	//SetAtlasStyle();
        //#ifdef __CINT__
        //  gROOT->LoadMacro("AtlasUtils.C");
        //#endif

	//get input histograms

	TFile* infile = new TFile(inputfile);
	TH1D* h_B0 = (TH1D*)infile->Get("Base_Bkg_l1pt0");
	TH1D* h_B1 = (TH1D*)infile->Get("Base_Bkg_l1pt1");
	TH1D* h_B2 = (TH1D*)infile->Get("Base_Bkg_l1pt2");
	TH1D* h_B3 = (TH1D*)infile->Get("Base_Bkg_l1pt3");
	TH1D* h_B4 = (TH1D*)infile->Get("Base_Bkg_l1pt4");
	TH1D* h_B5 = (TH1D*)infile->Get("Base_Bkg_l1pt5");
	TH1D* h_S0 = (TH1D*)infile->Get("Mcoll_signal_ME_l1pt0_rebin");
	TH1D* h_S1 = (TH1D*)infile->Get("Mcoll_signal_ME_l1pt1_rebin");
	TH1D* h_S2 = (TH1D*)infile->Get("Mcoll_signal_ME_l1pt2_rebin");
	TH1D* h_S3 = (TH1D*)infile->Get("Mcoll_signal_ME_l1pt3_rebin");
	TH1D* h_S4 = (TH1D*)infile->Get("Mcoll_signal_ME_l1pt4_rebin");
	TH1D* h_S5 = (TH1D*)infile->Get("Mcoll_signal_ME_l1pt5_rebin");
	TH1D* h_Scont0 = (TH1D*)infile->Get("Mcoll_wrong_signal_ME_l1pt0_rebin");
	TH1D* h_Scont1 = (TH1D*)infile->Get("Mcoll_wrong_signal_ME_l1pt1_rebin");
	TH1D* h_Scont2 = (TH1D*)infile->Get("Mcoll_wrong_signal_ME_l1pt2_rebin");
	TH1D* h_Scont3 = (TH1D*)infile->Get("Mcoll_wrong_signal_ME_l1pt3_rebin");
	TH1D* h_Scont4 = (TH1D*)infile->Get("Mcoll_wrong_signal_ME_l1pt4_rebin");
	TH1D* h_Scont5 = (TH1D*)infile->Get("Mcoll_wrong_signal_ME_l1pt5_rebin");
	TH1D* h_Fake0_ME = (TH1D*)infile->Get("Mcoll_Fakes_ME_l1pt0_rebin");
	TH1D* h_Fake0_EM = (TH1D*)infile->Get("Mcoll_Fakes_EM_l1pt0_rebin");
	TH1D* h_Fake1_ME = (TH1D*)infile->Get("Mcoll_Fakes_ME_l1pt1_rebin");
        TH1D* h_Fake1_EM = (TH1D*)infile->Get("Mcoll_Fakes_EM_l1pt1_rebin");
	TH1D* h_Fake2_ME = (TH1D*)infile->Get("Mcoll_Fakes_ME_l1pt2_rebin");
        TH1D* h_Fake2_EM = (TH1D*)infile->Get("Mcoll_Fakes_EM_l1pt2_rebin");
	TH1D* h_Fake3_ME = (TH1D*)infile->Get("Mcoll_Fakes_ME_l1pt3_rebin");
        TH1D* h_Fake3_EM = (TH1D*)infile->Get("Mcoll_Fakes_EM_l1pt3_rebin");
	TH1D* h_Fake4_ME = (TH1D*)infile->Get("Mcoll_Fakes_ME_l1pt4_rebin");
        TH1D* h_Fake4_EM = (TH1D*)infile->Get("Mcoll_Fakes_EM_l1pt4_rebin");
	TH1D* h_Fake5_ME = (TH1D*)infile->Get("Mcoll_Fakes_ME_l1pt5_rebin");
        TH1D* h_Fake5_EM = (TH1D*)infile->Get("Mcoll_Fakes_EM_l1pt5_rebin");




	//get the stuff from the workspace:
	
	TFile* file=TFile::Open(wsname);
	RooWorkspace* ws = (RooWorkspace*)file->Get("combined");
	ModelConfig  *mc = (ModelConfig*)ws->obj("ModelConfig");
	RooAbsData   *data = ws->data("obsData");
	RooSimultaneous* simPdf=(RooSimultaneous*)(mc->GetPdf());
	RooAbsReal* nll=simPdf->createNLL(*data);

	// FPT 0 **************************************	
	// EM channel
	
	RooCategory* chanCat = (RooCategory*) (&simPdf->indexCat());
        TIterator* iterat = chanCat->typeIterator() ;
        RooCatType* ttype = (RooCatType*)iterat->Next();

	RooAbsPdf  *pdf_stateEM  = simPdf->getPdf(ttype->GetName()) ;
	RooArgSet  *obstmpEM  = pdf_stateEM->getObservables( *mc->GetObservables() ) ;
	
	// get EM data
       	RooAbsData *dataEM = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));
		
	RooRealVar *obsEM     = ((RooRealVar*) obstmpEM->first());
	TString chanName1(ttype->GetName());

	// create data histogram
	TH1* hdataEM = dataEM->createHistogram("Data "+chanName1,*obsEM);
	// set errors to gaussian
        for (int ib=0 ; ib<hdataEM->GetNbinsX()+1 ; ib++) hdataEM->SetBinError(ib, sqrt(hdataEM->GetBinContent(ib)));

	double EMnorm = pdf_stateEM->expectedEvents(*obsEM);
	
	//****************************
	// ME channel
	ttype = (RooCatType*)iterat->Next();
	RooAbsPdf* pdf_stateME  = simPdf->getPdf(ttype->GetName()) ;
        RooArgSet* obstmpME  = pdf_stateME->getObservables( *mc->GetObservables() ) ;

	// get ME data
	RooAbsData *dataME = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));	
	RooRealVar* obsME = ((RooRealVar*) obstmpME->first());
	TString chanName2(ttype->GetName());

        // create data histogram
        TH1* hdataME = dataME->createHistogram("Data "+chanName2,*obsME);
        // set errors to gaussian
        for (int ib=0 ; ib<hdataME->GetNbinsX()+1 ; ib++) hdataME->SetBinError(ib, sqrt(hdataME->GetBinContent(ib)));
        
	
	double MEnorm = pdf_stateME->expectedEvents(*obsME);
	cout << "EM expected events = " << EMnorm << ", ME expected events = " << MEnorm << "." << endl;

	int nbins = hdataEM->GetNbinsX();

	// DO THE GLOBAL FIT
	
	minimize(nll);	
       
	// get gammas after fit
	double FinalGamma[nbins];
	TH1D* h_gammaB0 = (TH1D*)h_B0->Clone("gammaB0");
	for (int i=0; i<nbins; i++)
       	{
               	TString varname = "gamma_B0_l1pt0_bin_"+NumberToString(i);
               	FinalGamma[i] = ws->var(varname)->getVal();
               	cout << "Final gamma in bin "+NumberToString(i)+" = " << FinalGamma[i] << endl;
		h_gammaB0->SetBinContent(i+1,h_B0->GetBinContent(i+1)*FinalGamma[i]);
	}
	

	//FOR UNCONSTRAINED FPT - get final fpts
	double FinalFpt[5];
	double FinalFptErr[5];
	for (int k=0; k<5; k++){
		TString varname = "fl1pt_l1pt"+NumberToString(k);
		FinalFpt[k] = ws->var(varname)->getVal();
		FinalFptErr[k] =  ws->var(varname)->getError();
		cout << varname << " = "  << FinalFpt[k] << " +- " << FinalFptErr[k] << endl;
	}
	
	// get POI value
	double mu = ws->var("mu_BR_htm")->getVal();
	double muErr = ws->var("mu_BR_htm")->getError();

	TH1D* h_BG0 = (TH1D*)h_B0->Clone("h_BG0");
	TH1D* h_PositiveHalfMuSig0 = (TH1D*)h_S0->Clone("h_PosHalfMuSig0");
	h_PositiveHalfMuSig0->Scale(0.5*abs(mu));
	h_BG0->Add(h_PositiveHalfMuSig0,-1);


	TH1D* h_BG1 = (TH1D*)h_B1->Clone("h_BG1");
        TH1D* h_PositiveHalfMuSig1 = (TH1D*)h_S1->Clone("h_PosHalfMuSig1");
        h_PositiveHalfMuSig1->Scale(0.5*abs(mu));
        h_BG1->Add(h_PositiveHalfMuSig1,-1);
	
	TH1D* h_BG2 = (TH1D*)h_B2->Clone("h_BG2");
        TH1D* h_PositiveHalfMuSig2 = (TH1D*)h_S2->Clone("h_PosHalfMuSig2");
        h_PositiveHalfMuSig2->Scale(0.5*abs(mu));
        h_BG2->Add(h_PositiveHalfMuSig2,-1);

	TH1D* h_BG3 = (TH1D*)h_B3->Clone("h_BG3");
        TH1D* h_PositiveHalfMuSig3 = (TH1D*)h_S3->Clone("h_PosHalfMuSig3");
        h_PositiveHalfMuSig3->Scale(0.5*abs(mu));
        h_BG3->Add(h_PositiveHalfMuSig3,-1);

	TH1D* h_BG4 = (TH1D*)h_B4->Clone("h_BG4");
        TH1D* h_PositiveHalfMuSig4 = (TH1D*)h_S4->Clone("h_PosHalfMuSig4");
        h_PositiveHalfMuSig4->Scale(0.5*abs(mu));
        h_BG4->Add(h_PositiveHalfMuSig4,-1);

	TH1D* h_BG5 = (TH1D*)h_B5->Clone("h_BG5");
        TH1D* h_PositiveHalfMuSig5 = (TH1D*)h_S5->Clone("h_PosHalfMuSig5");
        h_PositiveHalfMuSig5->Scale(0.5*abs(mu));
        h_BG5->Add(h_PositiveHalfMuSig5,-1);
	
	TFile* outf = new TFile(outputfile,"RECREATE");

	h_BG0->Write();
	h_BG1->Write();
	h_BG2->Write();
	h_BG3->Write();
	h_BG4->Write();
	h_BG5->Write();
	h_Fake0_ME->Write();
	h_Fake0_EM->Write();
	h_Fake1_ME->Write();
        h_Fake1_EM->Write();
	h_Fake2_ME->Write();
        h_Fake2_EM->Write();
	h_Fake3_ME->Write();
        h_Fake3_EM->Write();
	h_Fake4_ME->Write();
        h_Fake4_EM->Write();
	h_Fake5_ME->Write();
        h_Fake5_EM->Write();

	double B = h_BG0->Integral()+h_BG1->Integral()+h_BG2->Integral()+h_BG3->Integral()+h_BG4->Integral()+h_BG5->Integral();
	
	cout << "B = " << B << endl;

}
