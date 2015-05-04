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







void PlotAll(TString wsname)
{
	char* binLabels[19] = {"60","70","80","90","100","110","120","130","140","150","160","170","180","190","200","250","300","400","1000"};	

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
        
	
	// get initial BG histogram
	//TH1* h_initial_BG_EM = pdf_stateEM->createHistogram("initial_BG_EM",*obsEM);
	//TH1* h_initial_BG_ME = pdf_stateME->createHistogram("initial_BG_ME",*obsME);
	
	double MEnorm = pdf_stateME->expectedEvents(*obsME);
	cout << "EM expected events = " << EMnorm << ", ME expected events = " << MEnorm << "." << endl;
	//h_initial_BG_EM->Scale(EMnorm);
	//h_initial_BG_ME->Scale(MEnorm);	

	// get initial gammas
	int nbins = hdataEM->GetNbinsX();
        double InitGamma[nbins];
        for (int i=0; i<nbins; i++)
        {
               	TString varname = "gamma_B0_l1pt0_bin_"+NumberToString(i);
               	InitGamma[i] = ws->var(varname)->getVal();
               	cout << "initial gamma"+NumberToString(i)+" = " << InitGamma[i] << endl;
        }
        double InitFpt = ws->var("fl1pt_l1pt0")->getVal();
        cout << "initial fpt_l1pt0 = " << InitFpt <<  endl;


	// DO THE GLOBAL FIT
	
	minimize(nll);	
       
	// get final BG histograms
	TH1* h_final_BG_EM = pdf_stateEM->createHistogram("final_BG_EM",*obsEM);
        TH1* h_final_BG_ME = pdf_stateME->createHistogram("final_BG_ME",*obsME); 
	h_final_BG_EM->Scale(EMnorm);
	h_final_BG_ME->Scale(MEnorm);
	
	// uncertainty bands
	TH1D* BuncertaintyEM = new TH1D("BuncertaintyEM","BuncertaintyEM",nbins,0,nbins);
	TH1D* BuncertaintyME = new TH1D("BuncertaintyME","BuncertaintyME",nbins,0,nbins);
	for (int i=1; i<=nbins; i++){
		double sigbEM = h_final_BG_EM->GetBinError(i);
		double bEM = h_final_BG_EM->GetBinContent(i);
		BuncertaintyEM->SetBinError(i,sigbEM); BuncertaintyEM->SetBinContent(i,bEM);
		double sigbME = h_final_BG_ME->GetBinError(i);
                double bME = h_final_BG_ME->GetBinContent(i);
                BuncertaintyME->SetBinError(i,sigbME); BuncertaintyME->SetBinContent(i,bME);
	}
	//BuncertaintyEM->SetFillStyle(3004); 
	BuncertaintyEM->SetFillColor(kGreen-9);
	BuncertaintyEM->SetLineColor(kBlack); BuncertaintyEM->SetLineStyle(2);
	//BuncertaintyME->SetFillStyle(3004); 
	BuncertaintyME->SetFillColor(kBlue-9);
        BuncertaintyME->SetLineColor(kBlack); BuncertaintyME->SetLineStyle(2);

	// get gammas after fit
	double FinalGamma[nbins];
	//TH1* h_initBG_times_gamma = (TH1*)h_initial_BG_EM->Clone("initBGEM_times_gamma");
	for (int i=0; i<nbins; i++)
       	{
               	TString varname = "gamma_B0_l1pt0_bin_"+NumberToString(i);
               	FinalGamma[i] = ws->var(varname)->getVal();
               	cout << "Final gamma in bin "+NumberToString(i)+" = " << FinalGamma[i] << endl;
       	//	h_initBG_times_gamma->SetBinContent(i+1,h_initial_BG_EM->GetBinContent(i+1)*FinalGamma[i]);
	}
	//double FinalFpt = ws->var("fl1pt_l1pt0")->getVal();
	
	// get final alpha (pull)
	RooRealVar* alphaVar = ws->var("alpha_l1ptsys_l1pt0");
	double alpha, alphaErr;
	if (alphaVar != NULL) {
		alpha = ws->var("alpha_l1ptsys_l1pt0")->getVal();
		alphaErr = ws->var("alpha_l1ptsys_l1pt0")->getError();
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
	
	// Draw
	TCanvas* c1 = new TCanvas("BG and Data "+chanName1+" "+chanName2,"BG and Data "+chanName1+" "+chanName2,600,600);
	BuncertaintyEM->Draw("E3 sames"); BuncertaintyME->Draw("E3 sames");
	//h_initial_BG_EM->SetLineColor(kGreen+2); h_initial_BG_EM->SetLineStyle(2); h_initial_BG_EM->Draw("sames");
	hdataEM->SetLineColor(kGreen+2); hdataEM->SetMarkerStyle(20); hdataEM->SetMarkerColor(kGreen+2);
	hdataEM->Draw("e1 sames");
	//h_initial_BG_ME->SetLineColor(kBlue); h_initial_BG_ME->SetLineStyle(2); h_initial_BG_ME->Draw("sames");
        hdataME->SetLineColor(kBlue); hdataME->SetMarkerStyle(20);  hdataME->SetMarkerColor(kBlue);
	hdataME->Draw("e1 sames");

	h_final_BG_EM->SetLineColor(kGreen+2); h_final_BG_EM->SetLineWidth(2); h_final_BG_EM->Draw("sames");
	h_final_BG_ME->SetLineColor(kBlue); h_final_BG_ME->SetLineWidth(2); h_final_BG_ME->Draw("sames");

	TLegend* leg = new TLegend(0.5,0.45,0.85,0.65);
        leg->SetFillColor(kWhite); leg->SetBorderSize(1); leg->SetLineColor(0); //leg->SetTextFont(14);
        leg->SetTextSize(.03);

	leg->AddEntry(hdataME,"DATA #mue","lep");
	leg->AddEntry(hdataEM,"DATA e#mu","lep");
	//leg->AddEntry(h_initial_BG_ME,"Initial #mue PDF","l");
	//leg->AddEntry(h_initial_BG_EM,"Initial e#mu PDF","l");
	leg->AddEntry(h_final_BG_ME,"#mue PDF = #gamma_{i}B_{i} + #muS_{i}","l");
	leg->AddEntry(h_final_BG_EM,"e#mu PDF = f(1+#alpha#sigma)(#gamma_{i}B_{i}+#muW_{i})","l");
	leg->Draw();

	cout << " ********************* Fit Values **************************** " <<  endl;
	if (alphaVar != NULL){cout << "alpha = " << alpha << " +- " << alphaErr << endl;}
	cout << "mu    = " << mu << " +- " << muErr << endl;

	TString WriteDownAlphaValue;
	TString WriteDownMuValue;
	WriteDownAlphaValue = "Fpt0 = ";
	WriteDownMuValue = "#mu = ";
	WriteDownAlphaValue += Form("%4.4f",FinalFpt[0]);
	WriteDownAlphaValue += "#pm";
	WriteDownAlphaValue += Form("%4.4f",FinalFptErr[0]);
	WriteDownMuValue += Form("%4.4f",mu);
        WriteDownMuValue += "#pm";
        WriteDownMuValue += Form("%4.4f",muErr);

	TLatex *texl = new TLatex(12,25,WriteDownAlphaValue);
   	texl->SetTextAlign(22); texl->SetTextSize(0.03); 
   	TLatex *texl2 = new TLatex(12,23,WriteDownMuValue);
        texl2->SetTextAlign(22); texl2->SetTextSize(0.03);
	texl->Draw(); 
	texl2->Draw();



	//FPT 1 ***********************************
	ttype = (RooCatType*)iterat->Next();

        RooAbsPdf  *pdf_stateEM1  = simPdf->getPdf(ttype->GetName()) ;
        RooArgSet  *obstmpEM1  = pdf_stateEM1->getObservables( *mc->GetObservables() ) ;
	RooAbsData *dataEM1 = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));

        RooRealVar *obsEM1     = ((RooRealVar*) obstmpEM1->first());
        TString chanName11(ttype->GetName());	
	TH1* hdataEM1 = dataEM1->createHistogram("Data "+chanName11,*obsEM1);
	for (int ib=0 ; ib<hdataEM1->GetNbinsX()+1 ; ib++) hdataEM1->SetBinError(ib, sqrt(hdataEM1->GetBinContent(ib)));

        double EMnorm1 = pdf_stateEM1->expectedEvents(*obsEM1);
	ttype = (RooCatType*)iterat->Next();
        RooAbsPdf* pdf_stateME1  = simPdf->getPdf(ttype->GetName()) ;
        RooArgSet* obstmpME1  = pdf_stateME1->getObservables( *mc->GetObservables() ) ;
	RooAbsData *dataME1 = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));
        RooRealVar* obsME1 = ((RooRealVar*) obstmpME1->first());
        TString chanName21(ttype->GetName());
	TH1* hdataME1 = dataME1->createHistogram("Data "+chanName21,*obsME1);

	for (int ib=0 ; ib<hdataME1->GetNbinsX()+1 ; ib++) hdataME1->SetBinError(ib, sqrt(hdataME1->GetBinContent(ib)));
	double MEnorm1 = pdf_stateME1->expectedEvents(*obsME1);
	TH1* h_final_BG_EM1 = pdf_stateEM1->createHistogram("final_BG_EM1",*obsEM1);
        TH1* h_final_BG_ME1 = pdf_stateME1->createHistogram("final_BG_ME1",*obsME1);
        h_final_BG_EM1->Scale(EMnorm1);
        h_final_BG_ME1->Scale(MEnorm1);
	TH1D* BuncertaintyEM1 = new TH1D("BuncertaintyEM1","BuncertaintyEM1",nbins,0,nbins);
        TH1D* BuncertaintyME1 = new TH1D("BuncertaintyME1","BuncertaintyME1",nbins,0,nbins);
        for (int i=1; i<=nbins; i++){
                double sigbEM = h_final_BG_EM1->GetBinError(i);
                double bEM = h_final_BG_EM1->GetBinContent(i);
                BuncertaintyEM1->SetBinError(i,sigbEM); BuncertaintyEM1->SetBinContent(i,bEM);
                double sigbME = h_final_BG_ME1->GetBinError(i);
                double bME = h_final_BG_ME1->GetBinContent(i);
                BuncertaintyME1->SetBinError(i,sigbME); BuncertaintyME1->SetBinContent(i,bME);
        }
	BuncertaintyEM1->SetFillColor(kGreen-9);
        BuncertaintyEM1->SetLineColor(kBlack); BuncertaintyEM1->SetLineStyle(2);
	BuncertaintyME1->SetFillColor(kBlue-9);
        BuncertaintyME1->SetLineColor(kBlack); BuncertaintyME1->SetLineStyle(2);
	double FinalGamma1[nbins];
        for (int i=0; i<nbins; i++)
        {
                TString varname = "gamma_B0_l1pt1_bin_"+NumberToString(i);
                FinalGamma1[i] = ws->var(varname)->getVal();
                cout << "Final gamma in bin "+NumberToString(i)+" = " << FinalGamma1[i] << endl;
        }
	TCanvas* c2 = new TCanvas("BG and Data "+chanName11+" "+chanName21,"BG and Data "+chanName11+" "+chanName21,600,600);
        BuncertaintyEM1->Draw("E3 sames"); BuncertaintyME1->Draw("E3 sames");
        hdataEM1->SetLineColor(kGreen+2); hdataEM1->SetMarkerStyle(20); hdataEM1->SetMarkerColor(kGreen+2);
        hdataEM1->Draw("e1 sames");
        hdataME1->SetLineColor(kBlue); hdataME1->SetMarkerStyle(20);  hdataME1->SetMarkerColor(kBlue);
        hdataME1->Draw("e1 sames");

        h_final_BG_EM1->SetLineColor(kGreen+2); h_final_BG_EM1->SetLineWidth(2); h_final_BG_EM1->Draw("sames");
        h_final_BG_ME1->SetLineColor(kBlue); h_final_BG_ME1->SetLineWidth(2); h_final_BG_ME1->Draw("sames");

        leg->Draw();

        cout << " ********************* Fit Values **************************** " <<  endl;
        cout << "mu    = " << mu << " +- " << muErr << endl;
	TString WriteDownAlphaValue1;
        WriteDownAlphaValue1 = "Fpt1 = ";
        WriteDownAlphaValue1 += Form("%4.4f",FinalFpt[1]);
        WriteDownAlphaValue1 += "#pm";
        WriteDownAlphaValue1 += Form("%4.4f",FinalFptErr[1]);

        TLatex *texl11 = new TLatex(12,25,WriteDownAlphaValue1);
        texl11->SetTextAlign(22); texl11->SetTextSize(0.03);
        texl11->Draw(); 
        texl2->Draw();

}
