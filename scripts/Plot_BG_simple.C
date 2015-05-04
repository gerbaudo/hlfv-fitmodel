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


template <typename T>
  string NumberToString ( T Number )
{
     ostringstream ss;
     ss << Number;
     return ss.str();
}




void Plot_BG_simple(TString wsname)
{
	double Fpt0Sigma = 0.13433;
	//get the stuff from the workspace:
	
	TFile* file=TFile::Open(wsname);
	RooWorkspace* ws = (RooWorkspace*)file->Get("combined");
	ModelConfig  *mc = (ModelConfig*)ws->obj("ModelConfig");
	RooAbsData   *data = ws->data("obsData");
	RooSimultaneous* simPdf=(RooSimultaneous*)(mc->GetPdf());
	RooAbsReal* nll=simPdf->createNLL(*data);
	
	//run on channels
	
	RooCategory* chanCat = (RooCategory*) (&simPdf->indexCat());
        TIterator* iterat = chanCat->typeIterator() ;
        RooCatType* ttype;
	bool stop = kFALSE;
	iterat->Next();
	while ((ttype = (RooCatType*) iterat->Next())&&!stop)
	{
		// run on one channel or all	
		stop = kTRUE;
		
		RooAbsPdf  *pdf_state  = simPdf->getPdf(ttype->GetName()) ;
		RooArgSet  *obstmp  = pdf_state->getObservables( *mc->GetObservables() ) ;
	
		// get data corresponding to chanel
        	RooAbsData *datatmp = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));
		
		RooRealVar *obs     = ((RooRealVar*) obstmp->first());
		TString chanName(ttype->GetName());

		// create data histogram
		TH1* hdata = datatmp->createHistogram("Data "+chanName,*obs);
		// set errors to gaussian
        	for (int ib=0 ; ib<hdata->GetNbinsX()+1 ; ib++) hdata->SetBinError(ib, sqrt(hdata->GetBinContent(ib)));
		
		// get initial BG histogram
		TH1* h_initial_BG = pdf_state->createHistogram("initial_BG_"+chanName,*obs);
		double BGnumEvents = pdf_state->expectedEvents(*obstmp);
		h_initial_BG->Scale(BGnumEvents);
			
		// get initial gammas
		int nbins = h_initial_BG->GetNbinsX();
        	double InitGamma[nbins];
        	for (int i=0; i<nbins; i++)
        	{
                	TString varname = "gamma_B0_0j_l1pt0_bin_"+NumberToString(i);
                	InitGamma[i] = ws->var(varname)->getVal();
                	cout << "initial gamma"+NumberToString(i)+" = " << InitGamma[i] << endl;
        	}
        	double InitFpt = ws->var("fl1pt_l1pt0")->getVal();
        	cout << "initial fpt_l1pt0 = " << InitFpt <<  endl;

		TCanvas* c1 = new TCanvas("BG and Data "+chanName,"BG and Data "+chanName,600,600);
		h_initial_BG->Draw();
		//hdata->DrawNormalized("sames E1");

		// DO THE GLOBAL FIT
		
		RooMinimizer minim(*nll);
        	//set some options:
        	minim.setPrintLevel(1);
        	minim.optimizeConst(1);
        	minim.setOffsetting(true);
        	minim.setMinimizerType("Minuit2");
        	minim.minimize("Minuit2");
        	minim.setStrategy(3); //0-3 where 0 is the fastest
        	minim.migrad();
        
		double FinalFpt = ws->var("fl1pt_l1pt0")->getVal();		
		double FinalAlpha = ws->var("alpha_l1ptsys_l1pt0")->getVal();

		// get gammas after fit
		double FinalGamma[nbins];
		TH1* h_initBG_times_gammaFpt = (TH1*)h_initial_BG->Clone("initBG_times_gamma");
		for (int i=0; i<nbins; i++)
        	{
                	TString varname = "gamma_B0_0j_l1pt0_bin_"+NumberToString(i);
                	FinalGamma[i] = ws->var(varname)->getVal();
                	cout << "Final gamma in bin "+NumberToString(i)+" = " << FinalGamma[i] << endl;
        		h_initBG_times_gammaFpt->SetBinContent(i+1,h_initial_BG->GetBinContent(i+1)*FinalGamma[i]*FinalFpt*(1+FinalAlpha*Fpt0Sigma));
		}
		cout << "initial fpt_l1pt0 = " << InitFpt <<  endl;
		cout << "final fpt_l1pt0 = " << FinalFpt <<  endl;
	
		TH1* h_final_BG = pdf_state->createHistogram("final_BG_"+chanName,*obs);
		h_final_BG->Scale(BGnumEvents);
        	//TCanvas* cf = new TCanvas("final BG","final BG",600,600);
		h_final_BG->Draw("sames");
		h_initBG_times_gammaFpt->Draw("sames");
		TH1* h_ratio = (TH1*)h_initial_BG->Clone("h_ratio");
		h_ratio->Divide(h_final_BG);
		//h_ratio->Draw();
		cout << "channel name = " << chanName << endl;	
		for ( int j=1; j<=nbins; j++)
		{
			double init = h_initial_BG->GetBinContent(j);
			double fina = h_final_BG->GetBinContent(j);
			double r = (fina)/init;
			cout << "in bin " << j << ", initial B = " << init << ", final B = " << fina << ", ratio = " << r << ", Gamma = " << FinalGamma[j-1] << ", gamma x fpt = " << FinalGamma[j-1]*FinalFpt*(1+FinalAlpha*Fpt0Sigma) << endl;
		}
		//TCanvas* c2 = new TCanvas(chanName+" & BG",chanName+" & BG",600,600);
		//TH1* h_initBG_scaled = (TH1*)h_initial_BG->Clone("initBG_scaled");
		//h_initBG_scaled->Scale(100);
		//h_initBG_scaled->Draw();
		hdata->Draw("E1 sames");	
	}


}
