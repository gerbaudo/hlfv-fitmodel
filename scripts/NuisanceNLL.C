
#include "TMath.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "TFile.h"
#include "TLine.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TString.h"
#include "TCanvas.h"

#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooNLLVar.h"
#include "RooStats/ModelConfig.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "Math/MinimizerOptions.h"
#include "TStopwatch.h"
#include "RooMinimizerFcn.h"
#include "RooMinimizer.h"
#include "RooCategory.h"
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


int minimize(RooAbsReal* fcn)
{
  static int nrItr = 0;
int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  int strat = 2;//2;//ROOT::Math::MinimizerOptions::DefaultStrategy();
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


void NuisanceNLL(TString filename,TString varname,double min, double max)
{
	//get the stuff from the workspace:
	TFile* file=TFile::Open(filename);
	RooWorkspace* ws = (RooWorkspace*)file->Get("combined");
	ModelConfig* mc = (ModelConfig*)ws->obj("ModelConfig");
	RooAbsData* data = ws->data("obsData");
	RooSimultaneous* simPdf=(RooSimultaneous*)(mc->GetPdf());
	RooArgSet *obs = simPdf->getObservables(ws->data("obsData"));
//	RooRealVar *mu =
//	dynamic_cast<RooRealVar*>(mc->GetParametersOfInterest()->first());
	
	RooAbsReal* nll=simPdf->createNLL(*data);
	
	minimize(nll);
	
	//get the value of the minimized likelihood:
	double unconditional_NLL = nll->getVal();
	double bestFitPval = ws->var(varname)->getVal();
	cout<<"1 nll = " <<nll->getVal();

	RooRealVar* var = ws->var(varname);


	int size = 10; 
	double Nll[size+1];
	double Pval[size+1];
	double nllAt2 = 0;
	double nllAt0 = 0;
	double nllAt1 = 0;
	for (int i=0; i<size+1; i++){
		Pval[i] =  min+i*(max-min)/size; //(i+1)*10./size;
		ws->var(varname)->setVal(Pval[i]);
		ws->var(varname)->setConstant(1);
		//ws->var(varname)->removeError();
		cout << "******************Parameter value = " << ws->var(varname)->getVal() << " *********************" << endl;
		RooAbsData* data1 = ws->data("obsData");
        	RooSimultaneous* simPdf1=(RooSimultaneous*)(mc->GetPdf());
        	RooArgSet *obs1 = simPdf1->getObservables(ws->data("obsData"));

        	RooAbsReal* nll2=simPdf1->createNLL(*data1);

		minimize(nll2);
		cout << "*********Parameter value after minim = " << ws->var(varname)->getVal() << " *********************" << endl;
		double Nllval = nll2->getVal();
		Nll[i]=2*(Nllval-unconditional_NLL);
		if (Nll[i] != Nll[i]) { Nll[i] = 0;}
		if (ws->var(varname)->getVal()==1) {nllAt1 = Nllval;}
		if (ws->var(varname)->getVal()==0) {nllAt0 = Nllval;}
		if (ws->var(varname)->getVal()==2) {nllAt2 = Nllval;}
	}	
	
	TGraph* g = new TGraph(size+1,Pval,Nll);
	g->GetXaxis()->SetTitle(varname);
	g->GetYaxis()->SetTitle("-2#DeltaLog(L)");
	g->Draw("AC*");

	TLine* line = new TLine(min,4,max,4);
	line->SetLineStyle(2); line->SetLineColor(kRed);
	line->Draw();

	TLine* line2 = new TLine(min,9,max,9);
        line2->SetLineStyle(2); line2->SetLineColor(kRed);
        line2->Draw();

	cout << "best fit value = " << bestFitPval << endl;
	cout << "unconditional nll = " << unconditional_NLL << endl;
	cout << "nll at 0 = " << nllAt0 << endl;
	cout << "nll at 1 = " << nllAt1 << endl;
	cout << "nll at 2 = " << nllAt2 << endl;


}
