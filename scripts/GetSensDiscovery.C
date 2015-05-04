
#include "TMath.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "TFile.h"
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

  int strat = 2;
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


void GetSensDiscovery(TString wsname,double muVal)
{
	//get the stuff from the workspace:
	TFile* file=TFile::Open(wsname);
	RooWorkspace* ws = (RooWorkspace*)file->Get("combined");
	ModelConfig* mc = (ModelConfig*)ws->obj("ModelConfig");
	RooAbsData* data = ws->data("obsData");
	RooSimultaneous* simPdf=(RooSimultaneous*)(mc->GetPdf());
	RooArgSet *obs = simPdf->getObservables(ws->data("obsData"));
	RooRealVar *mu =
	dynamic_cast<RooRealVar*>(mc->GetParametersOfInterest()->first());
	
	RooAbsReal* nll=simPdf->createNLL(*data);
	
	minimize(nll);
	//RooMinimizer minim(*nll);
	//set some options:
	//minim.setPrintLevel(-1);
	//minim.optimizeConst(1);
	//minim.setOffsetting(true);
	//minim.setMinimizerType("Minuit2");
	//minim.minimize("Minuit2");
	//minim.setStrategy(0); //0-2 where 0 is the fastest
	//minim.migrad();
	
	//get the value of the minimized likelihood:
	double uncond_NLL = nll->getVal();
	double muHat = mu->getVal();
	cout<<"unconditional nll = " <<uncond_NLL << endl;
	cout << "muHat = " << muHat << endl;
	
	mu->setVal(muVal);
	mu->setConstant(1);
	
	minimize(nll);
	//RooMinimizer minim2(*nll);
	//minim2.setPrintLevel(-1);
       // minim2.optimizeConst(1);
        //minim2.setOffsetting(true);
        //minim2.setMinimizerType("Minuit2");
        //minim2.minimize("Minuit2");
        //minim2.setStrategy(0); //0-2 where 0 is the fastest
        //minim2.migrad();
	
	double cond_NLL = nll->getVal();
	double qmu = 2*(cond_NLL-uncond_NLL);
	double signi = std::sqrt(abs(qmu));

	cout << "significance for mu = " << mu->getVal() << " is " << signi << endl;
	//cout << "conditional nll = " << cond_NLL << endl;

}
