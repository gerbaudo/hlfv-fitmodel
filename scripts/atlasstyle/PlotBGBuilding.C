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

char* binLabels[19] = {"60","70","80","90","100","110","120","130","140","150","160","170","180","190","200","250","300","400","1000"};

//char* binLabels[36] = {"60", "70","75", "80", "85", "90", "95", "100", "105", "110", "115", "120", "125", "130", "135", "140", "145", "150", "150", "160", "165", "170", "175", "180", "190", "200", "210", "220", "230", "240", "250", "300", "350", "400", "450", "1000"}


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







void PlotBGBuilding(TString wsname,TString inputfile)
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
	//double FinalFpt = ws->var("fl1pt_l1pt0")->getVal();
	
	// get final alpha (pull)
//	RooRealVar* alphaVar = ws->var("alpha_l1ptsys_l1pt0");
//	double alpha, alphaErr;
//	if (alphaVar != NULL) {
//		alpha = ws->var("alpha_l1ptsys_l1pt0")->getVal();
//		alphaErr = ws->var("alpha_l1ptsys_l1pt0")->getError();
//	}

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


	// Get 1/2mu(S-Scont)
	TH1D* h_halfMuSig0 = (TH1D*)h_S0->Clone("halfmuSig0");
	h_halfMuSig0->Add(h_Scont0,-1);
	h_halfMuSig0->Scale(0.5*mu);
	
	// Get gammaB+1/2muSig
	TH1D* h_gammaBplusHalfMuSig0 = (TH1D*)h_gammaB0->Clone("gammaBplusHalfMuSig0");
	h_gammaBplusHalfMuSig0->Add(h_halfMuSig0);
	TH1D* h_gammaBminusHalfMuSig0 = (TH1D*)h_gammaB0->Clone("gammaBminusHalfMuSig0");
        h_gammaBminusHalfMuSig0->Add(h_halfMuSig0,-1);

	TH1D* h_gammaBminusHalfMuSigTimesF0 = (TH1D*)h_gammaBminusHalfMuSig0->Clone("gammaB0minusHalfMuSigTimesF");
        h_gammaBminusHalfMuSigTimesF0->Scale(FinalFpt[0]); 
	
	h_gammaBminusHalfMuSigTimesF0->SetLineColor(kRed);
	h_gammaBminusHalfMuSigTimesF0->SetLineStyle(2);
	TH1D* h_gammaBTimesF0 = (TH1D*)h_gammaB0->Clone("gammaBTimesF0");       
        h_gammaBTimesF0->Scale(FinalFpt[0]);

	// uncertainty bands	
	TH1D* BuncertaintyEM = new TH1D("BuncertaintyEM","BuncertaintyEM",nbins,0,nbins);
        TH1D* BuncertaintyME = new TH1D("BuncertaintyME","BuncertaintyME",nbins,0,nbins);
        for (int i=1; i<=nbins; i++){
                double sigbEM = h_final_BG_EM->GetBinError(i);
                double bEM = h_gammaBminusHalfMuSigTimesF0->GetBinContent(i);
                BuncertaintyEM->SetBinError(i,sigbEM); BuncertaintyEM->SetBinContent(i,bEM);
                double sigbME = h_final_BG_ME->GetBinError(i);
                double bME = h_gammaBplusHalfMuSig0->GetBinContent(i);
                BuncertaintyME->SetBinError(i,sigbME); BuncertaintyME->SetBinContent(i,bME);
        }
        BuncertaintyEM->SetFillColor(kGreen-9);
        BuncertaintyEM->SetLineColor(kBlack); BuncertaintyEM->SetLineStyle(2);
        BuncertaintyME->SetFillColor(kBlue-9);
        BuncertaintyME->SetLineColor(kBlack); BuncertaintyME->SetLineStyle(2);


	// Draw
	TCanvas* c1 = new TCanvas("BG and Data "+chanName1,"BG and Data "+chanName1,600,600);
	SetLabels(BuncertaintyEM);
	BuncertaintyEM->Draw("E3 sames"); 
	hdataEM->SetLineColor(kGreen+2); hdataEM->SetMarkerStyle(20); hdataEM->SetMarkerColor(kGreen+2);
	hdataEM->Draw("e1 sames");
	//h_final_BG_EM->SetLineColor(kGreen+2); h_final_BG_EM->SetLineWidth(2); h_final_BG_EM->Draw("sames");
	
	h_gammaB0->SetLineStyle(2); h_gammaB0->Draw("hist sames");
	
	h_gammaBminusHalfMuSigTimesF0->Draw("hist sames");

	h_gammaBTimesF0->SetLineColor(kGreen+2); h_gammaBTimesF0->SetLineWidth(2);
        h_gammaBTimesF0->Draw("hist sames");

	TLegend* legEM = new TLegend(0.5,0.45,0.85,0.65);
        legEM->SetFillColor(kWhite); legEM->SetBorderSize(1); legEM->SetLineColor(0); //leg->SetTextFont(14);
        legEM->SetTextSize(.03);
	legEM->AddEntry(hdataEM,"MC e#mu","lep");
	legEM->AddEntry(h_gammaB0, "#gamma_{i}B_{i}","l");
	legEM->AddEntry(h_gammaBminusHalfMuSigTimesF0,"e#mu PDF = f(#gamma_{i}B_{i} - #frac{1}{2}#mu(S_{i}-S^{cont.}_{i}))","l");
	legEM->Draw();

	TCanvas* c2 = new TCanvas("BG and Data "+chanName2,"BG and Data "+chanName2,600,600);
	SetLabels(BuncertaintyME);
	BuncertaintyME->Draw("E3 sames");
	hdataME->SetLineColor(kBlue); hdataME->SetMarkerStyle(20);  hdataME->SetMarkerColor(kBlue);
        hdataME->Draw("e1 sames");
	//h_final_BG_ME->SetLineColor(kBlue); h_final_BG_ME->SetLineWidth(2); h_final_BG_ME->Draw("sames");
	//h_B0->Draw("hist sames");	
	h_gammaB0->Draw("hist sames");
	h_gammaBplusHalfMuSig0->SetLineColor(kBlue); h_gammaBplusHalfMuSig0->SetLineWidth(2);
	h_gammaBplusHalfMuSig0->Draw("hist sames");
	TLegend* legME = new TLegend(0.5,0.45,0.85,0.65);
	legME->SetFillColor(kWhite); legME->SetBorderSize(1); legME->SetLineColor(0); //leg->SetTextFont(14);
        legME->SetTextSize(.03);
        legME->AddEntry(hdataME,"MC #mue","lep");
	legME->AddEntry(h_gammaB0, "#gamma_{i}B_{i}","l");
	legME->AddEntry(h_gammaBplusHalfMuSig0,"#mue PDF = #gamma_{i}B_{i} + #frac{1}{2}#mu(S_{i}-S^{cont.}_{i})","l");
	legME->Draw();

	cout << " ********************* Fit Values **************************** " <<  endl;
	//if (alphaVar != NULL){cout << "alpha = " << alpha << " +- " << alphaErr << endl;}
	cout << "mu    = " << mu << " +- " << muErr << endl;

	TString WriteDownAlphaValue;
	TString WriteDownMuValue;
	WriteDownAlphaValue = "f0 = ";
	WriteDownMuValue = "#mu = ";
	WriteDownAlphaValue += Form("%4.4f",FinalFpt[0]);
	WriteDownAlphaValue += "#pm";
	WriteDownAlphaValue += Form("%4.4f",FinalFptErr[0]);
	WriteDownMuValue += Form("%4.4f",mu);
        WriteDownMuValue += "#pm";
        WriteDownMuValue += Form("%4.4f",muErr);

	TLatex *texl = new TLatex(7,35,WriteDownAlphaValue);
   	texl->SetTextAlign(22); texl->SetTextSize(0.03); 
   	TLatex *texl2 = new TLatex(3,35,WriteDownMuValue);
        texl2->SetTextAlign(22); texl2->SetTextSize(0.03);
	TLatex *tex0 = new TLatex(15,35,"l1pT0");
	texl2->Draw();
	tex0->Draw();
	c1->cd();
	texl->Draw();
        texl2->Draw();
	tex0->Draw();


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

	TH1D* h_gammaB1 = (TH1D*)h_B1->Clone("gammaB1");	
	double FinalGamma1[nbins];
        for (int i=0; i<nbins; i++)
        {
                TString varname = "gamma_B0_l1pt1_bin_"+NumberToString(i);
                FinalGamma1[i] = ws->var(varname)->getVal();
                cout << "l1pt1: Final gamma in bin "+NumberToString(i)+" = " << FinalGamma1[i] << endl;
		h_gammaB1->SetBinContent(i+1,h_B1->GetBinContent(i+1)*FinalGamma1[i]);
        }

	TH1D* h_halfMuSig1 = (TH1D*)h_S1->Clone("halfmuSig1");
        h_halfMuSig1->Add(h_Scont1,-1);
        h_halfMuSig1->Scale(0.5*mu);

	TH1D* h_gammaBplusHalfMuSig1 = (TH1D*)h_gammaB1->Clone("gammaBplusHalfMuSig1");
        h_gammaBplusHalfMuSig1->Add(h_halfMuSig1);
        TH1D* h_gammaBminusHalfMuSig1 = (TH1D*)h_gammaB1->Clone("gammaBminusHalfMuSig1");
        h_gammaBminusHalfMuSig1->Add(h_halfMuSig1,-1);

        TH1D* h_gammaBminusHalfMuSigTimesF1 = (TH1D*)h_gammaBminusHalfMuSig1->Clone("gammaB0minusHalfMuSigTimesF1");
        h_gammaBminusHalfMuSigTimesF1->Scale(FinalFpt[1]);

	TH1D* BuncertaintyEM1 = new TH1D("BuncertaintyEM1","BuncertaintyEM1",nbins,0,nbins);
        TH1D* BuncertaintyME1 = new TH1D("BuncertaintyME1","BuncertaintyME1",nbins,0,nbins);
        for (int i=1; i<=nbins; i++){
                double sigbEM = h_final_BG_EM1->GetBinError(i);
                double bEM = h_gammaBminusHalfMuSigTimesF1->GetBinContent(i);
                BuncertaintyEM1->SetBinError(i,sigbEM); BuncertaintyEM1->SetBinContent(i,bEM);
                double sigbME = h_final_BG_ME1->GetBinError(i);
                double bME = h_gammaBplusHalfMuSig1->GetBinContent(i);
                BuncertaintyME1->SetBinError(i,sigbME); BuncertaintyME1->SetBinContent(i,bME);
        }
        BuncertaintyEM1->SetFillColor(kGreen-9);
        BuncertaintyEM1->SetLineColor(kBlack); BuncertaintyEM1->SetLineStyle(2);
        BuncertaintyME1->SetFillColor(kBlue-9);
        BuncertaintyME1->SetLineColor(kBlack); BuncertaintyME1->SetLineStyle(2);
	
	TCanvas* c11 = new TCanvas("BG and Data "+chanName11,"BG and Data "+chanName11,600,600);
        SetLabels(BuncertaintyEM1);
	BuncertaintyEM1->Draw("E3 sames");
        hdataEM1->SetLineColor(kGreen+2); hdataEM1->SetMarkerStyle(20); hdataEM1->SetMarkerColor(kGreen+2);
        hdataEM1->Draw("e1 sames");
        //h_final_BG_EM1->SetLineColor(kGreen+2); h_final_BG_EM1->SetLineWidth(2); h_final_BG_EM1->Draw("sames");
        h_gammaB1->SetLineStyle(2); h_gammaB1->Draw("hist sames");
        h_gammaBminusHalfMuSigTimesF1->SetLineColor(kGreen+2);  h_gammaBminusHalfMuSigTimesF1->SetLineWidth(2);
        h_gammaBminusHalfMuSigTimesF1->Draw("hist sames");

        legEM->Draw();
	texl2->Draw();
	TLatex *tex1 = new TLatex(15,35,"l1pT1");
	tex1->Draw();

        TCanvas* c21 = new TCanvas("BG and Data "+chanName21,"BG and Data "+chanName21,600,600);
	SetLabels(BuncertaintyME1);
	BuncertaintyME1->Draw("E3 sames");
        hdataME1->SetLineColor(kBlue); hdataME1->SetMarkerStyle(20);  hdataME1->SetMarkerColor(kBlue);
        hdataME1->Draw("e1 sames");
        //h_final_BG_ME1->SetLineColor(kBlue); h_final_BG_ME1->SetLineWidth(2); h_final_BG_ME1->Draw("sames");
	h_gammaB1->Draw("hist sames");
        h_gammaBplusHalfMuSig1->SetLineColor(kBlue); h_gammaBplusHalfMuSig1->SetLineWidth(2);
        h_gammaBplusHalfMuSig1->Draw("hist sames");
	legME->Draw();
        
	cout << " ********************* Fit Values **************************** " <<  endl;
        cout << "mu    = " << mu << " +- " << muErr << endl;
	TString WriteDownAlphaValue1;
        WriteDownAlphaValue1 = "f1 = ";
        WriteDownAlphaValue1 += Form("%4.4f",FinalFpt[1]);
        WriteDownAlphaValue1 += "#pm";
        WriteDownAlphaValue1 += Form("%4.4f",FinalFptErr[1]);

        TLatex *texl11 = new TLatex(7,30,WriteDownAlphaValue1);
        texl11->SetTextAlign(22); texl11->SetTextSize(0.03);
        //texl11->Draw(); 
        texl2->Draw();
	tex1->Draw();
	c11->cd();
        texl11->Draw();


	// FPT2
	 ttype = (RooCatType*)iterat->Next();

        RooAbsPdf  *pdf_stateEM2  = simPdf->getPdf(ttype->GetName()) ;
        RooArgSet  *obstmpEM2  = pdf_stateEM2->getObservables( *mc->GetObservables() ) ;
        RooAbsData *dataEM2 = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));

        RooRealVar *obsEM2     = ((RooRealVar*) obstmpEM2->first());
        TString chanName12(ttype->GetName());
        TH1* hdataEM2 = dataEM2->createHistogram("Data "+chanName12,*obsEM2);
        for (int ib=0 ; ib<hdataEM2->GetNbinsX()+1 ; ib++) hdataEM2->SetBinError(ib, sqrt(hdataEM2->GetBinContent(ib)));

        double EMnorm2 = pdf_stateEM2->expectedEvents(*obsEM2);
        ttype = (RooCatType*)iterat->Next();
        RooAbsPdf* pdf_stateME2  = simPdf->getPdf(ttype->GetName()) ;
        RooArgSet* obstmpME2  = pdf_stateME2->getObservables( *mc->GetObservables() ) ;
        RooAbsData *dataME2 = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));
        RooRealVar* obsME2 = ((RooRealVar*) obstmpME2->first());
        TString chanName22(ttype->GetName());
        TH1* hdataME2 = dataME2->createHistogram("Data "+chanName22,*obsME2);

        for (int ib=0 ; ib<hdataME2->GetNbinsX()+1 ; ib++) hdataME2->SetBinError(ib, sqrt(hdataME2->GetBinContent(ib)));
        double MEnorm2 = pdf_stateME2->expectedEvents(*obsME2);
        TH1* h_final_BG_EM2 = pdf_stateEM2->createHistogram("final_BG_EM2",*obsEM2);
        TH1* h_final_BG_ME2 = pdf_stateME2->createHistogram("final_BG_ME2",*obsME2);
        h_final_BG_EM2->Scale(EMnorm2);
        h_final_BG_ME2->Scale(MEnorm2);
        
	TH1D* h_gammaB2 = (TH1D*)h_B2->Clone("gammaB2");
	double FinalGamma2[nbins];
        for (int i=0; i<nbins; i++)
        {
                TString varname = "gamma_B0_l1pt2_bin_"+NumberToString(i);
                FinalGamma2[i] = ws->var(varname)->getVal();
                cout << "l1pt2 : Final gamma in bin "+NumberToString(i)+" = " << FinalGamma2[i] << endl;
		h_gammaB2->SetBinContent(i+1,h_B2->GetBinContent(i+1)*FinalGamma2[i]);
        }

	TH1D* h_halfMuSig2 = (TH1D*)h_S2->Clone("halfmuSig2");
        h_halfMuSig2->Add(h_Scont2,-1);
        h_halfMuSig2->Scale(0.5*mu);

        TH1D* h_gammaBplusHalfMuSig2 = (TH1D*)h_gammaB2->Clone("gammaBplusHalfMuSig2");
        h_gammaBplusHalfMuSig2->Add(h_halfMuSig2);
        TH1D* h_gammaBminusHalfMuSig2 = (TH1D*)h_gammaB2->Clone("gammaBminusHalfMuSig2");
        h_gammaBminusHalfMuSig2->Add(h_halfMuSig2,-1);

        TH1D* h_gammaBminusHalfMuSigTimesF2 = (TH1D*)h_gammaBminusHalfMuSig2->Clone("gammaB0minusHalfMuSigTimesF2");
        h_gammaBminusHalfMuSigTimesF2->Scale(FinalFpt[2]);


	TH1D* BuncertaintyEM2 = new TH1D("BuncertaintyEM2","BuncertaintyEM2",nbins,0,nbins);
        TH1D* BuncertaintyME2 = new TH1D("BuncertaintyME2","BuncertaintyME2",nbins,0,nbins);
        for (int i=1; i<=nbins; i++){
                double sigbEM = h_final_BG_EM2->GetBinError(i);
                double bEM = h_gammaBminusHalfMuSigTimesF2->GetBinContent(i);
                BuncertaintyEM2->SetBinError(i,sigbEM); BuncertaintyEM2->SetBinContent(i,bEM);
                double sigbME = h_final_BG_ME2->GetBinError(i);
                double bME = h_gammaBplusHalfMuSig2->GetBinContent(i);
                BuncertaintyME2->SetBinError(i,sigbME); BuncertaintyME2->SetBinContent(i,bME);
        }

        BuncertaintyEM2->SetFillColor(kGreen-9);
        BuncertaintyEM2->SetLineColor(kBlack); BuncertaintyEM2->SetLineStyle(2);
        BuncertaintyME2->SetFillColor(kBlue-9);
        BuncertaintyME2->SetLineColor(kBlack); BuncertaintyME2->SetLineStyle(2);

        TCanvas* c12 = new TCanvas("BG and Data "+chanName12,"BG and Data "+chanName12,600,600);
        SetLabels(BuncertaintyEM2);
	BuncertaintyEM2->Draw("E3 sames");
        hdataEM2->SetLineColor(kGreen+2); hdataEM2->SetMarkerStyle(20); hdataEM2->SetMarkerColor(kGreen+2);
        hdataEM2->Draw("e1 sames");
//        h_final_BG_EM2->SetLineColor(kGreen+2); h_final_BG_EM2->SetLineWidth(2); h_final_BG_EM2->Draw("sames");
	h_gammaB2->SetLineStyle(2); h_gammaB2->Draw("hist sames");
        h_gammaBminusHalfMuSigTimesF2->SetLineColor(kGreen+2);  h_gammaBminusHalfMuSigTimesF2->SetLineWidth(2);
        h_gammaBminusHalfMuSigTimesF2->Draw("hist sames");
 	legEM->Draw();
        texl2->Draw();
	TLatex *tex2 = new TLatex(15,35,"l1pT2");
	tex2->Draw();
	TCanvas* c22 = new TCanvas("BG and Data "+chanName22,"BG and Data "+chanName22,600,600);
	SetLabels(BuncertaintyME2);
        BuncertaintyME2->Draw("E3 sames");
	hdataME2->SetLineColor(kBlue); hdataME2->SetMarkerStyle(20);  hdataME2->SetMarkerColor(kBlue);
        hdataME2->Draw("e1 sames");
  //      h_final_BG_ME2->SetLineColor(kBlue); h_final_BG_ME2->SetLineWidth(2); h_final_BG_ME2->Draw("sames");
    	h_gammaB2->Draw("hist sames");
        h_gammaBplusHalfMuSig2->SetLineColor(kBlue); h_gammaBplusHalfMuSig2->SetLineWidth(2);
        h_gammaBplusHalfMuSig2->Draw("hist sames");
	legME->Draw();

        cout << " ********************* Fit Values **************************** " <<  endl;
        cout << "mu    = " << mu << " +- " << muErr << endl;
        TString WriteDownAlphaValue2;
        WriteDownAlphaValue2 = "f2 = ";
        WriteDownAlphaValue2 += Form("%4.4f",FinalFpt[2]);
        WriteDownAlphaValue2 += "#pm";
        WriteDownAlphaValue2 += Form("%4.4f",FinalFptErr[2]);

        TLatex *texl12 = new TLatex(7,30,WriteDownAlphaValue2);
        texl12->SetTextAlign(22); texl12->SetTextSize(0.03);
        //texl12->Draw();
        texl2->Draw();
	tex2->Draw();
	c12->cd();
	texl12->Draw();

	// FPT 3
	ttype = (RooCatType*)iterat->Next();

        RooAbsPdf  *pdf_stateEM3  = simPdf->getPdf(ttype->GetName()) ;
        RooArgSet  *obstmpEM3  = pdf_stateEM3->getObservables( *mc->GetObservables() ) ;
        RooAbsData *dataEM3 = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));

        RooRealVar *obsEM3     = ((RooRealVar*) obstmpEM3->first());
        TString chanName13(ttype->GetName());
        TH1* hdataEM3 = dataEM3->createHistogram("Data "+chanName13,*obsEM3);
        for (int ib=0 ; ib<hdataEM3->GetNbinsX()+1 ; ib++) hdataEM3->SetBinError(ib, sqrt(hdataEM3->GetBinContent(ib)));

        double EMnorm3 = pdf_stateEM3->expectedEvents(*obsEM3);
        ttype = (RooCatType*)iterat->Next();
        RooAbsPdf* pdf_stateME3  = simPdf->getPdf(ttype->GetName()) ;
        RooArgSet* obstmpME3  = pdf_stateME3->getObservables( *mc->GetObservables() ) ;
        RooAbsData *dataME3 = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));
        RooRealVar* obsME3 = ((RooRealVar*) obstmpME3->first());
        TString chanName23(ttype->GetName());
        TH1* hdataME3 = dataME3->createHistogram("Data "+chanName23,*obsME3);

        for (int ib=0 ; ib<hdataME3->GetNbinsX()+1 ; ib++) hdataME3->SetBinError(ib, sqrt(hdataME3->GetBinContent(ib)));
        double MEnorm3 = pdf_stateME3->expectedEvents(*obsME3);
        TH1* h_final_BG_EM3 = pdf_stateEM3->createHistogram("final_BG_EM3",*obsEM3);
        TH1* h_final_BG_ME3 = pdf_stateME3->createHistogram("final_BG_ME3",*obsME3);
        h_final_BG_EM3->Scale(EMnorm3);
        h_final_BG_ME3->Scale(MEnorm3);
        
	TH1D* h_gammaB3 = (TH1D*)h_B3->Clone("gammaB3");
	double FinalGamma3[nbins];
        for (int i=0; i<nbins; i++)
        {
                TString varname = "gamma_B0_l1pt3_bin_"+NumberToString(i);
                FinalGamma3[i] = ws->var(varname)->getVal();
                cout << "l1pt3 : Final gamma in bin "+NumberToString(i)+" = " << FinalGamma3[i] << endl;
		h_gammaB3->SetBinContent(i+1,h_B3->GetBinContent(i+1)*FinalGamma3[i]);
        }


	TH1D* h_halfMuSig3 = (TH1D*)h_S3->Clone("halfmuSig3");
        h_halfMuSig3->Add(h_Scont3,-1);
        h_halfMuSig3->Scale(0.5*mu);

        TH1D* h_gammaBplusHalfMuSig3 = (TH1D*)h_gammaB3->Clone("gammaBplusHalfMuSig3");
        h_gammaBplusHalfMuSig3->Add(h_halfMuSig3);
        TH1D* h_gammaBminusHalfMuSig3 = (TH1D*)h_gammaB3->Clone("gammaBminusHalfMuSig3");
        h_gammaBminusHalfMuSig3->Add(h_halfMuSig3,-1);

        TH1D* h_gammaBminusHalfMuSigTimesF3 = (TH1D*)h_gammaBminusHalfMuSig3->Clone("gammaB0minusHalfMuSigTimesF3");
        h_gammaBminusHalfMuSigTimesF3->Scale(FinalFpt[3]);

	TH1D* BuncertaintyEM3 = new TH1D("BuncertaintyEM3","BuncertaintyEM3",nbins,0,nbins);
        TH1D* BuncertaintyME3 = new TH1D("BuncertaintyME3","BuncertaintyME3",nbins,0,nbins);
        for (int i=1; i<=nbins; i++){
                double sigbEM = h_final_BG_EM3->GetBinError(i);
                double bEM = h_gammaBminusHalfMuSigTimesF3->GetBinContent(i);
                BuncertaintyEM3->SetBinError(i,sigbEM); BuncertaintyEM3->SetBinContent(i,bEM);
                double sigbME = h_final_BG_ME3->GetBinError(i);
                double bME = h_gammaBplusHalfMuSig3->GetBinContent(i);
                BuncertaintyME3->SetBinError(i,sigbME); BuncertaintyME3->SetBinContent(i,bME);
        }
        BuncertaintyEM3->SetFillColor(kGreen-9);
        BuncertaintyEM3->SetLineColor(kBlack); BuncertaintyEM3->SetLineStyle(2);
        BuncertaintyME3->SetFillColor(kBlue-9);
        BuncertaintyME3->SetLineColor(kBlack); BuncertaintyME3->SetLineStyle(2);


        TCanvas* c13 = new TCanvas("BG and Data "+chanName13,"BG and Data "+chanName13,600,600);
	SetLabels(BuncertaintyEM3);
	BuncertaintyEM3->Draw("E3 sames");
        hdataEM3->SetLineColor(kGreen+2); hdataEM3->SetMarkerStyle(20); hdataEM3->SetMarkerColor(kGreen+2);
        hdataEM3->Draw("e1 sames");
        //h_final_BG_EM3->SetLineColor(kGreen+2); h_final_BG_EM3->SetLineWidth(2); h_final_BG_EM3->Draw("sames");
        h_gammaB3->SetLineStyle(2); h_gammaB3->Draw("hist sames");
        h_gammaBminusHalfMuSigTimesF3->SetLineColor(kGreen+2);  h_gammaBminusHalfMuSigTimesF2->SetLineWidth(2);
        h_gammaBminusHalfMuSigTimesF3->Draw("hist sames");
	legEM->Draw();
        texl2->Draw();
	TLatex *tex3 = new TLatex(15,35,"l1pT3");
	tex3->Draw();
        TCanvas* c23 = new TCanvas("BG and Data "+chanName23,"BG and Data "+chanName23,600,600);
	SetLabels(BuncertaintyME3);
        BuncertaintyME3->Draw("E3 sames");
	hdataME3->SetLineColor(kBlue); hdataME3->SetMarkerStyle(20);  hdataME3->SetMarkerColor(kBlue);
        hdataME3->Draw("e1 sames");
//        h_final_BG_ME3->SetLineColor(kBlue); h_final_BG_ME3->SetLineWidth(2); h_final_BG_ME3->Draw("sames");
	h_gammaB3->Draw("hist sames");
        h_gammaBplusHalfMuSig3->SetLineColor(kBlue); h_gammaBplusHalfMuSig3->SetLineWidth(2);
        h_gammaBplusHalfMuSig3->Draw("hist sames");
	legME->Draw();

        cout << " ********************* Fit Values **************************** " <<  endl;
        cout << "mu    = " << mu << " +- " << muErr << endl;
        TString WriteDownAlphaValue3;
        WriteDownAlphaValue3 = "f3 = ";
        WriteDownAlphaValue3 += Form("%4.4f",FinalFpt[3]);
        WriteDownAlphaValue3 += "#pm";
        WriteDownAlphaValue3 += Form("%4.4f",FinalFptErr[3]);

        TLatex *texl13 = new TLatex(7,30,WriteDownAlphaValue3);
        texl13->SetTextAlign(22); texl13->SetTextSize(0.03);
        //texl13->Draw();
        texl2->Draw();
	tex3->Draw();
	c13->cd();
        texl13->Draw();


	// FPT 4
	ttype = (RooCatType*)iterat->Next();

        RooAbsPdf  *pdf_stateEM4  = simPdf->getPdf(ttype->GetName()) ;
        RooArgSet  *obstmpEM4  = pdf_stateEM4->getObservables( *mc->GetObservables() ) ;
        RooAbsData *dataEM4 = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));

        RooRealVar *obsEM4     = ((RooRealVar*) obstmpEM4->first());
        TString chanName14(ttype->GetName());
        TH1* hdataEM4 = dataEM4->createHistogram("Data "+chanName14,*obsEM4);
        for (int ib=0 ; ib<hdataEM4->GetNbinsX()+1 ; ib++) hdataEM4->SetBinError(ib, sqrt(hdataEM4->GetBinContent(ib)));

        double EMnorm4 = pdf_stateEM4->expectedEvents(*obsEM4);
        ttype = (RooCatType*)iterat->Next();
        RooAbsPdf* pdf_stateME4  = simPdf->getPdf(ttype->GetName()) ;
        RooArgSet* obstmpME4  = pdf_stateME4->getObservables( *mc->GetObservables() ) ;
        RooAbsData *dataME4 = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));
        RooRealVar* obsME4 = ((RooRealVar*) obstmpME4->first());
        TString chanName24(ttype->GetName());
        TH1* hdataME4 = dataME4->createHistogram("Data "+chanName24,*obsME4);

        for (int ib=0 ; ib<hdataME4->GetNbinsX()+1 ; ib++) hdataME4->SetBinError(ib, sqrt(hdataME4->GetBinContent(ib)));
        double MEnorm4 = pdf_stateME4->expectedEvents(*obsME4);
        TH1* h_final_BG_EM4 = pdf_stateEM4->createHistogram("final_BG_EM4",*obsEM4);
        TH1* h_final_BG_ME4 = pdf_stateME4->createHistogram("final_BG_ME4",*obsME4);
        h_final_BG_EM4->Scale(EMnorm4);
        h_final_BG_ME4->Scale(MEnorm4);
        
	TH1D* h_gammaB4 = (TH1D*)h_B4->Clone("gammaB4");
	double FinalGamma4[nbins];
	for (int i=0; i<nbins; i++)
        {
                TString varname = "gamma_B0_l1pt4_bin_"+NumberToString(i);
                FinalGamma4[i] = ws->var(varname)->getVal();
                cout << "Final gamma in bin "+NumberToString(i)+" = " << FinalGamma4[i] << endl;
		h_gammaB4->SetBinContent(i+1,h_B4->GetBinContent(i+1)*FinalGamma4[i]);
        }


	TH1D* h_halfMuSig4 = (TH1D*)h_S4->Clone("halfmuSig4");
        h_halfMuSig4->Add(h_Scont4,-1);
        h_halfMuSig4->Scale(0.5*mu);

        TH1D* h_gammaBplusHalfMuSig4 = (TH1D*)h_gammaB4->Clone("gammaBplusHalfMuSig4");
        h_gammaBplusHalfMuSig4->Add(h_halfMuSig4);
        TH1D* h_gammaBminusHalfMuSig4 = (TH1D*)h_gammaB4->Clone("gammaBminusHalfMuSig4");
        h_gammaBminusHalfMuSig4->Add(h_halfMuSig4,-1);

        TH1D* h_gammaBminusHalfMuSigTimesF4 = (TH1D*)h_gammaBminusHalfMuSig4->Clone("gammaB0minusHalfMuSigTimesF4");
        h_gammaBminusHalfMuSigTimesF4->Scale(FinalFpt[4]);

	TH1D* BuncertaintyEM4 = new TH1D("BuncertaintyEM4","BuncertaintyEM4",nbins,0,nbins);
        TH1D* BuncertaintyME4 = new TH1D("BuncertaintyME4","BuncertaintyME4",nbins,0,nbins);
        for (int i=1; i<=nbins; i++){
                double sigbEM = h_final_BG_EM4->GetBinError(i);
                double bEM = h_gammaBminusHalfMuSigTimesF4->GetBinContent(i);
                BuncertaintyEM4->SetBinError(i,sigbEM); BuncertaintyEM4->SetBinContent(i,bEM);
                double sigbME = h_final_BG_ME4->GetBinError(i);
                double bME = h_gammaBplusHalfMuSig4->GetBinContent(i);
                BuncertaintyME4->SetBinError(i,sigbME); BuncertaintyME4->SetBinContent(i,bME);
        }
        BuncertaintyEM4->SetFillColor(kGreen-9);
        BuncertaintyEM4->SetLineColor(kBlack); BuncertaintyEM4->SetLineStyle(2);
        BuncertaintyME4->SetFillColor(kBlue-9);
        BuncertaintyME4->SetLineColor(kBlack); BuncertaintyME4->SetLineStyle(2);


        TCanvas* c14 = new TCanvas("BG and Data "+chanName14,"BG and Data "+chanName14,600,600);
	SetLabels(BuncertaintyEM4);
        BuncertaintyEM4->Draw("E3 sames");
        hdataEM4->SetLineColor(kGreen+2); hdataEM4->SetMarkerStyle(20); hdataEM4->SetMarkerColor(kGreen+2);
        hdataEM4->Draw("e1 sames");
//        h_final_BG_EM4->SetLineColor(kGreen+2); h_final_BG_EM4->SetLineWidth(2); h_final_BG_EM4->Draw("sames");
  	 h_gammaB4->SetLineStyle(2); h_gammaB4->Draw("hist sames");
        h_gammaBminusHalfMuSigTimesF4->SetLineColor(kGreen+2);  h_gammaBminusHalfMuSigTimesF2->SetLineWidth(2);
        h_gammaBminusHalfMuSigTimesF4->Draw("hist sames");

	legEM->Draw();
        texl2->Draw();
	TLatex *tex4 = new TLatex(15,35,"l1pT4");
	tex4->Draw();
        TCanvas* c24 = new TCanvas("BG and Data "+chanName24,"BG and Data "+chanName24,600,600);
	SetLabels(BuncertaintyME4);
        BuncertaintyME4->Draw("E3 sames");
        hdataME4->SetLineColor(kBlue); hdataME4->SetMarkerStyle(20);  hdataME4->SetMarkerColor(kBlue);
        hdataME4->Draw("e1 sames");
//        h_final_BG_ME4->SetLineColor(kBlue); h_final_BG_ME4->SetLineWidth(2); h_final_BG_ME4->Draw("sames");
  	h_gammaB4->Draw("hist sames");
        h_gammaBplusHalfMuSig4->SetLineColor(kBlue); h_gammaBplusHalfMuSig4->SetLineWidth(2);
        h_gammaBplusHalfMuSig4->Draw("hist sames");
	legME->Draw();

        cout << " ********************* Fit Values **************************** " <<  endl;
        cout << "mu    = " << mu << " +- " << muErr << endl;
        TString WriteDownAlphaValue4;
        WriteDownAlphaValue4 = "f4 = ";
        WriteDownAlphaValue4 += Form("%4.4f",FinalFpt[4]);
        WriteDownAlphaValue4 += "#pm";
        WriteDownAlphaValue4 += Form("%4.4f",FinalFptErr[4]);

        TLatex *texl14 = new TLatex(7,30,WriteDownAlphaValue4);
        texl14->SetTextAlign(22); texl14->SetTextSize(0.03);
        //texl14->Draw();
        texl2->Draw();
	tex4->Draw();
	c14->cd();
        texl14->Draw();


	// FPT 5
	ttype = (RooCatType*)iterat->Next();

        RooAbsPdf  *pdf_stateEM5  = simPdf->getPdf(ttype->GetName()) ;
        RooArgSet  *obstmpEM5  = pdf_stateEM5->getObservables( *mc->GetObservables() ) ;
        RooAbsData *dataEM5 = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));

        RooRealVar *obsEM5     = ((RooRealVar*) obstmpEM5->first());
        TString chanName15(ttype->GetName());
        TH1* hdataEM5 = dataEM5->createHistogram("Data "+chanName15,*obsEM5);
        for (int ib=0 ; ib<hdataEM5->GetNbinsX()+1 ; ib++) hdataEM5->SetBinError(ib, sqrt(hdataEM5->GetBinContent(ib)));

        double EMnorm5 = pdf_stateEM5->expectedEvents(*obsEM5);
        ttype = (RooCatType*)iterat->Next();
        RooAbsPdf* pdf_stateME5  = simPdf->getPdf(ttype->GetName()) ;
        RooArgSet* obstmpME5  = pdf_stateME5->getObservables( *mc->GetObservables() ) ;
        RooAbsData *dataME5 = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));
        RooRealVar* obsME5 = ((RooRealVar*) obstmpME5->first());
        TString chanName25(ttype->GetName());
        TH1* hdataME5 = dataME5->createHistogram("Data "+chanName25,*obsME5);

        for (int ib=0 ; ib<hdataME5->GetNbinsX()+1 ; ib++) hdataME5->SetBinError(ib, sqrt(hdataME5->GetBinContent(ib)));
        double MEnorm5 = pdf_stateME5->expectedEvents(*obsME5);
        TH1* h_final_BG_EM5 = pdf_stateEM5->createHistogram("final_BG_EM5",*obsEM5);
        TH1* h_final_BG_ME5 = pdf_stateME5->createHistogram("final_BG_ME5",*obsME5);
        h_final_BG_EM5->Scale(EMnorm5);
        h_final_BG_ME5->Scale(MEnorm5);
        
	 TH1D* h_gammaB5 = (TH1D*)h_B5->Clone("gammaB5");
	double FinalGamma5[nbins];
	for (int i=0; i<nbins; i++)
        {
                TString varname = "gamma_B0_l1pt5_bin_"+NumberToString(i);
                FinalGamma5[i] = ws->var(varname)->getVal();
                cout << "Final gamma in bin "+NumberToString(i)+" = " << FinalGamma5[i] << endl;
		h_gammaB5->SetBinContent(i+1,h_B5->GetBinContent(i+1)*FinalGamma5[i]);
        }

	TH1D* h_halfMuSig5 = (TH1D*)h_S5->Clone("halfmuSig5");
        h_halfMuSig5->Add(h_Scont5,-1);
        h_halfMuSig5->Scale(0.5*mu);

        TH1D* h_gammaBplusHalfMuSig5 = (TH1D*)h_gammaB5->Clone("gammaBplusHalfMuSig5");
        h_gammaBplusHalfMuSig5->Add(h_halfMuSig5);
        TH1D* h_gammaBminusHalfMuSig5 = (TH1D*)h_gammaB5->Clone("gammaBminusHalfMuSig5");
        h_gammaBminusHalfMuSig5->Add(h_halfMuSig5,-1);

        TH1D* h_gammaBminusHalfMuSigTimesF5 = (TH1D*)h_gammaBminusHalfMuSig5->Clone("gammaB0minusHalfMuSigTimesF5");

	TH1D* BuncertaintyEM5 = new TH1D("BuncertaintyEM5","BuncertaintyEM5",nbins,0,nbins);
        TH1D* BuncertaintyME5 = new TH1D("BuncertaintyME5","BuncertaintyME5",nbins,0,nbins);
        for (int i=1; i<=nbins; i++){
                double sigbEM = h_final_BG_EM5->GetBinError(i);
                double bEM = h_gammaBminusHalfMuSigTimesF5->GetBinContent(i);
                BuncertaintyEM5->SetBinError(i,sigbEM); BuncertaintyEM5->SetBinContent(i,bEM);
                double sigbME = h_final_BG_ME5->GetBinError(i);
                double bME = h_gammaBplusHalfMuSig5->GetBinContent(i);
                BuncertaintyME5->SetBinError(i,sigbME); BuncertaintyME5->SetBinContent(i,bME);
        }
        BuncertaintyEM5->SetFillColor(kGreen-9);
        BuncertaintyEM5->SetLineColor(kBlack); BuncertaintyEM5->SetLineStyle(2);
        BuncertaintyME5->SetFillColor(kBlue-9);
        BuncertaintyME5->SetLineColor(kBlack); BuncertaintyME5->SetLineStyle(2);


        TCanvas* c15 = new TCanvas("BG and Data "+chanName15,"BG and Data "+chanName15,600,600);
	SetLabels(BuncertaintyEM5);
        BuncertaintyEM5->Draw("E3 sames");
        hdataEM5->SetLineColor(kGreen+2); hdataEM5->SetMarkerStyle(20); hdataEM5->SetMarkerColor(kGreen+2);
        hdataEM5->Draw("e1 sames");
//        h_final_BG_EM5->SetLineColor(kGreen+2); h_final_BG_EM5->SetLineWidth(2); h_final_BG_EM5->Draw("sames");
	h_gammaB5->SetLineStyle(2); h_gammaB5->Draw("hist sames");
        h_gammaBminusHalfMuSigTimesF5->SetLineColor(kGreen+2);  h_gammaBminusHalfMuSigTimesF2->SetLineWidth(2);
        h_gammaBminusHalfMuSigTimesF5->Draw("hist sames");
        legEM->Draw();
        texl2->Draw();
	TLatex *tex5 = new TLatex(15,35,"l1pT5");
	tex5->Draw();

        TCanvas* c25 = new TCanvas("BG and Data "+chanName25,"BG and Data "+chanName25,600,600);
	SetLabels(BuncertaintyME5);
        BuncertaintyME5->Draw("E3 sames");
        hdataME5->SetLineColor(kBlue); hdataME5->SetMarkerStyle(20);  hdataME5->SetMarkerColor(kBlue);
        hdataME5->Draw("e1 sames");
//        h_final_BG_ME5->SetLineColor(kBlue); h_final_BG_ME5->SetLineWidth(2); h_final_BG_ME5->Draw("sames");
	h_gammaB5->Draw("hist sames");
        h_gammaBplusHalfMuSig5->SetLineColor(kBlue); h_gammaBplusHalfMuSig5->SetLineWidth(2);
        h_gammaBplusHalfMuSig5->Draw("hist sames");
  
	legME->Draw();

        cout << " ********************* Fit Values **************************** " <<  endl;
        cout << "mu    = " << mu << " +- " << muErr << endl;
        TString WriteDownAlphaValue5;
        WriteDownAlphaValue5 = "f5 = 1 ";

        TLatex *texl15 = new TLatex(7,30,WriteDownAlphaValue5);
        texl15->SetTextAlign(22); texl15->SetTextSize(0.03);
        //texl15->Draw();
        texl2->Draw();
	tex5->Draw();
	c15->cd();
        texl15->Draw();

}

#ifndef __CINT__
int main() {
//  BasicExample();
//  gPad->Print("basic.png");
    return 0;
    }
#endif

