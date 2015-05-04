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

RooWorkspace *w         ;
ModelConfig  *mc        ;
RooAbsData   *data      ;

void Initialize(const char* infile , const char* workspaceName, const char* modelConfigName, const char* ObsDataName) {

	// Load workspace, model and data
	TFile *file = TFile::Open(infile);
	if (!file) {
		cout << "The file " << infile << " is not found/created, will stop here." << endl;
		return;
   	}	
    	if(!(RooWorkspace*) file->Get(workspaceName)){
      		cout <<"workspace not found" << endl;
      		return;
    	}

    	w      = (RooWorkspace*) file->Get(workspaceName);
    	mc     = (ModelConfig*) w->obj(modelConfigName);
    	data   = w->data(ObsDataName);


    	w->SetName("w");
    	w->SetTitle("w");
    	// save snapshot before any fit has been done
	RooSimultaneous* pdf = (RooSimultaneous*) w->pdf("simPdf");
	RooArgSet* params = (RooArgSet*) pdf->getParameters(*data) ;
    	if(!w->loadSnapshot("snapshot_paramsVals_initial"))  w->saveSnapshot("snapshot_paramsVals_initial",*params);
    	else cout << endl << " Snapshot 'snapshot_paramsVals_initial' already exists in  workspace, will not overwrite it" << endl;
   	if(!data || !mc){
      		w->Print();
      		cout << "data or ModelConfig was not found" <<endl;
      		return;
    	}


}

 TCanvas* DrawShift(TString channel, TString var, TString comp, double mu, TH1* d, TH1* n, TH1* p1s, TH1* m1s) {
    cout << " " << comp << endl;
    cout << "N(-sigma) = " << m1s->Integral() << endl;
    cout << "N(+sigma) = " << p1s->Integral() << endl;
    cout << "N(nominal) = " << n->Integral() << endl;
    if(d) { cout << "N(Observed) = " << d->Integral() << endl; }

    var.ReplaceAll("alpha_Sys","");
    var.ReplaceAll("alpha_","");

    TString cname = "can_" + channel + "_" + comp + "_" + var + "_mu";
    cname += mu;
    cname.ReplaceAll("#","");
    cname.ReplaceAll("(","");
    cname.ReplaceAll(")","");
    cname.ReplaceAll("=","");
    TCanvas *canvas = new TCanvas(cname,cname,700,550);
    canvas->cd();
    TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
    pad1->SetBottomMargin(0.009);
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.25);
    pad2->SetTopMargin(0.009);
    pad2->SetBottomMargin(0.5);
    pad2->Draw();

    // style
    if(d) {
      d->SetLineColor(1);
      d->SetLineWidth(1);
      d->SetMarkerColor(1);
      d->SetMarkerSize(0.9);
      d->SetMarkerStyle(20);
    }
    n->SetLineWidth(2);
    p1s->SetLineColor(kRed);
    p1s->SetLineWidth(2);
    p1s->SetLineStyle(2);
    m1s->SetLineColor(kGreen);
    m1s->SetLineWidth(2);
    m1s->SetLineStyle(2);

    float max(0);
    float min(0);

    // put averages on the plot
    float avgUp = ( p1s->Integral() - n->Integral() ) / n->Integral();
    float avgDn = ( m1s->Integral() - n->Integral() ) / n->Integral();

    // Distribution in the upper pad
    pad1->cd();
    n->SetTitle(channel);
    max = p1s->GetMaximum();
    if(m1s->GetMaximum() > max) { max = m1s->GetMaximum(); }
    if(n->GetMaximum() > max) { max = n->GetMaximum(); }
    if(d) { if(d->GetMaximum() > max) { max = d->GetMaximum(); } }
    n->SetMaximum( 1.2*max );
    n->Draw("hist");
    if(d) { d->Draw("E1 same"); }
    p1s->Draw("hist same");
    m1s->Draw("hist same");

    // Distribution of the ratio in %
    pad2->cd();
    TH1F *p1s_ratio = (TH1F*) p1s->Clone();
    p1s_ratio->Add(n,-1); p1s_ratio->Divide(n); p1s_ratio->Scale(100);
    p1s_ratio->SetLineStyle(1);
    TH1F *m1s_ratio = (TH1F*) m1s->Clone();
    m1s_ratio->Add(n,-1); m1s_ratio->Divide(n); m1s_ratio->Scale(100);
    m1s_ratio->SetLineStyle(1);
    max = p1s_ratio->GetMaximum();
    if(m1s_ratio->GetMaximum() > max) { max = m1s_ratio->GetMaximum(); }
    min = p1s_ratio->GetMinimum();
    if(m1s_ratio->GetMinimum() < min) { min = m1s_ratio->GetMinimum(); }
    p1s_ratio->SetMaximum( 1.5*max );
    p1s_ratio->SetMinimum( min - 0.5*fabs(min) );
    p1s_ratio->GetYaxis()->SetNdivisions(004);
    p1s_ratio->GetXaxis()->SetTitleFont(43);
    p1s_ratio->GetXaxis()->SetTitleSize(16);
    p1s_ratio->GetXaxis()->SetTitleOffset(4);
    p1s_ratio->GetYaxis()->SetTitleOffset(1.1);
    p1s_ratio->GetYaxis()->SetTitleFont(43);
    p1s_ratio->GetYaxis()->SetTitleSize(13);
    p1s_ratio->GetXaxis()->SetLabelFont(43);
    p1s_ratio->GetXaxis()->SetLabelSize(13);
    p1s_ratio->GetYaxis()->SetLabelFont(43);
    p1s_ratio->GetYaxis()->SetLabelSize(13);
    p1s_ratio->SetTitle("");
    p1s_ratio->GetYaxis()->SetTitle("Rel. unc. (%)");
    p1s_ratio->Draw("hist");
    m1s_ratio->Draw("hist same");
    if (d){
      TH1F *d_ratio = (TH1F*) d->Clone();
      d_ratio->Add(n,-1); d_ratio->Divide(n); d_ratio->Scale(100);
      d_ratio->Draw("E1 same");
    }


    // write average shift on canvas
    pad1->cd();
    TString info(var+" "+comp);
    info.Append(Form(" %5.2f, %5.2f",avgUp*100,avgDn*100));
    info.Append('%');
    TLatex *niceinfo = new TLatex(0.12, 0.85, info);
    niceinfo->SetNDC();
    niceinfo->SetTextSize(0.045);
    niceinfo->Draw("same");

    // legend
     TLegend *leg = new TLegend(0.67, 0.64, 0.87, 0.86);
    TString varLegName(var);
    varLegName.ReplaceAll("alpha_Sys","");
    varLegName.ReplaceAll("alpha_","");
    if(d) { leg->AddEntry( d, "Data", "p" ); }
    leg->AddEntry( n, comp, "l" );
    leg->AddEntry( p1s, "+#sigma", "l" );
    leg->AddEntry( m1s, "-#sigma", "l" );
    leg->Draw();

    return canvas;
  } // DrawShift




template <typename T>
  string NumberToString ( T Number )
{
     ostringstream ss;
     ss << Number;
     return ss.str();
}

map <string,double> MapNuisanceParamNom;

void GetNominalValueNuisancePara(){
    TIterator *it = mc->GetNuisanceParameters()->createIterator();
    RooRealVar *var = NULL;
    if (MapNuisanceParamNom.size() > 0) MapNuisanceParamNom.clear();
    std::cout << "Nuisance parameter names and values" << std::endl;
    while ((var = (RooRealVar*)it->Next()) != NULL){
      const double val = var->getVal();
      MapNuisanceParamNom[(string)var->GetName()] = val;
    }
    return;
 }

bool IsAnormFactor(RooRealVar *var){
    bool result=false;
    string varname = (string) var->GetName();
    const double val =  MapNuisanceParamNom[varname];
    if (!((TString)varname).Contains("_stat_") && val==1.0) result=true;
    return result;
}

void SetPOI(double mu){
    RooRealVar * firstPOI = dynamic_cast<RooRealVar*>(mc->GetParametersOfInterest()->first());
    firstPOI->setVal(mu);
    return  ;
}

void SetNuisanceParaToSigma(RooRealVar *var, double Nsigma){

    string varname = (string) var->GetName();
    if ( varname.find("gamma_stat")!=string::npos ) return;

    if(strcmp(var->GetName(),"Lumi")==0){
      var->setVal(w->var("nominalLumi")->getVal()*(1+Nsigma*0.037));
    }
    else if (IsAnormFactor(var)) return;
    else var->setVal(Nsigma);

    return;
  }



void Plot_BG(TString wsname)
{
	//get the stuff from the workspace:
	
	TFile* file=TFile::Open(wsname);
	RooWorkspace* ws = (RooWorkspace*)file->Get("combined");
	mc = (ModelConfig*)ws->obj("ModelConfig");
	data = ws->data("obsData");
	RooSimultaneous* simPdf=(RooSimultaneous*)(mc->GetPdf());
	RooAbsReal* nll=simPdf->createNLL(*data);
	
	//run on channels
	
	RooCategory* chanCat = (RooCategory*) (&simPdf->indexCat());
        TIterator* iterat = chanCat->typeIterator() ;
        RooCatType* ttype;
	bool stop = kFALSE;
	while ((ttype = (RooCatType*) iterat->Next())&&!stop)
	{
		// bool toggle to run on one channel or all	
		stop = kTRUE;
		RooAbsPdf  *pdf_state  = simPdf->getPdf(ttype->GetName()) ;
		RooArgSet  *obstmp  = pdf_state->getObservables( *mc->GetObservables() ) ;
        	RooAbsData *datatmp = data->reduce(Form("%s==%s::%s",chanCat->GetName(),chanCat->GetName(),ttype->GetName()));
		RooRealVar *obs     = ((RooRealVar*) obstmp->first());
		TString chanName(ttype->GetName());

		// get data
		TH1* hdata = datatmp->createHistogram("Data "+chanName,*obs);
		// set errors to gaussian
        	for (int ib=0 ; ib<hdata->GetNbinsX()+1 ; ib++) hdata->SetBinError(ib, sqrt(hdata->GetBinContent(ib)));
		
		// get initial BG
		TH1* h_initial_BG = pdf_state->createHistogram("initial_BG_"+chanName,*obs);
	
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
        	minim.setPrintLevel(0);
        	minim.optimizeConst(1);
        	minim.setOffsetting(true);
        	minim.setMinimizerType("Minuit2");
        	minim.minimize("Minuit2");
        	minim.setStrategy(3); //0-3 where 0 is the fastest
        	minim.migrad();
        
		// get gammas after fit
		double FinalGamma[nbins];
		TH1* h_initBG_times_gamma = (TH1*)h_initial_BG->Clone("initBG_times_gamma");
		for (int i=0; i<nbins; i++)
        	{
                	TString varname = "gamma_B0_0j_l1pt0_bin_"+NumberToString(i);
                	FinalGamma[i] = ws->var(varname)->getVal();
                	cout << "Final gamma in bin "+NumberToString(i)+" = " << FinalGamma[i] << endl;
        		h_initBG_times_gamma->SetBinContent(i+1,h_initial_BG->GetBinContent(i+1)*FinalGamma[i]);
		}
		double FinalFpt = ws->var("fl1pt_l1pt0")->getVal();
		cout << "initial fpt_l1pt0 = " << InitFpt <<  endl;
		cout << "final fpt_l1pt0 = " << FinalFpt <<  endl;
	
		TH1* h_final_BG = pdf_state->createHistogram("final_BG_"+chanName,*obs);
        	//TCanvas* cf = new TCanvas("final BG","final BG",600,600);
		h_final_BG->Draw("sames");
		h_initBG_times_gamma->Draw("sames");
		TH1* h_ratio = (TH1*)h_initial_BG->Clone("h_ratio");
		h_ratio->Divide(h_final_BG);
		//h_ratio->Draw();
		cout << "channel name = " << chanName << endl;	
		for ( int j=1; j<=nbins; j++)
		{
			double init = h_initial_BG->GetBinContent(j);
			double fina = h_final_BG->GetBinContent(j);
			double r = (fina)/init;
			cout << "in bin " << j << ", initial B = " << init << ", final B = " << fina << ", ratio = " << r << ", Gamma = " << FinalGamma[j-1] << endl;
		}	
	}


}
