
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
	Initialize(wsname ,"combined", "ModelConfig", "obsData");
	TFile* file=TFile::Open(wsname);
	RooWorkspace* ws = (RooWorkspace*)file->Get("combined");
	//mc = (ModelConfig*)ws->obj("ModelConfig");
	//data = ws->data("obsData");
	RooSimultaneous* simPdf=(RooSimultaneous*)(mc->GetPdf());
	RooArgSet *obs = simPdf->getObservables(ws->data("obsData"));
	//obs->writeToFile("tempFile.txt");
	//RooRealVar *mu =
	//dynamic_cast<RooRealVar*>(mc->GetParametersOfInterest()->first());

	RooAbsReal* nll=simPdf->createNLL(*data);
	
	//get initial BG values
	int nbins = 15;
	double mu = 0;
	double nSigmaToVary = 1.0;
	double InitGamma[nbins];
	for (int i=0; i<nbins; i++)
	{
		TString varname = "gamma_B0_0j_l1pt0_bin_"+NumberToString(i);
		InitGamma[i] = ws->var(varname)->getVal();
		cout << "initial gamma"+NumberToString(i)+" = " << InitGamma[i] << endl;
	}
	double InitFpt = ws->var("fl1pt_l1pt0")->getVal();
	cout << "initial fpt_l1pt0 = " << InitFpt <<  endl;
	//ws->Print();
	//cout << "shape factor value = " << ws->var("emu_0j_l1pt0_B0_0j_l1pt0_shapeFactor[3]")->getVal() << endl;
	
	GetNominalValueNuisancePara();

    	RooCategory* channelCat = (RooCategory*) (&simPdf->indexCat());
    	TIterator* iter = channelCat->typeIterator() ;
    	RooCatType* tt = NULL;
    	TString dirName("trial");
    	while((tt=(RooCatType*) iter->Next()) ){

		cout << endl;
      		cout << endl;
	      	cout << " -- On category " << tt->GetName() << " " << endl;
      		ostringstream SubdirName;
	      	SubdirName << tt->GetName();
	      	gROOT->cd();
	
	      	// Get pdf associated with state from simpdf
	      	RooAbsPdf  *pdftmp  = simPdf->getPdf(tt->GetName()) ;
      		RooArgSet  *obstmp  = pdftmp->getObservables( *mc->GetObservables() ) ;
	      	RooAbsData *datatmp = data->reduce(Form("%s==%s::%s",channelCat->GetName(),channelCat->GetName(),tt->GetName()));
      		RooRealVar *obs     = ((RooRealVar*) obstmp->first());
      
		// Get the bin width
     		RooRealVar* binWidth = ((RooRealVar*) pdftmp->getVariables()->find(Form("binWidth_obs_x_%s_0",tt->GetName()))) ;
	      	if(!binWidth) { cout << "No bin width!" << tt->GetName() << endl; return; }
     		cout << "    Bin Width : " << binWidth->getVal() << endl;
      
      
		// Look at each component
	        cout << "    Contains the following components : " << endl;
	        TString modelName(tt->GetName());
	        modelName.Append("_model");
		RooRealSumPdf *pdfmodel = (RooRealSumPdf*) (pdftmp->getComponents())->find(modelName);
      		RooArgList funcList =  pdfmodel->funcList();
		RooLinkedListIter funcIter = funcList.iterator() ;
		RooProduct* comp = 0;
		float total(0);
		// want to see signal
		RooRealVar * firstPOI = dynamic_cast<RooRealVar*>(mc->GetParametersOfInterest()->first());
    		firstPOI->setVal(1);
		map<TString, TH1*> nominals;
		while( (comp = (RooProduct*) funcIter.Next())) {
        		cout << " Component : " << comp->GetName() << endl;
		        cout << "\t" << ( comp->createIntegral(*obs) )->getVal() * binWidth->getVal() << endl;
		        total += ( comp->createIntegral(*obs) )->getVal() * binWidth->getVal();
		}

	      	// Loop over nuisance params
      		TIterator* it = mc->GetNuisanceParameters()->createIterator();
      		RooRealVar* var = NULL;
      		bool IsAllStatDone = false;
                TString chanName(tt->GetName());
      	        while( (var = (RooRealVar*) it->Next()) ){
                	string varname = (string) var->GetName();
                        if ( varname.find("gamma_stat")!=string::npos ){
      				continue;
      			}
     		 
      		// one sigma not defined for floating parameters 
      		if (IsAnormFactor(var)) continue;
      		// user firendly label / name
      		TString varName(var->GetName());
      		varName.ReplaceAll("alpha_Sys","");
      		varName.ReplaceAll("alpha_","");
		
		// Not consider nuisance parameter being not assocaited to systematics
		if (MapNuisanceParamNom[varname]!=0.0 &&
		    MapNuisanceParamNom[varname]!=1.0 ) continue;
		
		cout << endl;
	        cout << "  -- On nuisance parameter : " << var->GetName() << endl;
                          		
		TString histName("");
		
		//-1 sigma
		SetNuisanceParaToSigma(var,-nSigmaToVary);
	        SetPOI(mu);
        	histName = chanName+"_"+varName+"_"+"_m1sigma";
	        TH1* hm1sigma = pdftmp->createHistogram(histName,*obs);	
	
		//+1 sigma
		SetNuisanceParaToSigma(var,+nSigmaToVary);
	        SetPOI(mu);
        	histName.ReplaceAll("m1sigma","p1sigma");
        	TH1* hp1sigma = pdftmp->createHistogram(histName,*obs);

		//nominal
		SetNuisanceParaToSigma(var,0.0);
        	SetPOI(mu);
	        histName.ReplaceAll("p1sigma","nominal");
        	TH1* hnominal = pdftmp->createHistogram(histName,*obs);			
 		
		//data
		histName.ReplaceAll("nominal","data");
	        TH1* hdata = datatmp->createHistogram(histName,*obs);
        	for (int ib=0 ; ib<hdata->GetNbinsX()+1 ; ib++) hdata->SetBinError(ib, sqrt(hdata->GetBinContent(ib)));


	        TString expName("AllBkg(#mu=");
	        expName += mu;
        	expName.Append(")");
	        TCanvas* c2 = DrawShift(chanName,(TString)var->GetName(),expName,mu,hdata,hnominal,hp1sigma,hm1sigma);
    		
	        c2->Write();
	        system(TString("mkdir -vp "+dirName));
        	c2->Print(dirName+"/totalExpected.eps");
	        c2->Print(dirName+"/totalExpected.png");
        	
	        c2->Close();
        	gROOT->cd();
		}
     	}	

	return;
	
	RooMinimizer minim(*nll);
	//set some options:
	minim.setPrintLevel(0);
	minim.optimizeConst(1);
	minim.setOffsetting(true);
	minim.setMinimizerType("Minuit2");
	minim.minimize("Minuit2");
	minim.setStrategy(3); //0-3 where 0 is the fastest
	//do the fit with various subroutines:
	minim.migrad();
	//minim.hesse();
	//minim.minos(); //this one takes long time
	
	//get 
	double unconditional_NLL = nll->getVal();
	double bestFitPval = ws->var("gamma_B0_0j_l1pt_bin5")->getVal();
	
	//RooArgSet* NP = mc->GetNuisanceParameters()
	
	int size = 20; 
	double Nll[size];
	double Pval[size];
	for (int i=0; i<size; i++){
		Pval[i] = -5 + (i+1)*10./size;
		//ws->var(varname)->setConstant(kTRUE);
		//ws->var(varname)->setVal(Pval[i]);
		//ws->var(varname)->removeError();
		//cout << "******************Parameter value = " << ws->var(varname)->getVal() << " *********************" << endl;
		//RooAbsData* data = ws->data("obsData");
        	RooSimultaneous* simPdf=(RooSimultaneous*)(mc->GetPdf());
        	RooArgSet *obs = simPdf->getObservables(ws->data("obsData"));
        	//RooRealVar *mu =
        	//dynamic_cast<RooRealVar*>(mc->GetParametersOfInterest()->first());

        	RooAbsReal* nll2=simPdf->createNLL(*data);

	}	
	
}
