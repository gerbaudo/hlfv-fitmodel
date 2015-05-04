// Author: Aaron Armbruster
// Date:   2011-11-16
//
// Description:
// Write XML, run histfactory, make asimov data, and edit and write workspace to file

#include "TFile.h"

#include "RooPoisson.h"

#include "macros/makeAsimovData.C"
#include "macros/setup.C"
#include "macros/printRand.C"
#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooRealSumPdf.h"
#include "RooNumIntConfig.h"
#include "RooProdPdf.h"
#include "RooHistFunc.h"
#include "RooDataHist.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <stdlib.h>
#include <sstream>

using namespace std;
using namespace RooFit;
using namespace RooStats;
// using namespace HistFactory;

struct WSInfo
{
  string name;
  vector<string> regionNames;
};

void joinSystematics(vector<vector<Response> >& vecResponse, vector<Response>& joined);
void writeXML2(string outFileName, double mass = 130, string version = "test", bool alt = false);
void writeXML(double mass = 130, string version = "test", bool alt = false);
RooArgList getMCMSList(RooWorkspace* w, ModelConfig* mc, RooArgSet& nuis, RooArgSet& globs, RooArgSet& obs, string passID, string failID, RooRealSumPdf* tbFail = NULL, RooRealSumPdf* tbFail_old = NULL);


void writeXML(double mass, string version, bool alt)
{
  setup(mass, alt);

  if (useLowPt && !useHighPt && (useHighMass || useHighMass2)) return;

  stringstream massStr;
  massStr << mass;
  string smass = massStr.str();

  stringstream altStr;
  if (alt) altStr << "_alt";
  string salt = altStr.str();

  RooArgSet obs;
  RooArgSet globs;
  RooArgSet nuis;
  RooArgSet poi;
  map<string, RooAbsPdf*> pdf_map;
  map<string, RooDataSet*> data_map;
  RooCategory* merged_cat = new RooCategory("merged_cat","merged_cat");


  //just do the whole damned thing in one call and be done with it
  string subfolder = "lvlv";
  writeXML2(subfolder, mass, version, alt);



  //grab the pdf, data, and sets
  TFile* file          = new TFile(("rev/"+version+"/xml/"+smass+salt+"/"+subfolder+"/output_combined_lvlv_model.root").c_str());
  RooWorkspace* w      = (RooWorkspace*)file->Get("combined");
  
  // Activate binned likelihood calculation for binned models
  RooFIter iter = w->components().fwdIterator() ;
  RooAbsArg* arg ;
  while((arg = iter.next())) {
    if (arg->IsA() == RooRealSumPdf::Class()) {
      arg->setAttribute("BinnedLikelihood");
    }
  }
  
  ModelConfig* mc_pre  = (ModelConfig*)w->obj("ModelConfig");
  RooDataSet* data     = (RooDataSet*)w->data("obsData");
  obs.add(*mc_pre->GetObservables());
      globs.add(*mc_pre->GetGlobalObservables());
      nuis.add(*mc_pre->GetNuisanceParameters());
      
      if (!useLumiAsPOI)
	{
	  // Use mu as nuisance parameter and spin0/spin2 fraction as poi
	  if (doSpin)
	    {
	      LOG(logWARNING) << "Replacing POI for spin";
	      RooRealVar* mu = (RooRealVar*)mc_pre->GetParametersOfInterest()->first();
	      mu->setVal(1);
	      mu->setConstant(0);
	      nuis.add(*mu);
	      poi.remove(*mu);
	      
          RooRealVar* epsilon = w->var("ATLAS_epsilon");
          epsilon->setRange(0, 1);
          poi.add(*w->var("ATLAS_epsilon"));
          nuis.remove(*w->var("ATLAS_epsilon"));
          nuis.remove(*w->var("ATLAS_epsilon_rejected"));
	    }
	  else if (doCRVRcheck)
	    {
	      LOG(logWARNING) << "Replacing POI for CR/VR test";
	      RooRealVar* mu = (RooRealVar*)mc_pre->GetParametersOfInterest()->first();
	      mu->setVal(1);
	      mu->setConstant(1);
	      nuis.add(*mu);
	      poi.remove(*mu);
	      
	      if (w->var("CRVR_HWW"))
		{
		  RooRealVar* alpha = w->var("CRVR_HWW");
		  alpha->setVal(1);
		  alpha->setConstant(0);
		  nuis.remove(*alpha);
		  poi.add(*alpha);
		}
	      else
		{
		  LOG(logWARNING) << "CRVR_HWW does not exist";
		}
	    } 
	  else
	    {
	      stringstream ggfMass;
	      ggfMass << 125;
	      string fixmass = ggfMass.str();
	      poi.add(*mc_pre->GetParametersOfInterest());
	      //if (doggf && !doVBF2j) poi.add(*w->var(("ATLAS_sampleNorm_ggf"+smass).c_str()));
	      //if (doggf && doVBF2j) poi.add(*w->var(("ATLAS_sampleNorm_ggf125"+fixmass).c_str()));
	      if (doggf) poi.add(*w->var(("ATLAS_sampleNorm_ggf"+smass).c_str()));
	      if (dovbf) poi.add(*w->var(("ATLAS_sampleNorm_vbf"+smass).c_str()));
	      if (dowh)  poi.add(*w->var(("ATLAS_sampleNorm_wh"+smass).c_str()));
	      if (dozh)  poi.add(*w->var(("ATLAS_sampleNorm_zh"+smass).c_str()));
	    }
	}
      else
	{
	  RooRealVar* mu = (RooRealVar*)mc_pre->GetParametersOfInterest()->first();
	  mu->setVal(1);
	  mu->setConstant(1);
	  
	  RooRealVar* lumi = w->var("Lumi");
	  lumi->setRange(0.1, 20);
	  poi.add(*w->var("Lumi"));
	}
      
      RooSimultaneous* simPdf = (RooSimultaneous*)mc_pre->GetPdf();
      RooCategory* cat = (RooCategory*)&simPdf->indexCat();
      RooCatType* tt = NULL;
      TIterator* itr = cat->typeIterator();
      
      simPdf->Print();


      if (w->obj("f_0j")) {
        RooRealVar* f_0j = w->var("f_0j");
        f_0j->setVal(ratio_S_NDY_SR_0j);
      }

      if (w->obj("f_1j")) {
        RooRealVar* f_1j = w->var("f_1j");
        f_1j->setVal(ratio_S_NDY_SR_1j);
      }

      if (w->obj("epsilon0_0j")) {
        RooRealVar* epsilon0_0j = w->var("epsilon0_0j");
        epsilon0_0j->setVal(epsilon0_data_0j);
      }

      if (w->obj("epsilon0_1j")) {
        RooRealVar* epsilon0_1j = w->var("epsilon0_1j");
        epsilon0_1j->setVal(0.62);
      }

      if (w->obj("PM_EFF_f_recoil_DY0j")) {
        RooRealVar* PM_EFF_f_recoil_DY0j = w->var("PM_EFF_f_recoil_DY0j");
        PM_EFF_f_recoil_DY0j->setVal(f_DY_all_0j);
      }

      if (w->obj("PM_EFF_f_recoil_DY1j")) {
        RooRealVar* PM_EFF_f_recoil_DY1j = w->var("PM_EFF_f_recoil_DY1j");
        PM_EFF_f_recoil_DY1j->setVal(f_DY_all_1j);
      }

      if(splitEfficiencies){
        if (w->obj("PM_EFF_f_recoil_NDY_WW0j")) {
          RooRealVar* PM_EFF_f_recoil_NDY_WW0j = w->var("PM_EFF_f_recoil_NDY_WW0j");
          PM_EFF_f_recoil_NDY_WW0j->setVal(f_NDY_WWCR_0j);
        }

        if (w->obj("PM_EFF_f_recoil_NDY_WW1j")){
          RooRealVar* PM_EFF_f_recoil_NDY_WW1j = w->var("PM_EFF_f_recoil_NDY_WW1j");
          PM_EFF_f_recoil_NDY_WW1j->setVal(f_NDY_WWCR_1j);
        }

        if (w->obj("PM_EFF_f_recoil_NDY_ZP0j")){
          RooRealVar* PM_EFF_f_recoil_NDY_ZP0j = w->var("PM_EFF_f_recoil_NDY_ZP0j");
          PM_EFF_f_recoil_NDY_ZP0j->setVal(f_NDY_ZP_0j);
        }

        if (w->obj("PM_EFF_f_recoil_NDY_ZP1j")){
          RooRealVar* PM_EFF_f_recoil_NDY_ZP1j = w->var("PM_EFF_f_recoil_NDY_ZP1j");
          PM_EFF_f_recoil_NDY_ZP1j->setVal(f_NDY_ZP_1j);
        }

        if (w->obj("PM_EFF_f_recoil_NDY_SR0j")){
          RooRealVar* PM_EFF_f_recoil_NDY_SR0j = w->var("PM_EFF_f_recoil_NDY_SR0j");
          PM_EFF_f_recoil_NDY_SR0j->setVal(f_NDY_SR_0j);
        }

        if (w->obj("PM_EFF_f_recoil_NDY_SR1j")){
          RooRealVar* PM_EFF_f_recoil_NDY_SR1j = w->var("PM_EFF_f_recoil_NDY_SR1j");
          PM_EFF_f_recoil_NDY_SR1j->setVal(f_NDY_SR_1j);
        }

      }else{

        if (w->obj("PM_EFF_f_recoil_NDY0j")){
          RooRealVar* PM_EFF_f_recoil_NDY0j = w->var("PM_EFF_f_recoil_NDY0j");
          PM_EFF_f_recoil_NDY0j->setVal(f_NDY_all);
        }

        if (w->obj("PM_EFF_f_recoil_NDY1j")){
          RooRealVar* PM_EFF_f_recoil_NDY1j = w->var("PM_EFF_f_recoil_NDY1j");
          PM_EFF_f_recoil_NDY1j->setVal(f_NDY_all);
        }
      }
      /*      
      RooRealVar* alpha_PM_theta_SR0j = NULL;
      bool has_alpha_PM_theta_SR0j = false;
      if (w->obj("alpha_PM_theta_SR0j")){
	alpha_PM_theta_SR0j = w->var("alpha_PM_theta_SR0j");
	has_alpha_PM_theta_SR0j=true;
	nuis.add(*w->var(("alpha_PM_theta_SR0j")));
      }
      
      RooRealVar* alpha_PM_theta_SR1j = NULL;
      bool has_alpha_PM_theta_SR1j = false;
      if (w->obj("alpha_PM_theta_SR1j")){
        alpha_PM_theta_SR1j = w->var("alpha_PM_theta_SR1j");
	has_alpha_PM_theta_SR1j = true;
	nuis.add(*w->var(("alpha_PM_theta_SR1j")));
      }
      
      */

      //get terms used in modified CMS method
      vector<RooArgList> mcmsLists;
      if (doModifiedCMS)
	{

          if (do1j || (do2j && cancel2j))
            {
              //get MCMS CR
              mcmsLists.push_back(getMCMSList(w, mc_pre, nuis, globs, obs, "OF_tbPass_2j", "OF_tbFail_2j"));
	      
              //get other mCMS regions
              for (int ir=0;ir<(int)regions->size();ir++)
                {
                  Region* r = &(*regions)[ir];
                  int nrChannels = r->channels.size();
                  for (int ic=0;ic<nrChannels;ic++)
                    {

                      string rname = r->name;
                      string clonedName = r->clonedName;
                      if (clonedName == "") continue; // skip over non-mcms regions
                      if (rname.find("tbPass") != string::npos || rname.find("tbFail") != string::npos) continue; // these have already been done
                      rname = r->channels[ic].name +"_" +rname + "_"+ r->channels[ic].jetName;
                      clonedName = r->channels[ic].name + "_"+clonedName +"_"+ r->channels[ic].jetName;
		      mcmsLists.push_back(getMCMSList(w, mc_pre, nuis, globs, obs, rname,  clonedName));
                    }
                }
            }

	  /*
	  //NEED TO CHECK IF THIS STILL WORKS PROPERLY
	  if (useAgnosticSignalRegionNaming)
	    for (unsigned int i_mll = 0; i_mll < parsedmllbounds.size()-1; i_mll++) 
	      for (unsigned int i_subleadpt = 0; i_subleadpt < parsedsubleadptbounds.size()-1; i_subleadpt++) {
		TString srName = TString::Format("em_signalLike%i%c_1j", i_mll, (char)(97+i_subleadpt));
		TString srTagName = TString::Format("em_sr%i%cTag_1j", i_mll, (char)(97+i_subleadpt));
		mcmsLists.push_back(getMCMSList(w, mc_pre, nuis, globs, obs, srTagName.Data(), srName.Data()));
		srName = TString::Format("me_signalLike%i%c_1j", i_mll, (char)(97+i_subleadpt));
		srTagName = TString::Format("me_sr%i%cTag_1j", i_mll, (char)(97+i_subleadpt));
		mcmsLists.push_back(getMCMSList(w, mc_pre, nuis, globs, obs, srTagName.Data(), srName.Data()));
	      }
	  */
   	  
	}//end of modified cms
      
      //Edit ATLAS_norm_btag to include nuis param ATLAS_btag21j_extrap
      RooFormulaVar* kappa_norm_btag  = NULL;
      bool has_ATLAS_norm_btag=false;
      if(w->var("ATLAS_norm_btag")){
        has_ATLAS_norm_btag=true;
        double uncertainty = 0.06; // get this from a config file maybe?
        stringstream formula;
        formula << "@0*pow(1+" << uncertainty << ",@1)";
	
        RooRealVar* the_nuis_param = w->var("ATLAS_btag21j_extrap");
        RooRealVar* norm_btag_temp = w->var("ATLAS_norm_btag");
	
        RooArgList formList;
        formList.add(*norm_btag_temp);
        formList.add(*the_nuis_param);

        //if this is only one uncertainty value, then no need to specify givePDFName and it can go outside of loop and import there
        kappa_norm_btag = new RooFormulaVar("kappa_ATLAS_norm_btag","kappa_ATLAS_norm_btag",formula.str().c_str(),formList);
	w->import(*kappa_norm_btag);
      }
      
      obs.remove(*cat);
     
      while ((tt = (RooCatType*)itr->Next()))//categories are em_SignalLike_1j for example, contain PDFs for each category
	{
	  RooAbsPdf* pdf = simPdf->getPdf(tt->GetName());
	  
	  LOG(logDEBUG) << "tt is " << tt->GetName();
	  
	  

	  if (w->obj("alpha_PM_theta_SR0j")|| w->obj("alpha_PM_theta_SR1j"))
            {
              stringstream editStr;
              editStr << "EDIT::" << pdf->GetName() << "_edit(" << pdf->GetName();

              if(w->obj("alpha_PM_theta_SR0j") && pdf->dependsOn(*w->var("alpha_PM_theta_temp_SR0j"))) editStr << ",alpha_PM_theta_temp_SR0j=alpha_PM_theta_SR0j";
              if(w->obj("alpha_PM_theta_SR1j") && pdf->dependsOn(*w->var("alpha_PM_theta_temp_SR1j"))) editStr << ",alpha_PM_theta_temp_SR1j=alpha_PM_theta_SR1j";

              editStr << ")";
              cout<< " NINA EDIT alpha" << editStr.str().c_str() <<endl;
              LOG(logINFO) << "editStr is " << editStr.str().c_str();
              if (editStr.str().find("=") != string::npos)
                {
                  w->factory(editStr.str().c_str());
                  pdf = w->pdf((string(pdf->GetName())+"_edit").c_str());
                  if (!pdf)
                    {
                      LOG(logERROR) << "Something went wrong in variable edit for alpha";
                      LOG(logERROR) << "editStr was " << editStr.str().c_str();
                      exit(1);
                    }
                }
            }







        // Do the actual workspace editing
	  if(mcmsLists.size())
	    {
	      stringstream editStr;
	      editStr << "EDIT::" << pdf->GetName() << "_edit(" << pdf->GetName();

	      for (int il=0;il<(int)mcmsLists.size();il++)
		{
		  if (mcmsLists[il].getSize()    && pdf->dependsOn(*mcmsLists[il].at(0))) editStr << "," << mcmsLists[il].at(0)->GetName() << "=" << mcmsLists[il].at(1)->GetName();
		}
 
	      //if(has_ATLAS_norm_btag && pdf->dependsOn(*w->var("ATLAS_norm_btag"))) editStr << ",ATLAS_norm_btag=kappa_ATLAS_norm_btag" ;
	      editStr << ")";
	      cout<< " NINA EDIT mcms" << editStr.str().c_str() <<endl;	      
	      LOG(logINFO) << "editStr is " << editStr.str().c_str();
	      if (editStr.str().find("=") != string::npos)
		{
		  w->factory(editStr.str().c_str());
		  pdf = w->pdf((string(pdf->GetName())+"_edit").c_str());
		  if (!pdf)
		    {
		      LOG(logERROR) << "Something went wrong in variable edit for mcms";
		      LOG(logERROR) << "editStr was " << editStr.str().c_str();
		      exit(1);
		    }
		}
	    }

	  //skip the inclusion of the region in the final pdf it's just a temporary MCMS CR
	  bool skip = false;
	  for (int i=0;i<(int)regions->size();i++)
	    {
	      if (!(*regions)[i].isMCMSCR) continue;
	      if (string(pdf->GetName()).find((*regions)[i].name) != string::npos) 
		{
		  skip = true;
		  break;
		}
	    }
	  if (skip) continue;
	  
	  merged_cat->defineType(tt->GetName());
	  pdf_map[tt->GetName()] = pdf;
	  
	  TList* data_list = doData ? data->split(*cat) : NULL;
	  int nrData = doData ? data_list->GetEntries() : 0;

	  bool found = false;
	  for (int id = 0; id < nrData; id++)
	    {
	      RooDataSet* thisData = (RooDataSet*)data_list->At(id);
	      if (string(thisData->GetName()).find(string(tt->GetName())) != string::npos)
		{
		  data_map[tt->GetName()] = thisData;
		  found = true;
		  break;
		}
	    }
	  if (!found && doData)
	    {
	      LOG(logERROR) << "Couldn't find data for cat " << tt->GetName();
	      exit(1);
	    }
	}
      
      
      //make sure the proper constraints exist in each prodpdf after inserting mcms histos
      cout << "getting observables" << endl;
      RooArgSet obsSet = *mc_pre->GetObservables();
      cout << "getting constraints" << endl;
      RooArgSet constraints = *simPdf->getAllConstraints(obsSet,nuis,false);
      cout << "resetting itr" << endl;
      delete itr;
      itr = merged_cat->typeIterator();
      while ((tt = (RooCatType*)itr->Next()))
	{
	RooProdPdf* prod = (RooProdPdf*)pdf_map[tt->GetName()];
	RooArgSet firstConstraints = *prod->getAllConstraints(obsSet,nuis,false);
	cout << "prod = " << prod << ", tt = " << tt->GetName() << endl;
	//get the data pdf
	RooArgList pdfs = prod->pdfList();
	pdfs.remove(constraints);
	if (pdfs.getSize() != 1)
	  {
	    cout << "ERROR in removing constraints in prod " << prod->GetName() << endl;
	    exit(1);
	  }
	RooRealSumPdf* data_pdf = (RooRealSumPdf*)pdfs.at(0);
	cout << "data_pdf = " << data_pdf << endl;
	
	RooArgSet theseNuis;
	TIterator* nitr = nuis.createIterator();
	RooRealVar* nui;
	while ((nui = (RooRealVar*)nitr->Next()))
	  {
	    if (data_pdf->dependsOn(*nui)) theseNuis.add(*nui);
	  }
	
	RooArgSet theseConstraints;
	cout << "starting citr loop" << endl;
	TIterator* citr = constraints.createIterator();
	RooAbsArg* arg;
	while ((arg = (RooAbsArg*)citr->Next()))
	{	
	  delete nitr;
	  nitr = theseNuis.createIterator();
	  //cout << "starting nitr loop" << endl;
	  while ((nui = (RooRealVar*)nitr->Next()))
	  {
	    if (arg->dependsOn(*nui))
	    {
	      theseConstraints.add(*arg);
	      break;
	    }
	  }
	  //cout << "end nitr loop" << endl;
	}

	RooArgSet theseConstraints2 = theseConstraints;
	theseConstraints2.remove(firstConstraints);
	if (theseConstraints2.getSize())
	{
	  cout << "Adding " << theseConstraints2.getSize() << " constraints to pdf " << prod->GetName() << ": " << endl;
	  TIterator* citr2 = theseConstraints2.createIterator();
	  RooAbsArg* arg2;
	  while ((arg2 = (RooAbsArg*)citr2->Next()))
	  {
	    cout << "-->" << arg2->GetName() << endl;
	  }
	  cout << endl;
	}

	cout << "adding data_pdf" << endl;
	theseConstraints.add(*data_pdf);
	cout << "creating newProd" << endl;
	RooProdPdf* newProd = new RooProdPdf((string(prod->GetName())+"_mod").c_str(),(string(prod->GetName())+"_mod").c_str(),theseConstraints);
	pdf_map[tt->GetName()] = newProd;
      }
      

  RooRealVar* weightVar = (RooRealVar*)obs.find("weightVar");
  obs.add(*merged_cat);
  if (!weightVar)
  {
    LOG(logERROR) << "Couldn't find weightVar in any of the input files.";
    exit(1);
  }

  // Merge models
  LOG(logINFO) << "Making pdf";
  RooSimultaneous* merged_pdf = new RooSimultaneous("mergedPdf","mergedPdf", pdf_map, *merged_cat);
  merged_pdf->Print();

  LOG(logINFO) << "Making data";
  RooDataSet* merged_data = doData ? new RooDataSet("obsData","obsData",obs,Index(*merged_cat), Import(data_map), WeightVar(*weightVar)) : NULL;

  // Get rid of constant parameters
  vector<string> constParams;
  constParams.push_back("Z_scaleF_ee0j");
  constParams.push_back("Z_scaleF_mm0j");
  constParams.push_back("Z_scaleF_em0j");
  constParams.push_back("Z_scaleF_ee1j");
  constParams.push_back("Z_scaleF_mm1j");
  constParams.push_back("Z_scaleF_ee2j");
  constParams.push_back("Z_scaleF_mm2j");
  constParams.push_back("Top_scaleF_0j");
  constParams.push_back("Lumi");
  constParams.push_back("ATLAS_norm_btag_squared");
  constParams.push_back("twice_ATLAS_norm_btag");
  constParams.push_back("f_0j");
  constParams.push_back("f_1j");
  constParams.push_back("epsilon0_0j");
  constParams.push_back("epsilon0_1j");

  int nrSamples = samples->size();
  for (int i = 0; i < nrSamples; i++)
  {
    constParams.push_back("ATLAS_sampleNorm_"+(*samples)[i].name);
  }
  constParams.push_back("mu_BR_WW");
  constParams.push_back("PM_EFF_f_recoil_rejected_DY0j_func");
  constParams.push_back("PM_EFF_f_recoil_rejected_DY1j_func");
  if (!splitEfficiencies)
  {
    constParams.push_back("PM_EFF_f_recoil_rejected_NDY0j_func");
    constParams.push_back("PM_EFF_f_recoil_rejected_NDY1j_func");
  }
  else
  {
    constParams.push_back("PM_EFF_f_recoil_rejected_NDY_SR0j_func");
    constParams.push_back("PM_EFF_f_recoil_rejected_NDY_ZP0j_func");
    constParams.push_back("PM_EFF_f_recoil_rejected_NDY_WW0j_func");
    constParams.push_back("PM_EFF_f_recoil_rejected_NDY_SR1j_func");
    constParams.push_back("PM_EFF_f_recoil_rejected_NDY_ZP1j_func");
    constParams.push_back("PM_EFF_f_recoil_rejected_NDY_WW1j_func");
    constParams.push_back("PM_EFF_f_recoil_rejected_NDY_SR0j");
    constParams.push_back("PM_EFF_f_recoil_rejected_NDY_SR1j");
  }
  constParams.push_back("alpha_PM_theta_temp_SR0j");
  constParams.push_back("alpha_PM_theta_temp_SR1j");
  //constParams.push_back("PM_EFF_f_recoil_NDY_S0j");
  //constParams.push_back("PM_EFF_f_recoil_rejected_NDY_S0j");
  //constParams.push_back("PM_EFF_f_recoil_NDY_S1j");
  //constParams.push_back("PM_EFF_f_recoil_rejected_NDY_S1j");
  // constParams.push_back("ATLAS_epsilon_rejected");
  // constParams.push_back("ATLAS_epsilon");
  // constParams.push_back("ATLAS_norm_WW0j");
  // constParams.push_back("ATLAS_norm_WW1j");
  // constParams.push_back("ATLAS_norm_Top1j");
  // constParams.push_back(("alpha_PM_theta_SR"+s_jj+"j").c_str());
  int nrConst = constParams.size();

  LOG(logWARNING) << "Removing constant nuisance parameters.";
  for (int i=0;i<nrConst;i++)
  {
    RooRealVar* param = (RooRealVar*)nuis.find(constParams[i].c_str());
    if (param)
    {
      param->setConstant(1);
      LOG(logWARNING) << "Setting param to const: " << constParams[i];

      nuis.remove(*param);
    }
  }




  globs.remove(*globs.find("nominalLumi"));
  nuis.sort();
  globs.sort();

  LOG(logINFO) << "Building workspaces and model config";
  RooWorkspace merged_ws("combined");
  ModelConfig mc("ModelConfig",&merged_ws);
  mc.SetPdf(*merged_pdf);
  mc.SetNuisanceParameters(nuis);
  mc.SetGlobalObservables(globs);
  mc.SetObservables(obs);
  mc.SetParametersOfInterest(poi);

  merged_ws.import(mc);
  if (doData) merged_ws.import(*merged_data);

  // For some reason these aren't being set to const in the above lines
  for (int i = 0; i < nrConst; i++)
  {
    RooRealVar* param = (RooRealVar*)merged_ws.var(constParams[i].c_str());
    if (param)
    {
      param->setConstant(1);
      LOG(logWARNING) << "Setting param to const again: " << constParams[i];
    }
  }

  ModelConfig* mcInWs = (ModelConfig*)merged_ws.obj("ModelConfig");

  // Bugfix for w->factory noRounding not working in histfactory
  TIterator* nItr = mcInWs->GetNuisanceParameters()->createIterator();
  RooRealVar* var;
  while ((var = (RooRealVar*)nItr->Next()))
  {
    string poisName = string(var->GetName()) + "_constraint";
    if (poisName.find("gamma_stat") == string::npos) continue;
    RooPoisson* pois = (RooPoisson*)merged_ws.pdf(poisName.c_str());
    if (!pois)
    {
      LOG(logERROR) << "ERROR::Couldn't find corresponding poisson for var " << var->GetName();
      exit(1);
    }
    pois->setNoRounding(true);
  }

  // Make asimov data
  RooArgSet funcs = merged_ws.allFunctions();
  TIterator* it = funcs.createIterator();
  TObject* tempObj = 0;
  while((tempObj=it->Next()))
  {
    HistFactory::FlexibleInterpVar* flex = dynamic_cast<HistFactory::FlexibleInterpVar*>(tempObj);
    if(flex)
    {
      flex->setAllInterpCodes(flatInterpCode);
    }
    PiecewiseInterpolation* piece = dynamic_cast<PiecewiseInterpolation*>(tempObj);
    if(piece)
    {
      piece->setAllInterpCodes(shapeInterpCode);
    }
  }

  LOG(logINFO) << "Making asimov data";
  RooDataSet* dataInWs = (RooDataSet*)merged_ws.data("obsData");
  makeAsimovData(mcInWs, conditionalAsimov && doData, &merged_ws, mcInWs->GetPdf(), dataInWs, 0);
  makeAsimovData(mcInWs, conditionalAsimov && doData, &merged_ws, mcInWs->GetPdf(), dataInWs, 1);
  makeAsimovData(mcInWs, conditionalAsimov && doData, &merged_ws, mcInWs->GetPdf(), dataInWs, 2);

  // merged_ws.Print();

  LOG(logINFO) << "Writing to file";
  system(("mkdir -vp workspaces/"+version).c_str());
  merged_ws.writeToFile(("workspaces/"+version+"/"+smass+salt+".root").c_str(), true); // write

  merged_ws.Print("t");

  printRand();
}




//
RooArgList getMCMSList(RooWorkspace* w, ModelConfig* mc, RooArgSet& nuis, RooArgSet& globs, RooArgSet& obs, string passID, string failID, RooRealSumPdf* tbFail, RooRealSumPdf* tbFail_old)
{
  RooArgList list;
  if (!doModifiedCMS) return list;

  RooSimultaneous* simPdf = (RooSimultaneous*)mc->GetPdf();
  RooCategory* cat = (RooCategory*)&simPdf->indexCat();
  RooCatType* tt = NULL;
  TIterator* itr = cat->typeIterator();


  RooProdPdf* tbFail_prod = NULL;
  //RooRealSumPdf* tbFail = NULL;
  RooRealSumPdf* tbFail_sum_new = NULL;

  cout << "Starting modification of terms for mcms method: " << passID << " -- " << failID << endl;
  //Insert tbPass top histograms into tbFail RooRealSumPdf (needs to be done before editing, since editing has to be done differently for tbPass and tbFail)
  RooProdPdf* tbPass_prod = NULL;

  while ((tt = (RooCatType*)itr->Next()))
  {
    RooAbsPdf* pdf = simPdf->getPdf(tt->GetName());
    cout << "tt = " << tt->GetName() << ", pdf = " << pdf << endl;
    if (string(pdf->GetName()).find(passID) != string::npos) tbPass_prod = (RooProdPdf*)pdf;
    if (string(pdf->GetName()).find(failID) != string::npos) tbFail_prod = (RooProdPdf*)pdf;
  }
  itr->Reset();

  if (tbPass_prod && tbFail_prod)
  {
    RooRealSumPdf* tbPass = NULL;
    RooArgList passList = tbPass_prod->pdfList();
    for (int i=0;i<(int)passList.getSize();i++)
    {
      RooAbsPdf* pdf = (RooAbsPdf*)passList.at(i);
      if (string(pdf->GetName()).find(passID) != string::npos) tbPass = (RooRealSumPdf*)pdf;
    }

    if (!tbFail)
    {
      RooArgList failList = tbFail_prod->pdfList();
      for (int i=0;i<(int)failList.getSize();i++)
      {
	RooAbsPdf* pdf = (RooAbsPdf*)failList.at(i);
	if (string(pdf->GetName()).find(failID) != string::npos) tbFail = (RooRealSumPdf*)pdf;
      }
    }


    tbPass->Print();
    RooArgList funcList = tbPass->funcList();
    RooArgList coefList = tbPass->coefList();

    RooAbsReal* ttbar_pass = NULL;
    RooAbsReal* ttbar_coef = NULL;
    RooAbsReal* st_pass = NULL;
    RooAbsReal* st_coef = NULL;
	
    cout << "Grabbing hist & coef" << endl;
    for (int i=0;i<(int)funcList.getSize();i++)
    {
      if (string(funcList.at(i)->GetName()).find("L_x_ttbar") != string::npos) 
      {
	ttbar_pass = (RooAbsReal*)funcList.at(i);
	ttbar_coef = (RooAbsReal*)coefList.at(i);
      }
      if (string(funcList.at(i)->GetName()).find("L_x_st") != string::npos) 
      {
	st_pass = (RooAbsReal*)funcList.at(i);
	st_coef = (RooAbsReal*)coefList.at(i);
      }
    }

    if (!ttbar_pass || !st_pass)
    {
      LOG(logERROR) << "ERROR::top hist can't be found in tb pass region (" << ttbar_pass << " " << st_pass << ")";
      exit(1);
    }

    //inject the histogram into the fail region, multiplying by (1-mu_eps) / mu_eps. (1/mu_eps cancels, leaving a 1-mu_eps coef)

    cout << "creating coefs" << endl;

    RooRealVar* norm_btag = w->var("ATLAS_norm_btag");
    if (!norm_btag)
    {
      LOG(logERROR) << "ERROR::ATLAS_norm_btag not in workspace";
      exit(1);
    }
    //create the coefficients
// 	  RooFormulaVar* ttbar_coef_fail = (RooFormulaVar*)w->function("ttbar_coef_fail");
// 	  RooFormulaVar* st_coef_fail = (RooFormulaVar*)w->function("st_coef_fail");
// 	  if (!ttbar_coef_fail)
// 	  {

// 	    RooArgList stList;
// 	    stList.add(*norm_btag);
// 	    stList.add(*st_coef);
// 	    st_coef_fail = new RooFormulaVar("st_coef_fail","st_coef_fail","(1-@0)/@0*@1",stList);
// // 	  w->import(*st_coef_fail);
// // 	  st_coef_fail = (RooFormulaVar*)w->function("st_coef_fail");

// 	    RooArgList ttbarList;
// 	    ttbarList.add(*norm_btag);
// 	    ttbarList.add(*ttbar_coef);
// 	    ttbar_coef_fail = new RooFormulaVar("ttbar_coef_fail","ttbar_coef_fail","(1-@0)/@0*@1",ttbarList);
// // 	  w->import(*ttbar_coef_fail);
// // 	  ttbar_coef_fail = (RooFormulaVar*)w->function("ttbar_coef_fail");
// 	  }

    cout << "creating roorealsumpdf" << endl;
    //create the new RooRealSumPdf for the fail term
    RooArgList funcList_fail = tbFail->funcList();
    RooArgList coefList_fail = tbFail->coefList();
    funcList_fail.add(*ttbar_pass);
    coefList_fail.add(*ttbar_coef);

    funcList_fail.add(*st_pass);
    coefList_fail.add(*st_coef);

    cout << "calling ctor" << endl;
    //use same name; will overwrite in workspace (no use for previous one). tbFail is pointer to obj in workspace, while new_sum is pointer to different obj in memory
    RooRealSumPdf* new_sum = new RooRealSumPdf((string(tbFail->GetName())+"_new").c_str(), (string(tbFail->GetName())+"_new").c_str(), funcList_fail, coefList_fail, kTRUE); 

// 	  //optimizations (stolen from histfactory code)
    new_sum->specialIntegratorConfig(kTRUE)->method1D().setLabel("RooBinIntegrator")  ;
    new_sum->specialIntegratorConfig(kTRUE)->method2D().setLabel("RooBinIntegrator")  ;
    new_sum->specialIntegratorConfig(kTRUE)->methodND().setLabel("RooBinIntegrator")  ;
    new_sum->forceNumInt();

// 	  // for mixed generation in RooSimultaneous
    new_sum->setAttribute("GenerateBinned"); // for use with RooSimultaneous::generate in mixed mode

    RooArgList formList(*norm_btag);
    RooFormulaVar one_minus("one_minus_ATLAS_norm_btag","one_minus_ATLAS_norm_btag","1-@0",formList);//the @ means the names it expects are not hard coded
    RooFormulaVar twice_one_minus("twice_one_minus_ATLAS_norm_btag","twice_one_minus_ATLAS_norm_btag","2*(1-@0)",formList);
    RooFormulaVar one_minus_squared("one_minus_ATLAS_norm_btag_squared","one_minus_ATLAS_norm_btag_squared","(1-@0)*(1-@0)",formList);

    //import & replace
    cout << "importing" << endl;
    w->import(*new_sum);//RecycleConflictNodes());
    w->import(one_minus);
    w->import(twice_one_minus);
    w->import(one_minus_squared);

    RooRealVar One("One","One",1);
    w->import(One);


    RooRealVar* pass_obs = (RooRealVar*)obs.find(string("obs_x_"+passID+"_"+(do2012 ? "2012":"2011")).c_str());
    RooRealVar* fail_obs = (RooRealVar*)obs.find(string("obs_x_"+failID+"_"+(do2012 ? "2012":"2011")).c_str());
    stringstream editStr;
    //editStr << "EDIT::" << new_sum->GetName() << "_edit(" << new_sum->GetName() << "," << norm_btag->GetName() << "=" << one_minus.GetName() << ")";
    editStr << "EDIT::" << new_sum->GetName() << "_edit(" << new_sum->GetName() << "," << pass_obs->GetName() << "=" << fail_obs->GetName();
    if (passID.find("tbTag") != string::npos)
    {
      editStr << "," << /*norm_btag->GetName()*/"twice_ATLAS_norm_btag" << "=" << twice_one_minus.GetName();
    }
    else if (useFullTagSample && failID.find("tbFail") != string::npos)
    {
      editStr << "," << norm_btag->GetName() << "=" << twice_one_minus.GetName();
    }
    else if (passID.find("x2Tag") != string::npos)
    {
      editStr << "," << /*norm_btag_squared->GetName()*/"ATLAS_norm_btag_squared" << "=" << one_minus_squared.GetName();
    }
    else
    {
      editStr << "," << norm_btag->GetName() << "=" << one_minus.GetName();
    }

    if (failID.find("tbFail") != string::npos) editStr << ")";
    else 
    {
      if (useStatSys)
      {
	editStr << ",mc_stat_" << passID + "_" + (do2012 ? "2012":"2011")<<"=One";
	int counter = 0;
	RooRealVar* nui;
	nuis.Print("v");
	while ((nui = (RooRealVar*)w->var(string("gamma_stat_"+passID+"_"+(do2012 ? "2012":"2011")+"_bin_"+itos(counter)).c_str())))
	{
	  cout << "Removing nui: " << nui->GetName() << endl;
	  nuis.remove(*nui);
	  globs.remove(*w->var(string("nom_gamma_stat_"+passID+"_"+(do2012 ? "2012":"2011")+"_bin_"+itos(counter)).c_str()));
	  counter++;
	}
      }
      if (pass_obs)
      {
	obs.remove(*pass_obs);
      }
      else
      {
	LOG(logERROR) << "Can't find observable to remove: " << passID;
	exit(1);
      }
      editStr << ")";
    }
    cout << "Edit string is " << endl;
    cout << editStr.str() << endl;
    w->factory(editStr.str().c_str());


    //cout << "failed ? " << failed << endl;
    new_sum->Print();
    tbFail_sum_new = (RooRealSumPdf*)w->pdf((string(new_sum->GetName())+"_edit").c_str());
    tbFail_sum_new->Print();
    //exit(1);
//make sure the observable in the dhist of the histfunc is correct
    RooArgSet* components = tbFail_sum_new->getComponents();
    TIterator* compItr = components->createIterator();
    RooAbsArg* arg;
    while ((arg = (RooAbsArg*)compItr->Next()))
    {
      if (string(arg->ClassName()) != "RooHistFunc") continue;
      RooHistFunc* hist = (RooHistFunc*)arg;
      hist->dataHist().changeObservableName(pass_obs->GetName(), fail_obs->GetName());
    }


    if (!tbFail_sum_new)
    {
      LOG(logERROR) << "ERROR::Problem retrieving new tbFail sum";
      exit(1);
    }
    if (tbFail_old) list.add(*tbFail_old);
    else list.add(*tbFail);
    list.add(*tbFail_sum_new);

    //return;
  }
  else
  {
    LOG(logERROR) << "ERROR::Can't find pass or fail prod! " << passID << ": " << tbPass_prod << ", " << failID << ": " << tbFail_prod;
    exit(1);
  }
  return list;
}

// writeXML2
void writeXML2(string subfolder, double mass, string version, bool alt)
{
  LOG(logINFO) << "Writing XML: " << mass << " " << version << " alt?" << alt;

  setup(mass, alt);

  stringstream massStr;
  massStr << mass;
  string smass = massStr.str();

  stringstream altStr;
  if (alt) altStr << "_alt";
  string salt = altStr.str();

  // if (version.find("Pacman") != string::npos && (!doee || !domm)) return;

  if (rewriteXML)
  {
    if (!overriden) // baseline config
    {
      useDetSys    = 1;
      useShape     = 1 && useDetSys;
      useThSys     = 1;
      useStatSys   = 1;
      useTheoryWW  = 0;
      useTheoryTop = 0;
      scaleZ       = 1;
      scaleTop     = 1;
    }

//     vector<string> channelNames;
//     if (doee) channelNames.push_back("ee");
//     if (doem) channelNames.push_back("em");
//     if (splitem && dome) channelNames.push_back("me");
//     if (domm) channelNames.push_back("mm");
//     if (combineCRs) channelNames.push_back("OF");
//     if (combineSFCRs) channelNames.push_back("SF");
//     int nrChannels = channelNames.size();

//     vector<string> jetNames;
//     if (do0j) jetNames.push_back("0j");
//     if (do1j) jetNames.push_back("1j");
//     if (do2j) jetNames.push_back("2j");
//     int nrJets = jetNames.size();

    int nrRegions = regions->size();
    int nrSamples = samples->size();
    int nrSys = fileNames->size();

    // load shape variations
    vector<Response> detR;
    if (useShape)
    {
      for (int isys = 0; isys < nrSys; isys++)
      {
        Sys* s = &(*fileNames)[isys];
        if (s->folder == "Nominal") continue;
        if (!s->isShape) continue;

        for (int isam = 0; isam < nrSamples; isam++)
        {
          // if ((*samples)[isam].name == "wjets") continue;

          if (s->sampleNames.find((*samples)[isam].name) == s->sampleNames.end()) continue;

          string hiFileName = "rev/"+version+"/normHists/" + s->folder + "/" + s->fileUp   + "/" + (*samples)[isam].name + ".root";
          string loFileName = "rev/"+version+"/normHists/" + s->folder + "/" + s->fileDown + "/" + (*samples)[isam].name + ".root";

	  for (int ir = 0; ir < nrRegions; ir++)
	  {
	    Region* r = &(*regions)[ir];
	    int nrChannels = r->channels.size();
	    
	    for (int ichan = 0; ichan < nrChannels; ichan++)
	    {
	      Channel* c = &r->channels[ichan];
	      string jetName = c->jetName;
	      string channelName = c->name;
	      
	      if (s->veto.find((*samples)[isam].name+"_"+channelName+"_"+(*regions)[ir].name+"_"+jetName) != s->veto.end()) continue;
	      
	      string name = s->folder;
	      
	      Response res;
	      res.name = name;
	      res.channel = channelName;
	      res.jet = jetName;
	      res.sample = (*samples)[isam].name;
	      res.region = (*regions)[ir].name;
	      res.hiFileName = hiFileName;
	      res.loFileName = loFileName;
	      res.hiHistName = channelName+"_"+(*regions)[ir].name+"_"+jetName;
	      res.loHistName = channelName+"_"+(*regions)[ir].name+"_"+jetName;
	      detR.push_back(res);
              
            }
          }
        }
      }
    }
    
    // Load flat systematics
    vector<Response> normR;
    if (useDetSys)
      {
	readInUncerts("rev/"+version+"/normHists/norms_"+smass+salt+".txt", normR, mass, 0, 0, bool(ZMode!=0?1:0));
      }
    
    // Load flat systematics for interpolated signal (was done at later stage)
    bool found = false;
    
    int nrBasePoints = 43;
    double* baseMassPoints = new double[nrBasePoints];
    double thisMass = 90;
    double step = 5;
    for (int i=0;i<nrBasePoints;i++)
    {
      baseMassPoints[i] = thisMass;
      if (baseMassPoints[i] == mass)
      {
        found = true;
        break;
      }
      if (thisMass >= 200) step = 20;

      thisMass += step;
    }

    if (!found) readInUncerts("rev/"+version+"/normHists/norms_"+smass+salt+"_sig.txt", normR, mass);

    // Load theory uncertainties
    vector<Response> thR;
    if (useThSys)
    {
      if (mode < 7) {
        if (!doWWCR2j) {
          readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"theory_constraints_other.txt", thR, mass, 0, 0, bool(ZMode!=0?1:0));
        } else if (doWWCR2j) {
          readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"theory_constraints_other_noww.txt", thR, mass, 0, 0, bool(ZMode!=0?1:0));
        }
      } else if (mode == 8) {
        readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"theory_constraints_other_3D.txt", thR, mass, 0, 0, bool(ZMode!=0?1:0));
      }

      if (mode == 0 && !doABCD2j) // FIXME: do we still need SF method?
      {
        readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"zjets_2j_cuts.txt", thR, mass, 0, 0);
      }
      else if(!doABCD2j)
      {
        readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"zjets_2j.txt", thR, mass, 0, 0);
      }

      if (useTheoryWW) readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"theory_constraints_ww.txt", thR, mass, 0, 0); // FIXME
      else if (mode != 7) 
      {
        if (useAltCRs) // FIXME: default CRs now, get rid of old definition
        {
          if (!doCRVRcheck)
          {
            if (skipSFWWCR) {

	      if (mode == 8){
                if(!useAgnosticSignalRegionNaming){
                  readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"alpha_uncerts_ww_alt_3D_skipSF.txt", thR, mass, 0, 0);
                }else{
                  readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"alpha_uncerts_ww_alt_3D_skipSF_AgnosticNaming.txt", thR, mass, 0, 0);
                }
              }else{
                readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"alpha_uncerts_ww_alt_skipSF.txt", thR, mass, 0, 0);
              }
            }
            else {
              readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"alpha_uncerts_ww_alt.txt", thR, mass, 0, 0);
            }
          }
          else readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"alpha_uncerts_ww_alt_CRVR.txt", thR, mass, 0, 0);
        }
        else readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"alpha_uncerts_ww.txt", thR, mass, 0, 0);
      }
      else if (mode == 7)
      {
        readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"mll_ratios.txt", thR, mass, 0, 0);
      }

      if (useTheoryTop) readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"theory_constraints_top.txt", thR, mass, 0, 0);
      else
      {
        // make sure that 2j numbers in both files
        if (skipSFtopCR) {
          if (mode == 8) readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"alpha_uncerts_top_3D_skipSF.txt", thR, mass, 0, 0);
          else readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"alpha_uncerts_top_skipSF.txt", thR, mass, 0, 0);
        }
        else readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"alpha_uncerts_top.txt", thR, mass, 0, 0);
      }
    }

    // Load signal cross section uncertainties
    stringstream xsFileName;
    xsFileName << "config"+string(do2012?"_2012/":"_2011/")+"xs_files/xs_" << mass << ".txt";
    if (useThSys) readInUncerts(xsFileName.str(), thR, mass);

    // Load signal branching ratio uncertainties
    stringstream brFileName;
    brFileName << "config"+string(do2012?"_2012/":"_2011/")+"br_files/br_" << mass << ".txt";
    if (useThSys) readInUncerts(brFileName.str(), thR, mass);

    // Load top0j and zjets uncerts
    if (ZMode == 2 && (doee || domm))
    {
      readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"zjets_pacman.txt", normR, mass);
    }
    if(!doZtautauCR){
      if (splitzjets)
	{
	  if(!doABCD2j)
	    { 
	      if(!NoMETCutDF){
		readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"ztautau_uncerts.txt", normR, mass);
	      }else{
		readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"ztautau_uncerts_NoMETCutDF.txt", normR, mass);
	      }
	      readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"zjets_uncerts_2j.txt", normR, mass);
	    }
	  else
	    {
	      readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"ztautau_uncerts_withabcd.txt", normR, mass); 
	      if (useDetSys) readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"zjets_abcd_uncerts_2j.txt", normR, mass);
	    }
	  //if(doABCD2j && do2012) readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"zjets_abcd_uncerts_2j.txt", normR, mass);
	  //else if(!doABCD2j && do2012) readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"zjets_uncerts_2j.txt", normR, mass);
	}
    }else{
      
      if (mode == 8) readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"alpha_uncerts_ztautau_3D.txt", thR, mass, 0, 0);
      else readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"alpha_uncerts_ztautau.txt", thR, mass, 0, 0);
    }
    
    if (useDetSys)
      {
	if (mode != 8) readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"top0j_uncerts.txt", normR, mass);
	else if (mode == 8){
	  if(!useAgnosticSignalRegionNaming)readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"top0j_uncerts_3D.txt", normR, mass);
	  else readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"top0j_uncerts_3D_AgnosticNaming.txt", normR, mass);
	}
	
	
	if (ZMode == 1)
	  {
	    if (useHighPt)
	      {
		stringstream fileName;
		fileName << "config"+string(do2012?"_2012/":"_2011/")+"zjets_abcdef";
		if (!useHighMass && !useHighMass2) fileName << "_lowm";
		else fileName << "_highm";
		fileName << ".txt";
		readInUncerts(fileName.str(), normR, mass);
	      }
	    
	    if (useLowPt)
	      {
		stringstream fileName;
		fileName << "config"+string(do2012?"_2012/":"_2011/")+"zjets_lowpt_abcdef.txt";
		
		readInUncerts(fileName.str(), normR, mass);
	      }
	  }
	else if (ZMode == 0)
	  {
	    readInUncerts("config"+string(do2012?"_2012/":"_2011/")+"zjets_scalefs.txt", normR, mass);
	  }
      }
    
    // Merge flat systemaitcs into one vector
    vector<vector<Response> > vecFlatSys;
    vecFlatSys.push_back(normR);
    // vecFlatSys.push_back(WjR); // now read from norms.txt
    vecFlatSys.push_back(thR);
    
    vector<Response> flatSys;
    joinSystematics(vecFlatSys, flatSys);

    int nrDetR = detR.size();
    int nrFlat = flatSys.size();

    // Prepare directory structures
    system(("mkdir -vp rev/"+version+"/xml/"+smass+salt+"/"+subfolder).c_str());
    // system(("cp config"+string(do2012?"_2012/":"_2011/")+"HistFactorySchema.dtd rev/"+version+"/xml/"+smass+salt+"/"+subfolder).c_str());
    
    RooStats::HistFactory::Measurement meas("lvlv", "lvlv");
    meas.SetOutputFilePrefix("./rev/"+version+"/xml/"+smass+salt+"/"+subfolder+"/output");
    meas.SetExportOnly(1);
    meas.SetPOI("SigXsecOverSM_HWW");
    meas.SetLumi(1.0);
    meas.SetLumiRelErr(0.037);



    //Put pacman equations into workspaces
    
    //grab the pdf, data, and sets
    
    for (int ireg = 0; ireg < nrRegions; ireg++)
      {
	Region* reg = &(*regions)[ireg];
	int nrChannels = reg->channels.size();
	
	for (int ichan = 0; ichan < nrChannels; ichan++)
	  {
	    Channel* c = &reg->channels[ichan];
	    string jetName = c->jetName;
	    string channelName = c->name;
	    
	    
	    if (ZMode == 2 && (jetName == "0j" || jetName == "1j") && (channelName == "ee" || channelName == "mm" || channelName == "SF")){
	      if ((*regions)[ireg].name == "ASR"||(*regions)[ireg].name.find("CZpeak") != string::npos ||(*regions)[ireg].name == "EWWCR")
		{
		  meas.AddPreprocessFunction(Form("PM_EFF_f_recoil_rejected_DY%s_func",jetName.c_str()),Form("1 - PM_EFF_f_recoil_DY%s",jetName.c_str()), Form("PM_EFF_f_recoil_DY%s[0,1]",jetName.c_str()));
		}	      
	      
	      //if(splitEfficiencies){
	      if ((*regions)[ireg].name.find("EfrecWWCR") != string::npos){
		  meas.AddPreprocessFunction(Form("PM_EFF_f_recoil_rejected_NDY_WW%s_func",jetName.c_str()),Form("1 - PM_EFF_f_recoil_NDY_WW%s",jetName.c_str()), Form("PM_EFF_f_recoil_NDY_WW%s[0,1]",jetName.c_str()));
		}
	      if ((*regions)[ireg].name.find("CZpeak") != string::npos)
		  {
		    meas.AddPreprocessFunction(Form("PM_EFF_f_recoil_rejected_NDY_ZP%s_func",jetName.c_str()),Form("1 - PM_EFF_f_recoil_NDY_ZP%s",jetName.c_str()), Form("PM_EFF_f_recoil_NDY_ZP%s[0,1]",jetName.c_str()));
		  }
	      for (int isam = 0; isam < nrSamples; isam++){
		
		if ((*regions)[ireg].name.find("AfrecSR") != string::npos){
		    if (((*samples)[isam].type == "signal" || ((*samples)[isam].name.find("ggf") != string::npos && doVBF2j)) && (channelName == "ee" || channelName == "mm" || channelName == "SF")){
		      meas.AddPreprocessFunction("PM_EFF_f_recoil_NDY_SR0j_func","(1-((1-@1)*((@0*(1-@0))/(@2*(1-@2)))))*pow((1+((@0*(1-@0))/(@2*(1-@2)))*0.1),@3)*@0", "PM_EFF_f_recoil_NDY_SR0j[1,0,1], f_0j[0,2], epsilon0_0j[0,1], alpha_PM_theta_temp_SR0j[0,-5,5]");
		      meas.AddPreprocessFunction("PM_EFF_f_recoil_NDY_SR1j_func","(1-((1-@1)*((@0*(1-@0))/(@2*(1-@2)))))*pow((1+((@0*(1-@0))/(@2*(1-@2)))*0.1),@3)*@0"," PM_EFF_f_recoil_NDY_SR1j[1,0,1], f_1j[0,2], epsilon0_1j[0,1], alpha_PM_theta_temp_SR1j[0,-5,5]");
		    }
		}
		if ((*regions)[ireg].name.find("ASR") != string::npos){
		  if (((*samples)[isam].type == "signal" || ((*samples)[isam].name.find("ggf") != string::npos && doVBF2j)) && (channelName == "ee" || channelName == "mm" || channelName == "SF")){
		    meas.AddPreprocessFunction("PM_EFF_f_recoil_rejected_NDY_SR0j_func","1-((1-((1-@1)*((@0*(1-@0))/(@2*(1-@2)))))*pow((1+((@0*(1-@0))/(@2*(1-@2)))*0.1),@3)*@0)", "PM_EFF_f_recoil_NDY_SR0j[1,0,1], f_0j[0,2], epsilon0_0j[0,1], alpha_PM_theta_temp_SR0j[0,-5,5]");
		    meas.AddPreprocessFunction("PM_EFF_f_recoil_rejected_NDY_SR1j_func","1-((1-((1-@1)*((@0*(1-@0))/(@2*(1-@2)))))*pow((1+((@0*(1-@0))/(@2*(1-@2)))*0.1),@3)*@0)", "PM_EFF_f_recoil_NDY_SR1j[1,0,1], f_1j[0,2], epsilon0_1j[0,1], alpha_PM_theta_temp_SR1j[0,-5,5]");
		  }else{
		    if(splitEfficiencies) meas.AddPreprocessFunction(Form("PM_EFF_f_recoil_rejected_NDY_SR%s",jetName.c_str()),Form("1 - PM_EFF_f_recoil_NDY_SR%s",jetName.c_str()), Form("PM_EFF_f_recoil_NDY_SR%s[1,0,1]",jetName.c_str()));
		    else meas.AddPreprocessFunction(Form("PM_EFF_f_recoil_rejected_NDY%s_func",jetName.c_str()),Form("1 - PM_EFF_f_recoil_NDY%s",jetName.c_str()), Form("PM_EFF_f_recoil_NDY%s[0,1]",jetName.c_str()));
		  }
		}
	      }
	      //}else{
		//}
	    }
	  }//end of chan loop
      }//end of reg loop
    
    // Main loop to make XML for each channel
    
    
    for (int ireg = 0; ireg < nrRegions; ireg++)
      {
	Region* reg = &(*regions)[ireg];
	int nrChannels = reg->channels.size();
	
	for (int ichan = 0; ichan < nrChannels; ichan++)
	  {
	    Channel* c = &reg->channels[ichan];
	    string jetName = c->jetName;
	    string channelName = c->name;

	    
	    bool skipLoop = skipRegion(reg, channelName, jetName);
	    if (skipLoop) continue;
	    
	    if (!doWWCR_MVA && reg->name.find("mainControl") != string::npos && ((jetName == "2j" && !doWWCR2j) || (doMVA == 1 && jetName == "0j"))) continue;
	    
	    string channelName_noyear = channelName+"_"+(*regions)[ireg].name+"_"+jetName;
	    string channelName_plusyear = channelName_noyear+(do2012?"_2012":"");
	    string folder = "rev/"+version+"/normHists/Nominal/Normal/";
	    
	    RooStats::HistFactory::Channel chan(channelName_plusyear);
	
	    string constraintType;
	    if (statMode == 0)
	      {
		constraintType = "Poisson";
	      }
	else if (statMode == 1)
          {
            constraintType = "Guassian";
          }
	else
          {
            LOG(logERROR) << "Undefined stat mode: " << statMode;
            exit(1);
          }
	    if (useStatSys) 
	      {
		chan.SetStatErrorConfig(stat_cutoff, constraintType);            
		//if (reg->isMCMSCR) chan.GetStatErrorConfig().SetRelErrorThreshold(10e9); // effectively remove
	      }
	    
	    if (doData) chan.SetData(channelName_noyear, folder+"data.root");
	    
	    for (int isam = 0; isam < nrSamples; isam++)
	      {
		
		if (jetName != "2j" && ((*samples)[isam] == "wwew" || (*samples)[isam] == "wzzzew" || (*samples)[isam] == "zjetsew" || (*samples)[isam] == "zleplepew" || (*samples)[isam] == "ztautauew")) continue;
		
		if ((*samples)[isam].type == "signal" || ((*samples)[isam].name.find("ggf") != string::npos && doVBF2j))
		  {
		    if ((channelName == "em" && (*regions)[ireg].name.find("AfrecSR")  != string::npos) ||
			(channelName == "em" && (*regions)[ireg].name.find("ASR")        != string::npos) ||
			(channelName == "em" && (*regions)[ireg].name.find("EfrecWWCR")  != string::npos) ||
			(channelName == "em" && (*regions)[ireg].name.find("EWWCR")      != string::npos) ||
			(channelName == "em" && (*regions)[ireg].name.find("CfrecZpeak") != string::npos) ||
			(channelName == "em" && (*regions)[ireg].name.find("CZpeak")     != string::npos) ||
			(channelName == "me" && (*regions)[ireg].name.find("AfrecSR")    != string::npos) ||
			(channelName == "me" && (*regions)[ireg].name.find("ASR")        != string::npos) ||
			(channelName == "me" && (*regions)[ireg].name.find("EfrecWWCR")  != string::npos) ||
			(channelName == "me" && (*regions)[ireg].name.find("EWWCR")      != string::npos) ||
			(channelName == "me" && (*regions)[ireg].name.find("CfrecZpeak") != string::npos) ||
			(channelName == "me" && (*regions)[ireg].name.find("CZpeak")     != string::npos) ||
			(channelName == "OF" && (*regions)[ireg].name.find("AfrecSR")    != string::npos) ||
			(channelName == "OF" && (*regions)[ireg].name.find("ASR")        != string::npos) ||
			(channelName == "OF" && (*regions)[ireg].name.find("EfrecWWCR")  != string::npos) ||
			(channelName == "OF" && (*regions)[ireg].name.find("EWWCR")      != string::npos) ||
			(channelName == "OF" && (*regions)[ireg].name.find("CfrecZpeak") != string::npos) ||
			(channelName == "OF" && (*regions)[ireg].name.find("CZpeak")     != string::npos))
		      {
			continue;
		      }
		  }
		
		string normByTheory = "True"; // apply lumi?
		if (((!useTheoryWW && ((*samples)[isam] == "ww" || (*samples)[isam] == "ggww" || (*samples)[isam] == "qqww") && (jetName == "1j" || jetName == "0j")) ||
		     (!useTheoryTop && ((*samples)[isam] == "ttbar" || (*samples)[isam] == "st") && (jetName == "1j" || jetName == "2j")) ||
		     (*samples)[isam] == "wjets" ||
		     (ZMode != 0 && ((*samples)[isam] == "zjets" || (*samples)[isam] == "zleplep" || (*samples)[isam] == "ztautau")) ||
		     (((*samples)[isam] == "ttbar" || (*samples)[isam] == "st") && (jetName == "0j")))) normByTheory = "False";
		if(doABCD2j && ((*samples)[isam] == "zleplep" || (*samples)[isam] == "ztautau" || (*samples)[isam] == "zjets") && (*regions)[ireg].name.find("signalLike") != string::npos && jetName == "2j") normByTheory = "False";
		if(doSameSignCR && ((*samples)[isam] == "wzzz" || (*samples)[isam] == "wg" || (*samples)[isam] == "wgs")) normByTheory = "False"; 
		if(doZtautauCR && ((*samples)[isam] == "zjets" || (*samples)[isam] == "zleplep" || (*samples)[isam] == "ztautau")) normByTheory = "False";
		if(doOSmSS && ((*samples)[isam] == "wzzz" || (*samples)[isam] == "wg" || (*samples)[isam] == "wgs")) normByTheory = "False"; // Don't apply Lumi systs if normalizing non WW dibosons from data
		
		// Use normByTheory when applying lumi sys later. here keep hf's lumi for possible scale tests
		
		RooStats::HistFactory::Sample sample((*samples)[isam].name, channelName_noyear, folder+(*samples)[isam].name+".root");
		
		if (!(channelName == "em" && (*regions)[ireg].name.find("AfrecSR")    != string::npos) &&
		    !(channelName == "em" && (*regions)[ireg].name.find("ASR")        != string::npos) &&
		    !(channelName == "em" && (*regions)[ireg].name.find("EfrecWWCR")  != string::npos) &&
		    !(channelName == "em" && (*regions)[ireg].name.find("EWWCR")      != string::npos) &&
		    !(channelName == "em" && (*regions)[ireg].name.find("CfrecZpeak") != string::npos) &&
		    !(channelName == "em" && (*regions)[ireg].name.find("CZpeak")     != string::npos) &&
		    !(channelName == "me" && (*regions)[ireg].name.find("AfrecSR")    != string::npos) &&
		    !(channelName == "me" && (*regions)[ireg].name.find("ASR")        != string::npos) &&
		    !(channelName == "me" && (*regions)[ireg].name.find("EfrecWWCR")  != string::npos) &&
		    !(channelName == "me" && (*regions)[ireg].name.find("EWWCR")      != string::npos) &&
		    !(channelName == "me" && (*regions)[ireg].name.find("CfrecZpeak") != string::npos) &&
		    !(channelName == "me" && (*regions)[ireg].name.find("CZpeak")     != string::npos) &&
		    !(channelName == "OF" && (*regions)[ireg].name.find("AfrecSR")    != string::npos) &&
		    !(channelName == "OF" && (*regions)[ireg].name.find("ASR")        != string::npos) &&
		    !(channelName == "OF" && (*regions)[ireg].name.find("EfrecWWCR")  != string::npos) &&
		    !(channelName == "OF" && (*regions)[ireg].name.find("EWWCR")      != string::npos) &&
		    !(channelName == "OF" && (*regions)[ireg].name.find("CfrecZpeak") != string::npos) &&
		    !(channelName == "OF" && (*regions)[ireg].name.find("CZpeak")     != string::npos))
		  {
		    for (int ir = 0; ir < nrDetR; ir++) // write shape systematics
		      {
			Response* r = &detR[ir];
			if (r->match(channelName, jetName, (*samples)[isam].name, (*regions)[ireg].name))
			  {
			    sample.AddHistoSys(r->name, r->loHistName, r->loFileName, "", r->hiHistName, r->hiFileName, "");
			  }
		      }
		    
              for (int ir = 0; ir < nrFlat; ir++) // write flat systematics
		{
		  Response* r = &flatSys[ir];
		  // 		if (r->name.find("TOP_THEO_1j") == string::npos &&
		  // 		    r->name.find("BTag_BEFF") == string::npos) continue;
		  if (r->match(channelName, jetName, (*samples)[isam].name, (*regions)[ireg].name))
		    {
		      sample.AddOverallSys(r->name, 1 + r->lo, 1 + r->hi);
		    }
		}

              if (useDetSys && normByTheory == "True")
		{
		  double lumi_val = 0.018;
		  if (do2012) lumi_val = 0.036; // <- new
		  sample.AddOverallSys("LUMI"+string(do2012?"_2012":"_2011"), 1./(1 + lumi_val), 1 + lumi_val);
		}
	      
	      
	      //need to add the 1j->2j btag extrapolation uncertainty somewhere in the xml so that the variables are in the workspace
	      double uncertainty = 0.06;// put in config file?
	      if (doModifiedCMS && (*regions)[ireg].name.find("topbox") != string::npos && jetName == "1j") // may need to edit later if cancel2j is used
		{
		  sample.AddOverallSys("ATLAS_btag21j_extrap",1./(1+uncertainty),1+uncertainty);// OverallSys is a contrained norm variation, include low and high uncertainty 
		}
	      
	      
              if ((*samples)[isam].type == "signal")
		{
		  sample.AddNormFactor("SigXsecOverSM_HWW", 1.0, 0.0, 50, 1);
		}
              else if (include125BG && (*samples)[isam].name.find("125bg") != string::npos)
              {
                sample.AddNormFactor("SigXsecOverSM_125_HWW", 1.0, 0.0, 50);
              }
	      
	      
	      if (!useTheoryWW && ((*samples)[isam] == "ww" || (*samples)[isam] == "ggww" || (*samples)[isam] == "qqww") && (jetName != "2j" || (doWWCR2j && jetName == "2j")))
		{
		  sample.AddNormFactor("ATLAS_norm_WW"+(splitNFs?channelName:"")+jetName, 1.0, 0.0, 10.0);
		}
	      
	      if(doSameSignCR && ((*samples)[isam] == "wwew" || (*samples)[isam] == "wzzz" || (*samples)[isam] == "wzzzew" || (*samples)[isam] == "wgs" || (*samples)[isam] == "wg") && (jetName != "2j"))
		{
		  sample.AddNormFactor("ATLAS_norm_Diboson"+(splitNFs?channelName:"")+jetName, 1.0, 0.0, 10.0); 
		}
	      if(doZtautauCR && ((*samples)[isam] == "zjets" || (*samples)[isam] == "ztautau" || (*samples)[isam] == "zleplep") && (jetName != "2j"))
		{
		  sample.AddNormFactor("ATLAS_norm_Ztautau"+(splitNFs?channelName:"")+jetName, 1.0, 0.0, 10.0);
		}
	      
	      
	      
	      //we are in 1j top cr
	      if (reg->name.find("tbPass") == string::npos && reg->name.find("tbFail") == string::npos)
		{
		  if (!useTheoryTop && ((*samples)[isam] == "ttbar" || (*samples)[isam] == "st") && jetName != "0j")
		    {
		      if (((reg->name.find("sr") != string::npos && reg->name.find("Tag") != string::npos) || 
			   reg->name.find("signalLike") != string::npos || reg->name.find("topboxLow") != string::npos) && (splitTopCR && jetName != "2j"))
			{
			  sample.AddNormFactor("ATLAS_norm_TopLow"+(splitNFs?channelName:"")+jetName, 1.0, 0.0, 10.0);
			}
		      else
			{
			  sample.AddNormFactor("ATLAS_norm_Top"+(splitNFs?channelName:"")+jetName, 1.0, 0.0, 10.0);
			}
		    }
		  
		  if (profileggf && doVBF2j && (*samples)[isam] == ("ggf"+smass).c_str())
		    {
		      sample.AddNormFactor("ATLAS_norm_ggf", 1.0, 0.0, 10.0);
		    }
		  
		  if (doCRVRcheck && (*regions)[ireg].name.find("signalLike") != string::npos && ((*samples)[isam] == "ww" || (*samples)[isam] == "ggww" || (*samples)[isam] == "qqww"))
		  {
		    sample.AddNormFactor("CRVR_HWW", 1.0, 0.0, 50, 1);
		  }
		
		if (doModifiedCMS && (reg->isMCMSCR || (reg->name.find("topbox") != string::npos && (jetName != "2j" || cancel2j))) && (jetName == "1j" || jetName == "2j" && cancel2j) && ((*samples)[isam] == "ttbar" || (*samples)[isam] == "st"))
		  {
		    if (reg->name.find("x2Tag") != string::npos) sample.AddNormFactor("ATLAS_norm_btag_squared", 1.0, 0.0, 1.5);
		    else sample.AddNormFactor("ATLAS_norm_btag", 1.0, 0.0, 1.5);
		    if (reg->name.find("tbTag") != string::npos || (useFullTagSample && reg->name.find("tbFail") != string::npos)) sample.AddNormFactor("twice_ATLAS_norm_btag", 1.0, 0.0, 1.5);
		  }
		}
	      else//we are in 2j top cr 
	      {
		if (((*samples)[isam] == "ttbar" || (*samples)[isam] == "st"))
		  {
		    sample.AddNormFactor("ATLAS_norm_TopPF"+jetName, 1.0, 0.0, 10.0);
		    if (reg->name.find("tbPass") != string::npos) sample.AddNormFactor("ATLAS_norm_btag", 1.0, 0.0, 1.5);
		  }
	      }
	      
              // Spin0/Spin2 fraction epsilon
              if (doSpin)
		{
		  if (((*samples)[isam].name.find("spin0p") != string::npos && useJHUspin0) || ((*samples)[isam].name.find("ggf") != string::npos && !useJHUspin0))
		    {
                  sample.AddNormFactor("ATLAS_epsilon", 0.5, -5.0, 5.0, 1);
		    }
		  else if ((*samples)[isam].name.find("spin2p") != string::npos)
		    {
		      sample.AddNormFactor("ATLAS_epsilon_rejected", 0.5, -5.0, 5.0, 1);
		    }
		}
	      
              if (doSpin && doMVA && ((*samples)[isam] == "zjets" || (*samples)[isam] == "zleplep" || (*samples)[isam] == "ztautau"))
		{
		  sample.AddNormFactor("ATLAS_norm_Z"+jetName, 1.0, 0.0, 10.0);
		}
            }
	    
            if (((*regions)[ireg].name.find("AfrecSR") != string::npos || (*regions)[ireg].name.find("ASR") != string::npos) && (channelName == "ee" || channelName == "mm" || channelName == "SF"))
	      {
		if ((*samples)[isam].type == "signal" || ((*samples)[isam].name.find("ggf") != string::npos && doVBF2j))
		  {
		    sample.AddOverallSys("PM_theta_SR"+jetName, 1.0, 1.0);
		  }
	      }
	    
            // Pacman control regions and nuisance parameters
            if (ZMode == 2 && (jetName == "0j" || jetName == "1j"))
	      {
		if ((*samples)[isam] == "zjets" || (*samples)[isam] == "zleplep" || (*samples)[isam] == "ztautau")
		  {
		    if (channelName == "ee" || channelName == "mm" || channelName == "SF")
		      {
			if ((*regions)[ireg].name == "AfrecSR")
			  {
			    sample.AddNormFactor("ATLAS_norm_SF_MUSR_DY"+string(splitPacmanNFs?"_"+channelName:"")+jetName, 1.0, 0.0, 50.0);
			    sample.AddNormFactor("PM_EFF_f_recoil_DY"+jetName, jetName=="0j"?f_DY_all_0j:f_DY_all_1j, 0.0, 1.0);
                  }
			
                  if ((*regions)[ireg].name == "ASR")
		    {
		      sample.AddNormFactor("ATLAS_norm_SF_MUSR_DY"+string(splitPacmanNFs?"_"+channelName:"")+jetName, 1.0, 0.0, 50.0);
		      sample.AddNormFactor("PM_EFF_f_recoil_rejected_DY"+jetName+"_func", jetName=="0j"?1-f_DY_all_0j:1-f_DY_all_1j, 0.0, 1.0);
		    }
		  
                  if ((*regions)[ireg].name.find("CfrecZpeak") != string::npos)
		    {
		      sample.AddNormFactor("ATLAS_norm_SF_MU_DY"+jetName, 1.0, 0.0, 50.0);
		      sample.AddNormFactor("PM_EFF_f_recoil_DY"+jetName, jetName=="0j"?f_DY_all_0j:f_DY_all_1j, 0.0, 1.0);
		    }
		  
                  if ((*regions)[ireg].name.find("CZpeak") != string::npos)
		    {
		      sample.AddNormFactor("ATLAS_norm_SF_MU_DY"+jetName, 1.0, 0.0, 50.0);
		      sample.AddNormFactor("PM_EFF_f_recoil_rejected_DY"+jetName+"_func", jetName=="0j"?1-f_DY_all_0j:1-f_DY_all_1j, 0.0, 1.0);
		    }
		  
                  if ((*regions)[ireg].name == "EfrecWWCR")
		    {
		      sample.AddNormFactor("ATLAS_norm_SF_MU_DY"+jetName, 1.0, 0.0, 50.0);
		      sample.AddNormFactor("PM_EFF_f_recoil_DY"+jetName, jetName=="0j"?f_DY_all_0j:f_DY_all_1j, 0.0, 1.0);
		    }
		  
                  if ((*regions)[ireg].name == "EWWCR")
		    {
		      sample.AddNormFactor("ATLAS_norm_SF_MU_DY"+jetName, 1.0, 0.0, 50.0);
		      sample.AddNormFactor("PM_EFF_f_recoil_rejected_DY"+jetName+"_func", jetName=="0j"?1-f_DY_all_0j:1-f_DY_all_1j, 0.0, 1.0);
		    }
		      }
		  }
		else
		  {
		    if ((*regions)[ireg].name.find("EfrecWWCR") != string::npos)
		      {
			if (splitEfficiencies) sample.AddNormFactor("PM_EFF_f_recoil_NDY_WW"+jetName, jetName=="0j"?f_NDY_WWCR_0j:f_NDY_WWCR_1j, 0.0, 1.0);
			else                   sample.AddNormFactor("PM_EFF_f_recoil_NDY"+jetName, f_NDY_all, 0.0, 1.0);
                }
		    
		    if ((*regions)[ireg].name.find("EWWCR") != string::npos)
		      {
			if (splitEfficiencies) sample.AddNormFactor("PM_EFF_f_recoil_rejected_NDY_WW"+jetName+"_func", jetName=="0j"?1-f_NDY_WWCR_0j:1-f_NDY_WWCR_1j, 0.0, 1.0);
                  else                   sample.AddNormFactor("PM_EFF_f_recoil_rejected_NDY"+jetName+"_func", 1-f_NDY_all, 0.0, 1.0);
		      }
		    
		    if ((*regions)[ireg].name.find("CfrecZpeak") != string::npos)
		      {
			if (splitEfficiencies) sample.AddNormFactor("PM_EFF_f_recoil_NDY_ZP"+jetName, jetName=="0j"?f_NDY_ZP_0j:f_NDY_ZP_1j, 0.0, 1.0);
			else                   sample.AddNormFactor("PM_EFF_f_recoil_NDY"+jetName, f_NDY_all, 0.0, 1.0);
		      }
		    
		    if ((*regions)[ireg].name.find("CZpeak") != string::npos)
		      {
			if (splitEfficiencies) sample.AddNormFactor("PM_EFF_f_recoil_rejected_NDY_ZP"+jetName+"_func", jetName=="0j"?1-f_NDY_ZP_0j:1-f_NDY_ZP_1j, 0.0, 1.0);
                  else                   sample.AddNormFactor("PM_EFF_f_recoil_rejected_NDY"+jetName+"_func", 1-f_NDY_all, 0.0, 1.0);
		      }
		    
		    if ((*regions)[ireg].name.find("AfrecSR") != string::npos)
		      {
			if (((*samples)[isam].type == "signal" || ((*samples)[isam].name.find("ggf") != string::npos && doVBF2j)) && (channelName == "ee" || channelName == "mm" || channelName == "SF"))
                  {
                    sample.AddNormFactor("PM_EFF_f_recoil_NDY_SR"+jetName+"_func", 1.0, 0.0, 1.0);
                  }
			else
			  {
			    if (splitEfficiencies) sample.AddNormFactor("PM_EFF_f_recoil_NDY_SR"+jetName, jetName=="0j"?f_NDY_SR_0j:f_NDY_SR_1j, 0.0, 1.0);
			    else                   sample.AddNormFactor("PM_EFF_f_recoil_NDY"+jetName, f_NDY_all, 0.0, 1.0);
			  }
		      }
		    
                if ((*regions)[ireg].name.find("ASR") != string::npos)
		  {
		    if (((*samples)[isam].type == "signal" || ((*samples)[isam].name.find("ggf") != string::npos && doVBF2j)) && (channelName == "ee" || channelName == "mm" || channelName == "SF"))
		      {
			sample.AddNormFactor("PM_EFF_f_recoil_rejected_NDY_SR"+jetName+"_func", 1.0, 0.0, 1.0);
		      }
		    else
		      {
			if (splitEfficiencies) sample.AddNormFactor("PM_EFF_f_recoil_rejected_NDY_SR"+jetName, jetName=="0j"?1-f_NDY_SR_0j:1-f_NDY_SR_1j, 0.0, 1.0);
			else                   sample.AddNormFactor("PM_EFF_f_recoil_rejected_NDY"+jetName+"_func", 1-f_NDY_all, 0.0, 1.0);
		      }
		  }
		  }
	      }
	    
            if ((channelName == "ee" || channelName == "mm" || channelName == "SF") && (ZMode == 2 && (jetName == "0j" || jetName == "1j")))
	      {
		if ((*samples)[isam] == "zjets" || (*samples)[isam] == "zleplep" || (*samples)[isam] == "ztautau") // DY
		  {
		    if ((*regions)[ireg].name == "topbox")
		      {
			sample.AddNormFactor("ATLAS_norm_SF_MU_DY"+jetName, 1.0, 0.0, 50.0);
		      }
		  }
	      }
	    
            // Normalization parameter for individual samples (various uses)
            sample.AddNormFactor("ATLAS_sampleNorm_"+(*samples)[isam].name, 1.0, 0.0, 10.0, 1);
	    
            // Normalization parameter for BR(HWW)
            if ((*samples)[isam].type == "signal" || ((*samples)[isam].name.find("ggf") != string::npos && doVBF2j))
            {
              sample.AddNormFactor("mu_BR_WW", 1.0, 0.0, 10.0, 1);
            }
	    
            bool activate = 1;
            if ((*samples)[isam].type == "signal" || ((*samples)[isam].name.find("ggf") != string::npos && doVBF2j)) activate = 0;
            if (((*samples)[isam].name == "zjets" || (*samples)[isam].name == "zleplep" || (*samples)[isam].name == "ztautau") && (ZMode == 1 || (ZMode == 2 && doABCD2j && jetName == "2j")) && channelName != "em") activate = 0;
            if (useStatSys && activate) sample.ActivateStatError();

            chan.AddSample(sample);

          }

          meas.AddChannel(chan);

        
      }
    }

    meas.CollectHistograms();
    meas.PrintTree();
    meas.PrintXML("./rev/"+version+"/xml/"+smass+salt+"/"+subfolder, meas.GetOutputFilePrefix());
    MakeModelAndMeasurementFast(meas);

  }

  LOG(logINFO) << "Done making workspace: " << subfolder;
}

// join vectors, removing duplicates
void joinSystematics(vector<vector<Response> >& vecResponse, vector<Response>& joined)
{
  int nrVec = vecResponse.size();
  for (int i=0;i<nrVec;i++)
  {
    vector<Response>* vecR = &vecResponse[i];
    int nrVec2 = vecR->size();
    for (int j=0;j<nrVec2;j++)
    {
      Response* thisR = &vecResponse[i][j];
      Response* r = NULL;
      int nrJoined = joined.size();
      for (int k=0;k<nrJoined;k++)
      {
        Response* res = &joined[k];
        if (res->match(thisR->channel, thisR->jet, thisR->sample, thisR->region) && res->name == thisR->name)
        {
          r = res;
          break;
        }
      }

      if (!r)
      {
        Response newR = *thisR;
        newR.hi += 1;
        newR.lo += 1;
        joined.push_back(newR);
      }
      else
      {
        r->hi *= 1+thisR->hi;
        r->lo *= 1+thisR->lo;
      }
    }
  }

  int nrJoined = joined.size();
  for (int i=0;i<nrJoined;i++)
  {
    joined[i].hi -= 1;
    joined[i].lo -= 1;
  }
}
