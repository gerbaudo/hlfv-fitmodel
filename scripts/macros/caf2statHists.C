// Macro to convert analysis root files to stat format hists

/*
Note: histograms containing "0j" will receive 0-jet ztt and top NFs.  Others will get the incl ztt NF.

*/

#include <utility>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "TSystem.h"
#include "TString.h"
#include "TObjArray.h"
#include "TFile.h"

//#include "macros/setup.C"


bool do2011 = false;
bool doCPS = false;


void setupSampleHist(std::vector<std::string> &samples, std::vector<std::string> &paths, std::vector<std::string> &histograms, std::string rebinningFile="")
{

  // Set samples and paths (X is a placeholder for flavors)

  samples.push_back(std::string("data")); paths.push_back(std::string("data/X"));

  samples.push_back(std::string("ggww")); paths.push_back(std::string("bkg/X/diboson/WW/ggWW"));
  samples.push_back(std::string("qqww")); paths.push_back(std::string("bkg/X/diboson/WW/qqWW"));
  samples.push_back(std::string("ttbar")); paths.push_back(std::string("bkg/X/top/ttbar"));
  samples.push_back(std::string("st")); paths.push_back(std::string("bkg/X/top/singletop"));
  samples.push_back(std::string("wg")); paths.push_back(std::string("bkg/X/diboson/NonWW/Wgamma/gammaWincl/Wgamma"));
  samples.push_back(std::string("wgs")); paths.push_back(std::string("bkg/X/diboson/NonWW/Wgamma/gammaWincl/WgammaStar"));
  samples.push_back(std::string("wzzz")); paths.push_back(std::string("bkg/X/diboson/NonWW/ZW + bkg/X/diboson/NonWW/ZZ"));
  samples.push_back(std::string("ztautau")); paths.push_back(std::string("bkg/X/Zjets/Nom/?/tt"));
  samples.push_back(std::string("zleplep")); paths.push_back(std::string("bkg/X/Zjets/Nom/?/ee + bkg/X/Zjets/Nom/?/mm"));
  
  samples.push_back(std::string("wwew")); paths.push_back(std::string("bkg/X/diboson/WW/161985"));
//  samples.push_back(std::string("wzzzew")); paths.push_back(std::string("bkg/X/diboson/NonWW/"));
  samples.push_back(std::string("zleplepew")); paths.push_back(std::string("bkg/X/ZjetsEW/ee + bkg/X/ZjetsEW/mm"));
  samples.push_back(std::string("ztautauew")); paths.push_back(std::string("bkg/X/ZjetsEW/tt"));

  samples.push_back(std::string("wjets")); paths.push_back(std::string("bkg/X/Wjets"));

  samples.push_back(std::string("ggfMASS")); paths.push_back(std::string("sig/X/mhMASS/ggf"));
  samples.push_back(std::string("vbfMASS")); paths.push_back(std::string("sig/X/mhMASS/vbf"));
  samples.push_back(std::string("zhMASS")); paths.push_back(std::string("sig/X/mhMASS/ZH"));
  samples.push_back(std::string("whMASS")); paths.push_back(std::string("sig/X/mhMASS/WH"));

  // Histograms
  if (rebinningFile!="") { //if a rebinning file is provided then we can read the list of histograms to dump from there
    std::cout << "Reading in rebinning configuration file: " << rebinningFile.c_str() << std::endl;
    ifstream inFile(rebinningFile.c_str());
    if (inFile.fail())  {
      std::cout << "ERROR::Couldn't open file: " << rebinningFile.c_str() << std::endl;
      exit(1);
    }
    string channels;
    string region;
    unsigned int nbins;
    while (!inFile.eof())  {
      inFile >> channels;
      inFile >> region;
      inFile >> nbins;
      std::string histname = channels + "_" + region;
      histograms.push_back(histname);
      std::cout << "Will take care of dumping " << histname.c_str() << std::endl;
    }
    inFile.close();
  } 
  else { //if no rebinning file is provided then we'll try to read a hard-coded list of hists
    histograms.push_back("em_signalLike1_0j");
    histograms.push_back("me_signalLike1_0j");
    histograms.push_back("em_signalLike1_1j");
    histograms.push_back("me_signalLike1_1j");
    histograms.push_back("em_signalLike2_0j");
    histograms.push_back("me_signalLike2_0j");
    histograms.push_back("em_signalLike2_1j");
    histograms.push_back("me_signalLike2_1j");
    histograms.push_back("OF_ASR_0j");
    histograms.push_back("OF_ASR_1j");
    histograms.push_back("SF_ASR_0j");
    histograms.push_back("SF_ASR_1j");
    histograms.push_back("OF_AfrecSR_0j");
    histograms.push_back("OF_AfrecSR_1j");
    histograms.push_back("SF_AfrecSR_0j");
    histograms.push_back("SF_AfrecSR_1j");
    histograms.push_back("OF_mainControl_0j");
    histograms.push_back("OF_mainControl_1j");
    histograms.push_back("OF_topbox_1j");
    histograms.push_back("OF_CZpeak_0j");
    histograms.push_back("OF_CZpeak_1j");
    histograms.push_back("SF_CZpeak_0j");
    histograms.push_back("SF_CZpeak_1j");
    histograms.push_back("OF_CfrecZpeak_0j");
    histograms.push_back("OF_CfrecZpeak_1j");
    histograms.push_back("SF_CfrecZpeak_0j");
    histograms.push_back("SF_CfrecZpeak_1j");

    // VBF
    histograms.push_back("em_signalLike1_2j");
    histograms.push_back("me_signalLike1_2j");


    /*  histograms.push_back("em_signalLike0301015_0j");
	histograms.push_back("me_signalLike0301015_0j");
	histograms.push_back("em_signalLike0301015_1j");
	histograms.push_back("me_signalLike0301015_1j");
	histograms.push_back("em_signalLike0301520_0j");
	histograms.push_back("me_signalLike0301520_0j");
	histograms.push_back("em_signalLike0301520_1j");
	histograms.push_back("me_signalLike0301520_1j");
	histograms.push_back("em_signalLike030201000_0j");
	histograms.push_back("me_signalLike030201000_0j");
	histograms.push_back("em_signalLike030201000_1j");
	histograms.push_back("me_signalLike030201000_1j");
	histograms.push_back("em_signalLike30501015_0j");
	histograms.push_back("me_signalLike30501015_0j");
	histograms.push_back("em_signalLike30501015_1j");
	histograms.push_back("me_signalLike30501015_1j");
	histograms.push_back("em_signalLike30501520_0j");
	histograms.push_back("me_signalLike30501520_0j");
	histograms.push_back("em_signalLike30501520_1j");
	histograms.push_back("me_signalLike30501520_1j");
	histograms.push_back("em_signalLike3050201000_0j");
	histograms.push_back("me_signalLike3050201000_0j");
	histograms.push_back("em_signalLike3050201000_1j");
	histograms.push_back("me_signalLike3050201000_1j");*/
 } //end of if no rebinning file is provided


}

int caf2statHists(std::string anaOutputDir = "", // Where are the output .root files
		  std::string version = "test", 
		  std::string syst="Normal",
		  bool useRebinningFile=false,
		  std::string prefix = "", // prefixes need a trailing underscore, and no other
		  std::string mass = "125",
		  std::string pathToLibrary = "HWWAnalysisCode/lib/libQFramework.so", 
		  std::string statDir = "./")
{

  if (anaOutputDir == "")
  {
    std::cout<<"Please specifiy anaOutputDir!"<<std::endl;
    return 1;
  }
    
  std::cout<<"Loading: "<<pathToLibrary.c_str()<<std::endl;
  gSystem->Load(pathToLibrary.c_str()); // load tqlibrary

  std::vector<std::string> samples;
  std::vector<std::string> paths;

  std::vector<std::string> histograms;
  std::string rebinningFileName="";
  if (useRebinningFile) rebinningFileName = "rev/" + version + "/rebinningList.txt";
  setupSampleHist(samples, paths, histograms, rebinningFileName);

  if (samples.size() != paths.size() )
  {
    std::cout<< "samples and paths are not the same size!" << std::endl;
    return 1;
  }
  

  // Load analysis file
  std::string anaName = syst;
  if (syst == "Normal") anaName = "Nominal";
  std::string inputfile = anaOutputDir + "samples_hww_analysis_" + version + "_" + anaName + ".root:samples_" + anaName;
  std::cout << "Will get histograms from " << inputfile.c_str() << std::endl;
  TQSampleFolder *samplesFolder = TQSampleFolder::loadSampleFolder(inputfile.c_str() );


  if (!samplesFolder)      {
    std::cout<<"Couldn't load: " << inputfile + "  Continuing" << std::endl;
    return 1;
  }

  //Deal with special cases where CAF name is different from stat variation name
  if (syst=="_BJetWeightUp")
    syst = "BtagUp";
  else if (syst=="_BJetWeightDown")
    syst = "BtagDown";
  else if (syst=="_CTJetWeightUp")
    syst = "CtagUp";
  else if (syst=="_CTJetWeightDown")
    syst = "CtagDown";
  else if (syst=="_MisTagJetWeightUp")
    syst = "MtagUp";
  else if (syst=="_MisTagJetWeightDown")
    syst = "MtagDown";
  else if (syst=="_lepIsoUp")
    syst = "IsoUp";
  else if (syst=="_lepIsoDown")
    syst = "IsoDown";
  else if (syst=="_lepTriggerSFup")
    syst = "TriggerUp";
  else if (syst=="_lepTriggerSFdown")
    syst = "TriggerDown";
  else if (syst=="_MTShapeUp")
    syst = "PowhegShapeUp";
  else if (syst=="_MTShapeDown")
    syst = "PowhegShapeDown";


  //Renaming the syst in case it is weight-related
  if (syst.at(0)=='_')
    syst = syst.substr(1);

  std::string statsystname = syst;

  for (unsigned int iSample = 0; iSample < samples.size(); iSample++)      
  {
    // Get parent syst folder from systsList.txt
    std::string parentfolder = "";
    std::fstream systsList;
    std::string line, one, two;
    std::stringstream ss;
    if (syst=="Nominal") {
      parentfolder="Nominal";
      statsystname="Normal";
    }
    else if (syst=="ElecFakeWeightUp") {
      parentfolder="FakeRate_EL_HWW";
      statsystname="SysFakeUp";
    }
    else if (syst=="ElecFakeWeightDown") {
      parentfolder="FakeRate_EL_HWW";
      statsystname="SysFakeDown";
    }
    else if (syst=="MuonFakeWeightUp") {
      parentfolder="FakeRate_MU_HWW";
      statsystname="SysFakeUp";
    }
    else if (syst=="MuonFakeWeightDown") {
      parentfolder="FakeRate_MU_HWW";
      statsystname="SysFakeDown";
    }
    else {
      systsList.open( ("rev/" + version + "/systsList.txt").c_str(), std::ios::in );
      if (systsList.is_open())
	{
	  while(getline(systsList, line))
	    {
	      ss.str("");
	      ss.clear();
	      ss << line;
	      ss >> one >> two;

	      if (two == statsystname)
		{
		  if (statsystname.find("FlavComp") != std::string::npos)
		    {
		      if (samples[iSample] == "qqww" || samples[iSample] == "ggww")
			parentfolder = "ATLAS_JES_FlavComp_HWW_WW";
		      else if (samples[iSample] == "ttbar" || samples[iSample] == "st")
			parentfolder = "ATLAS_JES_FlavComp_HWW_tt";
		      else
			parentfolder = "ATLAS_JES_FlavComp_HWW_other";
		      break;
		    }
		  else{
		    parentfolder = one;
		    break;
		  }
		}
	    }
	  if (parentfolder == "")
	    {
	      std::cout <<"Failed to find parent folder"<<std::endl;
	      return 1;
	    }
	}
      else
	{
	  std::cout<< "Failed to open " << "rev/" << version << "/systsList.txt" <<std::endl;
	  return 1;
	}
      systsList.close();
    }//end of looking up for parentfolder

    std::cout << "Will write histograms into " << parentfolder.c_str() << "/" << statsystname.c_str() << " for sample " << samples[iSample].c_str() << std::endl;

    if (samples[iSample] == "data" && parentfolder != "Nominal") continue; //only nominal for data

    system( ("mkdir -vp " + statDir + "/rev/" + version + "/hists/" + parentfolder + "/" +statsystname ).c_str() );

    string sampleName = TString(samples[iSample]).ReplaceAll("MASS", mass);
    TFile * tempOutput = new TFile(std::string(statDir + "/rev/" + version + "/hists/" + parentfolder + "/" +statsystname + "/" + sampleName + ".root").c_str(), "RECREATE");

    // Get Ztt NFs
    ifstream ztautau_norms_file(std::string("bsub/" + version + "/ztautau_norms.txt").c_str());
    if (ztautau_norms_file.bad() || ztautau_norms_file.fail())
    {
      std::cout<< << "ERROR::Couldn't open ztautau norms file: " << std::string("bsub/" + version + "/ztautau_norms.txt").c_str() << std::endl;
      exit(1);
    }
    double sf_Ztautau_0jet, sf_Ztautau_incl;
    while (!ztautau_norms_file.eof())
    {
      string term;
      ztautau_norms_file >> term;
      if (ztautau_norms_file.eof()) break;
      double val;
      ztautau_norms_file >> val;
      
      if (term == "em_signalLike_0j") sf_Ztautau_0jet = val;
      else if (term == "em_signalLike_1j") sf_Ztautau_incl = val;

    }
    ztautau_norms_file.close();
    std::cout << "Will apply following NFs:" << std::endl;
    std::cout << "  sf_Ztautau_0jet = " << sf_Ztautau_0jet << std::endl;
    std::cout << "  sf_Ztautau_incl = " << sf_Ztautau_incl << std::endl;


/*    samplesFolder->setScaleFactor("signalLike1_0j", sf_Ztautau_0jet, "bkg/?/Zjets/Z/tt"); // Ztt 0j
    samplesFolder->setScaleFactor("signalLike2_0j", sf_Ztautau_0jet, "bkg/?/Zjets/Z/tt"); // Ztt 0j
    samplesFolder->setScaleFactor("mainControl_0j", sf_Ztautau_0jet, "bkg/?/Zjets/Z/tt"); // Ztt 0j
    samplesFolder->setScaleFactor("signalLike1_0j", sf_Ztautau_0jet, "bkg/?/Zjets/DY/tt"); // Ztt 0j
    samplesFolder->setScaleFactor("signalLike2_0j", sf_Ztautau_0jet, "bkg/?/Zjets/DY/tt"); // Ztt 0j
    samplesFolder->setScaleFactor("mainControl_0j", sf_Ztautau_0jet, "bkg/?/Zjets/DY/tt"); // Ztt 0j



    samplesFolder->setScaleFactor("signalLike1_1j", sf_Ztautau_incl, "bkg/?/Zjets/Z/tt"); // Ztt incl
    samplesFolder->setScaleFactor("signalLike2_1j", sf_Ztautau_incl, "bkg/?/Zjets/Z/tt"); // Ztt incl
    samplesFolder->setScaleFactor("mainControl_1j", sf_Ztautau_incl, "bkg/?/Zjets/Z/tt"); // Ztt incl
    samplesFolder->setScaleFactor("topbox_1j",      sf_Ztautau_incl, "bkg/?/Zjets/Z/tt"); // Ztt incl
    samplesFolder->setScaleFactor("signalLike1_1j", sf_Ztautau_incl, "bkg/?/Zjets/DY/tt"); // Ztt incl
    samplesFolder->setScaleFactor("signalLike2_1j", sf_Ztautau_incl, "bkg/?/Zjets/DY/tt"); // Ztt incl
    samplesFolder->setScaleFactor("mainControl_1j", sf_Ztautau_incl, "bkg/?/Zjets/DY/tt"); // Ztt incl
    samplesFolder->setScaleFactor("topbox_1j",      sf_Ztautau_incl, "bkg/?/Zjets/DY/tt"); // Ztt incl */

    // Get Top 0j NF 
    ifstream top_norms_file(std::string("bsub/" + version + "/top_norms.txt").c_str());
    if (top_norms_file.bad() || top_norms_file.fail())    {
      std::cout<< << "ERROR::Couldn't open top norms file: " << std::string("bsub/" + version + "/top_norms.txt").c_str() << std::endl;
      exit(1);
    }
    double NFTop0jHiPT = 1.07;
    while (!top_norms_file.eof())    {
      string term;
      top_norms_file >> term;
      if (top_norms_file.eof()) break;
      double val;
      top_norms_file >> val;
      
      if (term == "allregions0j") NFTop0jHiPT = val;
    }
    top_norms_file.close();

    std::cout << "  NFTop0jHiPT = " << NFTop0jHiPT << std::endl;
/*    samplesFolder->setScaleFactor("signalLike1_0j", NFTop0jHiPT, "bkg/?/top"); // Top 0J
    samplesFolder->setScaleFactor("signalLike2_0j", NFTop0jHiPT, "bkg/?/top"); // Top 0J
    samplesFolder->setScaleFactor("mainControl_0j", NFTop0jHiPT, "bkg/?/top"); // Top 0J*/


    // Grab the histograms
    for (unsigned int iHist=0; iHist < histograms.size(); iHist++)	
    {

      TObjArray * histTokens = TString(prefix + histograms[iHist]).Tokenize("_");
      TString name="";

      int flavors = 0;
      if (prefix != "" ) flavors = 1;
      for (int a = flavors + 1 ; a < histTokens->GetEntriesFast(); a++)	  {
	  name += histTokens->At(a)->GetName();
	  if (a != histTokens->GetEntriesFast() -1) name += "_";	 
      }
      
      // Apply NFs
      if (name.Contains("0j"))
      {
	samplesFolder->setScaleFactor(name, NFTop0jHiPT, "bkg/?/top"); // Top 0J

	samplesFolder->setScaleFactor(name, sf_Ztautau_0jet, "bkg/?/Zjets/Z/tt"); // Ztt 0j
	samplesFolder->setScaleFactor(name, sf_Ztautau_0jet, "bkg/?/Zjets/DY/tt"); // Ztt 0j	
      }
      else
      {
	samplesFolder->setScaleFactor(name, sf_Ztautau_incl, "bkg/?/Zjets/Z/tt"); // Ztt incl
        samplesFolder->setScaleFactor(name, sf_Ztautau_incl, "bkg/?/Zjets/DY/tt"); // Ztt incl
      }


      TString path = paths[iSample];
      if (TString(histTokens->At(flavors)->GetName()).Contains("em"))
	path.ReplaceAll("X", "em");
      else if (TString(histTokens->At(flavors)->GetName()) == "me")
	path.ReplaceAll("X", "me");
      else if (TString(histTokens->At(flavors)->GetName()) == "OF")
	path = path.ReplaceAll("X", "em") + " + " + TString(paths[iSample]).ReplaceAll("X", "me");
      else if (TString(histTokens->At(flavors)->GetName()) == "SF")
	path = path.ReplaceAll("X", "ee") + " + " + TString(paths[iSample]).ReplaceAll("X", "mm");

      path.ReplaceAll("MASS", mass);

      //If we are dealing with a "pass fRecoil" histogram, we need to recover the "inefficiency" for MC rates to match the stat-code definition
      if (histograms[iHist].find("frec") != std::string::npos) {
	//We "cancel" the efficiency but only for MC, not data
	if (samples[iSample] != "data") {
	  //get the two counters corresponding to the efficiency
	  TQCounter * countAftCut = samplesFolder->getCounter(path.Data(), name.Data());
	  if (!countAftCut) {
	    std::cout << "FAILED to get samplesFolder->getCounter(" << path.Data() << ", " << name.Data() << ")" << std::endl;
	    continue;
	  }
	  float aftCut = countAftCut->getCounter();
	  TString nameBefCut = name;
	  nameBefCut.ReplaceAll("frec","");
	  TQCounter * countBefCut = samplesFolder->getCounter(path.Data(), nameBefCut.Data());
	  if (!countBefCut) {
	    std::cout << "FAILED to get samplesFolder->getCounter(" << path.Data() << ", " << nameBefCut.Data() << ")" << std::endl;
	    continue;
	  }
	  float befCut = countBefCut->getCounter();
	  //computing the efficiency
	  float eff = ((aftCut) ? befCut/aftCut : 1.);
	  //applying the 1/eff as a scale factor to the "pass fRecoil" region to emulate the "before fRecoil region
	  samplesFolder->setScaleFactor(name.Data(), 1./eff, path.Data());
	}//end of if !data
      }//end of fRecoil eff cancellation

      TString compname = name + "/MT_Remapped";
      //TString compname = name + "/MT_TrackHWW_Clj_Remapped"; // ecfa change -- to be deleted later
//      std::cout << " Getting samplesFolder->getHistogram(" << path.Data() << ", " << compname.Data() << ")" << std::endl;
      TH1 * tempHist = samplesFolder->getHistogram(path.Data(), compname.Data());

      //If we are dealing with a "before fRecoil" histogram, we need to get the "fail fRecoil" data to match the stat-code definition
      if (histograms[iHist].find("ASR") != std::string::npos || histograms[iHist].find("CZpeak") != std::string::npos) {
	//We create the "fail fRecoil" histogram only for data
	if (samples[iSample] == "data") {
	  //get the two "bef" and "pass" histograms to do the subtraction
	  TH1 * befCutHist = tempHist;
	  if (!befCutHist) continue;
	  TString compnameAftCut = compname;
	  compnameAftCut.ReplaceAll("ASR","AfrecSR");
	  compnameAftCut.ReplaceAll("CZpeak","CfrecZpeak");
	  TH1 * aftCutHist = samplesFolder->getHistogram(path.Data(), compnameAftCut.Data());
	  if (!aftCutHist) continue;
	  //doing the subtraction to get the "fail fRecoil" histogram for data
	  befCutHist->SetBinContent(1, befCutHist->GetBinContent(1) - aftCutHist->Integral());
	}//end of if data
      }//end of fRecoil eff cancellation

      if (tempHist) {
	TString newname = histTokens->At(flavors)->GetName();
	newname += "_"; newname += name;
	tempHist->SetName(newname);
	std::cout << "   written " << newname.Data() << std::endl;
	tempHist->SetDirectory(tempOutput);
	tempHist->Write();
      }
      else {
	std::cout << "FAILED: Getting samplesFolder->getHistogram(" << path.Data() << ", " << compname.Data() << ")" << std::endl;
      }
    }


    tempOutput->Write();
    tempOutput->Close();
  }//end of loop over samples
  std::cout << "Finished caf2statHists treatment" << std::endl;
  return 0;
}


