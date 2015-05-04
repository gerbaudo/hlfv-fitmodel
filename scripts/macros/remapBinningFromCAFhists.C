#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>

#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

#include "macros/optimizeBinning.C"

#include "macros/Enums.h"

#include "HWWAnalysisCode/include/TQSampleFolder.h"

#include "TROOT.h"
#include "TInterpreter.h"

#pragma link C++ class multiset<SampleEvent>;

bool doUseSysForRemapping=0;

//Gets SampleEvent objects from hist and put them into rates
void fillSampleEventsFromHist(TH1 * hist, SampleEvent * sampleType, multiset<SampleEvent> & rates) {
  if (!hist) {
    std::cout << "ERROR: histogram not found during rebinning procedure for " << sampleType->name << std::endl;
    return;
  }
  SampleEvent s(sampleType);
  for (int ibin=1; ibin<=hist->GetNbinsX(); ibin++) {
    s.mt = hist->GetBinCenter(ibin);
    s.w = hist->GetBinContent(ibin);
    s.w_err = hist->GetBinError(ibin)/s.w;
    if (s.w) 
      rates.insert(s);
  }
}

//Returns the full list of CAF paths to select for the relevant channels
std::string getFullListOfPaths(std::vector<std::string> paths, std::string channels) {
  std::string channelsCAFname="";
  for (unsigned int ipath=0; ipath<paths.size(); ipath++) {
    TString path = paths[ipath];
    if (channels=="em")
      path.ReplaceAll("X", "em");
    else if (channels=="me")
      path.ReplaceAll("X", "me");
    else if (channels=="OF")
      path = path.ReplaceAll("X", "em") + " + " + TString(paths[ipath]).ReplaceAll("X", "me");
    else if (channels=="SF")
      path = path.ReplaceAll("X", "ee") + " + " + TString(paths[ipath]).ReplaceAll("X", "mm");
    channelsCAFname += ((ipath) ? " + ": "") + path;
  }
  return channelsCAFname;
}


//Fills rates_s and rates_b with the content of the histograms histname for channels channels stored in CAF object samples
void getSampleEventRatesFromCAF(TQSampleFolder * samples, std::string channels, std::string histname, multiset<SampleEvent> & rates_s, multiset<SampleEvent> & rates_b) {
    std::vector<std::string> paths;
    paths.push_back(std::string("sig/X"));
    std::string channelsCAFname = getFullListOfPaths(paths, channels);
    SampleEvent sampleSig;
    sampleSig.name = "";
    sampleSig.eps = 0;
    sampleSig.sysClass = -999;
    TH1 * hist = samples->getHistogram(channelsCAFname,histname);
    fillSampleEventsFromHist(hist, &sampleSig, rates_s);

    paths.clear();
    paths.push_back(std::string("bkg/X/Zjets"));
    channelsCAFname = getFullListOfPaths(paths, channels);
    SampleEvent sampleZDY;
    sampleZDY.name = "zdy";
    sampleZDY.eps = (doUseSysForRemapping) ? 0.35 : 0.;
    sampleZDY.sysClass = 3;
    hist = samples->getHistogram(channelsCAFname,histname);
    fillSampleEventsFromHist(hist, &sampleZDY, rates_b);

    paths.clear();
    paths.push_back(std::string("bkg/X/diboson/WW"));
    channelsCAFname = getFullListOfPaths(paths, channels);
    SampleEvent sampleWW;
    sampleWW.name = "ww";
    sampleWW.eps = (doUseSysForRemapping) ? 0.05 : 0.;
    sampleWW.sysClass = 0;
    hist = samples->getHistogram(channelsCAFname,histname);
    fillSampleEventsFromHist(hist, &sampleWW, rates_b);

    paths.clear();
    paths.push_back(std::string("bkg/X/top"));
    channelsCAFname = getFullListOfPaths(paths, channels);
    SampleEvent sampleTop;
    sampleTop.name = "top";
    sampleTop.eps = (doUseSysForRemapping) ? 0.05 : 0.;
    sampleTop.sysClass = 0;
    hist = samples->getHistogram(channelsCAFname,histname);
    fillSampleEventsFromHist(hist, &sampleTop, rates_b);

    paths.clear();
    paths.push_back(std::string("bkg/X/Wjets"));
    channelsCAFname = getFullListOfPaths(paths, channels);
    SampleEvent sampleWjets;
    sampleWjets.name = "wjets";
    sampleWjets.eps = (doUseSysForRemapping) ? 0.1 : 0.;
    sampleWjets.sysClass = 1;
    hist = samples->getHistogram(channelsCAFname,histname);
    fillSampleEventsFromHist(hist, &sampleWjets, rates_b);

    paths.clear();
    paths.push_back(std::string("bkg/X/diboson/NonWW"));
    channelsCAFname = getFullListOfPaths(paths, channels);
    SampleEvent sampleDib;
    sampleDib.name = "nonww";
    sampleDib.eps = (doUseSysForRemapping) ? 0.1 : 0.;
    sampleDib.sysClass = 2;
    hist = samples->getHistogram(channelsCAFname,histname);
    fillSampleEventsFromHist(hist, &sampleDib, rates_b);

}


//Returns the histogram definition corresponding to the passed parameters
std::string writeHistDef(std::string region, std::string chan, int nbins, std::string formula) {
  std::string histDef ="# ================================================================================= \n"
                       "	% cuts = '" + region + "*', channel = '" + chan + "' \n"
                       "# ================================================================================= \n"
                       "TH1F('MT_Remapped',''," + TString::Format("%i",nbins).Data() + ",0.,1.) << ( ( " + formula + " ) : 'Mapped m_{T} [arbitrary]' ) \n"
                       "\n\n";
  return histDef;
}


//Remaps the binning of the histograms listed in rebinningListFile from "workspace" version
//Remaps with flat bkg if mappingMode==RemappingFlatBkg
//Remaps to optimize significance if mappingMode==RemappingOptimizedSignificance
void remapBinningFromCAFhists(std::string rebinningListFile, std::string version, std::string anaPath = "HWWAnalysisCode/analysis/HWWlvlv_2012") {
  //Opening the file containing the "unbinned" histograms
  gSystem->Load("HWWAnalysisCode/lib/libQFramework.so");
  std::string unbinnedhistfile=anaPath + "/samples_hww_analysis_Nominal_unbinnedhists_" + version + ".root";
  TFile f(unbinnedhistfile.c_str());

  //Opening the "histograms.txt" file 
  std::string rebinnedHistFile=std::string(anaPath + "/histograms_remappedbinning_" + version + ".txt");
  ifstream iRebinnedHistFile;
  iRebinnedHistFile.open (rebinnedHistFile.c_str(), ifstream::in);
  ofstream oRebinnedHistFile(rebinnedHistFile.c_str());

  //Reading the binning configuration file and deal with binning remapping
  std::cout << "Reading in rebinning configuration file: " << rebinningListFile << std::endl;
  ifstream inFile(rebinningListFile.c_str());
  if (inFile.fail())  {
    std::cout << "ERROR::Couldn't open file: " << rebinningListFile.c_str() << std::endl;
    exit(1);
  }

  //Will store regions for which we need to postpone the rebinning after the remapping of another hist has been done
  std::vector<std::string> ppRegions;
  std::vector<std::string> ppChans;
  std::vector<std::string> ppTakeFrom;
  std::map<std::string, unsigned int> ppNbins;
  std::map<std::string, std::string> ppFormulas;

  //For each remapping to do
  string channels;
  string region;
  string str_nbins;
  unsigned int remappingMode;
  unsigned int mTdef;
  while (!inFile.eof())  {
    inFile >> channels;
    inFile >> region;
    inFile >> str_nbins;
    inFile >> remappingMode;
    inFile >> mTdef;
    TString mtVariable = TransvMassDef_CAFnames[static_cast<TransvMassDef>(mTdef)];

    MappingModes mappingMode= ((remappingMode==2) ? RemappingFlatBkg : RemappingOptimizedSignificance);

    std::string channelsCAFname;
    std::string chan;
    if (channels=="OF") {
      channelsCAFname="bkg/em+bkg/me";
      chan="em,me";
    }
    else if (channels=="SF") {
      channelsCAFname="bkg/ee+bkg/mm";
      chan="ee,mm";
    }
    else {
      channelsCAFname="bkg/"+channels;
      chan=channels;
    }

    char* c; int nbins = strtol(str_nbins.c_str(), &c, 10); 
    if (*c) { //If the binning of this region has to be taken from the remapping of another hist, we postpone the treatment of it
      ppRegions.push_back(region); 
      ppChans.push_back(chan); 
      ppTakeFrom.push_back(str_nbins);
    }
    else { //Otherwise we do the remapping right now
      std::string formula="";
      if (nbins>1) {
	std::vector<Int_t> * borders = new std::vector<Int_t> ;
	TQSampleFolder * samples_Nominal = (TQSampleFolder*)f.Get("samples_Nominal");
	if (!samples_Nominal) {
	  std::cout << "ERROR: samples_Nominal could not be opened from " << unbinnedhistfile.c_str() << std::endl;
	  return;
	}
	TH1 * hist = NULL;
	if (mappingMode==RemappingFlatBkg) {
	  std::cout << "Flat Remapping" << std::endl;
	  int res = samples_Nominal->createRemappedHistogramFlat(region+"/" + mtVariable.Data(), nbins, channelsCAFname, region+"/"+mtVariable.Data()+"_Remapped", borders);
	  hist = samples_Nominal->getHistogram(channelsCAFname,region+"/"+mtVariable.Data()+"_Remapped");
	}
	else if (mappingMode==RemappingOptimizedSignificance) {
	  std::cout << "Optimized Remapping" << std::endl;
	  hist = samples_Nominal->getHistogram(channelsCAFname,region+"/"+mtVariable.Data());
          if (!hist) 
	    std::cout << "ERROR: Remapping of histogram will fail because " << region.c_str() << " does not exist !" << std::endl;
	  else {
	    if (hist->Integral()==0.) std::cout << "ERROR: Remapping of histogram will fail because " << region.c_str() << " is empty !" << std::endl;
	    //Transforms the unbinned histograms into "SampleEvent" objects
	    multiset<SampleEvent> rates_s, rates_b;
	    getSampleEventRatesFromCAF(samples_Nominal, channels, region+"/"+mtVariable.Data(), rates_s, rates_b);

	    //Obtains bins corresponding to the remapping
	    double Zinit[3]={0,0,0}, Zfin[3]={0,0,0};
	    std::vector<double> remappedBins[3];
	    for (unsigned short iMode=0; iMode<1; iMode++)
		optimizeBinning(iMode, rates_s, rates_b, Zinit[iMode], Zfin[iMode], remappedBins[iMode], nbins);
	    unsigned short bestMode=-1;
	    std::cout << "Zfin = {" << Zfin[0] << ", " << Zfin[1] << ", " << Zfin[2] << "}" << std::endl;
	    if (Zfin[2]>Zfin[1])
	      if (Zfin[2]>Zfin[0])
		bestMode=2;
	      else
		bestMode=0;
	    else
	      if (Zfin[1]>Zfin[0])
		bestMode=1;
	      else
		bestMode=0;
	    std::cout << "Best remapping corresponds to remappingMode==" << bestMode << std::endl;
	    Double_t xbins[nbins+1];
	    xbins[0]=0;
	    for (int ibin=0;ibin<(int)remappedBins[bestMode].size();ibin++) {
	      std::cout << "Got bin boundary at " << remappedBins[bestMode][ibin] << std::endl;
	      xbins[ibin+1] = remappedBins[bestMode][ibin];
	      borders->push_back(hist->FindBin(remappedBins[bestMode][ibin]));
	    }
	    xbins[nbins]=1e9;
	    //hist = new TH1F(std::string(region+"_"+mtVariable.Data()+"_Remapped").c_str(), std::string(region+"/"+mtVariable.Data()+"_Remapped").c_str(), (Int_t)nbins, xbins);
	  }
	}
	if (!hist || hist->Integral()==0.) 
	  std::cout << "Remapping of " << region+"/"+ mtVariable.Data()+" for " << channelsCAFname << " failed !" << std::endl;
	else { //printing the binning
	  hist = samples_Nominal->getHistogram(channelsCAFname,region+"/"+mtVariable.Data());
	  std::cout << "Remapping binning of " << region+"/"+ mtVariable.Data()+" for " << channelsCAFname << " is :" << std::endl;
	  formula = "-0.099 ";
	  for (unsigned int ibin=0; ibin<borders->size(); ibin++) {
	    float lowbound = ((ibin>0) ? hist->GetBinLowEdge(borders->at(ibin-1))*1000. : 0.);
	    float upbound = hist->GetBinLowEdge(borders->at(ibin))*1000.;
	    formula += "+";
	    formula += TString::Format("%.3f*(%s>%.0f)*(%s<%.0f)",(ibin+1)*(1./nbins),mtVariable.Data(),lowbound,mtVariable.Data(),upbound).Data();
	  }
	  formula += TString::Format(" +1.0*(%s>%.0f)",mtVariable.Data(),hist->GetBinLowEdge(borders->at(borders->size()-1))*1000.).Data();
	}
	borders->clear();
	delete borders;
      }//end of if remapping necessary (not a single bin)
      else
	formula="0.5";
   
      if (formula!="") {
	std::cout << formula << std::endl;
	//and dumping it into the histogram file      
	oRebinnedHistFile << writeHistDef(region, chan, nbins, formula);
      }

      //We keep the information about this remapping in case another region needs it
      std::string key = channels+"_"+region;
      ppNbins[key] = nbins;
      ppFormulas[key] = formula;
    }//end of if (direct) remapping to be done for region
  }//end of loop of rebinningList file entries

  //Now dealing with the regions which have been postponed to reuse the remapping from another hist
  for (int iR=0; iR<ppRegions.size(); iR++) {
    std::cout << "Reusing mapping from " << ppTakeFrom[iR].c_str() << " for " << ppChans[iR] << " " << ppRegions[iR] << std::endl;
    oRebinnedHistFile << writeHistDef(ppRegions[iR], ppChans[iR], ppNbins[ppTakeFrom[iR]], ppFormulas[ppTakeFrom[iR]]);
  }

  std::cout << "Remapped binning written in " << rebinnedHistFile.c_str() << std::endl;
  //closing the histogram.txt file
  oRebinnedHistFile.close();
  iRebinnedHistFile.close();
}
