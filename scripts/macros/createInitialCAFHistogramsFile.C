#ifndef _CREATE_INITIAL_CAF_HISTFILE_C_
#define _CREATE_INITIAL_CAF_HISTFILE_C_

// Author: O.Arnaez & J.Long
// Date:   2013-10-16
//
// Description:
// Creates the initial histograms_StatPlots.txt file to feed into CAF in order to obtain the "unbinned" discriminant plots

#include <iostream>
#include <fstream>
#include <string>

#include "TString.h"

#include "macros/Enums.h"


void createInitialCAFHistogramsFile(std::string rebinningListFile, std::string version) {
  //Opening the "histograms.txt" file 
  std::string initialHistFile=std::string("bsub/" + version + "/histograms_statPlots_" + version + ".txt");
  ifstream iInitialHistFile;
  iInitialHistFile.open (initialHistFile.c_str(), ifstream::in);
  ofstream oInitialHistFile(initialHistFile.c_str());

  oInitialHistFile <<
"# ================================================================================= \n"
"# Author: O.Arnaez & J.Long                                                         \n"
"#										     \n"
"# This files defines the histograms/plots to be produced by the anaylsis code when  \n"
"# producing histograms to feed into the stat code. The syntax follows the usual     \n"
"# HWWAnalysisCode rules.							     \n"
"#										     \n"
"# Lines starting with \"#\" will be ignored (comments).			     \n"
"# Lines starting with \"%\" are control lines to define parameters:		     \n"
"#   cuts = ...     defines the list of cuts the following histograms will be	     \n"
"#                  filled at (e.g. cuts = 'Cut_0jet*' for all cuts Cut_0jet	     \n"
"#                  downstream).						     \n"
"#   channel = ...  defines the lepton flavour channels the histograms will be	     \n"
"#                  filled for (e.g. channel = 'ee,mm' for ee and mm channel).	     \n"
"# ================================================================================= \n" << std::endl;

  //Reading the binning configuration file and deal with binning remapping
  std::cout << "Reading in rebinning configuration file: " << rebinningListFile << std::endl;
  ifstream inFile(rebinningListFile.c_str());
  if (inFile.fail())  {
    std::cout << "ERROR::Couldn't open file: " << rebinningListFile.c_str() << std::endl;
    exit(1);
  }
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
    char* c; int nbins = strtol(str_nbins.c_str(), &c, 10); if (*c) nbins=1;

    std::string chan;
    if (channels=="OF") {
      chan="em,me";
    }
    else if (channels=="SF") {
      chan="ee,mm";
    }
    else {
      chan=channels;
    }

    std::string formula="";
    if (nbins>1) {
      nbins=5000;
      formula = TString::Format("%s /1000.",mtVariable.Data()).Data();
    }//end of if remapping necessary (not a single bin)
    else
      formula="0.5";
   
    std::cout << formula << std::endl;

    //and dumping it into the histogram file      
    oInitialHistFile << "# ================================================================================= \n";
    oInitialHistFile << "	% cuts = '" << region << "*', channel = '" << chan << "' \n";
    oInitialHistFile << "# ================================================================================= \n";
    oInitialHistFile << "TH1F('" << mtVariable.Data() << "',''," << nbins << ",0.," << ((nbins>1) ? "500." : "1.") << ") << ( ( " << formula << " ) : 'Original m_{T} [GeV]' ) \n";
    oInitialHistFile << "\n\n";

  }

  std::cout << "Initial histograms definitions written in " << initialHistFile.c_str() << std::endl;
  //closing the histogram.txt file
  oInitialHistFile.close();

}



#endif
