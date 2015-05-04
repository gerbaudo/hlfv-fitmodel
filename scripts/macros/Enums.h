#ifndef _ENUM_H_
#define _ENUM_H_

#include <string>


//-----------------------------------------------------------------------------------------
//Discriminant variable-related enums
enum MappingModes { RemappingFlatBkg=2, RemappingOptimizedSignificance=5};

enum TransvMassDef{ TransvMassDef_MT=0, TransvMassDef_MT_TrackHWW_Cl, TransvMassDef_MT_TrackHWW_Clj, n_TransvMassDefs };
static const std::string TransvMassDef_CAFnames[] = { "MT", "MT_TrackHWW_Cl", "MT_TrackHWW_Clj" };



//-----------------------------------------------------------------------------------------
//Contains definitions for Wjets split uncertainties
//Note that these can be split "per EL/MU" with the corresponding string appended to all names
enum WjetsFakeSysts {WjetsFakeSysts_SysFakeStat_10_15=0, WjetsFakeSysts_SysFakeStat_15_20, WjetsFakeSysts_SysFakeStat_20_25, WjetsFakeSysts_SysFakeStat_GT25, 
		     WjetsFakeSysts_SysFakeFlav_10_15, WjetsFakeSysts_SysFakeFlav_15_20, WjetsFakeSysts_SysFakeFlav_20_25, WjetsFakeSysts_SysFakeFlav_GT25, 
		     WjetsFakeSysts_SysFakeOther_10_15, WjetsFakeSysts_SysFakeOther_15_20, WjetsFakeSysts_SysFakeOther_20_25, WjetsFakeSysts_SysFakeOther_GT25, 
		     n_WjetsFakeSysts};
enum WjetsFakeSysts_OSSS {WjetsFakeSysts_SysFakeOSSS_10_15=0, WjetsFakeSysts_SysFakeOSSS_15_20, WjetsFakeSysts_SysFakeOSSS_20_25, WjetsFakeSysts_SysFakeOSSS_GT25, n_WjetsFakeSysts_OSSS};
static const std::string WjetsFakeSysts_CAFnames[] = {"SysFakeStat_10_15", "SysFakeStat_15_20", "SysFakeStat_20_25", "SysFakeStat_GT25",
						      "SysFakeFlav_10_15", "SysFakeFlav_15_20", "SysFakeFlav_20_25", "SysFakeFlav_GT25",
						      "SysFakeOther_10_15", "SysFakeOther_15_20", "SysFakeOther_20_25", "SysFakeOther_GT25"};
static const std::string WjetsFakeSysts_OSSS_CAFnames[] = {"SysFakeOSSS_10_15", "SysFakeOSSS_15_20", "SysFakeOSSS_20_25", "SysFakeOSSS_GT25"};
static const std::string WjetsFakeSysts_NPnames[] =  {"FakeRate_Stat_10_15_HWW",  "FakeRate_Stat_15_20_HWW",  "FakeRate_Stat_20_25_HWW",  "FakeRate_Stat_GT25_HWW",
						      "FakeRate_Flav_10_15_HWW",  "FakeRate_Flav_15_20_HWW",  "FakeRate_Flav_20_25_HWW",  "FakeRate_Flav_GT25_HWW",
						      "FakeRate_Other_10_15_HWW", "FakeRate_Other_15_20_HWW", "FakeRate_Other_20_25_HWW", "FakeRate_Other_GT25_HWW"};
static const std::string WjetsFakeSysts_OSSS_NPnames[] =  {"FakeRate_OSSS_10_15_HWW",  "FakeRate_OSSS_15_20_HWW",  "FakeRate_OSSS_20_25_HWW",  "FakeRate_OSSS_GT25_HWW"};



//-----------------------------------------------------------------------------------------
//Contains enums for the split of electron uncertainties
//Number of bin indices and indices to use as sources of "uncorrelated" electron uncertainties
static const unsigned int nElecIDSystsBins = 12;
static unsigned int elecIDBinsList[] ={90020, 90025, 100007, 100008, 100012, 100013, 100017, 100018, 100027, 100028, 100032, 100033};
static const unsigned int nElecRecoSystsBins = 21;
static unsigned int elecRecoBinsList[] ={90020, 90025, 100000, 100013, 100025, 100037, 100049, 100073, 100085, 100097, 100109, 100121, 100133, 100145, 100157, 100169, 100181, 100205, 100217, 100229, 100241};
static const unsigned int nElecRecoIDSystsBins = 2;
static unsigned int elecRecoIDBinsList[] ={80015, 80010};


#endif
