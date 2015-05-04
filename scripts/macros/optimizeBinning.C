#include <iostream>
#include <set>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TH1D.h"

#include "macros/setup.C"

#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"


using namespace std;


struct SampleEvent
{
  SampleEvent()
  {

  }

  SampleEvent(string _name, string _type, int _sysClass, double _eps)
  {
    name=_name;
    type=_type;
    sysClass=_sysClass;
    eps=_eps;
  }

  bool operator<(const SampleEvent& rhs) const
  {
    return mt < rhs.mt;
  }

  string name;
  string type;
  int sysClass;
  double eps;
  double mt;
  double w;
};

void mapFlatBackground(multiset<SampleEvent>& rates_s, multiset<SampleEvent>& rates_b, 
		       vector<multiset<SampleEvent> >& binned_rates_s, vector<multiset<SampleEvent> >& binned_rates_b,
		       vector<double>& bins, int nrBins);

void mapFlatSignal(multiset<SampleEvent>& rates_s, multiset<SampleEvent>& rates_b, 
		   vector<multiset<SampleEvent> >& binned_rates_s, vector<multiset<SampleEvent> >& binned_rates_b,
		   vector<double>& bins, int nrBins);

void mapGausSignal(multiset<SampleEvent>& rates_s, multiset<SampleEvent>& rates_b, 
		   vector<multiset<SampleEvent> >& binned_rates_s, vector<multiset<SampleEvent> >& binned_rates_b,
		   vector<double>& bins, int nrBins);

void fillBins(multiset<SampleEvent>& rates_s, multiset<SampleEvent>& rates_b, 
	      vector<multiset<SampleEvent> >& binned_rates_s, vector<multiset<SampleEvent> >& binned_rates_b,
	      vector<double>& bins);

void optimizeBinning(multiset<SampleEvent>& rates_s, multiset<SampleEvent>& rates_b, 
		     vector<multiset<SampleEvent> >& binned_rates_s, vector<multiset<SampleEvent> >& binned_rates_b,
		     vector<double>& bins, double& Zinit, double& Zfin, int iitr);
void optimizeBinning(int _mappingMode, multiset<SampleEvent>& rates_s, multiset<SampleEvent>& rates_b, double& Zinit, double& Zfin, vector<double>& bins, int nrBins);

double getZdiff(double totSi, double errSi_sq, double totBi, double errBi_sq, 
		double totSi1, double errSi1_sq, double totBi1, double errBi1_sq, 
		double totSi_orig, double errSi_sq_orig, double totBi_orig, double errBi_sq_orig, 
		double totSi1_orig, double errSi1_sq_orig, double totBi1_orig, double errBi1_sq_orig);
double getZ(double s, double b);
double getHalfZ2(double s1, double b1);
double getdHalfZ2dS(double s, double b);
double getdHalfZ2dB(double s, double b);
double getSigmaHalfZ2(double s, double b, double sigma_s, double sigma_b);
double getSum(multiset<SampleEvent>& rates, double& sum_err_sq);
void printBins(vector<double> bins);

double getEps(multiset<SampleEvent>& rates, map<int, double>& bg_map, map<int, double>& eps_map);
double getEps(multiset<SampleEvent>& rates);
TGraph* plotMetric(multiset<pair<double, double> >& metric);

string dtos(double d)
{
  stringstream s;
  s << d;
  return s.str();
}

TRandom3 rndm2;

bool checkMetric(double totSi, double totBi, double errSi_sq, double errBi_sq, 
		 double totSi1, double totBi1, double errSi1_sq, double errBi1_sq, 
		 double totSi_orig, double totBi_orig, double errSi_sq_orig, double errBi_sq_orig, 
		 double totSi1_orig, double totBi1_orig, double errSi1_sq_orig, double errBi1_sq_orig, 
		 double& best_metric1, double& best_metric2, double& best_metric, 
		 double& best_sigmaHalfZ2_1, double& best_sigmaHalfZ2_2, double& best_sigmaHalfZ2, 
		 double& best_x, double x, double Zsig,
		 multiset<pair<double, double> >& set_metric1, multiset<pair<double, double> >& set_metric2, multiset<pair<double, double> >& set_metric,
		 double* useThisMetric1, double* useThisMetric2, double* useThisMetric);

//void optimizeBinning(string fileName, double& Zinit, double& Zfin);

// int mappingMode = 0;
// int _nrBins = 5;

bool drawPlots = 1;
bool useSys = 1;
string __channelName = "";
string _version = "";
int prescale = 0;

vector<double> binned_Zinit;
vector<double> binned_Zfin;

void optimizeBinning(int _mappingMode, string _channelName, double& Zinit, double& Zfin, vector<double>& bins, int nrBins)
{
  __channelName = _channelName;

  stringstream fileName;
  fileName << "opt/" << __channelName << ".root";
  TFile f(fileName.str().c_str());
  TTree* tree = (TTree*)f.Get("Opt");

  double mt,w;
  bool isS,isB;
  string* sampleName = new string;
  tree->SetBranchAddress("mt",&mt);
  tree->SetBranchAddress("w",&w);
  tree->SetBranchAddress("isS",&isS);
  tree->SetBranchAddress("isB",&isB);
  tree->SetBranchAddress("sampleName",&sampleName);

  rndm2.SetSeed(1);

  multiset<SampleEvent> rates_s, rates_b;
  int nrEntries = tree->GetEntries();
  for (int entry=0;entry<nrEntries;entry++)
  {
    tree->GetEntry(entry);

    if (*sampleName == "zleplep" || *sampleName == "ztautau") continue;
    //if (*sampleName == "wjets") w *= 2;
    //if (sampleName == "wg")
    SampleEvent sample;
    sample.name = "";
    sample.eps = 0;
    sample.sysClass = -999;
    if (*sampleName == "ggww" || *sampleName == "qqww")
    {
      sample.name = "ww";
      sample.eps = 0.05;
      sample.sysClass = 0;
    }
    if (*sampleName == "st" || *sampleName == "ttbar")
    {
      sample.name = "top";
      sample.eps = 0.05;
      sample.sysClass = 0;
    }
    if (*sampleName == "wjets")
    {
      sample.name = "wjets";
      sample.eps = 0.1;
      sample.sysClass = 1;
    }
    if (*sampleName == "wg" || *sampleName == "wgs" || *sampleName == "wzzz")
    {
      sample.name = "nonww";
      sample.eps = 0.1;
      sample.sysClass = 2;
    }
    sample.mt = mt;
    sample.w = w;

    if (isS) rates_s.insert(sample);
    if (isB) 
    {
      rates_b.insert(sample);
    }
  }

  f.Close();

  optimizeBinning(_mappingMode, rates_s, rates_b, Zinit, Zfin, bins, nrBins);
}

void optimizeBinning(int _mappingMode, string fileName, int nrBins)
{
  double Zinit, Zfin;
  vector<double> bins;
  optimizeBinning(_mappingMode, fileName, Zinit, Zfin, bins, nrBins);
}



void optimizeBinning(int _mappingMode = 2, int nrBins = 5, string version = "")// Nina changed to 1 for testing
{
//   mappingMode = _mappingMode;
//   _nrBins = nrBins;
  _version=version;

  vector<string> _fileNames;
  _fileNames.push_back("em_signalLike1_0j");
  _fileNames.push_back("em_signalLike2_0j");
  _fileNames.push_back("me_signalLike1_0j");
  _fileNames.push_back("me_signalLike2_0j");
  _fileNames.push_back("em_signalLike1_1j");
  _fileNames.push_back("em_signalLike2_1j");
  _fileNames.push_back("me_signalLike1_1j");
  _fileNames.push_back("me_signalLike2_1j");
  int nrFiles = _fileNames.size();

//   vector<vector<double> > vec_bins;
  double Zinit=0, Zfin=0;
  for (int i=0;i<nrFiles;i++)
  {
    double thisZinit=0,thisZfin=0;

    vector<double> bins;

    optimizeBinning(_mappingMode, _fileNames[i], thisZinit, thisZfin, bins, nrBins);

    TH1D* hist_init = new TH1D(string(_fileNames[i]+"_init").c_str(),string(_fileNames[i]+"_init").c_str(),nrBins,0,nrBins);
    TH1D* hist_fin = new TH1D(string(_fileNames[i]+"_fin").c_str(),string(_fileNames[i]+"_fin").c_str(),nrBins,0,nrBins);
    for (int ibin=0;ibin<nrBins;ibin++)
    {
      hist_init->SetBinContent(ibin+1,binned_Zinit[ibin]);
      hist_fin->SetBinContent(ibin+1,binned_Zfin[ibin]);
    }

    TCanvas* c1 = new TCanvas("c1","c1",1024,768);

    hist_init->SetMinimum(0);
    hist_init->SetMaximum(max(hist_init->GetMaximum(),hist_fin->GetMaximum())*1.2);

    hist_init->Draw();
    hist_fin->SetLineColor(kRed);
    hist_fin->Draw("same");
    
    c1->SaveAs(string("plots/"+_version+"/"+_fileNames[i]+".eps").c_str());

    Zinit += thisZinit*thisZinit;
    Zfin += thisZfin*thisZfin;
//     vec_bins.push_back(bins);
  }

  Zinit = sqrt(Zinit);
  Zfin = sqrt(Zfin);
  cout << "Zinit = " << Zinit << ", Zfin = " << Zfin << endl;
}


void optimizeBinning(int _mappingMode, multiset<SampleEvent>& rates_s, multiset<SampleEvent>& rates_b, double& Zinit, double& Zfin, vector<double>& bins, int nrBins)
{
//  int nrBins = _nrBins;
  int nrItr = 10;

  //if (fileName.find("1j") != string::npos) nrBins -= 2;

  //string regionName = fileName;
  //fileName = "opt/" + fileName + ".root";

  vector<multiset<SampleEvent> > binned_rates_s, binned_rates_b;

//initial mapping
  if (_mappingMode == 0) mapFlatSignal(rates_s, rates_b, binned_rates_s, binned_rates_b, bins, nrBins);
  else if (_mappingMode == 1) mapFlatBackground(rates_s, rates_b, binned_rates_s, binned_rates_b, bins, nrBins);
  else if (_mappingMode == 2) mapGausSignal(rates_s, rates_b, binned_rates_s, binned_rates_b, bins, nrBins);
  printBins(bins);

  for (int i=0;i<nrItr;i++)
  {
    double thisZinit=0,thisZfin=0;
    cout << "---->Starting iteration " << i << endl;
    optimizeBinning(rates_s, rates_b, binned_rates_s, binned_rates_b, bins, thisZinit, thisZfin, i);
    if (i == 0) Zinit = thisZinit;
    if (i == nrItr-1) Zfin = thisZfin;
//     fillBins(rates_s, rates_b, binned_rates_s, binned_rates_b, bins);
    printBins(bins);
    if (thisZinit == thisZfin) 
    {
      Zfin = thisZfin;
      break;
    }
//     optimizeBinning(binned_rates_s, binned_rates_b, bins, 1);
//     fillBins(rates_s, rates_b, binned_rates_s, binned_rates_b, bins);
//     printBins(bins);
  }

//try moving bins to different places 

//   vector<double> originalBins = bins;
//   for (int ibin=0;ibin<nrBins;ibin++)
//   {
//     cout << "Moving bin: " << ibin << endl;
//     double this_Z0 = Z0;
//     double this_Z = Z0;
//     for (int jbin=0;jbin<nrBins;jbin++)
//     {
      
//     }
//     for (int iitr=0;iitr<nrItr;iitr++)
//     {
//       cout << "---->Starting iteration " << iitr << endl;
//       double thisZinit=0,thisZfin=0;
//       optimizeBinning(rates_s, rates_b, binned_rates_s, binned_rates_b, bins, thisZinit, thisZfin, i);
//       printBins(bins);
//       if (i == 0) this_Z0 = thisZinit;
//       if (i == nrItr-1) this_Z = thisZfin;
//       if (thisZinit == thisZfin) 
//       {
// 	this_Z = thisZfin;
// 	break;
//       }
//     }

//     if (
//   }




//   ofstream outFile((string(fileName)+".txt").c_str());
//   outFile << regionName;
//   for (int i=0;i<(int)bins.size();i++)
//   {
//     outFile << " " << bins[i];
//   }
//   outFile << "\n";
//   outFile.close();


}

void optimizeBinning(multiset<SampleEvent>& rates_s, multiset<SampleEvent>& rates_b, 
		     vector<multiset<SampleEvent> >& binned_rates_s, vector<multiset<SampleEvent> >& binned_rates_b,
		     vector<double>& bins, double& Zinit, double& Zfin, int iitr)
{
  int nrBins = bins.size();
  double Zsig = 5;





  if (iitr == 0) binned_Zinit.clear();
  binned_Zfin.clear();

  double Ztot_init = 0;
  if (useSys)
  {
    double Vinv_00 = 0;
    double Vinv_01 = 0;
    double Vinv_11 = 1;
    for (int i=0;i<=nrBins;i++)
    {
      multiset<SampleEvent>& these_rates_s  = binned_rates_s[i];
      multiset<SampleEvent>& these_rates_b  = binned_rates_b[i];
    
      double eps_eff = getEps(these_rates_b);

      double errs, errb;
      double s = getSum(these_rates_s, errs);
      double b = getSum(these_rates_b, errb);
    
      Vinv_00 += s*s/(s+b);
      Vinv_11 += b*b*eps_eff*eps_eff/(s+b);
      Vinv_01 += s*b*eps_eff/(s+b);

      if (iitr == 0) binned_Zinit.push_back(getZ(s,b));
    }

    double Vinvinv_00 = Vinv_11 / (Vinv_00*Vinv_11-Vinv_01*Vinv_01);
    Ztot_init = 1./Vinvinv_00;
  }
  else
  {
    for (int i=0;i<=nrBins;i++)
    {
      multiset<SampleEvent>& these_rates_s  = binned_rates_s[i];
      multiset<SampleEvent>& these_rates_b  = binned_rates_b[i];

      double errs, errb;
      double s = getSum(these_rates_s, errs);
      double b = getSum(these_rates_b, errb);

      Ztot_init += pow(getZ(s, b), 2);

      if (iitr == 0) binned_Zinit.push_back(getZ(s,b));
    }
  }
  Ztot_init = sqrt(Ztot_init);
  Zinit = Ztot_init;










  for (int odd = 0; odd <= 1; odd++)
  {
    vector<double> new_bins;
    for (int ibin=0;ibin<nrBins;ibin++)
    {
      if (odd && double(ibin / 2) == ibin / 2.)
      {
	new_bins.push_back(bins[ibin]);
	continue;
      }
      if (!odd && !(double(ibin / 2) == ibin / 2.))
      {
	new_bins.push_back(bins[ibin]);
	continue;
      }

      multiset<SampleEvent>& rates_s_i  = binned_rates_s[ibin];
      multiset<SampleEvent>& rates_s_i1 = binned_rates_s[ibin+1];

      multiset<SampleEvent>& rates_b_i  = binned_rates_b[ibin];
      multiset<SampleEvent>& rates_b_i1 = binned_rates_b[ibin+1];

      multiset<pair<double, double> > set_metric,set_metric1,set_metric2;

      double xi = bins[ibin];

      double errSi_sq, errBi_sq, errSi1_sq, errBi1_sq;
      double totSi  = getSum(rates_s_i, errSi_sq);
      double totSi1 = getSum(rates_s_i1, errSi1_sq);
      double totBi  = getSum(rates_b_i, errBi_sq);
      double totBi1 = getSum(rates_b_i1, errBi1_sq);

      map<int, double> Bi_map,Bi1_map;
      map<int, double> eps_map;
      double eps_eff = getEps(rates_b_i, Bi_map, eps_map);
      double eps_eff1 = getEps(rates_b_i1, Bi1_map, eps_map);

      double best_eps_eff = eps_eff;
      double best_eps_eff1 = eps_eff1;

      double Vinv_00 = totSi*totSi/(totSi+totBi) + totSi1*totSi1/(totSi1+totBi1);
      double Vinv_11 = eps_eff*eps_eff*totBi*totBi/(totSi+totBi) + eps_eff1*eps_eff1*totBi1*totBi1/(totSi1+totBi1)+1;
      double Vinv_01 = totSi*totBi*eps_eff/(totSi+totBi)+totSi1*totBi1*eps_eff1/(totSi1+totBi1);
      
      double Vinvinv_00 = Vinv_11 / (Vinv_00*Vinv_11-Vinv_01*Vinv_01);

      double best_metric1 = getZ(totSi, totBi);
      double best_metric2 = getZ(totSi1, totBi1);
      double best_metric = 0.5*(pow(best_metric1, 2) + pow(best_metric2, 2));

      if (useSys) best_metric = sqrt(1./Vinvinv_00);

      double best_x = xi;
      double best_sigmaHalfZ2_1 = getSigmaHalfZ2(totSi, totBi, sqrt(errSi_sq), sqrt(errBi_sq));
      double best_sigmaHalfZ2_2 = getSigmaHalfZ2(totSi1, totBi1, sqrt(errSi1_sq), sqrt(errBi1_sq));
      double best_sigmaHalfZ2 = sqrt(pow(best_sigmaHalfZ2_1, 2) + pow(best_sigmaHalfZ2_2, 2));
      cout << "First metric (" << ibin << ") = " << best_metric << " +/- " << best_sigmaHalfZ2 << ", m1 = " << 0.5*best_metric1*best_metric1 << " +/- " << best_sigmaHalfZ2_1 << ", m2 = " << 0.5*best_metric2*best_metric2 << " +/- " << best_sigmaHalfZ2_2 << ", x = " << best_x << ", rb = " << totBi1 / (totBi1 + totBi) << ", rs = " << totSi1 / (totSi1 + totSi) << endl;

//       Ztot_init += best_metric;

      double best_totSi=totSi, best_totSi1=totSi1, best_totBi=totBi, best_totBi1=totBi1;
      double best_errSi_sq=errSi_sq, best_errSi1_sq=errSi1_sq, best_errBi_sq=errBi_sq, best_errBi1_sq=errBi1_sq;

      multiset<SampleEvent>::iterator itr_si1  = rates_s_i1.begin();
      multiset<SampleEvent>::iterator itr_si1E = rates_s_i1.end();
      multiset<SampleEvent>::iterator itr_bi1  = rates_b_i1.begin();
      multiset<SampleEvent>::iterator itr_bi1E = rates_b_i1.end();

      map<int, double> Bi_map_orig = Bi_map;
      map<int, double> Bi1_map_orig = Bi1_map;
      double totSi_orig  = totSi;
      double totSi1_orig = totSi1;
      double totBi_orig  = totBi;
      double totBi1_orig = totBi1;
      double errSi_sq_orig  = errSi_sq;
      double errSi1_sq_orig = errSi1_sq;
      double errBi_sq_orig  = errBi_sq;
      double errBi1_sq_orig = errBi1_sq;
      while (itr_si1 != itr_si1E && itr_bi1 != itr_bi1E)
      {
	double x;
	if (itr_si1->mt < itr_bi1->mt)
	{
	  double w = itr_si1->w;
	  totSi += w;
	  totSi1 -= w;
	  errSi_sq += w*w;
	  errSi1_sq -= w*w;
	  x = itr_si1->mt;
	  itr_si1++;	  
	}
	else if (itr_si1->mt == itr_bi1->mt)
	{
	  double w = itr_si1->w;
	  totSi += w;
	  totSi1 -= w;
	  errSi_sq += w*w;
	  errSi1_sq -= w*w;
	  x = itr_si1->mt;
	  itr_si1++;

	  w = itr_bi1->w;
	  totBi += w;
	  totBi1 -= w;
	  errBi_sq += w*w;
	  errBi1_sq -= w*w;
	  x = itr_bi1->mt;
	  Bi_map[itr_bi1->sysClass] += w;
	  Bi1_map[itr_bi1->sysClass] -= w;
	  itr_bi1++;
	}
	else
	{
	  double w = itr_bi1->w;
	  totBi += w;
	  totBi1 -= w;
	  errBi_sq += w*w;
	  errBi1_sq -= w*w;
	  x = itr_bi1->mt;
	  Bi_map[itr_bi1->sysClass] += w;
	  Bi1_map[itr_bi1->sysClass] -= w;
	  itr_bi1++;
	}

	eps_eff = 0;
	for (map<int, double>::iterator itr=Bi_map.begin();itr!=Bi_map.end();itr++)
	{
	  eps_eff += itr->second*itr->second*eps_map[itr->first]*eps_map[itr->first];
	}
	eps_eff = sqrt(eps_eff)/totBi;

	eps_eff1 = 0;
	for (map<int, double>::iterator itr=Bi1_map.begin();itr!=Bi1_map.end();itr++)
	{
	  eps_eff1 += itr->second*itr->second*eps_map[itr->first]*eps_map[itr->first];
	}
	eps_eff1 = sqrt(eps_eff1)/totBi1;
	
	Vinv_00 = totSi*totSi/(totSi+totBi) + totSi1*totSi1/(totSi1+totBi1);
	Vinv_11 = eps_eff*eps_eff*totBi*totBi/(totSi+totBi) + eps_eff1*eps_eff1*totBi1*totBi1/(totSi1+totBi1)+1;
	Vinv_01 = totSi*totBi*eps_eff/(totSi+totBi)+totSi1*totBi1*eps_eff1/(totSi1+totBi1);
	
	Vinvinv_00 = Vinv_11 / (Vinv_00*Vinv_11-Vinv_01*Vinv_01);

	double useThisMetric = sqrt(1./Vinvinv_00);
	double useThisMetric1 = totSi/sqrt(totSi+totBi+eps_eff*eps_eff*totBi*totBi);
	double useThisMetric2 = totSi1/sqrt(totSi1+totBi1+eps_eff1*eps_eff1*totBi1*totBi1);

	
	double uni = rndm2.Uniform(0,prescale);
	bool pass = uni > 1 ? 0 : checkMetric(totSi, totBi, errSi_sq, errBi_sq, totSi1, totBi1, errSi1_sq, errBi1_sq,
					      totSi_orig, totBi_orig, errSi_sq_orig, errBi_sq_orig, totSi1_orig, totBi1_orig, errSi1_sq_orig, errBi1_sq_orig, 
					      best_metric1, best_metric2, best_metric, 
					      best_sigmaHalfZ2_1, best_sigmaHalfZ2_2, best_sigmaHalfZ2, best_x, x, Zsig,
					      set_metric1, set_metric2, set_metric, 
					      useSys ? &useThisMetric1 : NULL, useSys ? &useThisMetric2 : NULL, useSys ? &useThisMetric : NULL);

	if (pass)
	{
	  best_totSi = totSi;
	  best_totSi1 = totSi1;
	  best_totBi = totBi;
	  best_totBi1 = totBi1;

	  best_errSi_sq = errSi_sq;
	  best_errSi1_sq = errSi1_sq;
	  best_errBi_sq = errBi_sq;
	  best_errBi1_sq = errBi1_sq;

	  best_eps_eff = eps_eff;
	  best_eps_eff1 = eps_eff1;
	}

//       double metric = 0.5*(pow(getZ(totSi, totBi), 2) + pow(getZ(totSi1, totBi1), 2));
//       double sigmaHalfZ2 = sqrt(pow(getSigmaHalfZ2(totSi, totBi, sqrt(errSi_sq), sqrt(errBi_sq)), 2) + pow(getSigmaHalfZ2(totSi1, totBi1, sqrt(errSi1_sq), sqrt(errBi1_sq)), 2));

//       if (ibin == 1) cout << "metric = " << metric << " +/- " << sigmaHalfZ2 << ", x = " << x << endl;
//       //double rb = totBi1 / (totBi1 + totBi);
//       if ((metric > best_metric && metric / sigmaHalfZ2 > Zsig) || best_metric / best_sigmaHalfZ2 < Zsig)
//       {
// 	best_metric = metric;
// 	best_sigmaHalfZ2 = sigmaHalfZ2;
// // 	cout << "Si = " << totSi << ", Bi = " << totBi << ", Si1 = " << totSi1 << ", Bi1 = " << totBi1 << endl;
// // 	cout << "Si + Si1 = " << totSi+totSi1 << ", Bi + Bi1 = " << totBi + totBi1 << endl;
//  	//cout << "metric = " << metric << " +/- " << sigmaHalfZ2 << ", x = " << x << endl;
// 	best_x = x;
//       }
      }

    
    

      multiset<SampleEvent>::reverse_iterator itr_si  = rates_s_i.rbegin();
      multiset<SampleEvent>::reverse_iterator itr_siE = rates_s_i.rend();
      multiset<SampleEvent>::reverse_iterator itr_bi  = rates_b_i.rbegin();
      multiset<SampleEvent>::reverse_iterator itr_biE = rates_b_i.rend();
    
      Bi_map = Bi_map_orig;
      Bi1_map = Bi1_map_orig;
      totSi = totSi_orig;
      totSi1 = totSi1_orig;
      totBi = totBi_orig;
      totBi1 = totBi1_orig;
      errSi_sq = errSi_sq_orig;
      errSi1_sq = errSi1_sq_orig;
      errBi_sq = errBi_sq_orig;
      errBi1_sq = errBi1_sq_orig;
      while (itr_si != itr_siE && itr_bi != itr_biE)
      {
	double x;
	if (itr_si->mt > itr_bi->mt)
	{
	  double w = itr_si->w;
	  totSi -= w;
	  totSi1 += w;
	  x = itr_si->mt;
	  itr_si++;
	}
	else if (itr_si->mt == itr_bi->mt)
	{
	  double w = itr_si->w;
	  totSi -= w;
	  totSi1 += w;
	  x = itr_si->mt;
	  itr_si++;


	  w = itr_bi->w;
	  totBi -= w;
	  totBi1 += w;
	  x = itr_bi->mt;
	  Bi_map[itr_bi->sysClass] -= w;
	  Bi1_map[itr_bi->sysClass] += w;
	  itr_bi++;
	}
	else
	{
	  double w = itr_bi->w;
	  totBi -= w;
	  totBi1 += w;
	  x = itr_bi->mt;
	  Bi_map[itr_bi->sysClass] -= w;
	  Bi1_map[itr_bi->sysClass] += w;
	  itr_bi++;
	}



	eps_eff = 0;
	for (map<int, double>::iterator itr=Bi_map.begin();itr!=Bi_map.end();itr++)
	{
	  eps_eff += itr->second*itr->second*eps_map[itr->first]*eps_map[itr->first];
	}
	eps_eff = sqrt(eps_eff)/totBi;

	eps_eff1 = 0;
	for (map<int, double>::iterator itr=Bi1_map.begin();itr!=Bi1_map.end();itr++)
	{
	  eps_eff1 += itr->second*itr->second*eps_map[itr->first]*eps_map[itr->first];
	}
	eps_eff1 = sqrt(eps_eff1)/totBi1;
	
	Vinv_00 = totSi*totSi/(totSi+totBi) + totSi1*totSi1/(totSi1+totBi1);
	Vinv_11 = eps_eff*eps_eff*totBi*totBi/(totSi+totBi) + eps_eff1*eps_eff1*totBi1*totBi1/(totSi1+totBi1)+1;
	Vinv_01 = totSi*totBi*eps_eff/(totSi+totBi)+totSi1*totBi1*eps_eff1/(totSi1+totBi1);
	
	Vinvinv_00 = Vinv_11 / (Vinv_00*Vinv_11-Vinv_01*Vinv_01);

	double useThisMetric = sqrt(1./Vinvinv_00);
	double useThisMetric1 = totSi/sqrt(totSi+totBi+eps_eff*eps_eff*totBi*totBi);
	double useThisMetric2 = totSi1/sqrt(totSi1+totBi1+eps_eff1*eps_eff1*totBi1*totBi1);

	double uni = rndm2.Uniform(0,prescale);
	bool pass = uni > 1 ? 0 : checkMetric(totSi, totBi, errSi_sq, errBi_sq, totSi1, totBi1, errSi1_sq, errBi1_sq, 
					      totSi_orig, totBi_orig, errSi_sq_orig, errBi_sq_orig, totSi1_orig, totBi1_orig, errSi1_sq_orig, errBi1_sq_orig, 
					      best_metric1, best_metric2, best_metric, 
					      best_sigmaHalfZ2_1, best_sigmaHalfZ2_2, best_sigmaHalfZ2, best_x, x, Zsig,
					      set_metric1, set_metric2, set_metric, 
					      useSys ? &useThisMetric1 : NULL, useSys ? &useThisMetric2 : NULL, useSys ? &useThisMetric : NULL);

	if (pass)
	{
	  best_totSi = totSi;
	  best_totSi1 = totSi1;
	  best_totBi = totBi;
	  best_totBi1 = totBi1;

	  best_errSi_sq = errSi_sq;
	  best_errSi1_sq = errSi1_sq;
	  best_errBi_sq = errBi_sq;
	  best_errBi1_sq = errBi1_sq;

	  best_eps_eff = eps_eff;
	  best_eps_eff1 = eps_eff1;
	}

//       double metric = 0.5*(pow(getZ(totSi, totBi), 2) + pow(getZ(totSi1, totBi1), 2));
//       double sigmaHalfZ2 = sqrt(pow(getSigmaHalfZ2(totSi, totBi, sqrt(errSi_sq), sqrt(errBi_sq)), 2) + pow(getSigmaHalfZ2(totSi1, totBi1, sqrt(errSi1_sq), sqrt(errBi1_sq)), 2));

//       //cout << "metric = " << metric << " +/- " << sigmaHalfZ2 << endl;
//       if ((metric > best_metric && metric / sigmaHalfZ2 > Zsig) || best_metric / best_sigmaHalfZ2 < Zsig)
//       {
// 	best_metric = metric;
// 	best_sigmaHalfZ2 = sigmaHalfZ2;
// // 	cout << "Si = " << totSi << ", Bi = " << totBi << ", Si1 = " << totSi1 << ", Bi1 = " << totBi1 << endl;
// // 	cout << "Si + Si1 = " << totSi+totSi1 << ", Bi + Bi1 = " << totBi + totBi1 << endl;
// //  	cout << "metric = " << metric << " +/- " << sigmaHalfZ2 << ", x = " << x << endl;
// 	best_x = x;
//       }
      }


      double Zdiff = getZdiff(best_totSi, best_errSi_sq, best_totBi, best_errBi_sq, 
			      best_totSi1, best_errSi1_sq, best_totBi1, best_errBi1_sq, 
			      totSi_orig, errSi_sq_orig, totBi_orig, errBi_sq_orig, 
			      totSi1_orig, errSi1_sq_orig, totBi1_orig, errBi1_sq_orig);

      double Z = sqrt(pow(getZ(best_totSi, best_totBi), 2) + pow(getZ(best_totSi1, best_totBi1), 2));
      double Z0 = sqrt(pow(getZ(totSi_orig, totBi_orig), 2) + pow(getZ(totSi1_orig, totBi1_orig), 2));
      double dZ = sqrt(Z*Z-Z0*Z0);
      double dZ2 = Z-Z0;


      eps_eff = 0;
      for (map<int, double>::iterator itr=Bi_map.begin();itr!=Bi_map.end();itr++)
      {
	eps_eff += itr->second*itr->second*eps_map[itr->first]*eps_map[itr->first];
      }
      eps_eff = sqrt(eps_eff)/totBi;
      
      eps_eff1 = 0;
      for (map<int, double>::iterator itr=Bi1_map.begin();itr!=Bi1_map.end();itr++)
      {
	eps_eff1 += itr->second*itr->second*eps_map[itr->first]*eps_map[itr->first];
      }
      eps_eff1 = sqrt(eps_eff1)/totBi;
      
      Vinv_00 = best_totSi*best_totSi/(best_totSi+best_totBi) + best_totSi1*best_totSi1/(best_totSi1+best_totBi1);
      Vinv_11 = eps_eff*eps_eff*best_totBi*best_totBi/(best_totSi+best_totBi) + eps_eff1*eps_eff1*best_totBi1*best_totBi1/(best_totSi1+best_totBi1)+1;
      Vinv_01 = best_totSi*best_totBi*eps_eff/(best_totSi+best_totBi)+best_totSi1*best_totBi1*eps_eff1/(best_totSi1+best_totBi1);
      
      Vinvinv_00 = Vinv_11 / (Vinv_00*Vinv_11-Vinv_01*Vinv_01);
      
      double useThisMetric = sqrt(1./Vinvinv_00);
      double useThisMetric1 = best_totSi/sqrt(best_totSi+best_totBi+eps_eff*eps_eff*best_totBi*best_totBi);
      double useThisMetric2 = best_totSi1/sqrt(best_totSi1+best_totBi1+eps_eff1*eps_eff1*best_totBi1*best_totBi1);

      totSi  = getSum(rates_s_i, errSi_sq);
      totSi1 = getSum(rates_s_i1, errSi1_sq);
      totBi  = getSum(rates_b_i, errBi_sq);
      totBi1 = getSum(rates_b_i1, errBi1_sq);

      cout << "Last metric (" << ibin << ") = " << best_metric << " +/- " << best_sigmaHalfZ2 << ", m1 = " << 0.5*best_metric1*best_metric1 << " +/- " << best_sigmaHalfZ2_1 << ", m2 = " << 0.5*best_metric2*best_metric2 << " +/- " << best_sigmaHalfZ2_2 << ", x = " << best_x << ", rb = " << totBi1 / (totBi1 + totBi) << ", rs = " << totSi1 / (totSi1 + totSi) << endl;
  
      if (set_metric.size())
      {
	TGraph* plot1 = plotMetric(set_metric1);
	TGraph* plot2 = plotMetric(set_metric2);
	TGraph* plot = plotMetric(set_metric);

	double plot_max = 0;
	for (multiset<pair<double, double> >::iterator itr = set_metric.begin(); itr != set_metric.end(); itr++)
	{
	  plot_max = max(plot_max, itr->second);
	}

	double win_min = set_metric.size() ? set_metric.begin()->first : 0;
	double win_max = set_metric.size() ? set_metric.rbegin()->first : 0;

	TCanvas* c1 = new TCanvas("c1","c1",1024,768);


	stringstream plotName;
	plotName << "plots/" << _version << "/" << __channelName << "/";
	system((string("mkdir -vp ") + plotName.str()).c_str());
	plotName << "plot_bin" << ibin << "_itr" << iitr << ".eps";

	plot->SetTitle(plotName.str().c_str());      
	if (drawPlots) plot->Draw("al");

	double plot_min = min(plot1->GetMinimum(), plot2->GetMinimum());
	plot_min = min(plot_min, plot->GetMinimum());
	plot->SetMinimum(0);//plot_min);
	plot->SetMaximum(1.2*plot_max);//1.2*plot->GetMaximum());

      
	plot1->SetLineColor(kBlue);
	if (drawPlots) plot1->Draw("l");
      
	plot2->SetLineColor(kRed);
	if (drawPlots) plot2->Draw("l");

	TLine l;
	l.SetLineStyle(2);
	l.DrawLine(xi, 0, xi, 1.2*plot_max);

	l.SetLineColor(kGreen-2);
	l.DrawLine(best_x, 0, best_x, 1.2*plot_max);


	l.SetLineColor(kBlack);
	l.DrawLine(best_x, plot_max, xi, plot_max);

	if (Zdiff < 2) best_x = xi;
	else
	{
	  double dxrand = rndm2.Uniform(0.0*(win_max-win_min));
	  best_x = best_x + dxrand;
	  best_x = min(best_x, win_max);
	  best_x = max(best_x, win_min);
	}

	l.SetLineColor(kCyan);
	l.DrawLine(best_x, 0, best_x, 1.2*plot_max);


	stringstream tStr;
	tStr << setprecision(3);
	tStr << "#frac{Z^{2}-Z^{2}_{0}}{#sigma} = " << Zdiff << " (" << dZ << " +/- " << dZ/Zdiff << ", #Delta(Z) = " << dZ2 << ")";
	TLatex t;
	t.SetNDC();
	t.DrawLatex(0.175, 0.87, tStr.str().c_str());
      
	if (drawPlots) c1->SaveAs(plotName.str().c_str());
      }
//       Ztot_fin += best_metric;
      new_bins.push_back(best_x);
    }
    bins = new_bins;
    fillBins(rates_s, rates_b, binned_rates_s, binned_rates_b, bins);
  }





  double Ztot_fin = 0;
  if (useSys)
  {
    double Vinv_00 = 0;
    double Vinv_01 = 0;
    double Vinv_11 = 1;
    for (int i=0;i<=nrBins;i++)
    {
      multiset<SampleEvent>& these_rates_s  = binned_rates_s[i];
      multiset<SampleEvent>& these_rates_b  = binned_rates_b[i];
    
      double eps_eff = getEps(these_rates_b);

      double errs, errb;
      double s = getSum(these_rates_s, errs);
      double b = getSum(these_rates_b, errb);
    
      Vinv_00 += s*s/(s+b);
      Vinv_11 += b*b*eps_eff*eps_eff/(s+b);
      Vinv_01 += s*b*eps_eff/(s+b);

      binned_Zfin.push_back(getZ(s,b));
    }

    double Vinvinv_00 = Vinv_11 / (Vinv_00*Vinv_11-Vinv_01*Vinv_01);
  
    Ztot_fin = 1./Vinvinv_00;
  }
  else
  {
    for (int i=0;i<=nrBins;i++)
    {
      multiset<SampleEvent>& these_rates_s  = binned_rates_s[i];
      multiset<SampleEvent>& these_rates_b  = binned_rates_b[i];

      double errs, errb;
      double s = getSum(these_rates_s, errs);
      double b = getSum(these_rates_b, errb);

      Ztot_fin += pow(getZ(s, b), 2);

      binned_Zfin.push_back(getZ(s,b));
    }
  }
  Ztot_fin = sqrt(Ztot_fin);
  Zfin = Ztot_fin;

  cout << "Ztot_init = " << Ztot_init << ", Ztot_fin = " << Ztot_fin << endl;

}

double getZdiff(double totSi, double errSi_sq, double totBi, double errBi_sq, 
		double totSi1, double errSi1_sq, double totBi1, double errBi1_sq, 
		double totSi_orig, double errSi_sq_orig, double totBi_orig, double errBi_sq_orig, 
		double totSi1_orig, double errSi1_sq_orig, double totBi1_orig, double errBi1_sq_orig)
{
  double term1 = (pow(getdHalfZ2dS(totSi, totBi), 2) + pow(getdHalfZ2dS(totSi_orig, totBi_orig), 2))*errSi_sq_orig;
  double term2 = (pow(getdHalfZ2dS(totSi, totBi), 2) + pow(getdHalfZ2dS(totSi1_orig, totBi1_orig), 2))*fabs(errSi_sq-errSi_sq_orig);
  double term3 = (pow(getdHalfZ2dS(totSi1, totBi1), 2) + pow(getdHalfZ2dS(totSi1_orig, totBi1_orig), 2))*errSi1_sq;

  double term4 = (pow(getdHalfZ2dB(totSi, totBi), 2) + pow(getdHalfZ2dB(totSi_orig, totBi_orig), 2))*errBi_sq_orig;
  double term5 = (pow(getdHalfZ2dB(totSi, totBi), 2) + pow(getdHalfZ2dB(totSi1_orig, totBi1_orig), 2))*fabs(errBi_sq-errBi_sq_orig);
  double term6 = (pow(getdHalfZ2dB(totSi1, totBi1), 2) + pow(getdHalfZ2dB(totSi1_orig, totBi1_orig), 2))*errBi1_sq;

  double sigma = sqrt(term1*term1+term2*term2+term3*term3+term4*term4+term5*term5+term6*term6);
  //cout << "sigma = " << sigma << endl;

  double Z = sqrt(pow(getZ(totSi, totBi), 2) + pow(getZ(totSi1, totBi1), 2));
  double Z0 = sqrt(pow(getZ(totSi_orig, totBi_orig), 2) + pow(getZ(totSi1_orig, totBi1_orig), 2));
  //cout << "Z = " << Z << ", Z0 = " << Z0 << endl;
  double dZ = sqrt(Z*Z-Z0*Z0);
  double sigma_dZ = sigma/dZ;
  //cout << "dZ = " << dZ << " +/- " << sigma_dZ << endl;

  return sigma_dZ > 0 ? dZ/sigma_dZ : 0;
}

TGraph* plotMetric(multiset<pair<double, double> >& metric)
{
  int counter = 0;
  double* x = new double[metric.size()];
  double* y = new double[metric.size()];
  for (multiset<pair<double, double> >::iterator itr = metric.begin();itr!=metric.end();itr++)
  {
    x[counter] = itr->first;
    y[counter] = itr->second;
    counter++;
  }
  return new TGraph(metric.size(), x, y);
}


bool checkMetric(double totSi, double totBi, double errSi_sq, double errBi_sq,
		 double totSi1, double totBi1, double errSi1_sq, double errBi1_sq,
		 double totSi_orig, double totBi_orig, double errSi_sq_orig, double errBi_sq_orig,
		 double totSi1_orig, double totBi1_orig, double errSi1_sq_orig, double errBi1_sq_orig,
		 double& best_metric1, double& best_metric2, double& best_metric, 
		 double& best_sigmaHalfZ2_1, double& best_sigmaHalfZ2_2, double& best_sigmaHalfZ2, 
		 double& best_x, double x, double Zsig,
		 multiset<pair<double, double> >& set_metric1, multiset<pair<double, double> >& set_metric2, multiset<pair<double, double> >& set_metric,
		 double* useThisMetric1, double* useThisMetric2, double* useThisMetric)
{
  double z1 = getZ(totSi, totBi);
  double z2 = getZ(totSi1, totBi1);
  double halfZ2 = 0.5*(z1*z1+z2*z2);

  double Zdiff = getZdiff(totSi, errSi_sq, totBi, errBi_sq, 
			  totSi1, errSi1_sq, totBi1, errBi1_sq, 
			  totSi_orig, errSi_sq_orig, totBi_orig, errBi_sq_orig, 
			  totSi1_orig, errSi1_sq_orig, totBi1_orig, errBi1_sq_orig);
    
  double sigmaHalfZ2_1 = getSigmaHalfZ2(totSi, totBi, sqrt(errSi_sq), sqrt(errBi_sq));
  double sigmaHalfZ2_2 = getSigmaHalfZ2(totSi1, totBi1, sqrt(errSi1_sq), sqrt(errBi1_sq));
  double sigmaHalfZ2 = sqrt(pow(sigmaHalfZ2_1, 2) + pow(sigmaHalfZ2_2, 2));

  double sig1 = 0.5*z1*z1 / sigmaHalfZ2_1;
  double sig2 = 0.5*z2*z2 / sigmaHalfZ2_2;
  double sig  = halfZ2 / sigmaHalfZ2;

  double metric1 = useThisMetric1 ? *useThisMetric1 : z1;
  double metric2 = useThisMetric2 ? *useThisMetric2 : z2;
  double metric  = useThisMetric  ? *useThisMetric  : halfZ2;

  bool isSig = false;
  if (sig1 > Zsig && sig2 > Zsig && sig > Zsig)
  {
    isSig = true;
    set_metric1.insert(make_pair(x, metric1));
    set_metric2.insert(make_pair(x, metric2));
    set_metric.insert(make_pair(x, useThisMetric ? metric : sqrt(2*metric)));
  }

  //cout << "metric = " << metric << " +/- " << sigmaHalfZ2 << ", x = " << x << ", isSig ? " << isSig << ", best_metric = " << best_metric;// << endl;
  //if (ibin == 1) cout << "metric = " << metric << " +/- " << sigmaHalfZ2 << ", x = " << x << endl;
  //double rb = totBi1 / (totBi1 + totBi);
  if ((metric > best_metric && isSig) || 
      (/*best_metric / best_sigmaHalfZ2 < Zsig || */best_metric1 / best_sigmaHalfZ2_1 < Zsig || best_metric2 / best_sigmaHalfZ2_2 < Zsig))
  //if ((Zdiff > best_metric && metric / sigmaHalfZ2 > Zsig && 0.5*z1*z1 / sigmaHalfZ2_1 > Zsig && 0.5*z2*z2 / sigmaHalfZ2_2 > Zsig) || 
  //    (best_metric1 / best_sigmaHalfZ2_1 < Zsig || best_metric2 / best_sigmaHalfZ2_2 < Zsig))
  {
    best_metric = metric;
    best_sigmaHalfZ2 = sigmaHalfZ2;

    best_metric1 = 0.5*z1*z1;
    best_sigmaHalfZ2_1 = sigmaHalfZ2_1;

    best_metric2 = 0.5*z2*z2;
    best_sigmaHalfZ2_2 = sigmaHalfZ2_2;

// 	cout << "Si = " << totSi << ", Bi = " << totBi << ", Si1 = " << totSi1 << ", Bi1 = " << totBi1 << endl;
// 	cout << "Si + Si1 = " << totSi+totSi1 << ", Bi + Bi1 = " << totBi + totBi1 << endl;
    //cout << "metric = " << metric << " +/- " << sigmaHalfZ2 << ", x = " << x << endl;
    //cout << " passed" << endl;

    best_x = x;
    return true;
  }
  //cout << endl;

  return false;
}


double getZ(double s, double b)
{
  if (b <= 0) return 0;
  return sqrt(2*((s+b)*TMath::Log(1 + s/b) - s));
}

double getHalfZ2(double s, double b)
{
  if (b <= 0) return 0;
  return (s+b)*TMath::Log(1 + s/b) - s;
}

double getSigmaHalfZ2(double s, double b, double sigma_s, double sigma_b)
{
  if (b <= 0) return 10e9;
  double sigma = sqrt(pow(getdHalfZ2dS(s,b)*sigma_s, 2) + pow(getdHalfZ2dB(s,b)*sigma_b, 2));
  if (sigma != sigma) return 10e9;
  return sigma;
}

double getdHalfZ2dS(double s, double b)
{
  return TMath::Log(1+s/b);
}

double getdHalfZ2dB(double s, double b)
{
  return TMath::Log(1+s/b)-s/b;
}

void printBins(vector<double> bins)
{
  for (int i=0;i<(int)bins.size();i++)
  {
    cout << "Bin " << i << " = " << bins[i] << endl;
  }
  cout << endl;
}

double getSum(multiset<SampleEvent>& rates, double& sum_err)
{
  double tot = 0;
  sum_err = 0;
  for (multiset<SampleEvent>::iterator itr = rates.begin(); itr != rates.end();itr++)
  {
    tot += itr->w;
    sum_err += itr->w*itr->w;
  }

  return tot;
}

double getEps(multiset<SampleEvent>& rates, map<int, double>& bg_map, map<int, double>& eps_map)
{
  vector<int> sysClasses;
  for (multiset<SampleEvent>::iterator itr = rates.begin(); itr != rates.end();itr++)
  {
    if (bg_map.find(itr->sysClass) != bg_map.end())
    {
      eps_map[itr->sysClass] = itr->eps;
      bg_map[itr->sysClass] = 0;
      sysClasses.push_back(itr->sysClass);
    }
    bg_map[itr->sysClass] += itr->w;
  }

  double eps = 0;
  double bg = 0;
  for (int i=0;i<(int)sysClasses.size();i++)
  {
    int c = sysClasses[i];
    double b = bg_map[c];
    double e = eps_map[c];
    eps += e*e*b*b;
    bg += b;
  }

  return sqrt(eps)/bg;
}


double getEps(multiset<SampleEvent>& rates)
{
  map<int, double> bg_map, eps_map;
  return getEps(rates, bg_map, eps_map);
}


void mapFlatSignal(multiset<SampleEvent>& rates_s, multiset<SampleEvent>& rates_b, 
		   vector<multiset<SampleEvent> >& binned_rates_s, vector<multiset<SampleEvent> >& binned_rates_b,
		   vector<double>& bins, int nrBins)
{
  double tot = 0;

  for (multiset<SampleEvent>::iterator set_itr = rates_s.begin();set_itr!=rates_s.end();set_itr++)
  {
    tot += set_itr->w;
  }

  double sum_tot = 0;
  double nrPerBin = tot/nrBins;
  double totOrig = tot;
  tot = 0;
  //gaussian mapping stuff
  double cdf = 0;
  //bool isOdd = nrBins / 2 == nrBins / 2.;
//   double minZ = -2.5;//-(nrBins - isOdd - 1) / 2;
//   double maxZ = +2.5;
  //int maxZ = +(nrBins - isOdd - 1) / 2;

  multiset<SampleEvent> these_rates_s, these_rates_b;
  for (multiset<SampleEvent>::iterator set_itr = rates_s.begin();set_itr!=rates_s.end();set_itr++)
  {
    sum_tot += set_itr->w;
    tot += set_itr->w;
    cdf = sum_tot / totOrig;

//     if (mappingMode == 2 || mappingMode == 3)
//     {
//     cout << "mt = " << set_itr->first << ", w = " << set_itr->second << endl;
//     cout << "tot = " << tot << ", nrPerBin = " << nrPerBin << ", nrPerBin * nrBins = " << nrPerBin * nrBins << endl;
//     cout << "nrBins = " << nrBins << ", bins.size() = " << bins.size() << endl;
      if (tot > nrPerBin && int(bins.size()+1) < nrBins && set_itr->w >= 0 /* fu*%! */ && tot < nrPerBin*nrBins) // and find the mT points that satisfy that
      {
	//cout << "adding bin: " << set_itr->first << endl;
	bins.push_back(set_itr->mt);
	tot = 0;
      }
//     }
//     else if (mappingMode == 4)
//     {
//       double Z = ROOT::Math::gaussian_quantile(cdf, 1);
	      
	      
//       for (int ib=1;ib<nrBins;ib++)
//       {
// 	//if (Z - minZ + 1 >= ib && int(bins.size()) < ib)
// 	if (Z >= (maxZ - minZ) / nrBins * ib + minZ && int(bins.size()) < ib)
// 	{
// // 		  double boundary = double(ib)/nrBins;
// // 		  bins.push_back(boundary);
// 	  bins.insert(set_itr->first);
// // 		  cout << "boundary = " << boundary << endl;
// 	  break;
// 	}
//       }
//     }
  }

  fillBins(rates_s, rates_b, binned_rates_s, binned_rates_b, bins);


}

void mapGausSignal(multiset<SampleEvent>& rates_s, multiset<SampleEvent>& rates_b, 
		   vector<multiset<SampleEvent> >& binned_rates_s, vector<multiset<SampleEvent> >& binned_rates_b,
		   vector<double>& bins, int nrBins)
{
  double tot = 0;

  for (multiset<SampleEvent>::iterator set_itr = rates_s.begin();set_itr!=rates_s.end();set_itr++)
  {
    tot += set_itr->w;
  }

  double sum_tot = 0;
  //double nrPerBin = tot/nrBins;
  double totOrig = tot;
  tot = 0;
  //gaussian mapping stuff
  double cdf = 0;
  //bool isOdd = nrBins / 2 == nrBins / 2.;
//   double minZ = -2.5;//-(nrBins - isOdd - 1) / 2;
//   double maxZ = +2.5;
  //int maxZ = +(nrBins - isOdd - 1) / 2;

  multiset<SampleEvent> these_rates_s, these_rates_b;
  for (multiset<SampleEvent>::iterator set_itr = rates_s.begin();set_itr!=rates_s.end();set_itr++)
  {
    sum_tot += set_itr->w;
    tot += set_itr->w;
    cdf = sum_tot / totOrig;

//     if (mappingMode == 2 || mappingMode == 3)
//     {
//     cout << "mt = " << set_itr->first << ", w = " << set_itr->second << endl;
//     cout << "tot = " << tot << ", nrPerBin = " << nrPerBin << ", nrPerBin * nrBins = " << nrPerBin * nrBins << endl;
//     cout << "nrBins = " << nrBins << ", bins.size() = " << bins.size() << endl;

/*
      if (tot > nrPerBin && int(bins.size()+1) < nrBins && set_itr->second >= 0 && tot < nrPerBin*nrBins) // and find the mT points that satisfy that
      {
	cout << "adding bin: " << set_itr->first << endl;
	bins.push_back(set_itr->first);
	tot = 0;
      }
*/

//     }
//     else if (mappingMode == 4)
//     {
    double Z = ROOT::Math::gaussian_quantile(cdf, 1);
    
    double minZ = -2.5;
    double maxZ = +2.5;
    
    for (int ib=1;ib<nrBins;ib++)
    {
      //if (Z - minZ + 1 >= ib && int(bins.size()) < ib)
      if (Z >= (maxZ - minZ) / nrBins * ib + minZ && int(bins.size()) < ib)
      {
	// 		  double boundary = double(ib)/nrBins;
	// 		  bins.push_back(boundary);
	bins.push_back(set_itr->mt);
	// 		  cout << "boundary = " << boundary << endl;
	break;
      }
    }
  }


  fillBins(rates_s, rates_b, binned_rates_s, binned_rates_b, bins);


}

void mapFlatBackground(multiset<SampleEvent>& rates_s, multiset<SampleEvent>& rates_b, 
		   vector<multiset<SampleEvent> >& binned_rates_s, vector<multiset<SampleEvent> >& binned_rates_b,
		   vector<double>& bins, int nrBins)
{
  double tot = 0;

  for (multiset<SampleEvent>::iterator set_itr = rates_b.begin();set_itr!=rates_b.end();set_itr++)
  {
    tot += set_itr->w;
  }

  double sum_tot = 0;
  double nrPerBin = tot/nrBins;
  double totOrig = tot;
  tot = 0;
  //gaussian mapping stuff
  double cdf = 0;
  //bool isOdd = nrBins / 2 == nrBins / 2.;
//   double minZ = -2.5;//-(nrBins - isOdd - 1) / 2;
//   double maxZ = +2.5;
  //int maxZ = +(nrBins - isOdd - 1) / 2;

  multiset<SampleEvent> these_rates_s, these_rates_b;
  for (multiset<SampleEvent>::iterator set_itr = rates_b.begin();set_itr!=rates_b.end();set_itr++)
  {
    sum_tot += set_itr->w;
    tot += set_itr->w;
    cdf = sum_tot / totOrig;

//     if (mappingMode == 2 || mappingMode == 3)
//     {
//     cout << "mt = " << set_itr->first << ", w = " << set_itr->second << endl;
//     cout << "tot = " << tot << ", nrPerBin = " << nrPerBin << ", nrPerBin * nrBins = " << nrPerBin * nrBins << endl;
//     cout << "nrBins = " << nrBins << ", bins.size() = " << bins.size() << endl;
      if (tot > nrPerBin && int(bins.size()+1) < nrBins && set_itr->w >= 0 /* fu*%! */ && tot < nrPerBin*nrBins) // and find the mT points that satisfy that
      {
	//cout << "adding bin: " << set_itr->first << endl;
	bins.push_back(set_itr->mt);
	tot = 0;
      }
//     }
//     else if (mappingMode == 4)
//     {
//       double Z = ROOT::Math::gaussian_quantile(cdf, 1);
	      
	      
//       for (int ib=1;ib<nrBins;ib++)
//       {
// 	//if (Z - minZ + 1 >= ib && int(bins.size()) < ib)
// 	if (Z >= (maxZ - minZ) / nrBins * ib + minZ && int(bins.size()) < ib)
// 	{
// // 		  double boundary = double(ib)/nrBins;
// // 		  bins.push_back(boundary);
// 	  bins.insert(set_itr->first);
// // 		  cout << "boundary = " << boundary << endl;
// 	  break;
// 	}
//       }
//     }
  }

  fillBins(rates_s, rates_b, binned_rates_s, binned_rates_b, bins);


}


void fillBins(multiset<SampleEvent>& rates_s, multiset<SampleEvent>& rates_b, 
	      vector<multiset<SampleEvent> >& binned_rates_s, vector<multiset<SampleEvent> >& binned_rates_b,
	      vector<double>& bins)
{
  int nrBins = bins.size();

  binned_rates_s.clear();
  binned_rates_b.clear();

  multiset<SampleEvent> rates;
  for (int i=0;i<=(int)bins.size();i++)
  {
    binned_rates_s.push_back(rates);
    binned_rates_b.push_back(rates);
  }

  int bin = 0;
  for (multiset<SampleEvent>::iterator set_itr = rates_s.begin();set_itr!=rates_s.end();set_itr++)
  {
    if (bin < nrBins && set_itr->mt > bins[bin]) bin++;
    binned_rates_s[bin].insert(*set_itr);
  }

  bin=0;
  for (multiset<SampleEvent>::iterator set_itr = rates_b.begin();set_itr!=rates_b.end();set_itr++)
  {
    if (bin < nrBins && set_itr->mt > bins[bin]) bin++;
    binned_rates_b[bin].insert(*set_itr);
  }
}

