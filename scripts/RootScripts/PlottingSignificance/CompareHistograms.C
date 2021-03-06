/*
  Code for "Plotting the Differences Between Data and Expectation"
  by Georgios Choudalakis and Diego Casadei
  (http://arxiv.org/abs/1111.2062)

  -----------------------------------------------------------------

  Diego Casadei <casadei@cern.ch> 7 Dec 2011
 */





#include "pValuePoissonError.h"
#include "CompareHistograms.h"



/*
  Given two ROOT histograms (with the same binning!) containing the
  observed and expected counts, create and return a histogram showing
  the significance of their bin-to-bin discrepancies.

  If the histogram representing the expectation (second input
  parameter) has non-zero bin "errors", these are considered the
  standard deviations representing the full uncertainty and the
  significance is computed accordingly, unless this is disabled (third
  parameter).

  -----------------------------------------------------------------

  Diego Casadei <casadei@cern.ch> 7 Dec 2011
  (last update: 24 Feb 2012)
*/
TH1F* CompareHistograms(TH1* hObs, TH1* hExp,
			bool neglectUncertainty)
{
  if (hObs==0 || hExp==0) return 0;

  TString name=hObs->GetName();
  name+="_cmp_";
  name+=hExp->GetName();
  int Nbins = hObs->GetNbinsX();
  /*
  // this assumes constant binning
  TH1F* hOut = new TH1F(name, "",
			Nbins,
			hObs->GetXaxis()->GetXmin(),
			hObs->GetXaxis()->GetXmax());
  hOut->GetXaxis()->SetTitle( hObs->GetXaxis()->GetTitle() );
  hOut->GetYaxis()->SetTitle("significance");
  hOut->SetFillColor(2);
  */
  /*
  // OK for variable binning but ALL properties are cloned...
  TH1F* hOut = (TH1F*) hObs->Clone(name);
  hOut->Clear();
  hOut->GetYaxis()->SetTitle("significance");
  hOut->SetFillColor(2);
  */
  // variable binning
  TH1F* hOut = new TH1F(name, "",
			hObs->GetXaxis()->GetNbins(),
			hObs->GetXaxis()->GetXbins()->GetArray());
  hOut->GetXaxis()->SetTitle( hObs->GetXaxis()->GetTitle() );
  hOut->GetYaxis()->SetTitle("significance");
  hOut->SetFillColor(2);

  for (int i=1; i<=Nbins; ++i) { // SKIP UNDER- AND OVER-FLOWS
    if (hObs->GetBinContent(i)<=0) continue;

    unsigned nObs = (int) hObs->GetBinContent(i);
    float nExp = hExp->GetBinContent(i);
    float vrnc = hExp->GetBinError(i);
    vrnc *= vrnc; // variance
    float sig = 0;
    if (vrnc>0 && !neglectUncertainty) {
      // account for systematic uncertainty
      float pValue = pValuePoissonError(nObs, nExp, vrnc);
      if (pValue<0.5) sig = pValueToSignificance(pValue, (nObs>nExp));
    } else {
      // assume perfect knowledge of Poisson parameter
      float pValue = pValuePoisson(nObs,nExp);
      if (pValue<0.5) sig = pValueToSignificance(pValue, (nObs>nExp));
    }
    hOut->SetBinContent(i, sig);
  }

  return hOut;

}

