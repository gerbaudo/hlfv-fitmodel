


#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TString.h"
#include "TRandom.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMath.h"
#include "TGraph.h"

#include <iostream>
#include <fstream>
#include <string>


double H0 = 70.0;
double c =  9.7156034880*TMath::Power(10,-15);

double Eta(double a,double Lam)
{
	double s = TMath::Power((1-Lam)/Lam,(1/3));
	double ugly = ((1/(a*a*a*a))-0.1540*(s/(a*a*a))+0.19097*(s*s*s/a)+0.066941*s*s*s*s);
	return 2*TMath::Power((s*s*s+1),0.5)*TMath::Power(ugly,(-1/8));
}

double D_L(double z,double Lambda)
{
	double eta1 = Eta(1,Lambda);
	double eta2 = Eta(1/(1+z),Lambda);

	return (c/H0)*(1+z)*(eta1-eta2);

}


double Mu(double z,double L)
{
	double D = D_L(z,L);
	return 25 + 5*TMath::Log10(D);
}


void Plot_mu_z()
{
	double Lambda_m[] = {0.2,0.3,0.4,0.5};
	double mu[5][100];
	double z[100];

	for (int i=0;i<5;i++)
	{
		for (int j = 0; j<100; j++)
		{
			z[j] = 0 + 0.02*j;
			mu[i][j] = Mu(z,Lambda_m[i]);
		}

	}

	TGraph* g1 = new TGraph(100,z,mu[1][]);

	g1->Draw();

}

void ReadData(TString filename)
{

	string line;
	ifstream myfile ("example.txt");
	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{
			cout << line << '\n';
		}
		myfile.close();
	}

	else cout << "Unable to open file";

	return 0;

}



