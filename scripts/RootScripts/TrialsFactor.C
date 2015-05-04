void TrialsFactor(Int_t crossings=1, Double_t crossing_significance=0, Double_t local_significance=2.5){
  TF1 *f00 = new TF1("Trials factor","exp(-(pow(x,2)-pow(0.0,2))/2)/(ROOT::Math::chisquared_cdf_c(pow(x,2),1)/2)",2.0,8);
  TF1 *f05 = new TF1("Trials factor","exp(-(pow(x,2)-pow(0.5,2))/2)/(ROOT::Math::chisquared_cdf_c(pow(x,2),1)/2)",2.5,8);
  TF1 *f10 = new TF1("Trials factor","exp(-(pow(x,2)-pow(1.0,2))/2)/(ROOT::Math::chisquared_cdf_c(pow(x,2),1)/2)",3.0,8);
  TF1 *f15 = new TF1("Trials factor","exp(-(pow(x,2)-pow(1.5,2))/2)/(ROOT::Math::chisquared_cdf_c(pow(x,2),1)/2)",3.5,8);
  // TF1 *f20 = new TF1("Trials factor","exp(-(pow(x,2)-pow(2.0,2))/2)/(ROOT::Math::chisquared_cdf_c(pow(x,2),1)/2)",4.0,8);

  f00->SetMaximum(50);
  f00->SetMinimum(5);

  f00->SetLineColor(1);
  f05->SetLineColor(2);
  f10->SetLineColor(3);
  f15->SetLineColor(4);
  //f20->SetLineColor(6);

  f00->Draw();
  f05->Draw("same");
  f10->Draw("same");
  f15->Draw("same");
  // f20->Draw("same");

  c1->SetGrid();
  //c1->SetLogx();

  TLegend *legend = new TLegend(0.65,0.45,0.85,0.85);
  legend->AddEntry(f00,"Per 0.0 #\sigma crossing","l");
  legend->AddEntry(f05,"Per 0.5 #\sigma crossing","l");
  legend->AddEntry(f10,"Per 1.0 #\sigma crossing","l");
  legend->AddEntry(f15,"Per 1.5 #\sigma crossing","l");
  // legend->AddEntry(f20,"Per 2.0 #\sigma crossing","l");
  legend->Draw();

  //f1->SetTitle("Trials factor per crossing;Signficance;Trials factor/crossing");
  // f00->SetTitle("exp(-(x^2-S^2)/2)/(ROOT::Math::chisquared_cdf_c(x^2,1)/2);Local signficance (#sigma);(Trials factor-1)/crossing");
  f00->SetTitle("p_{0}^{global} = p_{0}^{local} + <N>*k, k=(Trials factor-1)/crossing;Local signficance (#sigma);(Trials factor-1)/crossing");
  
  Double_t local_p0 = ROOT::Math::chisquared_cdf_c(pow(local_significance,2),1)/2;
  Double_t global_p0 = local_p0 + crossings * exp(-(pow(local_significance,2)-pow(crossing_significance,2))/2);
  Double_t trials_factor = global_p0/local_p0;
  Double_t global_significance = sqrt(ROOT::Math::chisquared_quantile_c(global_p0*2,1));
  cout << "Local  p-value, significance: " <<  local_p0 << ", " << local_significance << " sigma" << endl;
  cout << "Trials factor: " << trials_factor << " (Trials factor minus 1: " << trials_factor-1 << ")" << endl;
  cout << "Global p-value, significance: " <<   global_p0 << ", " << global_significance << " sigma" << endl;
}
