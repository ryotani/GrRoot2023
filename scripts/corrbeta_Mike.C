void plot(char* name){
  //TFile *f = new TFile("hist/hist_all.root");
  //TFile *f = new TFile("/shared/storage/physnp/su562/Calibrations/Z_correction_6.1.2020/47Mn/Calfiles/Calhists/firsthalf/collated.root"); //46Cr
  //TFile *f = new TFile("/shared/storage/physnp/su562/Calibrations/46Mn/Calhists/collated.root");
TFile *f = new TFile("/shared/storage/physnp/mab503/GrROOT_Files/Calculate_out/47Ti/46Titest.root");

  TH1F *h[2];
  h[0] = (TH1F*)f->Get(Form("egamdc_AFW_%s",name));
  h[1] = (TH1F*)f->Get(Form("egamdc_ABW_%s",name));
  h[0]->Rebin(2);
  h[1]->Rebin(2);
  h[0]->Draw("");
  h[0]->SetMinimum(0);
  h[1]->Draw("same");
  h[0]->SetLineColor(2);
  h[1]->SetLineColor(4);
}

void correct(double e1, double e2, double beta){
  //double theta1 = 51.56;
  //double theta2 = 97.4;
  //double theta1 = 57.2958; //46V (b=0.29)
  //double theta2 = 96.6866279;
  //double theta1 = 57.2958; //47Mn (b=0.29)
  //double theta2 = 97.4028;
  //double theta1 = 57.86874; //46V (b=0.39)
  //double theta2 = 97.689304;
  //double theta1 = 57.86874; //47Mn (b=0.39)
  //double theta2 = 95.11099;
  //double theta1 = 56.436343; //46V (b=0.387)
  //double theta2 = 95.11099;
  //double theta1 = 58.155216; //47Mn (b=0.387)
  //double theta2 = 95.11099;
  //double theta1 = 57.86874; //44Ti (b=0.386)
  //double theta2 = 95.11099;
  //double theta1 = 57.2958; //47Ti (b=0.386)
  //double theta2 = 97.4028;
  double theta1 = 58.28; //45Sc/45Cr 
  double theta2 = 90.00;
  
  
  theta1 = cos(theta1*3.141592/180.);
  theta2 = cos(theta2*3.141592/180.);

  double gamma = 1./sqrt(1.-beta*beta);
  
  double EL1 = e1/(gamma*(1.-beta*theta1));
  double EL2 = e2/(gamma*(1.-beta*theta2));

  double d1 = beta*pow(gamma,3) * (1.-beta*theta1) - gamma *theta1;
  double d2 = beta*pow(gamma,3) * (1.-beta*theta2) - gamma *theta2;
cout << e1 <<"      " <<e2 << endl;
  cout << "d1 " << d1 << "\t" << "d2 " << d2 << endl;
  cout << "EL1 " << EL1 << "\t" << "EL2 " << EL2 << endl;
  cout << "t1 " << theta1 << "\t" << "t2 " << theta2 << endl;
  d1*=EL1;
  d2*=EL2;

  double diff = e1-e2;
  double dbeta = diff/(d2 - d1);
  cout << "--------------------------------------------------" << endl;
  cout << "dbeta " << dbeta << " new beta " << beta+dbeta << " corrected to " << e1+d1*dbeta << "  " << e2+d2*dbeta << endl;

}
