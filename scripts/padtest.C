char *mainch = "crdcpad0";
char *ddchar = "crdcmaxpad0";
Double_t fonegaus(Double_t *x, Double_t *par);
void PressEnterToContinue(){
  int c;
  printf( "Press ENTER to continue... " );
  fflush( stdout );
  do c = getchar(); 
  while ((c != '\n') && (c != EOF));
}
void doit(int cr){
  //TFile *f = new TFile("padcal.root");
  TFile *f = new TFile("padcal2.root");
  TGraph *g[224];
  TH1F *h[4][224];
  TF1 *fu[4][224];
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(3,2);
  for(int p=0;p<224;p++){
    if(p==70)
      continue;
    h[0][p] = (TH1F*)f->Get(Form("%s_%d_%d_in36Arout35Ar",mainch,cr,p));
    h[1][p] = (TH1F*)f->Get(Form("%s_%d_%d_in36Arout32S",mainch,cr,p));
    h[2][p] = (TH1F*)f->Get(Form("%s_%d_%d_in36Arout28Si",mainch,cr,p));
    fu[0][p] = (TF1*)f->Get(Form("f_%s_%d_%d_in36Arout35Ar",mainch,cr,p));
    fu[1][p] = (TF1*)f->Get(Form("f_%s_%d_%d_in36Arout32S",mainch,cr,p));
    fu[2][p] = (TF1*)f->Get(Form("f_%s_%d_%d_in36Arout28Si",mainch,cr,p));
    g[p] = (TGraph*)f->Get(Form("line_%d_%d",cr,p));
    for(int i=0;i<3;i++){
      c->cd(1+i);
      h[i][p]->Draw();
      fu[i][p]->Draw("same");
      h[i][p]->GetXaxis()->SetRangeUser(0,1000);
    }
    c->cd(5);
    g[p]->Draw("AP*");
    c->Modified();
    c->Update();
    PressEnterToContinue();

  }


}
void check(int cr, int p){
  TFile *f = new TFile("padcal2.root");
  //TFile *f = new TFile("padcal.root");
  TH1F *h[4];
  TF1 *fu[4];
  TGraph *g;
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(3,2);
  h[0] = (TH1F*)f->Get(Form("%s_%d_%d_in36Arout35Ar",mainch,cr,p));
  h[1] = (TH1F*)f->Get(Form("%s_%d_%d_in36Arout32S",mainch,cr,p));
  h[2] = (TH1F*)f->Get(Form("%s_%d_%d_in36Arout28Si",mainch,cr,p));
  fu[0] = (TF1*)f->Get(Form("f_%s_%d_%d_in36Arout35Ar",mainch,cr,p));
  fu[1] = (TF1*)f->Get(Form("f_%s_%d_%d_in36Arout32S",mainch,cr,p));
  fu[2] = (TF1*)f->Get(Form("f_%s_%d_%d_in36Arout28Si",mainch,cr,p));
  g = (TGraph*)f->Get(Form("line_%d_%d",cr,p));
  for(int i=0;i<3;i++){
    c->cd(1+i);
    h[i]->Draw();
    fu[i]->Draw("same");
    h[i]->GetXaxis()->SetRangeUser(0,1000);
  }
  c->cd(5);
  g->Draw("AP*");
}
void checkmatch(){
  TFile *f = new TFile("padcal2.root");
  //TFile *f = new TFile("padcal.root");
  TH1F *h[2][4];
  TF1 *fu[2][4];
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(4,2);
  for(int cr=0;cr<2;cr++){
    int p=70;
    h[cr][0] = (TH1F*)f->Get(Form("%s_%d_%d_in36Arout35Ar",mainch,cr,p));
    h[cr][1] = (TH1F*)f->Get(Form("%s_%d_%d_in36Arout32S",mainch,cr,p));
    h[cr][2] = (TH1F*)f->Get(Form("%s_%d_%d_in36Arout28Si",mainch,cr,p));
    fu[cr][0] = (TF1*)f->Get(Form("f_%s_%d_%d_in36Arout35Ar",mainch,cr,p));
    fu[cr][1] = (TF1*)f->Get(Form("f_%s_%d_%d_in36Arout32S",mainch,cr,p));
    fu[cr][2] = (TF1*)f->Get(Form("f_%s_%d_%d_in36Arout28Si",mainch,cr,p));
    for(int i=0;i<3;i++){
      c->cd(1+i+cr*4);
      h[cr][i]->Draw();
      fu[cr][i]->Draw("same");
     h[cr][i]->GetXaxis()->SetRangeUser(0,1000);
    }
  }
}
void checkall(){
  TFile *f = new TFile("his/his0056_new.root");
  TH2F *h[2][4];
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(4,2);
  for(int cr=0;cr<2;cr++){
    h[cr][0] = (TH2F*)f->Get(Form("%s_%d_in36Arout35Ar",ddchar,cr));
    h[cr][1] = (TH2F*)f->Get(Form("%s_%d_in36Arout32S",ddchar,cr));
    h[cr][2] = (TH2F*)f->Get(Form("%s_%d_in36Arout28Si",ddchar,cr));
    for(int i=0;i<3;i++){
      c->cd(1+i+cr*4);
      h[cr][i]->Draw("colz");
      h[cr][i]->GetYaxis()->SetRangeUser(0,1000);
    }
  }
}
void refit(int cr,int p){
  TFile *f = new TFile("padcal2.root");
  //TFile *f = new TFile("padcal.root");
  TF1* fit[5];
  //TH1F *h[4];
  TH1F *h[3];
  TGraph *g;
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(3,2);
  //double mean[4];
  double mean[3];
  h[0] = (TH1F*)f->Get(Form("%s_%d_%d_in36Arout35Ar",mainch,cr,p));
  h[1] = (TH1F*)f->Get(Form("%s_%d_%d_in36Arout32S",mainch,cr,p));
  h[2] = (TH1F*)f->Get(Form("%s_%d_%d_in36Arout28Si",mainch,cr,p));
  for(int i=0;i<3;i++){
    c->cd(1+i);
    //h[i]->Rebin(2);
    h[i]->Draw();
    h[i]->GetXaxis()->SetRangeUser(0,1000);


    fit[i] = new TF1(Form("fit_%d",i),fonegaus,0,1500,3);
    fit[i]->SetLineColor(3);
    fit[i]->SetLineWidth(1);
    
    fit[i]->SetParameter(0,h[i]->Integral());

    //good for crdc1 76 - 152 failures
    //fit[i]->SetParameter(1,130);
    //fit[i]->SetRange(100,200);

    fit[i]->SetParameter(1,100);
    fit[i]->SetParLimits(1,300,600);
    //fit[i]->SetRange(80,150);
    fit[i]->SetRange(150,700);
 

    // if(i==1)
    //   fit[i]->SetRange(70,140);
    //fit[i]->SetParLimits(1,center-25,center+25);
    
    fit[i]->SetParameter(2,30);
    fit[i]->SetParLimits(2,10,50);
    h[i]->Fit(fit[i],"Rq");
    mean[i] = fit[i]->GetParameter(1);
    
 }
  //Matching channel: 140.083 152.072 143.63 152.636
  //Matching channel: 130.752 142.59 136.227 144.149
  double match[3];
  //double match[4];
  if(cr==0){
    match[2] = 362.452;
    match[1] = 401.774;
    match[0] = 441.911;
    //match[0] = 140.083;
    //match[1] = 152.072;
    //match[2] = 143.63;
    //match[3] = 152.636;
  }
  if(cr==1){
    match[2] = 364.225;
    match[1] = 426.191;
    match[0] = 488.898;
    //match[0] = 130.752;
    //match[1] = 142.59;
    //match[2] = 136.227;
    //match[3] = 144.149;
  }
  g = new TGraph(3,mean,match);
  //g = new TGraph(4,mean,match);
  c->cd(5);
  g->Draw("AP*");
  TF1* line = new TF1("line","[0]+[1]*x");
  line->SetParameter(0,0.0);
  line->SetParameter(1,1.0);
  //line->SetParLimits(1,0.1,2.0);
  line->SetParLimits(1,0.1,3.0);
  g->Fit(line,"q");
  cout << line->GetParameter(0) << "\t" << line->GetParameter(1) << endl;
  TEnv* output = new TEnv("hist/padcal0.fixed.dat");
  output->SetValue(Form("Crdc.%d.Slope.%03d",cr,p),line->GetParameter(1));
  output->SetValue(Form("Crdc.%d.Offset.%03d",cr,p),line->GetParameter(0));
  output->SaveLevel(kEnvLocal);
}
Double_t fonegaus(Double_t *x, Double_t *par){
  static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
  Double_t arg;
  /*
  par[0]   gauss content
  par[1]   gauss mean
  par[2]   gauss width
  */
  Double_t result =0;
  arg = (x[0]-par[1])/(sqrt2*par[2]);
  result += 1/(sqrt2pi*par[2]) * par[0] * exp(-arg*arg);
  return result;
}
