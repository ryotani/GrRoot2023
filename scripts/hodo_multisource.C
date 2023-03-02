#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"

Double_t gBgConstant, gBgSlope, gContent, gMean, gContent_1, gMean_1, gContent_2, gMean_2, gSigma, gSigma_1, gSigma_2, gBinW, gChi2pNDF;

Double_t gaus_lbg(Double_t *x, Double_t *par);
void fit(){
  //TFile *f = new TFile("hodocal/hist/hist0039.root");
  TFile *f1 = new TFile("his/his0110.root");  //22Na
  TFile *f2 = new TFile("his/his0111.root");  //22Na

  TCanvas *c1 = new TCanvas("c1","c1",0,0,1800,1000);
  c1->cd();
  c1->Divide(8,4);
  TCanvas *c2 = new TCanvas("c2","c2",0,0,1800,1000);
  c2->cd();
  c2->Divide(8,4);
  TCanvas *c3 = new TCanvas("c3","c3",0,0,1800,1000);
  c3->cd();
  c3->Divide(8,4);
  TH1F *h[32][2];
  TF1 *fu[32][3];
  TGraph *g[32];

  //ofstream fout("hodocal.dat");

  for(int i=0;i<32;i++){
    f1->cd();
    h[i][0] = (TH1F*)f1->Get(Form("hodo_%d",i));
    c1->cd(1+i);
    h[i][0]->Rebin(8);
    h[i][0]->Draw();
    h[i][0]->GetXaxis()->SetRangeUser(100,4000);
    //cout << i <<"\t"<<h[i]->Integral(50,2000) << endl;
    
    double gLowX =350;
    double gUpX =750;
    fu[i][0] = new TF1(Form("fu_0_%d",i),gaus_lbg,gLowX,gUpX,5);
    gContent = h[i][0]->Integral(h[i][0]->FindBin(gLowX),h[i][0]->FindBin(gUpX)); 
    gMean    = 0.5 * ( gLowX + gUpX);  
    gSigma   = 0.3 * ( gUpX  - gLowX); 

    //gSigma   = 0.3 * ( gUpX  - gLowX); 
    gBinW = h[i][0]->GetBinWidth(1);
    fu[i][0]->SetParameters(0, 0, gSigma, gContent, gMean); 
    //fu[i][0]->SetParLimits(2,1,100);
    fu[i][0]->SetRange(gLowX, gUpX);
    fu[i][0]->SetLineColor(3);
    fu[i][0]->SetLineWidth(1);
    h[i][0]->Fit(fu[i][0],"Rq");
    gLowX =1100;
    gUpX =1600;

    //if(i==5){
    //  gLowX =3350;
    //  gUpX =4000;
    //}
    //if(i==25){
    //  gLowX =3200;
    //  gUpX =3800;
    //}
//     if(i==18||i==16){
//       gLowX =1300;
//       gUpX =1700;
//     }
     if(i==9||i==24||i==26){
       gLowX =1200;
       gUpX =1600;
     }

    fu[i][1] = new TF1(Form("fu_1_%d",i),gaus_lbg,gLowX,gUpX,5);
    gContent = h[i][0]->Integral(h[i][0]->FindBin(gLowX),h[i][0]->FindBin(gUpX)); 
    gMean    = 0.5 * ( gLowX + gUpX);  
    gSigma   = 0.3 * ( gUpX  - gLowX); 
    gBinW = h[i][0]->GetBinWidth(1);
    fu[i][1]->SetParameters(0, 0, gSigma, gContent, gMean); 
    //cout<<fu[i][1]->GetParameter(3)<<" "<<fu[i][1]->GetParameter(4)<<" "<<fu[i][1]->GetParameter(2)<<endl;
    //fu[i][1]->SetParLimits(3,gLowX*1.05,gUpX*0.95);
    //fu[i][1]->SetParLimits(3,1,1e6);
    //fu[i][1]->SetParLimits(2,1,100);
    fu[i][1]->SetRange(gLowX, gUpX);
    fu[i][1]->SetLineColor(3);
    fu[i][1]->SetLineWidth(1);
    h[i][0]->Fit(fu[i][1],"Rq");
    fu[i][1]->Draw("same");
    fu[i][0]->Draw("same");
    //cout<<fu[i][1]->GetParameter(3)<<" "<<fu[i][1]->GetParameter(4)<<" "<<fu[i][1]->GetParameter(2)<<" "<<fu[i][1]->GetParameter(0)<<" "<<fu[i][1]->GetParameter(1)<<endl;

    f2->cd();
    h[i][1] = (TH1F*)f2->Get(Form("hodo_%d",i));
    c2->cd(1+i);
    gLowX =30;
    gUpX =300;
    h[i][1]->Rebin(8);
    h[i][1]->Draw();
    h[i][1]->GetXaxis()->SetRangeUser(0,500);
    fu[i][2] = new TF1(Form("fu_2_%d",i),gaus_lbg,gLowX,gUpX,5);
    gContent = h[i][1]->Integral(h[i][1]->FindBin(gLowX),h[i][1]->FindBin(gUpX)); 
    gMean    = 0.5 * ( gLowX + gUpX);  
    gSigma   = 0.3 * ( gUpX  - gLowX); 
    gBinW = h[i][1]->GetBinWidth(1);
    fu[i][2]->SetParameters(0, 0, gSigma, gContent, gMean); 
    fu[i][2]->SetRange(gLowX, gUpX);
    fu[i][2]->SetLineColor(3);
    fu[i][2]->SetLineWidth(1);
    h[i][1]->Fit(fu[i][2],"Rq");
    fu[i][2]->Draw("same");

    double x[3],y[3];
    x[0] = fu[i][0]->GetParameter(4);
    x[1] = fu[i][1]->GetParameter(4);
    x[2] = fu[i][2]->GetParameter(4);
    //cout<<x0<<" "<<x1<<endl;
    y[0] = 511;
    y[1] = 1275;
    y[2] = 122;
 

    c3->cd(1+i);
    g[i] = new TGraph(3,x,y);
    TF1* line = new TF1("line","[0]+[1]*x");
    line->SetParameter(0,0.0);
    line->SetParameter(1,1.0);
    line->SetParLimits(1,0.1,2.0);
    g[i]->SetMarkerStyle(3);
    g[i]->SetMarkerSize(1);
    g[i]->Draw("AP");
    g[i]->Fit(line,"q");
    double a = line->GetParameter(1);
    double b = line->GetParameter(0);


    //double a = (y1-y0)/(x1-x0);
    //double b = (y1+y0)-a*(x1+x0);
    //b/=2.;

    //cout << a <<"\t" << b << endl;
    cout << Form("Hodoscope.Offset.%d\t",i) << b << endl;
    cout << Form("Hodoscope.Gain.%d\t",i) << a << endl;
    
  }
  c2->cd();
  c2->SaveAs("fit2.ps");
  c3->cd();
  c3->SaveAs("fit3.ps");
}
Double_t gaus_lbg(Double_t *x, Double_t *par)
{
/*
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss width
  par[3]   gauss0 constant
  par[4]   gauss0 mean
*/
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   Double_t arg;
   if (par[2] == 0) par[2]=1;                 //  force widths /= 0
   arg = (x[0] - par[4])/(sqrt2*par[2]);
   Double_t fitval = par[0] + x[0]*par[1]
              + gBinW/(sqrt2pi*par[2]) * par[3] * exp(-arg*arg);
   return fitval;
}
