#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"

Double_t gBgConstant, gBgSlope, gContent, gMean, gContent_1, gMean_1, gContent_2, gMean_2, gSigma, gSigma_1, gSigma_2, gBinW, gChi2pNDF;

Double_t gaus_lbg(Double_t *x, Double_t *par);
void fit(){
  //TFile *f = new TFile("hodocal/hist/hist0039.root");
  TFile *f = new TFile("his/his0110.root");

  TCanvas *c = new TCanvas("c","c",0,0,1800,1000);
  c->Divide(8,4);
  TH1F *h[32];
  TF1 *fu[32][2];

  //ofstream fout("hodocal.dat");

  for(int i=0;i<32;i++){
    h[i] = (TH1F*)f->Get(Form("hodo_%d",i));
    c->cd(1+i);
    h[i]->Rebin(8);
    h[i]->Draw();
    h[i]->GetXaxis()->SetRangeUser(100,4000);
    //cout << i <<"\t"<<h[i]->Integral(50,2000) << endl;
    
    double gLowX =350;
    double gUpX =750;
    fu[i][0] = new TF1(Form("fu_0_%d",i),gaus_lbg,gLowX,gUpX,5);
    gContent = h[i]->Integral(h[i]->FindBin(gLowX),h[i]->FindBin(gUpX)); 
    gMean    = 0.5 * ( gLowX + gUpX);  
    gSigma   = 0.3 * ( gUpX  - gLowX); 

    //gSigma   = 0.3 * ( gUpX  - gLowX); 
    gBinW = h[i]->GetBinWidth(1);
    fu[i][0]->SetParameters(0, 0, gSigma, gContent, gMean); 
    //fu[i][0]->SetParLimits(2,1,100);
    fu[i][0]->SetRange(gLowX, gUpX);
    fu[i][0]->SetLineColor(3);
    fu[i][0]->SetLineWidth(1);
    h[i]->Fit(fu[i][0],"Rq");
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
    gContent = h[i]->Integral(h[i]->FindBin(gLowX),h[i]->FindBin(gUpX)); 
    gMean    = 0.5 * ( gLowX + gUpX);  
    gSigma   = 0.3 * ( gUpX  - gLowX); 
    gBinW = h[i]->GetBinWidth(1);
    fu[i][1]->SetParameters(0, 0, gSigma, gContent, gMean); 
    //cout<<fu[i][1]->GetParameter(3)<<" "<<fu[i][1]->GetParameter(4)<<" "<<fu[i][1]->GetParameter(2)<<endl;
    //fu[i][1]->SetParLimits(3,gLowX*1.05,gUpX*0.95);
    //fu[i][1]->SetParLimits(3,1,1e6);
    //fu[i][1]->SetParLimits(2,1,100);
    fu[i][1]->SetRange(gLowX, gUpX);
    fu[i][1]->SetLineColor(3);
    fu[i][1]->SetLineWidth(1);
    h[i]->Fit(fu[i][1],"Rq");
    fu[i][1]->Draw("same");
    fu[i][0]->Draw("same");
    //cout<<fu[i][1]->GetParameter(3)<<" "<<fu[i][1]->GetParameter(4)<<" "<<fu[i][1]->GetParameter(2)<<" "<<fu[i][1]->GetParameter(0)<<" "<<fu[i][1]->GetParameter(1)<<endl;

    double x0,x1,y0,y1;
    x0 = fu[i][0]->GetParameter(4);
    x1 = fu[i][1]->GetParameter(4);
    //cout<<x0<<" "<<x1<<endl;
    y0 = 511;
    y1 = 1275;
      
    double a = (y1-y0)/(x1-x0);
    double b = (y1+y0)-a*(x1+x0);
    b/=2.;

    //cout << a <<"\t" << b << endl;
    cout << Form("Hodoscope.Offset.%d\t",i) << b << endl;
    cout << Form("Hodoscope.Gain.%d\t",i) << a << endl;
    
  }
  c->SaveAs("fit.ps");
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
