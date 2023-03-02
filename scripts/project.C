#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1S.h"
#include "TH2S.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
using namespace TMath;
using namespace std;

void project(TH2F *hist){
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(2,2);
  c->cd(1);
  hist->Draw("colz");
  c->cd(2);
  hist->ProjectionY()->Draw();
  c->cd(3);
  hist->ProjectionX()->Draw();
}
void project(TH2S *hist){
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(2,2);
  c->cd(1);
  hist->Draw("colz");
  c->cd(2);
  hist->ProjectionY()->Draw();
  c->cd(3);
  hist->ProjectionX()->Draw();
}
void projectbins(TH2F *hist, int widthx, int widthy =0, int startx =0, int starty =0){

  if(widthy ==0)
    widthy = widthx;


  int xstep = hist->GetNbinsX()/widthx;
  if(xstep>10){
    cout << "more than 10 steps, not plotting more than 10 " << endl;
    xstep = 10;
  }
  TH1F* px[10];
  int ystep = hist->GetNbinsY()/widthy;
  if(ystep>10){
    cout << "more than 10 steps, not plotting more than 10 " << endl;
    ystep = 10;
  }
  TH1F* py[10];
  TCanvas *c = new TCanvas("c","c",0,0,800,800);
  c->Divide(2,2);
  c->cd(1);
  hist->Draw("colz");
  c->cd(2);  
  double max =0;
  for(int i=0;i<xstep;i++){
    //cout << "bin " << i*widthx<<"\t"<<(i+1)*widthx<<endl;
    px[i] = NULL;
    px[i] = (TH1F*)hist->ProjectionY(Form("px_%d",i),i*widthx+startx,(i+1)*widthx+startx);
    //    cout <<"\t"<< i*widthx<<"\t"<<(i+1)*widthx<<endl;
    px[i]->GetYaxis()->UnZoom();
    if(i==0)
      px[0]->Draw();
    px[i]->SetLineColor(1+i);
    if(px[i]->GetMaximum()>max)
      max = px[i]->GetMaximum();
    px[i]->Draw("same");
  }
  px[0]->GetYaxis()->SetRangeUser(0,max*1.2);
  c->cd(3);  
  max =0;
  for(int i=0;i<ystep;i++){
    py[i] = NULL;
    py[i] = (TH1F*)hist->ProjectionX(Form("py_%d",i),i*widthy+starty,(i+1)*widthy+starty);
    py[i]->GetYaxis()->UnZoom();
    if(i==0)
      py[0]->Draw();
    py[i]->SetLineColor(1+i);
    if(py[i]->GetMaximum()>max)
      max = py[i]->GetMaximum();
    py[i]->Draw("same");
  }
  py[0]->GetYaxis()->SetRangeUser(0,max*1.2);
}
