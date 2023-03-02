//   This fitting program was written by Oleg Ivanov
//      How to use:
//   1. Load the program in root by typing .L FitPeaks.C - the code will be compiled and all the functions loaded to memory
//   2. Call GetFitPeaks() to see the list of fitting functions available
//   3. Call a function of choice with appropriate variables
//      In case of questions or comments please call me at 016 32 72 73

#include "TF1.h"
#include "TMath.h"
double pi = TMath::Pi();
double rad2deg = 180./pi;
double deg2rad = pi/180.;
 
Double_t gBgConstant, gBgSlope, gContent, gMean, gContent_1, gMean_1, gContent_2, gMean_2, gSigma, gSigma_1, gSigma_2, gBinW, gChi2pNDF;

//_____________________________________________________

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
//_____________________________________________________
Double_t gaus_nbg(Double_t *x, Double_t *par)
{
/*
  par[0]   gauss width
  par[1]   gauss0 constant
  par[2]   gauss0 mean
*/
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   Double_t arg;
   if (par[0] == 0) par[0]=1;                 //  force widths /= 0
   arg = (x[0] - par[2])/(sqrt2*par[0]);
   Double_t fitval =  gBinW/(sqrt2pi*par[0]) * par[1] * exp(-arg*arg);
   return fitval;
}
//_____________________________________________________

Double_t db_gaus_lbg(Double_t *x, Double_t *par)
{
/*
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss width
  par[3]   gauss0 constant
  par[4]   gauss0 mean
  par[5]   gauss1 constant
  par[6]   gauss1 mean
*/
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   Double_t arg_1, arg_2;
   if (par[2] == 0) par[2]=1;                 //  force widths /= 0
   arg_1 = (x[0] - par[4])/(sqrt2*par[2]);
   arg_2 = (x[0] - par[6])/(sqrt2*par[2]);
   Double_t fitval = par[0] + x[0]*par[1]
                   + gBinW/(sqrt2pi*par[2]) * par[3] * exp(-arg_1*arg_1)
                   + gBinW/(sqrt2pi*par[2]) * par[5] * exp(-arg_2*arg_2);
   return fitval;
}
//_____________________________________________________

Double_t db_gaus_lbg_diff(Double_t *x, Double_t *par)
{
/*
  par[0]   background constant
  par[1]   background slope
  par[2]   gauss0 width
  par[3]   gauss0 constant
  par[4]   gauss0 mean
  par[5]   gauss1 width
  par[6]   gauss1 constant
  par[7]   gauss1 mean
*/
   static Float_t sqrt2pi = TMath::Sqrt(2*TMath::Pi()), sqrt2 = TMath::Sqrt(2.);
   Double_t arg_1, arg_2;
   if (par[2] == 0) par[2]=1;                 //  force widths /= 0
   arg_1 = (x[0] - par[4])/(sqrt2*par[2]);
   arg_2 = (x[0] - par[7])/(sqrt2*par[5]);
   Double_t fitval = par[0] + x[0]*par[1]
                   + gBinW/(sqrt2pi*par[2]) * par[3] * exp(-arg_1*arg_1)
                   + gBinW/(sqrt2pi*par[5]) * par[6] * exp(-arg_2*arg_2);
   return fitval;
}
//_____________________________________________________

void GetFitHelp(void)
{
 printf("============================================================================\n");
 printf("   The following functions are available:\n");
 printf("----------------------------------------------------------------------------\n");
 printf("single peak without background\n");
 printf("          FindPeak(TH1F *HistogramName, Double_t LeftLimit, Double_t Energy)\n");
 printf("          FitPeak(TH1F *HistogramName, Double_t LeftLimit, Double_t RightLimit)\n");
 printf("----------------------------------------------------------------------------\n");
 printf("single peak with background\n");
 printf("          FindSinglePeak(TH1F *HistogramName, Double_t LeftLimit, Double_t Energy)\n");
 printf("          FitSinglePeak(TH1F *HistogramName, Double_t LeftLimit, Double_t RightLimit)\n");
 printf("----------------------------------------------------------------------------\n");

 printf("two peaks with background\n");

 printf("          FitDoublePeak(TH1F *HistogramName, Double_t LeftLimit_Peak_1, Double_t RightLimit_Peak_1,\n");
 printf("                                             Double_t LeftLimit_Peak_2, Double_t RightLimit_Peak_2)\n");
 printf("----------------------------------------------------------------------------\n");

 printf("two peaks with background and commond width\n");
 printf("   DeconvoluteDoublePeak(TH1F *HistogramName, Double_t LeftLimit_Peak_1, Double_t RightLimit_Peak_1,\n");
 printf("                                              Double_t LeftLimit_Peak_2, Double_t RightLimit_Peak_2)\n");
 printf("----------------------------------------------------------------------------\n");
 printf(" * Information\n");
 printf("----------------------------------------------------------------------------\n");
 printf(" *     TH1F : histograms with one float per channel. Maximum precision 7 digits.\n");
 printf(" * Double_t : the same as 'double' - Float 8 bytes\n");
 printf(" *  Float_t : the same as 'float'  - Float 4 bytes\n");
 printf("----------------------------------------------------------------------------\n");
 printf("   Have a Nice Work!\n");
 printf("============================================================================\n");
}

void FindSinglePeak(TH1F *hist, Double_t gEnergy)
{
 gBinW = hist->GetBinWidth(1);
 gLowX = gEnergy - gBinW*10.;
 gUpX  = gEnergy + gBinW*10.;
 hist->GetXaxis()->SetRange(gLowX,gUpX);
 FitSinglePeak(hist, gLowX, gUpX);
}
void FindPeak(TH1F *hist, Double_t gEnergy)
{
 gBinW = hist->GetBinWidth(1);
 gLowX = gEnergy - gBinW*10.;
 gUpX  = gEnergy + gBinW*10.;
 hist->GetXaxis()->SetRange(gLowX,gUpX);
 FitPeak(hist, gLowX, gUpX);
}

void FitSinglePeak(TH1F *hist, Double_t gLowX, Double_t gUpX)
{
 if(gLowX < gUpX)
  { 
// *** Creating the function of the form 'gaus_lbg' defined above ***
   TF1 fitfunc("gauss_linbg",gaus_lbg, 0, 1, 5);
// *** Obtaining and specifying the start values for the fit ***
   gContent = hist->Integral(hist->FindBin(gLowX),hist->FindBin(gUpX)); 
   gMean    = 0.5 * ( gLowX + gUpX);  
   gSigma   = 0.3 * ( gUpX  - gLowX); 
   gBinW = hist->GetBinWidth(1);
//   printf("__________________\n_The Start Values_\n  Bin Width: %d\n Mean Value: %d\n    Content: %d\n      Sigma: %d\n__________________\n",gBinW,gMean,gContent,gSigma);
   fitfunc.SetParameters(0, 0, gSigma, gContent, gMean); 
   fitfunc.SetRange(gLowX, gUpX);
   fitfunc.SetLineColor(3);
   fitfunc.SetLineWidth(1);
  
   fitfunc.SetParName(0,"BgConstant");
   fitfunc.SetParName(1,"BgSlope   ");
   fitfunc.SetParName(2,"Sigma     ");
   fitfunc.SetParName(3,"Content   ");
   fitfunc.SetParName(4,"Mean      ");

// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("gauss_linbg", "R", "SAME");

   gBgConstant = fitfunc.GetParameter(0);
   gBgSlope    = fitfunc.GetParameter(1);
   gSigma      = fitfunc.GetParameter(2);
   gContent    = fitfunc.GetParameter(3);
   gMean       = fitfunc.GetParameter(4);
   gChi2pNDF   = fitfunc.GetChisquare() / fitfunc.GetNDF();

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
   printf("            FWHM: %f +- %f\n",2*gSigma*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(2));
  } else cout << "Couldn't fit! Error: The Lower Limit is larger than the Upper Limit!" << endl;
}
//_____________________________________________________
void FitPeak(TH1F *hist, Double_t gLowX, Double_t gUpX)
{
 if(gLowX < gUpX)
  { 
// *** Creating the function of the form 'gaus_lbg' defined above ***
   TF1 fitfunc("gauss_linbg",gaus_nbg, 0, 1, 3);
// *** Obtaining and specifying the start values for the fit ***
   gContent = hist->Integral(hist->FindBin(gLowX),hist->FindBin(gUpX)); 
   gMean    = 0.5 * ( gLowX + gUpX);  
   gSigma   = 0.3 * ( gUpX  - gLowX); 
   gBinW = hist->GetBinWidth(1);
//   printf("__________________\n_The Start Values_\n  Bin Width: %d\n Mean Value: %d\n    Content: %d\n      Sigma: %d\n__________________\n",gBinW,gMean,gContent,gSigma);
   fitfunc.SetParameters(gSigma, gContent, gMean); 
   fitfunc.SetRange(gLowX, gUpX);
   fitfunc.SetLineColor(3);
   fitfunc.SetLineWidth(1);
  
   fitfunc.SetParName(0,"Sigma     ");
   fitfunc.SetParName(1,"Content   ");
   fitfunc.SetParName(2,"Mean      ");

// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("gauss_linbg", "R", "SAME");

   gSigma      = fitfunc.GetParameter(0);
   gContent    = fitfunc.GetParameter(1);
   gMean       = fitfunc.GetParameter(2);
   gChi2pNDF   = fitfunc.GetChisquare() / fitfunc.GetNDF();

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
   printf("            FWHM: %f +- %f\n",2*gSigma*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(2));
  } else cout << "Couldn't fit! Error: The Lower Limit is larger than the Upper Limit!" << endl;
}
//_____________________________________________________

void DeconvoluteDoublePeak(TH1F *hist, Double_t gLowX_1, Double_t gUpX_1, Double_t gLowX_2, Double_t gUpX_2)
{
 if(gLowX_1 < gUpX_1)
  {
 if(gLowX_2 < gUpX_2)
  {
// *** Creating the function of the form '2_gaus_lbg' defined above ***
   TF1 fitfunc("db_gauss_linbg",db_gaus_lbg, 0, 1, 7);
// *** Obtaining and specifying the start values for the fit ***
   gBinW      = hist->GetBinWidth(1);
   gContent_1 = gBinW*(hist->Integral(hist->FindBin(gLowX_1),hist->FindBin(gUpX_1)));
   gContent_2 = gBinW*(hist->Integral(hist->FindBin(gLowX_2),hist->FindBin(gUpX_2)));
// *** Searching for maximum Y value through the bins specified by limits
   int i, i_1, i_2;
   Double_t V, V_max;
   
   i_1 = int(hist->FindBin(gLowX_1));
   i_2 = int(hist->FindBin(gUpX_1));
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean_1 = double(i);
      }
    }
   gMean_1 = gBinW*gMean_1; 

   i_1 = int(hist->FindBin(gLowX_2));
   i_2 = int(hist->FindBin(gUpX_2));
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean_2 = double(i);
      }
    }
   gMean_2 = gBinW*gMean_2; 

   gSigma     = 0.5 * (0.3 * (gUpX_1 - gLowX_1) + 0.3 * (gUpX_2 - gLowX_2));
   printf("__________________\n_Peak 1: The Start Values_\n  Bin Width: %d\n Mean Value: %d\n    Content: %d\n      Sigma: %d\n__________________\n",gBinW,gMean_1,gContent_1,gSigma);
   printf("__________________\n_Peak 2: The Start Values_\n  Bin Width: %d\n Mean Value: %d\n    Content: %d\n      Sigma: %d\n__________________\n",gBinW,gMean_2,gContent_2,gSigma);
   fitfunc.SetParameters(0, 0, gSigma, gContent_1, gMean_1, gContent_2, gMean_2);
   fitfunc.SetRange(gLowX_1, gUpX_2);
   fitfunc.SetLineColor(3);
   fitfunc.SetLineWidth(1);
  
   fitfunc.SetParName(0,"BgConstant");
   fitfunc.SetParName(1,"BgSlope   ");
   fitfunc.SetParName(2,"Sigma     ");
   fitfunc.SetParName(3,"Content 1 ");
   fitfunc.SetParName(4,"Mean 1    ");
   fitfunc.SetParName(5,"Content 2 ");
   fitfunc.SetParName(6,"Mean 2    ");

// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("db_gauss_linbg", "R", "SAME");

   gBgConstant = fitfunc.GetParameter(0);
   gBgSlope    = fitfunc.GetParameter(1);
   gSigma      = fitfunc.GetParameter(2);
   gContent_1  = fitfunc.GetParameter(3);
   gMean_1     = fitfunc.GetParameter(4);
   gContent_2  = fitfunc.GetParameter(5);
   gMean_2     = fitfunc.GetParameter(6);
   gChi2pNDF   = fitfunc.GetChisquare() / fitfunc.GetNDF();

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
   printf("            FWHM: %f +- %f\n",2*gSigma*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(2));
  } else cout << "Couldn't fit! Error: Peak 2: The Lower Limit is larger than the Upper Limit!" << endl;
  } else cout << "Couldn't fit! Error: Peak 1: The Lower Limit is larger than the Upper Limit!" << endl;
}
//_____________________________________________________

void FitDoublePeak(TH1F *hist, Double_t gLowX_1, Double_t gUpX_1, Double_t gLowX_2, Double_t gUpX_2)
{
 if(gLowX_1 < gUpX_1)
  {
 if(gLowX_2 < gUpX_2)
  {
// *** Creating the function of the form '2_gaus_lbg' defined above ***
   TF1 fitfunc("db_gauss_linbg_diff",db_gaus_lbg_diff, 0, 1, 8);
// *** Obtaining and specifying the start values for the fit ***
   gBinW      = hist->GetBinWidth(1);
   gContent_1 = gBinW*(hist->Integral(hist->FindBin(gLowX_1),hist->FindBin(gUpX_1)));
   gContent_2 = gBinW*(hist->Integral(hist->FindBin(gLowX_2),hist->FindBin(gUpX_2)));
// *** Searching for maximum Y value through the bins specified by limits
   int i, i_1, i_2;
   Double_t V, V_max;
   
   i_1 = int(hist->FindBin(gLowX_1));
   i_2 = int(hist->FindBin(gUpX_1));
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean_1 = double(i);
      }
    }
   gMean_1 = gBinW*gMean_1; 

   i_1 = int(hist->FindBin(gLowX_2));
   i_2 = int(hist->FindBin(gUpX_2));
   V_max = hist->GetBinContent(i_1);
   for(i = i_1; i <= i_2; i++)
    {
     V = hist->GetBinContent(i);
     if(V > V_max)
      {
       V_max = V;
       gMean_2 = double(i);
      }
    }
   gMean_2 = gBinW*gMean_2; 

   gSigma_1 = 0.3 * (gUpX_1 - gLowX_1);
   gSigma_2 = 0.3 * (gUpX_2 - gLowX_2);
   
   printf("__________________\n_Peak 1: The Start Values_\n  Bin Width: %d\n Mean Value: %d\n    Content: %d\n      Sigma: %d\n__________________\n",gBinW,gMean_1,gContent_1,gSigma);
   printf("__________________\n_Peak 2: The Start Values_\n  Bin Width: %d\n Mean Value: %d\n    Content: %d\n      Sigma: %d\n__________________\n",gBinW,gMean_2,gContent_2,gSigma);
   fitfunc.SetParameters(0, 0, gSigma_1, gContent_1, gMean_1, gSigma_2, gContent_2, gMean_2);
   fitfunc.SetRange(gLowX_1, gUpX_2);
   fitfunc.SetLineColor(3);
   fitfunc.SetLineWidth(1);
  
   fitfunc.SetParName(0,"BgConstant");
   fitfunc.SetParName(1,"BgSlope   ");
   fitfunc.SetParName(2,"Sigma 1   ");
   fitfunc.SetParName(3,"Content 1 ");
   fitfunc.SetParName(4,"Mean 1    ");
   fitfunc.SetParName(5,"Sigma 2   ");
   fitfunc.SetParName(6,"Content 2 ");
   fitfunc.SetParName(7,"Mean 2    ");

// *** Fitting: 'R' means within the range specified above ***
   hist->Fit("db_gauss_linbg_diff", "R", "SAME");

   gBgConstant = fitfunc.GetParameter(0);
   gBgSlope    = fitfunc.GetParameter(1);
   gSigma_1    = fitfunc.GetParameter(2);
   gContent_1  = fitfunc.GetParameter(3);
   gMean_1     = fitfunc.GetParameter(4);
   gSigma_2    = fitfunc.GetParameter(5);
   gContent_2  = fitfunc.GetParameter(6);
   gMean_2     = fitfunc.GetParameter(7);
   gChi2pNDF   = fitfunc.GetChisquare() / fitfunc.GetNDF();

   printf("      Chi Square: %f\n",fitfunc.GetChisquare());
   printf("            FWHM 1: %f +- %f\n",2*gSigma_1*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(2));
   printf("            FWHM 2: %f +- %f\n",2*gSigma_2*sqrt(2*log(2)),2*sqrt(2*log(2))*fitfunc.GetParError(5));
  } else cout << "Couldn't fit! Error: Peak 2: The Lower Limit is larger than the Upper Limit!" << endl;
  } else cout << "Couldn't fit! Error: Peak 1: The Lower Limit is larger than the Upper Limit!" << endl;
}
