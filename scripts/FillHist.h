#ifndef FillHist_h
#define FillHist_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2F.h>
#include <TCutG.h>
#include <TList.h>

#include "/npdisks/data4/rdol500/e14027/GrROOTold/src/Gretina.hh"
#include "/npdisks/data4/rdol500/e14027/GrROOTold/src/GretinaCalc.hh"
#include "/npdisks/data4/rdol500/e14027/GrROOTold/src/S800.hh"
#include "/npdisks/data4/rdol500/e14027/GrROOTold/src/S800Calc.hh"
#include "/npdisks/data4/rdol500/e14027/GrROOTold/src/RunInfo.hh"
#include "/npdisks/data4/rdol500/e14027/GrROOTold/src/Scaler.hh"
#include "/npdisks/data4/rdol500/e14027/GrROOTold/src/Settings.hh"

#include <vector>

using namespace std;

class FillHist
{

public:

			//--- constructor and destructors
			FillHist(char *inFile, char *inCuts) { mInput = new TFile(inFile,"read"); mCutFile = new TFile(inCuts,"read"); mSet(); }
			~FillHist() {}
			
			//--- public member functions
			void			mSet();
			void			mSetCuts();
			void			mSetHistos();
			
			void 			Run();
			
			void 			mTerminate();
			
private:		

			//--- private member variables

			TFile							*mInput;
			TFile							*mCutFile;
			TFile							*mOutput;
			TTree							*mTree;
			
			TList							*mCuts;				
			
   		//--- List of pointers to data branches
   		GretinaCalc					*mgretina;
	   	S800Calc						*ms800;
	   	
	   	//--- List of the histograms to be populated
	   	TH2F							*mxfp_vs_objm;
	   	TH2F							*mICde_vs_objmc_in79Sr;
	   	TH2F							*mICde_vs_x0_in79Sr;
	   	TH2F							*mICde_vs_x0_in79SroutC;
	   	
			//--- private member functions
			
			TCutG							*mFindCut(const char *cutName)	{ return	(TCutG*)mCuts->FindObject(cutName); }			   	

};

#endif // #ifdef FillHist_cxx
