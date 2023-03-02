#define FillHist_cxx

#include "FillHist.h"

//--- Set the tree, branches and call the initializing functions
void FillHist::mSet()
{
	mTree = (TTree*)mInput->Get("ctr");
	
	mgretina = new GretinaCalc;
	ms800 = new S800Calc;	
		
	mTree->SetBranchAddress("gretinacalc", &mgretina);
	mTree->SetBranchAddress("s800calc", &ms800);
	
	mOutput = new TFile("BeamCuts.root","recreate");
	
	mSetCuts();
	mSetHistos();
}

//--- Set the cuts, at the moment the cuts must be called inBeam and inBeamoutA
void FillHist::mSetCuts()
{
	TCutG *cut;
	mCuts = new TList();
	
	cut = (TCutG*)mCutFile->FindObjectAny("in79Sr");
	mCuts->Add( cut );
	
	cut = (TCutG*)mCutFile->FindObjectAny("in79SroutC");
	mCuts->Add( cut );
}

//--- Set the histograms to be filled
void FillHist::mSetHistos()
{
	mxfp_vs_objm = new TH2F("xfp_vs_objm", "A1900 scintillator TOF (xfp) Vs S800 analysis line scintillator TOF (obj)", 800, -2000, -1200, 1200, -1800, -600 );
	mxfp_vs_objm->GetXaxis()->SetTitle("objm (time)");
	mxfp_vs_objm->GetYaxis()->SetTitle("xfpm (time)"); 

	mICde_vs_objmc_in79Sr = new TH2F("ICde_vs_objmc_in79Sr","beam Chamber", 450,-1250,-800,1000, 1000, 3000);
	mICde_vs_objmc_in79Sr->GetXaxis()->SetTitle("objmc (au(time))");
	mICde_vs_objmc_in79Sr->GetYaxis()->SetTitle("ICde (au (energy))");	
	
	mICde_vs_x0_in79Sr = new TH2F("ICde_vs_x0_in79Sr","beam Chamber", 600,-300,300,1500, 500, 3500);
	mICde_vs_x0_in79Sr->GetXaxis()->SetTitle("x0 [distance]");
	mICde_vs_x0_in79Sr->GetYaxis()->SetTitle("ICde [au (energy)]");	
	
	mICde_vs_x0_in79SroutC = new TH2F("ICde_vs_x0_in79SroutC","beam Chamber", 600,-300,300,1500, 500, 3500);
	mICde_vs_x0_in79SroutC->GetXaxis()->SetTitle("x0 [distance]");
	mICde_vs_x0_in79SroutC->GetYaxis()->SetTitle("ICde [au (energy)]");	
}

void FillHist::Run()
{
	for( int i = 0; i < mTree->GetEntries(); i++ )
	{
		mTree->GetEntry(i);
	
		mxfp_vs_objm->Fill( ms800->GetTOF()->GetMOBJ(), ms800->GetTOF()->GetMXFP() );
		
		if( (this->mFindCut("in79Sr"))->IsInside( ms800->GetTOF()->GetMOBJ(), ms800->GetTOF()->GetMXFP() ) )
		{
			mICde_vs_objmc_in79Sr->Fill( ms800->GetTOF()->GetMOBJC(), ms800->GetIC()->GetDE() );			
			
			mICde_vs_x0_in79Sr->Fill( ms800->GetPAD(0)->GetX(), ms800->GetIC()->GetDE() );			
		}
		
		if( (this->mFindCut("in79SroutC"))->IsInside( ms800->GetTOF()->GetMOBJC(), ms800->GetIC()->GetDE() ) )
		{
			mICde_vs_x0_in79SroutC->Fill( ms800->GetPAD(0)->GetX(), ms800->GetIC()->GetDE() );
		}		
	}
	
	//--- call the termination function
	mTerminate();
}

void FillHist::mTerminate()
{
	mOutput->Write();
	mOutput->Close();
}
