#include <cmath>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <ctime>
#include <iostream>

#include <gsl/span>

// #include <TSystem.h>
#include <Math/Vector4D.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <TStyle.h>

//using namespace o2;


//_________________________________________________________________________________________________
void ResPlot(std::string RefResFileName, std::string NewResFileName, std::string opt_direction = "x", int NumChamber = 0, double limit_abs = 5.0, int Nbin = 200){



	TFile *RefFile = TFile::Open(RefResFileName.c_str());
	TTree *RefTree = (TTree*)RefFile->Get("TreeE");

	TFile *NewFile = TFile::Open(NewResFileName.c_str());
	TTree *NewTree = (TTree*)NewFile->Get("TreeE");


	int RefClDetElem;
  	int RefClDetElemNumber;
  	double RefClusterX;
  	double RefClusterY;
  	double RefTrackX;
  	double RefTrackY;
  	double RefTrackSlopeX;
  	double RefTrackSlopeY;
  	RefTree->SetBranchAddress("fClusterX",&RefClusterX);
   	RefTree->SetBranchAddress("fClusterY",&RefClusterY);
    RefTree->SetBranchAddress("fTrackX",&RefTrackX);
    RefTree->SetBranchAddress("fTrackY",&RefTrackY);
    RefTree->SetBranchAddress("fTrackSlopeX",&RefTrackSlopeX);
    RefTree->SetBranchAddress("fTrackSlopeY",&RefTrackSlopeY);
    RefTree->SetBranchAddress("fClDetElem",&RefClDetElem);
    RefTree->SetBranchAddress("fClDetElemNumber",&RefClDetElemNumber);


	int NewClDetElem;
  	int NewClDetElemNumber;
  	double NewClusterX;
  	double NewClusterY;
  	double NewTrackX;
  	double NewTrackY;
  	double NewTrackSlopeX;
  	double NewTrackSlopeY;
  	NewTree->SetBranchAddress("fClusterX",&NewClusterX);
   	NewTree->SetBranchAddress("fClusterY",&NewClusterY);
    NewTree->SetBranchAddress("fTrackX",&NewTrackX);
    NewTree->SetBranchAddress("fTrackY",&NewTrackY);
    NewTree->SetBranchAddress("fTrackSlopeX",&NewTrackSlopeX);
    NewTree->SetBranchAddress("fTrackSlopeY",&NewTrackSlopeY);
    NewTree->SetBranchAddress("fClDetElem",&NewClDetElem);
    NewTree->SetBranchAddress("fClDetElemNumber",&NewClDetElemNumber);



    int Ref_NbEntries = RefTree->GetEntries();
    int New_NbEntries = NewTree->GetEntries();

    if(Ref_NbEntries != New_NbEntries){
    	cout << "Trees being processed are not compatible in their size!"<<endl;
    }

    TH1F *Hist_Ref = new TH1F("Hist_Ref","Residual", Nbin, (-1.0)*limit_abs, limit_abs);
    TH1F *Hist_New = new TH1F("Hist_New","Residual", Nbin, (-1.0)*limit_abs, limit_abs);

    TCanvas *c = new TCanvas("c", "Ref vs New");
    TCanvas *c1 = new TCanvas("c1", "Ref");
    TCanvas *c2 = new TCanvas("c2", "New");

    if(opt_direction=="x"){

        Hist_Ref->SetTitle(RefResFileName.c_str());
        Hist_New->SetTitle(NewResFileName.c_str());

        Hist_Ref->SetXTitle("ClusterX-TrackX [cm]");
        Hist_New->SetXTitle("ClusterX-TrackX [cm]");


        Hist_Ref->SetYTitle("Counts");
        Hist_New->SetYTitle("Counts");


        for(int i=0; i < Ref_NbEntries;i++){
            RefTree->GetEntry(i);
            double Res = RefClusterX - RefTrackX;
            if(NumChamber==0){
                Hist_Ref->Fill(Res);
            }else{
                if(NumChamber==int(RefClDetElem/100)) Hist_Ref->Fill(Res);
            }
        }

        for(int i=0; i < New_NbEntries;i++){
            NewTree->GetEntry(i);
            double Res = NewClusterX - NewTrackX;
            if(NumChamber==0){
                Hist_New->Fill(Res);
            }else{
                if(NumChamber==int(NewClDetElem/100)) Hist_New->Fill(Res);
            }
        }


    }else{

        Hist_Ref->SetTitle(RefResFileName.c_str());
        Hist_New->SetTitle(NewResFileName.c_str());

        Hist_Ref->SetXTitle("ClusterY-TrackY [cm]");
        Hist_New->SetXTitle("ClusterY-TrackY [cm]");

        Hist_Ref->SetYTitle("Counts");
        Hist_New->SetYTitle("Counts");


        for(int i=0; i < Ref_NbEntries;i++){
            RefTree->GetEntry(i);
            double Res = RefClusterY - RefTrackY;
            if(NumChamber==0){
                Hist_Ref->Fill(Res);
            }else{
                if(NumChamber==int(RefClDetElem/100)) Hist_Ref->Fill(Res);
            }
        }

        for(int i=0; i < New_NbEntries;i++){
            NewTree->GetEntry(i);
            double Res = NewClusterY - NewTrackY;
            if(NumChamber==0){
                Hist_New->Fill(Res);
            }else{
                if(NumChamber==int(NewClDetElem/100)) Hist_New->Fill(Res);
            }
        }        
    }

    //gStyle->SetOptStat(0);
    c1->cd();
    Hist_Ref->Draw();
    c2->cd();
    Hist_New->Draw();
    c->cd();
    Hist_New->SetLineColor(kRed);
    Hist_New->Draw();
    Hist_Ref->SetLineColor(kBlue);
    Hist_Ref->Draw("same");

    auto legend = new TLegend();
    legend->AddEntry("Hist_Ref",RefResFileName.c_str());
    legend->AddEntry("Hist_New",NewResFileName.c_str());
    legend->Draw();

}