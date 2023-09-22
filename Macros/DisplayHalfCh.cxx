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
#include <THStack.h>
#include <TColor.h>

//using namespace o2;
const Int_t NDetElemCh[10] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};

const Int_t colorsR[13] = {kRed+2,kRed+1,kOrange+10,kOrange+9,kOrange+8,kOrange+7,kOrange+6,kOrange+5,kOrange+4,kOrange+3,kOrange+2,kOrange+1,kOrange};
const Int_t colorsL[13] = {kAzure+2,kAzure+1,kViolet+10,kViolet+9,kViolet+8,kViolet+7,kViolet+6,kViolet+5,kViolet+4,kViolet+3,kViolet+2,kViolet+1,kViolet};

//________________________________________________________________________________
void DisplayHalfCh(string ResFileName, string Axe){
	
	TFile *f = TFile::Open(ResFileName.c_str());
	TTree *t = (TTree*)f->Get("TreeE");

	TFile *Res_File = TFile::Open(Form("ResidualPlot%s.root",Axe.c_str()),"RECREATE");

	for(int NumCh=1;NumCh<=10;NumCh++){
		int id_ch = NumCh-1;
		int NDetElem = NDetElemCh[id_ch];
		string cTitle = Form("Residual %s for Ch%d",Axe.c_str(),NumCh);
		auto c = new TCanvas(cTitle.c_str(),cTitle.c_str());
		gStyle->SetOptStat(0);
		auto legend = new TLegend();
		legend->SetHeader("Det Elem","C");

		TH1F *hist[NDetElem];
		THStack *hs = new THStack(Form("Residual%s_Ch%d",Axe.c_str(),NumCh), Form("Residual%s_Ch%d",Axe.c_str(),NumCh));

		for(int i = 0; i < NDetElem;i++){
			string DetTitle = Form("%s%d","DET",NumCh*100+i);
			hist[i] = new TH1F(DetTitle.c_str(),DetTitle.c_str(),150,-2.0,2.0);
			int DetElem = 100*NumCh+i;
			t->Draw(Form("fResidu%sGlobal>>%s",Axe.c_str(),hist[i]->GetName()),Form("fClDetElem==%d && fBendingMomentum < 80 && abs(fResidu%sGlobal) <2",DetElem,Axe.c_str()),"goff");
			
			hist[i]->Scale(1./hist[i]->Integral());

			if(i>NDetElem/4 && i<=NDetElem/4+NDetElem/2){
				hist[i]->SetLineColor(colorsL[i-(NDetElem/4+1)]);
			}else{
				if(i<=NDetElem/4){
					hist[i]->SetLineColor(colorsR[i]);	
				}else{
					hist[i]->SetLineColor(colorsR[i-NDetElem/2]);
				}
			}

			hist[i]->SetLineWidth(2);
			hs->Add(hist[i]);
			legend->AddEntry(hist[i],DetTitle.c_str(),"l");
		}
		hs->Draw("hist C nostack");
		legend->Draw();
		c->SetGrid(1,0);
		Res_File->WriteObjectAny(c,"TCanvas",hs->GetName());
	}

}








