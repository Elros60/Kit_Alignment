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
#include <TChain.h>
#include <TGeoMatrix.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TString.h>

#include "CommonConstants/LHCConstants.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"

#include "SimulationDataFormat/MCCompLabel.h"

//  Test with alignment codes
#include "DataFormatsMCH/Cluster.h"
#include "MCHAlign/Alignment.h"
#include "MCHTracking/Track.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackFitter.h"

#include "DetectorsCommonDataFormats/AlignParam.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"





void ClCompare(int DetElem, std::string Axis, std::string file1, std::string file2, int nBin=200){

	TFile *f1 = TFile::Open(file1.c_str());
	TFile *f2 = TFile::Open(file2.c_str());

	TTree *t1 = (TTree*)f1->Get("TreeE");
	TTree *t2 = (TTree*)f2->Get("TreeE");


	t1->Draw(Form("%s>>hist1(%d,-1.0,1.0)",Axis.c_str(),nBin),Form("fClDetElem==%d",DetElem),"goff");
	t2->Draw(Form("%s>>hist2(%d,-1.0,1.0)",Axis.c_str(),nBin),Form("fClDetElem==%d",DetElem),"goff");

	TH1F *hist1 = (TH1F*)gDirectory->Get("hist1");
	TH1F *hist2 = (TH1F*)gDirectory->Get("hist2");

	hist1->Scale(1./hist1->Integral());
	hist2->Scale(1./hist2->Integral());

	//auto c = new TCanvas();
	hist1->SetLineColor(kRed);
	hist1->Draw("H");
	hist2->SetLineColor(kBlue);
	hist2->Draw("Hsame");

	auto legend = new TLegend();
	legend->AddEntry("hist1", file1.c_str());
	legend->AddEntry("hist2", file2.c_str());
	legend->Draw();




}