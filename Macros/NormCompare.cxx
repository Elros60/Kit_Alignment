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

//  Test with alignment codes
#include "MCHAlign/Alignment.h"

using namespace o2;

const int fgNCh = 10;
const int fgNDetElemCh[fgNCh] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
const int fgSNDetElemCh[fgNCh + 1] = {0,  4,  8,   12,  16, 34,
                                      52, 78, 104, 130, 156};

const int fgNDetElemHalfCh[fgNCh*2] = {2, 2, 2, 2, 2, 2, 2, 2, 9,
									 9, 9, 9, 13, 13, 13, 13, 13, 13, 13, 13};

const int fgSNDetElemHalfCh[fgNCh*2+1] = {0, 3, 6, 9, 12, 15, 18, 21, 24, 34, 44, 54, 64, 78, 92, 106, 120, 134, 148, 162, 176};
	
const int fgNDetHalfChMax = 13;
const int fgDetElemHalfCh[fgNCh*2][fgNDetHalfChMax] =
  {
    {100, 103, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {101, 102, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

    {200, 203, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {201, 202, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

    {300, 303, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {301, 302, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

    {400, 403, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {401, 402, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},

    {500, 501, 502, 503, 504, 514, 515, 516, 517, 0, 0, 0, 0},
    {505, 506, 507, 508, 509, 510, 511, 512, 513, 0, 0, 0, 0},

    {600, 601, 602, 603, 604, 614, 615, 616, 617, 0, 0, 0, 0},
    {605, 606, 607, 608, 609, 610, 611, 612, 613, 0, 0, 0, 0},

    {700, 701, 702, 703, 704, 705, 706, 720, 721, 722, 723, 724, 725},
    {707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719},

    {800, 801, 802, 803, 804, 805, 806, 820, 821, 822, 823, 824, 825},
    {807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819},

    {900, 901, 902, 903, 904, 905, 906, 920, 921, 922, 923, 924, 925},
    {907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919},

    {1000, 1001, 1002, 1003, 1004, 1005, 1006, 1020, 1021, 1022, 1023, 1024, 1025},
    {1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019}

};

o2::mch::geo::TransformationCreator transformation;
std::map<int, o2::math_utils::Transform3D> transformRef;
std::map<int, o2::math_utils::Transform3D> transformRefFit;
std::map<int, o2::math_utils::Transform3D> transformTracks;

Int_t GetDetElemNumber(Int_t iDetElemId);
Int_t GetDetElemId(Int_t iDetElemNumber);
constexpr double pi() { return 3.14159265358979323846; }


//___________________________________________________________________________
void NormCompare(){

	auto &segmentation = mch::mapping::segmentation(300);
	if (segmentation.nofPads() != 27873) {
  	LOG(error) << "wrong mapping implementation";
  	LOG(error) << "do gSystem->Load(\"libO2MCHMappingImpl4\") before compiling "
                "this macro";
  	exit(-1);
	}

  const auto grp = parameters::GRPObject::loadFrom(base::NameConf::getGRPFileName());
  base::Propagator::initFieldFromGRP(grp);


  const int NbIter = 7;
  double NumIter[NbIter] = {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0};
  int nDetElem = 156;
  double normDiffs[6][NbIter] = {};

  for(int it = 0 ; it<NbIter;it++){

    std::string GeoFile_Ref = Form("%s%g%s","o2sim_geometry_RefAligned_It", NumIter[it], ".root");
    std::string GeoFile_RefFit = Form("%s%g%s","o2sim_geometry_RefFitAligned_It", NumIter[it], ".root");
    std::string GeoFile_Tracks = Form("%s%g%s","o2sim_geometry_TracksAligned_It", NumIter[it], ".root");

    std::string ParamFile_Ref = Form("%s%g%s","Ch7Fixed_AlignParam_Ref_It", NumIter[it], ".root");
    std::string ParamFile_RefFit = Form("%s%g%s","Ch7Fixed_AlignParam_RefFit_It", NumIter[it], ".root");
    std::string ParamFile_TrackRef = Form("%s%g%s","Ch7Fixed_AlignParam_TracksRef_It", NumIter[it], ".root");


    TFile *Params_Ref = TFile::Open(ParamFile_Ref.c_str());
    std::vector<o2::detectors::AlignParam> ParamsTrack_Ref = *(Params_Ref->Get<std::vector<o2::detectors::AlignParam>>("alignment"));

    TFile *Params_RefFit = TFile::Open(ParamFile_RefFit.c_str());
    std::vector<o2::detectors::AlignParam> ParamsTrack_RefFit = *(Params_RefFit->Get<std::vector<o2::detectors::AlignParam>>("alignment"));

    TFile *Params_TrackRef = TFile::Open(ParamFile_TrackRef.c_str());
    std::vector<o2::detectors::AlignParam> ParamsTrack_TrackRef = *(Params_TrackRef->Get<std::vector<o2::detectors::AlignParam>>("alignment"));


    base::GeometryManager::loadGeometry(GeoFile_Ref.c_str());
    transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);
    for (int i = 0; i < 156; i++) { 
      int iDEN = GetDetElemId(i);
      transformRef[iDEN] = transformation(iDEN);
    }

    base::GeometryManager::loadGeometry(GeoFile_RefFit.c_str());
    transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);
    for (int i = 0; i < 156; i++) { 
      int iDEN = GetDetElemId(i);
      transformRefFit[iDEN] = transformation(iDEN);
    }

    base::GeometryManager::loadGeometry(GeoFile_Tracks.c_str());
    transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);
    for (int i = 0; i < 156; i++) { 
      int iDEN = GetDetElemId(i);
      transformTracks[iDEN] = transformation(iDEN);
    }


    double E2NormSum = 0.0;
    double E2NormSumFit = 0.0;
    double E2NormSumRefFit = 0.0;

    for(int i = 0; i<156;i++){

      int iDEN = GetDetElemId(i);
      auto transform3DRef = transformRef[iDEN];
      auto transform3DRefFit = transformRefFit[iDEN];
      auto transform3DTracks = transformTracks[iDEN];

      TMatrixD MTransRef(3,4);
      TMatrixD MTransRefFit(3,4);
      TMatrixD MTransTracks(3,4);

      transform3DRef.GetTransformMatrix(MTransRef);
      transform3DRefFit.GetTransformMatrix(MTransRefFit);
      transform3DTracks.GetTransformMatrix(MTransTracks);

      //TMatrixD MTrans = MTransRef - MTransTracks;
      //TMatrixD MTrans_Fit = MTransRefFit - MTransTracks;
      //TMatrixD MTrans_RefFit = MTransRef - MTransRefFit;

      //E2NormSum += MTrans.E2Norm();
      //E2NormSumFit += MTrans_Fit.E2Norm();
      //E2NormSumRefFit += MTrans_RefFit.E2Norm();

      E2NormSum += MTransRef.E2Norm()/MTransTracks.E2Norm();
      E2NormSumFit += MTransRefFit.E2Norm()/MTransTracks.E2Norm();
      E2NormSumRefFit += MTransRef.E2Norm()/MTransRefFit.E2Norm();
    }

    E2NormSum = E2NormSum/nDetElem;
    E2NormSumFit = E2NormSumFit/nDetElem;
    E2NormSumRefFit = E2NormSumRefFit/nDetElem;

    normDiffs[0][it] = E2NormSum;
    normDiffs[1][it] = E2NormSumFit;
    normDiffs[2][it] = E2NormSumRefFit;




    
    double E2NormSumParam_RefTrack = 0.0;
    double E2NormSumParam_RefFitTrack = 0.0;
    double E2NormSumParam_RefFit = 0.0;

    for(int i=0; i<ParamsTrack_Ref.size();i++){


      o2::detectors::AlignParam param_Ref = ParamsTrack_Ref.at(i);
      o2::detectors::AlignParam param_RefFit = ParamsTrack_RefFit.at(i);
      o2::detectors::AlignParam param_TrackRef = ParamsTrack_TrackRef.at(i);
      double x,y,z,psi,theta,phi;


      x = param_Ref.getX();
      y = param_Ref.getY();
      z = param_Ref.getZ();
      psi = param_Ref.getPsi();
      theta = param_Ref.getTheta();
      phi = param_Ref.getPhi();

      TGeoHMatrix matrix_temp1;
      matrix_temp1.SetDx(x);
      matrix_temp1.SetDy(y);
      matrix_temp1.SetDz(z);
      matrix_temp1.RotateX(psi*180/pi());
      matrix_temp1.RotateY(theta*180/pi());
      matrix_temp1.RotateZ(phi*180/pi());

      o2::math_utils::Transform3D mat_Ref(matrix_temp1);
      TMatrixD trans_Ref(3,4);
      mat_Ref.GetTransformMatrix(trans_Ref);



      x = param_RefFit.getX();
      y = param_RefFit.getY();
      z = param_RefFit.getZ();
      psi = param_RefFit.getPsi();
      theta = param_RefFit.getTheta();
      phi = param_RefFit.getPhi();

      TGeoHMatrix matrix_temp2;
      matrix_temp2.SetDx(x);
      matrix_temp2.SetDy(y);
      matrix_temp2.SetDz(z);
      matrix_temp2.RotateX(psi*180/pi());
      matrix_temp2.RotateY(theta*180/pi());
      matrix_temp2.RotateZ(phi*180/pi());

      o2::math_utils::Transform3D mat_RefFit(matrix_temp2);
      TMatrixD trans_RefFit(3,4);
      mat_RefFit.GetTransformMatrix(trans_RefFit);



      x = param_TrackRef.getX();
      y = param_TrackRef.getY();
      z = param_TrackRef.getZ();
      psi = param_TrackRef.getPsi();
      theta = param_TrackRef.getTheta();
      phi = param_TrackRef.getPhi();

      TGeoHMatrix matrix_temp3;
      matrix_temp3.SetDx(x);
      matrix_temp3.SetDy(y);
      matrix_temp3.SetDz(z);
      matrix_temp3.RotateX(psi*180/pi());
      matrix_temp3.RotateY(theta*180/pi());
      matrix_temp3.RotateZ(phi*180/pi());

      o2::math_utils::Transform3D mat_TrackRef(matrix_temp3);
      TMatrixD trans_TrackRef(3,4);
      mat_TrackRef.GetTransformMatrix(trans_TrackRef);


      //TMatrixD diff_RefTracks = trans_Ref - trans_TrackRef;
      //TMatrixD diff_RefFitTracks = trans_RefFit - trans_TrackRef;
      //TMatrixD diff_RefFit = trans_Ref - trans_RefFit;

      //E2NormSumParam_RefTrack += diff_RefTracks.E2Norm();
      //E2NormSumParam_RefFitTrack += diff_RefFitTracks.E2Norm();
      //E2NormSumParam_RefFit += diff_RefFit.E2Norm();

      E2NormSumParam_RefTrack += trans_Ref.E2Norm()/trans_TrackRef.E2Norm();
      E2NormSumParam_RefFitTrack += trans_RefFit.E2Norm()/trans_TrackRef.E2Norm();
      E2NormSumParam_RefFit += trans_Ref.E2Norm()/trans_RefFit.E2Norm();


    }

    E2NormSumParam_RefTrack = E2NormSumParam_RefTrack/nDetElem;
    E2NormSumParam_RefFitTrack = E2NormSumParam_RefFitTrack/nDetElem;
    E2NormSumParam_RefFit = E2NormSumParam_RefFit/nDetElem;

    normDiffs[3][it] = E2NormSumParam_RefTrack;
    normDiffs[4][it] = E2NormSumParam_RefFitTrack;
    normDiffs[5][it] = E2NormSumParam_RefFit;
    
    // LOG(info) << Form("%s%g%s%f","Iterations: ", NumIter[it], " Norm difference: ", E2NormSum);

  }

  auto c1 = new TCanvas("c1", "Norm Compare Geometry");
  gStyle->SetOptStat(0);
  c1->cd();
  c1->SetGrid();
  auto mg1 = new TMultiGraph("mg1", "Norm compare | geometry;Number of iterations;E2 norm");

  TGraph *graphNorm = new TGraph(NbIter, NumIter, normDiffs[0]);
  graphNorm->SetName("gr1");
  //graphNorm->SetLineColor(kBlue);
  graphNorm->SetMarkerSize(1);
  graphNorm->SetLineWidth(1.5);
  graphNorm->SetMarkerColor(kBlue);
  mg1->Add(graphNorm, "C*");
  
  TGraph *graphNormFit = new TGraph(NbIter, NumIter, normDiffs[1]);
  graphNormFit->SetName("gr2");
  //graphNormFit->SetLineColor(kRed);
  graphNormFit->SetMarkerStyle(21);
  graphNormFit->SetMarkerSize(1);
  graphNormFit->SetLineWidth(1.5);
  graphNormFit->SetMarkerColor(kRed);
  mg1->Add(graphNormFit,"CP");  


  mg1->Draw("A");

  auto legend1 = new TLegend();
  legend1->AddEntry("gr1","Norm(Ref - Tracks)");
  legend1->AddEntry("gr2","Norm(RefFit - Tracks)");
  legend1->Draw();


  auto c2 = new TCanvas("c2", "Compare Ref vs Ref(Fit) | geometry");
  c2->cd();
  c2->SetGrid();
  TGraph *graphNormRefFit = new TGraph(NbIter, NumIter, normDiffs[2]);
  graphNormRefFit->SetName("gr3");
  graphNormRefFit->SetTitle("Ref vs Ref(Fit) | geometry;Number of iteration; E2 norm");
  //graphNormFit->SetLineColor(kRed);
  graphNormRefFit->SetMarkerSize(2);
  graphNormRefFit->SetLineWidth(1.5);
  graphNormRefFit->SetMarkerColor(kBlue);
  graphNormRefFit->Draw("APC*");



  auto c3 = new TCanvas("c3", "Norm Compare Params");
  gStyle->SetOptStat(0);
  c3->cd();
  c3->SetGrid();
  auto mg2 = new TMultiGraph("mg2", "Norm compare | params;Number of iterations;E2 norm");

  TGraph *graphNormParam = new TGraph(NbIter, NumIter, normDiffs[3]);
  graphNormParam->SetName("gr4");
  //graphNormFit->SetLineColor(kRed);
  graphNormParam->SetMarkerSize(1);
  graphNormParam->SetLineWidth(1.5);
  graphNormParam->SetMarkerColor(kBlue);
  mg2->Add(graphNormParam,"C*");

  TGraph *graphNormParamFit = new TGraph(NbIter, NumIter, normDiffs[4]);
  graphNormParamFit->SetName("gr5");
  //graphNormFit->SetLineColor(kRed);
  graphNormParamFit->SetMarkerStyle(21);
  graphNormParamFit->SetMarkerSize(1);
  graphNormParamFit->SetLineWidth(1.5);
  graphNormParamFit->SetMarkerColor(kRed);
  mg2->Add(graphNormParamFit, "CP");

  mg2->Draw("A");

  auto legend2 = new TLegend();
  legend2->AddEntry("gr4","NormParam(Ref - Tracks)");
  legend2->AddEntry("gr5","NormParam(RefFit - Tracks)");
  legend2->Draw();


  auto c4 = new TCanvas("c4", "Compare Ref vs Ref(Fit) | params");
  c4->cd();
  c4->SetGrid();
  TGraph *graphNormParamRefFit = new TGraph(NbIter, NumIter, normDiffs[5]);
  graphNormParamRefFit->SetName("gr6");
  graphNormParamRefFit->SetTitle("Ref vs Ref(Fit) | params;Number of iteration; E2 norm");
  //graphNormFit->SetLineColor(kRed);
  graphNormParamRefFit->SetMarkerSize(2);
  graphNormParamRefFit->SetLineWidth(1.5);
  graphNormParamRefFit->SetMarkerColor(kBlue);
  graphNormParamRefFit->Draw("APC*");

}


//_________________________________________________________________________________________________
Int_t GetDetElemNumber(Int_t iDetElemId) {
  /// get det element number from ID
  // get chamber and element number in chamber
  const Int_t iCh = iDetElemId / 100;
  const Int_t iDet = iDetElemId % 100;

  // make sure detector index is valid
  if (!(iCh > 0 && iCh <= fgNCh && iDet < fgNDetElemCh[iCh - 1])) {
    LOG(fatal) << "Invalid detector element id: " << iDetElemId;
  }

  // add number of detectors up to this chamber
  return iDet + fgSNDetElemCh[iCh - 1];
}


//_________________________________________________________________________________________________
Int_t GetDetElemId(Int_t iDetElemNumber) {
  // make sure detector number is valid
  if (!(iDetElemNumber >= fgSNDetElemCh[0] &&
        iDetElemNumber < fgSNDetElemCh[fgNCh])) {
    LOG(fatal) << "Invalid detector element number: " << iDetElemNumber;
  }
  /// get det element number from ID
  // get chamber and element number in chamber
  int iCh = 0;
  int iDet = 0;
  for (int i = 1; i <= fgNCh; i++) {
    if (iDetElemNumber < fgSNDetElemCh[i]) {
      iCh = i;
      iDet = iDetElemNumber - fgSNDetElemCh[i - 1];
      break;
    }
  }

  // make sure detector index is valid
  if (!(iCh > 0 && iCh <= fgNCh && iDet < fgNDetElemCh[iCh - 1])) {
    LOG(fatal) << "Invalid detector element id: " << 100 * iCh + iDet;
  }

  // add number of detectors up to this chamber
  return 100 * iCh + iDet;
}


