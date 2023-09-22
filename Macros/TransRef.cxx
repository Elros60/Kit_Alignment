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
std::map<int, o2::math_utils::Transform3D> transformTrack;

Int_t GetDetElemNumber(Int_t iDetElemId);
Int_t GetDetElemId(Int_t iDetElemNumber);
constexpr double pi() { return 3.14159265358979323846; }


//___________________________________________________________________________
void TransRef(std::string GeoFile_Ref, std::string GeoFile_Tracks, std::string ParamFile, std::string OutParams = "AlignmentParams@IdealGeo.root"
  , std::string OutGeo = "o2sim_geometry-alignedCCDB.root"){

  	auto &segmentation = mch::mapping::segmentation(300);
  	if (segmentation.nofPads() != 27873) {
    	LOG(error) << "wrong mapping implementation";
    	LOG(error) << "do gSystem->Load(\"libO2MCHMappingImpl4\") before compiling "
                  "this macro";
    	exit(-1);
  	}


  	time_t time_now = time(0);
  	tm* now = localtime(&time_now);
  	std::string time_date = Form("%d%d%d%d%d%d",now->tm_year+1900,now->tm_mon+1,now->tm_mday,now->tm_hour,now->tm_min,now->tm_sec);



  	const auto grp = parameters::GRPObject::loadFrom(base::NameConf::getGRPFileName());
  	base::Propagator::initFieldFromGRP(grp);


    base::GeometryManager::loadGeometry(GeoFile_Tracks.c_str());
    transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);
    for (int i = 0; i < 156; i++) { 
      int iDEN = GetDetElemId(i);
      transformTrack[iDEN] = transformation(iDEN);
    }  


	  base::GeometryManager::loadGeometry(GeoFile_Ref.c_str());
  	transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);
  	for (int i = 0; i < 156; i++) { 
    	int iDEN = GetDetElemId(i);
    	transformRef[iDEN] = transformation(iDEN);
  	}



  	TFile *fParams = TFile::Open(ParamFile.c_str());
  	std::vector<o2::detectors::AlignParam> ParamsTrack = *(fParams->Get<std::vector<o2::detectors::AlignParam>>("alignment"));


  	// std::string outName = Form("%s%s%s",OutParams.c_str(), time_date.c_str(),".root");
  	TFile *fOut = TFile::Open(OutParams.c_str(),"RECREATE");



  	std::vector<o2::detectors::AlignParam> ParamsRef;
    o2::detectors::AlignParam param_Ref;
  	for (int hc = 0; hc < 20; hc++) {

  		ParamsRef.emplace_back(ParamsTrack.at(fgSNDetElemHalfCh[hc]));
  		cout << hc <<endl;

	    for (int de = 0; de < fgNDetElemHalfCh[hc]; de++) {

	    	int iDEN = fgDetElemHalfCh[hc][de];
	    	o2::detectors::AlignParam param_Track = ParamsTrack.at(fgSNDetElemHalfCh[hc]+1+de);

	    	LOG(info) << Form("%s%s","Processing DET Elem: ", (param_Track.getSymName()).c_str());

	    	TGeoHMatrix delta_track;
	    	TGeoRotation r("Rotation/Track", param_Track.getPsi()/pi()*180, param_Track.getTheta()/pi()*180, param_Track.getPhi()/pi()*180);
	    	delta_track.SetRotation(r.GetRotationMatrix());
	    	delta_track.SetDx(param_Track.getX());
	    	delta_track.SetDy(param_Track.getY());
	    	delta_track.SetDz(param_Track.getZ());

	    	TGeoHMatrix transRef = transformRef[iDEN];
	    	TGeoHMatrix transTrack = transformTrack[iDEN];
	    	TGeoHMatrix transRefTrack = transTrack * transRef.Inverse();
	    	transRefTrack.Print();
	    	TGeoHMatrix delta_ref = delta_track * transRefTrack;
	    	delta_ref.Print();

	    	param_Ref.setSymName((param_Track.getSymName()).c_str()); 
	    	param_Ref.setGlobalParams(delta_ref);
        param_Ref.applyToGeometry();
	    	ParamsRef.emplace_back(param_Ref);

   		}
	}

	fOut->WriteObjectAny(&ParamsRef, "std::vector<o2::detectors::AlignParam>",
                           "alignment");
  fOut->Close();	
  gGeoManager->Export(OutGeo.c_str());

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







