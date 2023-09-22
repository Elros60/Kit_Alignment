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
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>

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


struct TrackInfo {
  TrackInfo(const mch::TrackMCH &mch) : mchTrack(mch) {}

  const mch::TrackMCH &mchTrack;
  mch::TrackParam paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
  std::vector<const mch::Digit *> digits{};
  double mchTime = -1.;
  double mchTimeRMS = 0.;
  double mchTimeSt12 = -1.;
  double mchTimeRMSSt12 = 0.;
  double mchTimeSt345 = -1.;
  double mchTimeRMSSt345 = 0.;
  std::vector<const mch::Digit *> digitsAtClusterPos{};
  double mchTimeAtClusterPos = -1.;
  double mchTimeRMSAtClusterPos = 0.;
  double mchTimeAtClusterPosSt12 = -1.;
  double mchTimeRMSAtClusterPosSt12 = 0.;
  double mchTimeAtClusterPosSt345 = -1.;
  double mchTimeRMSAtClusterPosSt345 = 0.;
  int midTime = -1;
};

const int fgNCh = 10;
const int fgNDetElemCh[fgNCh] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
const int fgSNDetElemCh[fgNCh + 1] = {0,  4,  8,   12,  16, 34,
                                      52, 78, 104, 130, 156};

///////////////////////////////////////////////////////////////////////////////
o2::mch::geo::TransformationCreator transformation;
o2::mch::Alignment *test_align = new o2::mch::Alignment();
o2::mch::TrackFitter *trackFitter = new o2::mch::TrackFitter();

std::map<int, o2::math_utils::Transform3D> transformRef; // reference geometry w.r.t track data
std::map<int, o2::math_utils::Transform3D> transformNew; // new geometry
///////////////////////////////////////////////////////////////////////////////

static const double muMass =
    TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
uint16_t minNSamplesSignal = 17;
double signalParam[4] = {80., 16., 12., 1.2};
uint16_t minNSamplesBackground = 14;
double backgroundParam[4] = {18., 24., -20., 7.};
int bcIntegrationRange = 6; // time window ([-range, range]) to integrate digits
int minNDigitsSignal = 10;  // minimum number of digits passing the signal cuts
                            // to select signal events

constexpr double pi() { return 3.14159265358979323846; }
std::tuple<TFile *, TTreeReader *> LoadData(const char *fileName,
                                            const char *treeName);
void LoadDigits(TrackInfo &trackInfo, const std::vector<mch::Cluster> &clusters,
                const std::vector<mch::Digit> &digits, bool selectSignal,
                bool rejectBackground);
void computeMCHTime(const std::vector<const mch::Digit *> &digits, double &mean,
                    double &rms, int deMin = 100, int deMax = 1025);
bool FindMuon(int iMCHTrack, std::vector<dataformats::TrackMCHMID> &muonTracks);
bool ExtrapToVertex(TrackInfo &trackInfo);
bool IsSelected(TrackInfo &trackInfo);
bool IsSignal(TrackInfo &trackInfo);
bool IsReconstructible(TrackInfo &trackInfo);
void CreateHistosAtVertex(std::vector<TH1 *> &histos, const char *extension);
void FillHistosAtVertex(const TrackInfo &trackInfo, std::vector<TH1 *> &histos);
void DrawHistosAtVertex(std::vector<TH1 *> histos[2]);
void CreateTimeHistos(std::vector<TH1 *> &histos, const char *extension);
void FillTimeHistos(const std::vector<const mch::Digit *> &digits,
                    double mchTime, double mchTimeRMS, int midTime,
                    gsl::span<TH1 *> histos, int deMin = 100, int deMax = 1025);
void DrawTimeHistos(gsl::span<TH1 *> histos, const char *extension);
void CreateChargeHistos(std::vector<TH1 *> &histos, const char *extension);
void FillChargeHistos(const std::vector<const mch::Digit *> &digits,
                      gsl::span<TH1 *> histos, int deMin = 100,
                      int deMax = 1025);
void DrawChargeHistos(gsl::span<TH1 *> histos, const char *extension);
void CreateCorrelationHistos(std::vector<TH1 *> &histos);
void FillCorrelationHistos(const std::vector<const mch::Digit *> &digits,
                           TH1 *hist, double timeRef, int deMin = 100,
                           int deMax = 1025);
void DrawCorrelationHistos(std::vector<TH1 *> &histos);
double signalCut(double *x, double *p);
double backgroundCut(double *x, double *p);
void WriteHistos(TFile *f, const char *dirName,
                 const std::vector<TH1 *> &histos);

mch::Track MCHFormatConvert(mch::TrackMCH &mchTrack,
                            std::vector<mch::Cluster> &mchClusters, bool doReAlign);
bool RemoveTrack(mch::Track &track, double ImproveCut);
void drawHisto(double *params, double *errors, double *pulls, TTree &Res_Tree, std::string outFileName);

Int_t GetDetElemNumber(Int_t iDetElemId);
Int_t GetDetElemId(Int_t iDetElemNumber);


// Load gSystem->Load("libO2MCHMappingImpl4"); in ROOT before compile the marco.

//_________________________________________________________________________________________________
void mchAlignment(std::string prefix, std::string mchFileName,
                     std::string muonFileName = "", std::string recDataFileName = "",
                     std::string recConsFileName = "",
                     std::string outFileName = "Alignment",
                     std::string RefGeoFileName = "",
                     std::string NewGeoFileName = "",
                     bool doAlign = false,
                     bool doReAlign = false,
                     std::string param_config = "pp",
                     Double_t weightRecord = 1) {



  // prefix defines the path for geometry file

  auto &segmentation = mch::mapping::segmentation(300);
  if (segmentation.nofPads() != 27873) {
    LOG(error) << "wrong mapping implementation";
    LOG(error) << "do gSystem->Load(\"libO2MCHMappingImpl4\") before compiling "
                  "this macro";
    exit(-1);
  }

  double Reso_X;
  double Reso_Y;
  double ImproveCut;

  if(param_config == "PbPb"){
    Reso_X = 0.2;
    Reso_Y = 0.2;
    ImproveCut = 4.0;
    LOG(info) << "Using PbPb parameter set.";
  }else if(param_config == "pp"){
    Reso_X = 0.4;
    Reso_Y = 0.4;
    ImproveCut = 6.0;
    LOG(info) << "Using pp parameter set.";
  }else{
    LOG(fatal) << "Please enter a correct parameter configuration option.";
    exit(-1);
  }
  // pp set: Reso 0.4   Sigma Improve 6.0
  // PbPb set: Reso 0.2   Sigma Improve 4.0

  ////////////////////////////////////////////////////////////////////////////
  // load magnetic field (from o2sim_grp.root) and geometry (from           //
  // o2sim_geometry.root) and prepare track extrapolation to vertex (0,0,0) //
  ////////////////////////////////////////////////////////////////////////////
  

  cout << "Loading magnetic field and geometry from: " << prefix
            << endl;

  const auto grp = parameters::GRPObject::loadFrom(base::NameConf::getGRPFileName());
  base::Propagator::initFieldFromGRP(grp);
  mch::TrackExtrap::setField(); // need to set field value in TrackExtrap.cxx
  mch::TrackExtrap::useExtrapV2();

  // Config for trackfitter(will be used in track conversion)
  trackFitter->initField(grp->getL3Current(), grp->getDipoleCurrent());
  trackFitter->smoothTracks(true);
  trackFitter->setChamberResolution(Reso_X, Reso_Y);
  trackFitter->useChamberResolution();


  // Store current geometry(w.r.t tracks data) into transformations
  base::GeometryManager::loadGeometry(RefGeoFileName.c_str());
  transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);
  for (int i = 0; i < 156; i++) { 
    int iDEN = GetDetElemId(i);
    transformRef[iDEN] = transformation(iDEN);
  }


  // Store new geometry transformations for evaluation
  if(doReAlign){
    base::GeometryManager::loadGeometry(NewGeoFileName.c_str());
    transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);
    for (int i = 0; i < 156; i++) {
      int iDEN = GetDetElemId(i);
      transformNew[iDEN] = transformation(iDEN);
    }
  }

  

  /////////////////
  // Load tracks //
  /////////////////

  
  // Input file: MCH -> mchtracks.root
  cout << "Reading tracks..." << endl;
  cout << "Loading MCH tracks..." <<endl;
  
  std::string DataFile = Form("%s/%s", prefix.c_str(), mchFileName.c_str());
  auto [fMCH, mchReader] = LoadData(DataFile.c_str(), "o2sim");

  TTreeReaderValue<std::vector<mch::ROFRecord>> mchROFs = {*mchReader,
                                                           "trackrofs"};
  TTreeReaderValue<std::vector<mch::TrackMCH>> mchTracks = {*mchReader,
                                                            "tracks"};
  TTreeReaderValue<std::vector<mch::Cluster>> mchClusters = {*mchReader,
                                                             "trackclusters"};


  
  
  /*                                                           
  // For tracking matching MCH-MID
  // Input file: MUON -> muontracks; contains MCH–MID matching
  LOG(info) << "matching mode";
  cout << "Loading MID muon tracks..." <<endl;
  auto [fMUON, muonReader] = LoadData(muonFileName.c_str(), "o2sim");
  TTreeReaderValue<std::vector<dataformats::TrackMCHMID>> muonTracks = {*muonReader, "tracks"};
  int nTF = muonReader->GetEntries(false);
  if (mchReader->GetEntries(false) != nTF) {
    LOG(error) << mchFileName << " and " << muonFileName
             << " do not contain the same number of TF";
    exit(-1);
  }     
  */                                             


  // Reading data
  /*                                                          
  TChain *mchChain = new TChain("o2sim");
  TChain *muonChain = new TChain("o2sim");
  int runNumbers[] = {520495, 520496};
  std::string DataFile;

  for (int iRun = 0; iRun < 1; iRun++) {
    DataFile = Form("%s/%d/%s", prefix.c_str(), runNumbers[iRun], mchFileName.c_str());
    mchChain->AddFile(DataFile.c_str());
    // DataFile = Form("%s/%d/%s", prefix.c_str(), runNumbers[iRun], muonFileName.c_str());
    // muonChain->AddFile(DataFile.c_str());
  }


  TTreeReader *mchReader = new TTreeReader(mchChain);
  TTreeReaderValue<std::vector<mch::ROFRecord>> mchROFs = {*mchReader, "trackrofs"};
  TTreeReaderValue<std::vector<mch::TrackMCH>> mchTracks = {*mchReader, "tracks"};
  TTreeReaderValue<std::vector<mch::Cluster>> mchClusters = {*mchReader, "trackclusters"};
  
  */

  ///////////////////////                                                           
  //  Start alignment  //
  ///////////////////////


  LOG(info) << "Start alignment...";

  //////////////////////////////////////////////////////////////
  // Output files: 1. alignment records 2. aligned parameters //
  //////////////////////////////////////////////////////////////


  //////////////////////////////////////////
  // Configurations for alignment process //
  //////////////////////////////////////////

  test_align->SetDoEvaluation(kTRUE);


  // Fix chambers
  const Int_t chambers[] = {1,6,0};
  for (Int_t i = 0; chambers[i] > 0; ++i) {
    std::cout << "Fixing chamber " << chambers[i] << std::endl;
    test_align->FixChamber(chambers[i]);
  }


  // Variation range for parameters
  test_align->SetAllowedVariation(0, 2.0);
  test_align->SetAllowedVariation(1, 0.3);
  test_align->SetAllowedVariation(2, 0.002);
  test_align->SetAllowedVariation(3, 2.0);


  // Initialize alignment algorithm
  test_align->init(recDataFileName, recConsFileName);
  test_align->SetBFieldOn(mch::TrackExtrap::isFieldON());

  /////////////////////////////////////////////////////////////////////////
  // create root file for saving results of alignment process(GlobalFit) //
  // 3 branches: params, errors, pulls                                   //
  /////////////////////////////////////////////////////////////////////////

  const Int_t kSplitlevel = 98;
  const Int_t kBufsize = 32000;

  int NGlobalPar = test_align->fNGlobal;

  double params[624];
  double errors[624];
  double pulls[624];


  // Open output file if any
  std::string Align_file = Form("%s%s",outFileName.c_str(),"_AlignParam.root");

  //////////////////////////////////////
  // Main loop for alignement process //
  //////////////////////////////////////

  cout << "Start processing..." << endl;
  cout <<endl;
  int tracksGood = 0;
  int tracksGoodwithoutFit = 0;
  int tracksAll = 0;

  // processing for each track
  while (mchReader->Next()) {
    int id_event = mchReader->GetCurrentEntry();
    cout << "=========================================================="
                "===="
                "=="
                "======="
             << endl;

    cout << "=========================================================="
                "===="
                "=="
                "======="
             << endl;         
    cout << "Start reading event: " << id_event <<endl;
    tracksAll += mchTracks->size();
    cout << "   "<< mchROFs->size() << " MCH ROF records loaded..." <<endl;
    cout << "   "<< mchTracks->size() << " MCH tracks loaded..." <<endl;
    cout << "=========================================================="
                "===="
                "=="
                "======="
             << endl;
    cout << "=========================================================="
                "===="
                "=="
                "======="
             << endl;
    for (const auto &mchROF : *mchROFs) {

      for (int iMCHTrack = mchROF.getFirstIdx();
        iMCHTrack <= mchROF.getLastIdx(); ++iMCHTrack) {        

        auto mchTrack = mchTracks->at(iMCHTrack);
        int id_track = iMCHTrack;
        int nb_clusters = mchTrack.getNClusters();

        // Track selection, saving only tracks having exactly 10 clusters
        if(nb_clusters <= 9) continue;
        tracksGoodwithoutFit += 1;

        cout << "=========================================================="
                "===="
                "=="
                "======="
                << endl;
        cout << "Start processing for track: " << id_track << " at event: "<< id_event <<endl;

        cout << nb_clusters << " clusters attached for track " << id_track << endl;


        // Format conversion from TrackMCH to Track(MCH internal use)
        mch::Track convertedTrack = MCHFormatConvert(mchTrack, *mchClusters, doReAlign);

        // Erase removable track
        if(RemoveTrack(convertedTrack, ImproveCut)){
          LOG(info) << "Current track has been excluded for alignment process.";
          continue;
        }else{
          LOG(info) << "Current track is being considered into alignment process.";
          tracksGood += 1;
        }

        //  Track processing, saving residuals
        AliMillePedeRecord *mchRecord = test_align->ProcessTrack(
            convertedTrack, transformation, doAlign, weightRecord);

        cout << "Processing done for track: " << iMCHTrack << endl;
        cout << "=========================================================="
                "===="
                "=="
                "======="
             << endl;
        
      }
    }
    cout << endl;
    cout << endl;
    cout << endl;

  }


  cout << "Start global fitting..." << endl;
  cout << "=============================================================="
          "======"
          "==="
       << endl;


  // Process global fit for each track
  if(doAlign) test_align->GlobalFit(params, errors, pulls);


  
  cout << "=============================================================="
          "======"
          "==="
       << endl;
  cout << "Global fitting done." << endl;
  cout << "=============================================================="
          "======"
          "==="
       << endl;

  cout <<endl;
  cout <<endl;


  // Evaluation for track removing and selection
  LOG(info) << Form("%s%d", "Number of good tracks used in alignment process: ",tracksGood);
  LOG(info) << Form("%s%d", "Number of good tracks without fit processing: ",tracksGoodwithoutFit);
  LOG(info) << Form("%s%d","Total number of tracks loaded: ", tracksAll);
  LOG(info) << Form("%s%f","Ratio before fit: ", double(tracksGoodwithoutFit)/tracksAll);    
  LOG(info) << Form("%s%f","Ratio after fit: ", double(tracksGood)/tracksAll);
  cout<<endl;
  cout<<endl;

  // Generate new geometry w.r.t alignment results
  if(doAlign){

    LOG(info) << "Generating new geometry using global parameters...";
    std::vector<o2::detectors::AlignParam> ParamAligned;
    test_align->ReAlign(ParamAligned, params);

    TFile *FileAlign = TFile::Open(Align_file.c_str(), "RECREATE");
    FileAlign->cd();
    FileAlign->WriteObjectAny(&ParamAligned, "std::vector<o2::detectors::AlignParam>", "alignment");
    FileAlign->Close();

    string Geo_file;

    if(doReAlign){
      Geo_file = Form("%s%s","o2sim_geometry_ReAlign",".root");
    }else{
      Geo_file = Form("%s%s","o2sim_geometry_Align",".root");
    }

    // Store aligned geometry
    gGeoManager->Export(Geo_file.c_str());

    // Store param plots
    drawHisto(params, errors, pulls,*(test_align->GetResTree()), outFileName);

  }

  // Close files and store all tracks' records
  test_align->terminate();
  LOG(info) << "Alignment finished";
  LOG(info) << "Test done!";

}



//_________________________________________________________________________________________________
mch::Track MCHFormatConvert(mch::TrackMCH &mchTrack,
                            std::vector<mch::Cluster> &mchClusters, bool doReAlign) {


  LOG(info) << "Start track format conversion...";

          
  mch::Track convertedTrack = mch::Track();
  auto Param0 = mchTrack.getParameters();
  double Z0 = mchTrack.getZ();
  mch::TrackParam extrapParam = mch::TrackParam(Z0, Param0);


  // Get clusters for current track
  int id_cluster_first = mchTrack.getFirstClusterIdx();
  int id_cluster_last = mchTrack.getLastClusterIdx();
  cout << "Current track's cluster starts at index: " << id_cluster_first
       << " stops at: " << id_cluster_last << endl;


  for (int id_cluster = id_cluster_first;
       id_cluster < id_cluster_last + 1; ++id_cluster) {


    mch::Cluster *cluster = &(mchClusters.at(id_cluster));
    const int DEId_cluster = cluster->getDEId();
    const int CId_cluster = cluster->getChamberId();
    const int ind_cluster = cluster->getClusterIndex();


    // Transformations to new geometry from reference geometry
    if(doReAlign){

      o2::math_utils::Point3D<double> local;
      o2::math_utils::Point3D<double> master;

      
      master.SetXYZ(cluster->getX(), cluster->getY(), cluster->getZ());

      transformRef[cluster->getDEId()].MasterToLocal(master, local);
      transformNew[cluster->getDEId()].LocalToMaster(local, master);
      
      cluster->x = master.x();
      cluster->y = master.y();
      cluster->z = master.z();

    }

    convertedTrack.createParamAtCluster(*cluster);

  }

  // Get trackparameters by calling trackFitter
  //convertedTrack.print();
  //trackFitter->fit(convertedTrack,false,false);
  //convertedTrack.print();

  return mch::Track(convertedTrack);

}


//_________________________________________________________________________________________________
bool RemoveTrack(mch::Track &track, double ImproveCut){

  const double maxChi2Cluster = 2*ImproveCut*ImproveCut;
  bool removeTrack = false;

  try{
    trackFitter->fit(track, false);
  }catch(exception const& e){
    removeTrack = true;
    return removeTrack;
  }

  auto itStartingParam = std::prev(track.rend());
  //if(track.isRemovable()) return true;
  while(true){

    try {
        trackFitter->fit(track, true, false, (itStartingParam == track.rbegin()) ? nullptr : &itStartingParam);
      } catch (exception const&) {
        removeTrack = true;
        break;
    }

    double worstLocalChi2 = -1.0;

    track.tagRemovableClusters(0x1F, false);

    auto itWorstParam = track.end();

    for(auto itParam = track.begin(); itParam != track.end(); ++itParam){
      if(itParam->getLocalChi2() > worstLocalChi2){
        worstLocalChi2 = itParam->getLocalChi2();
        itWorstParam = itParam;
      }
    }

    if(worstLocalChi2 < maxChi2Cluster) break;

    if(!itWorstParam->isRemovable()){
        removeTrack = true;
        track.removable();
        break;
    }


    auto itNextParam = track.removeParamAtCluster(itWorstParam);
    auto itNextToNextParam = (itNextParam == track.end()) ? itNextParam : std::next(itNextParam);
    itStartingParam = track.rbegin();

    if(track.getNClusters()<10){
      removeTrack = true;
      break;
    }else{
      while (itNextToNextParam != track.end()) {
        if (itNextToNextParam->getClusterPtr()->getChamberId() != itNextParam->getClusterPtr()->getChamberId()) {
          itStartingParam = std::make_reverse_iterator(++itNextParam);
          break;
        }
        ++itNextToNextParam;
      }
    }


  }

  if(!removeTrack){
    for (auto& param : track) {
        param.setParameters(param.getSmoothParameters());
        param.setCovariances(param.getSmoothCovariances());
    }
    LOG(info) << "Track has been finalised.";
  }

  return removeTrack;

}

//_________________________________________________________________________________________________
void drawHisto(double *params, double *errors, double *pulls, TTree &Res_Tree, std::string outFileName){

  TH1F *hPullX = new TH1F("hPullX", "hPullX", 201, -10, 10);
  TH1F *hPullY = new TH1F("hPullY", "hPullY", 201, -10, 10);
  TH1F *hPullZ = new TH1F("hPullZ", "hPullZ", 201, -10, 10);
  TH1F *hPullPhi = new TH1F("hPullPhi", "hPullPhi", 201, -10, 10);

  double deNumber[156];

  double alignX[156];
  double alignY[156];
  double alignZ[156];
  double alignPhi[156];
  double pullX[156];
  double pullY[156];
  double pullZ[156];
  double pullPhi[156];

  for (int iDEN = 0; iDEN < 156; iDEN++) {
    deNumber[iDEN] = iDEN + 0.5;
    alignX[iDEN] = params[iDEN * 4];
    alignY[iDEN] = params[iDEN * 4 + 1];
    alignZ[iDEN] = params[iDEN * 4 + 3];
    alignPhi[iDEN] = params[iDEN * 4 + 2];
    pullX[iDEN] = pulls[iDEN * 4];
    pullY[iDEN] = pulls[iDEN * 4 + 1];
    pullZ[iDEN] = pulls[iDEN * 4 + 3];
    pullPhi[iDEN] = pulls[iDEN * 4 + 2];
    if (params[iDEN * 4]) {

      hPullX->Fill(pulls[iDEN * 4]);
      hPullY->Fill(pulls[iDEN * 4 + 1]);
      hPullZ->Fill(pulls[iDEN * 4 + 3]);
      hPullPhi->Fill(pulls[iDEN * 4 + 2]);
    }
  }

  TGraph *graphAlignX = new TGraph(156, deNumber, alignX);
  TGraph *graphAlignY = new TGraph(156, deNumber, alignY);
  TGraph *graphAlignZ = new TGraph(156, deNumber, alignZ);
  TGraph *graphAlignPhi = new TGraph(156, deNumber, alignPhi);
  TGraph *graphAlignYZ = new TGraph(156, alignY, alignZ);

  TGraph *graphPullX = new TGraph(156, deNumber, pullX);
  TGraph *graphPullY = new TGraph(156, deNumber, pullY);
  TGraph *graphPullZ = new TGraph(156, deNumber, pullZ);
  TGraph *graphPullPhi = new TGraph(156, deNumber, pullPhi);


  graphAlignX->SetMarkerStyle(24);
  graphPullX->SetMarkerStyle(25);

  //   graphAlignX->Draw("AP");

  graphAlignY->SetMarkerStyle(24);
  graphPullY->SetMarkerStyle(25);

  // graphAlignY->Draw("Psame");

  graphAlignZ->SetMarkerStyle(24);
  graphPullZ->SetMarkerStyle(25);

  //   graphAlignZ->Draw("AP");
  graphAlignPhi->SetMarkerStyle(24);
  graphPullPhi->SetMarkerStyle(25);

  graphAlignYZ->SetMarkerStyle(24);
  // graphAlignYZ->Draw("AP");

  // Saving plots
  std::string PlotFiles_name = Form("%s%s",outFileName.c_str(),"_results.root");
  TFile *PlotFiles = TFile::Open(PlotFiles_name.c_str(),"RECREATE");
  PlotFiles->WriteObjectAny(hPullX,"TH1F","hPullX");
  PlotFiles->WriteObjectAny(hPullY,"TH1F","hPullY");
  PlotFiles->WriteObjectAny(hPullZ,"TH1F","hPullZ");
  PlotFiles->WriteObjectAny(hPullPhi,"TH1F","hPullPhi");
  PlotFiles->WriteObjectAny(graphAlignX,"TGraph","graphAlignX");
  PlotFiles->WriteObjectAny(graphAlignY,"TGraph","graphAlignY");
  PlotFiles->WriteObjectAny(graphAlignZ,"TGraph","graphAlignZ");
  PlotFiles->WriteObjectAny(graphAlignYZ,"TGraph","graphAlignYZ");
  

  TCanvas *cvn1 = new TCanvas("cvn1", "cvn1", 1200, 1600);
  //cvn1->Draw();
  cvn1->Divide(1, 4);
  TLine limLine(4, -5, 4, 5);
  TH1F *aHisto = new TH1F("aHisto", "AlignParam", 161, 0, 160);
  aHisto->SetXTitle("Det. Elem. Number");
  for (int i = 1; i < 5; i++) {
    cvn1->cd(i);
    double Range[4] = {5.0, 1.0, 5.0, 0.01};
    switch (i) {
    case 1:
      aHisto->SetYTitle("#delta_{#X} (cm)");
      aHisto->GetYaxis()->SetRangeUser(-5.0, 5.0);
      aHisto->DrawCopy();
      graphAlignX->Draw("Psame");
      limLine.DrawLine(4, -Range[i-1], 4, Range[i-1]);
      limLine.DrawLine(8, -Range[i-1], 8, Range[i-1]);
      limLine.DrawLine(12, -Range[i-1], 12, Range[i-1]);
      limLine.DrawLine(16, -Range[i-1], 16, Range[i-1]);
      limLine.DrawLine(16 + 18, -Range[i-1], 16 + 18, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18, -Range[i-1], 16 + 2 * 18, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18 + 26, -Range[i-1], 16 + 2 * 18 + 26, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18 + 2 * 26, -Range[i-1], 16 + 2 * 18 + 2 * 26, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18 + 3 * 26, -Range[i-1], 16 + 2 * 18 + 3 * 26, Range[i-1]);
      break;
    case 2:
      aHisto->SetYTitle("#delta_{#Y} (cm)");
      aHisto->GetYaxis()->SetRangeUser(-1.0, 1.0);
      aHisto->DrawCopy();
      graphAlignY->Draw("Psame");
      limLine.DrawLine(4, -Range[i-1], 4, Range[i-1]);
      limLine.DrawLine(8, -Range[i-1], 8, Range[i-1]);
      limLine.DrawLine(12, -Range[i-1], 12, Range[i-1]);
      limLine.DrawLine(16, -Range[i-1], 16, Range[i-1]);
      limLine.DrawLine(16 + 18, -Range[i-1], 16 + 18, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18, -Range[i-1], 16 + 2 * 18, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18 + 26, -Range[i-1], 16 + 2 * 18 + 26, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18 + 2 * 26, -Range[i-1], 16 + 2 * 18 + 2 * 26, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18 + 3 * 26, -Range[i-1], 16 + 2 * 18 + 3 * 26, Range[i-1]);      
      break;
    case 3:
      aHisto->SetYTitle("#delta_{#Z} (cm)");
      aHisto->GetYaxis()->SetRangeUser(-5.0, 5.0);
      aHisto->DrawCopy();
      graphAlignZ->Draw("Psame");
      limLine.DrawLine(4, -Range[i-1], 4, Range[i-1]);
      limLine.DrawLine(8, -Range[i-1], 8, Range[i-1]);
      limLine.DrawLine(12, -Range[i-1], 12, Range[i-1]);
      limLine.DrawLine(16, -Range[i-1], 16, Range[i-1]);
      limLine.DrawLine(16 + 18, -Range[i-1], 16 + 18, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18, -Range[i-1], 16 + 2 * 18, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18 + 26, -Range[i-1], 16 + 2 * 18 + 26, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18 + 2 * 26, -Range[i-1], 16 + 2 * 18 + 2 * 26, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18 + 3 * 26, -Range[i-1], 16 + 2 * 18 + 3 * 26, Range[i-1]);  
      break;
    case 4:
      aHisto->SetYTitle("#delta_{#varphi} (cm)");
      aHisto->GetYaxis()->SetRangeUser(-0.01, 0.01);
      aHisto->DrawCopy();
      graphAlignPhi->Draw("Psame");
      limLine.DrawLine(4, -Range[i-1], 4, Range[i-1]);
      limLine.DrawLine(8, -Range[i-1], 8, Range[i-1]);
      limLine.DrawLine(12, -Range[i-1], 12, Range[i-1]);
      limLine.DrawLine(16, -Range[i-1], 16, Range[i-1]);
      limLine.DrawLine(16 + 18, -Range[i-1], 16 + 18, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18, -Range[i-1], 16 + 2 * 18, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18 + 26, -Range[i-1], 16 + 2 * 18 + 26, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18 + 2 * 26, -Range[i-1], 16 + 2 * 18 + 2 * 26, Range[i-1]);
      limLine.DrawLine(16 + 2 * 18 + 3 * 26, -Range[i-1], 16 + 2 * 18 + 3 * 26, Range[i-1]);  
      break;
    }
  }

  Int_t RefClDetElem;
  Int_t RefClDetElemNumber;
  Float_t RefClusterX;
  Float_t RefClusterY;
  Float_t RefTrackX;
  Float_t RefTrackY;
  Float_t RefTrackSlopeX;
  Float_t RefTrackSlopeY;
  Res_Tree.SetBranchAddress("fClusterX",&RefClusterX);
  Res_Tree.SetBranchAddress("fClusterY",&RefClusterY);
  Res_Tree.SetBranchAddress("fTrackX",&RefTrackX);
  Res_Tree.SetBranchAddress("fTrackY",&RefTrackY);
  Res_Tree.SetBranchAddress("fTrackSlopeX",&RefTrackSlopeX);
  Res_Tree.SetBranchAddress("fTrackSlopeY",&RefTrackSlopeY);
  Res_Tree.SetBranchAddress("fClDetElem",&RefClDetElem);
  Res_Tree.SetBranchAddress("fClDetElemNumber",&RefClDetElemNumber);

  TH1F *Histos_Res[2][11];
  for(int i=0;i<2;i++){
    for(int j=0;j<11;j++){
      if(i==0){
        Histos_Res[i][j] = new TH1F(Form("%s%d","Residual_X_Ch",j),Form("%s%d","Residual_x_Ch",j),200,-5,5);
      }

      if(i==1){
        Histos_Res[i][j] = new TH1F(Form("%s%d","Residual_Y_Ch",j),Form("%s%d","Residual_y_Ch",j),200,-5,5);
      }
    }
  }


  TH1F *Histos_DETRes[2][156];
  for(int i=0;i<2;i++){
    for(int j=0;j<156;j++){
      if(i==0) Histos_DETRes[i][j] = new TH1F(Form("%s%d","Hist_x_DET",j+1),Form("%s%d","Hist_x_DET",j+1),200,-5,5);
      if(i==1) Histos_DETRes[i][j] = new TH1F(Form("%s%d","Hist_y_DET",j+1),Form("%s%d","Hist_y_DET",j+1),200,-5,5);
    }
  }

  int Ref_NbEntries = Res_Tree.GetEntries();
  for(int i=0; i < Ref_NbEntries;i++){

    Res_Tree.GetEntry(i);
    double Res_X = RefClusterX - RefTrackX;
    double Res_Y = RefClusterY - RefTrackY;

    for(int iCh=0; iCh<11;iCh++){
      if(iCh==0){
        Histos_Res[0][iCh]->Fill(Res_X);
        Histos_Res[1][iCh]->Fill(Res_Y);
      }else{
        if(iCh == int(RefClDetElem/100)){
          Histos_Res[0][iCh]->Fill(Res_X);
          Histos_Res[1][iCh]->Fill(Res_Y);
        }
      }
    }

    for(int iDEN = 0; iDEN < 156; iDEN++){
      if(RefClDetElemNumber == iDEN){
        Histos_DETRes[0][iDEN]->Fill(Res_X);
        Histos_DETRes[1][iDEN]->Fill(Res_Y);
      }

    }

  }

  for(int i=0;i<2;i++){
    for(int j=0;j<11;j++){
      if(i==0) PlotFiles->WriteObjectAny(Histos_Res[i][j],"TH1F",Form("%s%d","Residual_X_Ch",j));
      if(i==1) PlotFiles->WriteObjectAny(Histos_Res[i][j],"TH1F",Form("%s%d","Residual_Y_Ch",j));
    }
  }

  double ResX_mean[156]={};
  double ResX_err[156]={};
  double ResY_mean[156]={};
  double ResY_err[156]={};
  for(int iDEN = 0; iDEN < 156; iDEN++){

    ResX_mean[iDEN]=Histos_DETRes[0][iDEN]->GetMean();
    ResX_err[iDEN]=Histos_DETRes[0][iDEN]->GetRMS();

    ResY_mean[iDEN]=Histos_DETRes[1][iDEN]->GetMean();
    ResY_err[iDEN]=Histos_DETRes[1][iDEN]->GetRMS();
  }

  TGraphErrors *graphResX = new TGraphErrors(156, deNumber, ResX_mean, nullptr, ResX_err);
  TGraphErrors *graphResY = new TGraphErrors(156, deNumber, ResY_mean, nullptr, ResY_err);

  TCanvas *graphRes = new TCanvas("graphRes","graphRes",1200, 1600);
  TH1F *gHisto = new TH1F("gHisto", "Residuals", 161, 0, 160);
  gHisto->SetXTitle("Det. Elem. Number");

  graphRes->Divide(1,2);

  graphRes->cd(1);
  gHisto->SetYTitle("TrackX - ClusterX (cm)");
  gHisto->GetYaxis()->SetRangeUser(-5.0, 5.0);
  gHisto->DrawCopy();
  graphResX->SetMarkerStyle(8);
  graphResX->SetMarkerSize(0.7);
  graphResX->SetLineColor(kBlue);
  graphResX->Draw("PZsame");
  limLine.DrawLine(4, -5, 4, 5);
  limLine.DrawLine(8, -5, 8, 5);
  limLine.DrawLine(12, -5, 12, 5);
  limLine.DrawLine(16, -5, 16, 5);
  limLine.DrawLine(16 + 18, -5, 16 + 18, 5);
  limLine.DrawLine(16 + 2 * 18, -5, 16 + 2 * 18, 5);
  limLine.DrawLine(16 + 2 * 18 + 26, -5, 16 + 2 * 18 + 26, 5);
  limLine.DrawLine(16 + 2 * 18 + 2 * 26, -5, 16 + 2 * 18 + 2 * 26, 5);
  limLine.DrawLine(16 + 2 * 18 + 3 * 26, -5, 16 + 2 * 18 + 3 * 26, 5);

  graphRes->cd(2);
  gHisto->SetYTitle("TrackY - ClusterY (cm)");
  gHisto->GetYaxis()->SetRangeUser(-5.0, 5.0);
  gHisto->DrawCopy();
  graphResY->SetMarkerStyle(8);
  graphResY->SetMarkerSize(0.7);
  graphResY->SetLineColor(kBlue);
  graphResY->Draw("PZsame");
  limLine.DrawLine(4, -5, 4, 5);
  limLine.DrawLine(8, -5, 8, 5);
  limLine.DrawLine(12, -5, 12, 5);
  limLine.DrawLine(16, -5, 16, 5);
  limLine.DrawLine(16 + 18, -5, 16 + 18, 5);
  limLine.DrawLine(16 + 2 * 18, -5, 16 + 2 * 18, 5);
  limLine.DrawLine(16 + 2 * 18 + 26, -5, 16 + 2 * 18 + 26, 5);
  limLine.DrawLine(16 + 2 * 18 + 2 * 26, -5, 16 + 2 * 18 + 2 * 26, 5);
  limLine.DrawLine(16 + 2 * 18 + 3 * 26, -5, 16 + 2 * 18 + 3 * 26, 5);


  PlotFiles->WriteObjectAny(cvn1,"TCanvas","AlignParam");
  PlotFiles->WriteObjectAny(graphRes,"TCanvas","ResGraph");
  PlotFiles->Delete();

}


//_________________________________________________________________________________________________
std::tuple<TFile *, TTreeReader *> LoadData(const char *fileName,
                                            const char *treeName) {
  /// open the input file and get the intput tree

  TFile *f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    LOG(error) << "opening file " << fileName << " failed";
    exit(-1);
  }

  TTreeReader *r = new TTreeReader(treeName, f);
  if (r->IsZombie()) {
    LOG(error) << "tree " << treeName << " not found";
    exit(-1);
  }

  return std::make_tuple(f, r);
}

//_________________________________________________________________________________________________
void LoadDigits(TrackInfo &trackInfo, const std::vector<mch::Cluster> &clusters,
                const std::vector<mch::Digit> &digits, bool selectSignal,
                bool rejectBackground) {
  /// fill the lists of digits associated to the track

  int nClusterOnTopOfNoDigit(0);

  for (int iCl = trackInfo.mchTrack.getFirstClusterIdx();
       iCl <= trackInfo.mchTrack.getLastClusterIdx(); ++iCl) {

    const auto &cluster = clusters[iCl];

    // get the pads at the cluster position
    math_utils::Point3D<float> global{cluster.x, cluster.y, cluster.z};
    auto t = transformation(cluster.getDEId());
    auto local = t ^ (global);
    int padIDNB(-1), padIDB(-1);
    auto &segmentation = mch::mapping::segmentation(cluster.getDEId());
    bool padsFound = segmentation.findPadPairByPosition(local.x(), local.y(),
                                                        padIDB, padIDNB);
    bool padFoundNB = padsFound || segmentation.isValid(padIDNB);
    bool padFoundB = padsFound || segmentation.isValid(padIDB);
    if (!padFoundNB && !padFoundB) {
      LOG(warning) << "cluster on top of no pad";
    }

    bool digitFound(false);
    for (uint32_t iDig = 0; iDig < cluster.nDigits; ++iDig) {
      const auto &digit = digits[cluster.firstDigit + iDig];
      if (selectSignal) {
        double nSample = digit.getNofSamples();
        if (digit.getNofSamples() < minNSamplesSignal ||
            digit.getADC() < signalCut(&nSample, signalParam)) {
          continue;
        }
      }
      if (rejectBackground) {
        double nSample = digit.getNofSamples();
        if (digit.getNofSamples() < minNSamplesBackground ||
            digit.getADC() < backgroundCut(&nSample, backgroundParam)) {
          continue;
        }
      }
      trackInfo.digits.push_back(&digit);
      if ((padFoundNB && digit.getPadID() == padIDNB) ||
          (padFoundB && digit.getPadID() == padIDB)) {
        trackInfo.digitsAtClusterPos.push_back(&digit);
        digitFound = true;
      }
    }
    if (!digitFound) {
      ++nClusterOnTopOfNoDigit;
    }
  }

  if (nClusterOnTopOfNoDigit > 0 && trackInfo.midTime >= 0) {
    LOG(warning) << "matched track with " << nClusterOnTopOfNoDigit << "/"
                 << trackInfo.mchTrack.getNClusters()
                 << " clusters on top of no digit";
  }
}

//_________________________________________________________________________________________________
void computeMCHTime(const std::vector<const mch::Digit *> &digits, double &mean,
                    double &rms, int deMin, int deMax) {
  /// compute the average time and time dispersion of MCH digits

  if (digits.empty()) {
    LOG(error) << "cannot compute mch time";
    return;
  }

  mean = 0.;
  double t2 = 0.;
  double n = 0.;
  for (const auto digit : digits) {
    if (digit->getDetID() >= deMin && digit->getDetID() <= deMax) {
      mean += digit->getTime();
      t2 += static_cast<double>(digit->getTime()) * digit->getTime();
      n += 1.;
    }
  }

  if (n == 0.) {
    LOG(error) << "cannot compute mch time";
    mean = -1.;
    return;
  }

  mean /= n;
  t2 /= n;
  rms = sqrt(t2 - mean * mean);
}

//_________________________________________________________________________________________________
bool FindMuon(int iMCHTrack, std::vector<dataformats::TrackMCHMID> &muonTracks) {
  /// find the MCH-MID matched track corresponding to this MCH track
  for (const auto &muon : muonTracks) {
    // cout << "Muon track index: " << muon.getMCHRef().getIndex()<<endl;
    if (muon.getMCHRef().getIndex() == iMCHTrack) {
      return true;
    }
  }
  return false;
}

//_________________________________________________________________________________________________
bool ExtrapToVertex(TrackInfo &trackInfo) {
  /// compute the track parameters at vertex, at DCA and at the end of the
  /// absorber return false if the propagation fails

  // extrapolate to vertex
  trackInfo.paramAtVertex.setZ(trackInfo.mchTrack.getZ());
  trackInfo.paramAtVertex.setParameters(trackInfo.mchTrack.getParameters());
  if (!mch::TrackExtrap::extrapToVertex(trackInfo.paramAtVertex, 0., 0., 0., 0.,
                                        0.)) {
    return false;
  }

  // extrapolate to DCA
  mch::TrackParam trackParamAtDCA(trackInfo.mchTrack.getZ(),
                                  trackInfo.mchTrack.getParameters());
  if (!mch::TrackExtrap::extrapToVertexWithoutBranson(trackParamAtDCA, 0.)) {
    return false;
  }
  double dcaX = trackParamAtDCA.getNonBendingCoor();
  double dcaY = trackParamAtDCA.getBendingCoor();
  trackInfo.dca = sqrt(dcaX * dcaX + dcaY * dcaY);

  // extrapolate to the end of the absorber
  mch::TrackParam trackParamAtRAbs(trackInfo.mchTrack.getZ(),
                                   trackInfo.mchTrack.getParameters());
  if (!mch::TrackExtrap::extrapToZ(trackParamAtRAbs, -505.)) {
    return false;
  }
  double xAbs = trackParamAtRAbs.getNonBendingCoor();
  double yAbs = trackParamAtRAbs.getBendingCoor();
  trackInfo.rAbs = sqrt(xAbs * xAbs + yAbs * yAbs);

  return true;
}

//_________________________________________________________________________________________________
bool IsSelected(TrackInfo &trackInfo) {
  /// apply standard track selections + pDCA

  static const double sigmaPDCA23 = 80.;
  static const double sigmaPDCA310 = 54.;
  static const double nSigmaPDCA = 6.;
  static const double relPRes = 0.0004;
  static const double slopeRes = 0.0005;

  double thetaAbs = TMath::ATan(trackInfo.rAbs / 505.) * TMath::RadToDeg();
  if (thetaAbs < 2. || thetaAbs > 10.) {
    return false;
  }

  double p = trackInfo.paramAtVertex.p();
  double eta = 0.5 * log((p + trackInfo.paramAtVertex.pz()) /
                         (p - trackInfo.paramAtVertex.pz()));
  if (eta < -4. || eta > -2.5) {
    return false;
  }

  double pDCA = trackInfo.mchTrack.getP() * trackInfo.dca;
  double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
  double nrp = nSigmaPDCA * relPRes * p;
  double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
  double slopeResEffect = 535. * slopeRes * p;
  double sigmaPDCAWithRes =
      TMath::Sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
  if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
    return false;
  }

  return true;
}

//_________________________________________________________________________________________________
bool IsSignal(TrackInfo &trackInfo) {
  /// check if the track has still enough digits in time to pass the
  /// signal selection

  int nDigits(0);

  for (const auto digit : trackInfo.digits) {
    if (digit->getTime() >= trackInfo.mchTime - bcIntegrationRange &&
        digit->getTime() <= trackInfo.mchTime + bcIntegrationRange) {
      ++nDigits;
    }
  }

  return nDigits > minNDigitsSignal;
}

//_________________________________________________________________________________________________
bool IsReconstructible(TrackInfo &trackInfo) {
  /// check if the track has still enough digits to be reconstructible

  bool hasDigits[10] = {false, false, false, false, false,
                        false, false, false, false, false};
  for (const auto digit : trackInfo.digits) {
    hasDigits[digit->getDetID() / 100 - 1] = true;
  }

  int nFiredChambersSt45 = 0;
  for (int i = 6; i < 10; ++i) {
    if (hasDigits[i]) {
      ++nFiredChambersSt45;
    }
  }

  return (hasDigits[0] || hasDigits[1]) && (hasDigits[2] || hasDigits[3]) &&
         (hasDigits[4] || hasDigits[5]) && nFiredChambersSt45 >= 3;
}

//_________________________________________________________________________________________________
void CreateHistosAtVertex(std::vector<TH1 *> &histos, const char *extension) {
  /// create single muon histograms at vertex

  histos.emplace_back(
      new TH1F(Form("pT%s", extension), "pT;p_{T} (GeV/c)", 300, 0., 30.));
  histos.emplace_back(
      new TH1F(Form("eta%s", extension), "eta;eta", 200, -4.5, -2.));
  histos.emplace_back(
      new TH1F(Form("phi%s", extension), "phi;phi", 360, 0., 360.));
  histos.emplace_back(
      new TH1F(Form("dca%s", extension), "DCA;DCA (cm)", 500, 0., 500.));
  histos.emplace_back(new TH1F(Form("pDCA23%s", extension),
                               "pDCA for #theta_{abs} < 3#circ;pDCA (GeV.cm/c)",
                               2500, 0., 5000.));
  histos.emplace_back(new TH1F(Form("pDCA310%s", extension),
                               "pDCA for #theta_{abs} > 3#circ;pDCA (GeV.cm/c)",
                               2500, 0., 5000.));
  histos.emplace_back(
      new TH1F(Form("rAbs%s", extension), "rAbs;R_{abs} (cm)", 1000, 0., 100.));
  histos.emplace_back(new TH1F(Form("nClusters%s", extension),
                               "number of clusters per track;n_{clusters}", 20,
                               0., 20.));
  histos.emplace_back(new TH1F(Form("chi2%s", extension),
                               "normalized #chi^{2};#chi^{2} / ndf", 500, 0.,
                               50.));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillHistosAtVertex(const TrackInfo &trackInfo,
                        std::vector<TH1 *> &histos) {
  /// fill single muon histograms at vertex

  double thetaAbs = TMath::ATan(trackInfo.rAbs / 505.) * TMath::RadToDeg();
  double pDCA = trackInfo.mchTrack.getP() * trackInfo.dca;
  double pT = sqrt(trackInfo.paramAtVertex.px() * trackInfo.paramAtVertex.px() +
                   trackInfo.paramAtVertex.py() * trackInfo.paramAtVertex.py());
  double p = trackInfo.paramAtVertex.p();
  double eta = 0.5 * log((p + trackInfo.paramAtVertex.pz()) /
                         (p - trackInfo.paramAtVertex.pz()));
  double phi = 180. + atan2(-trackInfo.paramAtVertex.px(),
                            -trackInfo.paramAtVertex.py()) /
                          pi() * 180.;

  histos[0]->Fill(pT);
  histos[1]->Fill(eta);
  histos[2]->Fill(phi);
  histos[3]->Fill(trackInfo.dca);
  if (thetaAbs < 3) {
    histos[4]->Fill(pDCA);
  } else {
    histos[5]->Fill(pDCA);
  }
  histos[6]->Fill(trackInfo.rAbs);
  histos[7]->Fill(trackInfo.mchTrack.getNClusters());
  histos[8]->Fill(trackInfo.mchTrack.getChi2OverNDF());
}

//_________________________________________________________________________________________________
void DrawHistosAtVertex(std::vector<TH1 *> histos[2]) {
  /// draw histograms at vertex

  // find the optimal number of pads
  int nPadsx(1), nPadsy(1);
  while ((int)histos[0].size() > nPadsx * nPadsy) {
    if (nPadsx == nPadsy) {
      ++nPadsx;
    } else {
      ++nPadsy;
    }
  }

  // draw histograms
  TCanvas *cHist =
      new TCanvas("histos", "histos", 10, 10, TMath::Max(nPadsx * 300, 1200),
                  TMath::Max(nPadsy * 300, 900));
  cHist->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos[0].size(); ++i) {
    cHist->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
    gPad->SetLogy();
    histos[0][i]->SetStats(false);
    histos[0][i]->SetLineColor(4);
    histos[0][i]->Draw();
    histos[1][i]->SetLineColor(2);
    histos[1][i]->Draw("same");
  }

  // add a legend
  TLegend *lHist = new TLegend(0.5, 0.65, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0][0],
                  Form("%g mch tracks", histos[0][0]->GetEntries()), "l");
  lHist->AddEntry(histos[1][0],
                  Form("%g muon tracks", histos[1][0]->GetEntries()), "l");
  cHist->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void CreateTimeHistos(std::vector<TH1 *> &histos, const char *extension) {
  /// create track time histograms

  histos.emplace_back(new TH1F(Form("timeResVsMCH%s", extension),
                               "#Deltat vs <MCH time>;#Deltat (BC)", 8001,
                               -2000.25, 2000.25));
  histos.emplace_back(new TH1F(Form("timeResVsMCH%sMatch", extension),
                               "#Deltat vs <MCH time>;#Deltat (BC)", 8001,
                               -2000.25, 2000.25));
  histos.emplace_back(
      new TH1F(Form("timeResVsMID%s", extension),
               "#Deltat vs MID time (matched tracks);#Deltat (BC)", 4001,
               -2000.5, 2000.5));
  histos.emplace_back(new TH1F(Form("timeRMS%s", extension),
                               "MCH time dispersion;#sigmat (BC)", 4001, -0.25,
                               2000.25));
  histos.emplace_back(new TH1F(Form("timeRMS%sMatch", extension),
                               "MCH time dispersion;#sigmat (BC)", 4001, -0.25,
                               2000.25));
  histos.emplace_back(
      new TH1F(Form("timeDiffMCHMID%s", extension),
               "<MCH time> - MID time (matched tracks);#Deltat (BC)", 8001,
               -2000.25, 2000.25));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillTimeHistos(const std::vector<const mch::Digit *> &digits,
                    double mchTime, double mchTimeRMS, int midTime,
                    gsl::span<TH1 *> histos, int deMin, int deMax) {
  /// fill track time histograms

  for (const auto digit : digits) {
    if (digit->getDetID() >= deMin && digit->getDetID() <= deMax) {
      histos[0]->Fill(digit->getTime() - mchTime);
    }
  }
  histos[3]->Fill(mchTimeRMS);

  if (midTime >= 0) {
    for (const auto digit : digits) {
      if (digit->getDetID() >= deMin && digit->getDetID() <= deMax) {
        histos[1]->Fill(digit->getTime() - mchTime);
        histos[2]->Fill(digit->getTime() - midTime);
      }
    }
    histos[4]->Fill(mchTimeRMS);
    histos[5]->Fill(mchTime - midTime);
  }
}

//_________________________________________________________________________________________________
void DrawTimeHistos(gsl::span<TH1 *> histos, const char *extension) {
  /// draw track time histograms

  TCanvas *cHist = new TCanvas(Form("cTime%s", extension),
                               Form("cTime%s", extension), 10, 10, 800, 800);
  cHist->Divide(2, 2);
  cHist->cd(1);
  gPad->SetLogy();
  histos[0]->SetStats(false);
  histos[0]->SetLineColor(4);
  histos[0]->Draw();
  histos[1]->SetLineColor(2);
  histos[1]->Draw("same");
  cHist->cd(2);
  gPad->SetLogy();
  histos[2]->SetStats(false);
  histos[2]->SetLineColor(2);
  histos[2]->Draw();
  cHist->cd(3);
  gPad->SetLogy();
  histos[3]->SetStats(false);
  histos[3]->SetLineColor(4);
  histos[3]->Draw();
  histos[4]->SetLineColor(2);
  histos[4]->Draw("same");
  cHist->cd(4);
  gPad->SetLogy();
  histos[5]->SetStats(false);
  histos[5]->SetLineColor(2);
  histos[5]->Draw();

  TLegend *lHist = new TLegend(0.1, 0.8, 0.5, 0.9);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0], "all tracks", "l");
  lHist->AddEntry(histos[1], "matched tracks", "l");
  cHist->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void CreateChargeHistos(std::vector<TH1 *> &histos, const char *extension) {
  /// create track charge histograms

  histos.emplace_back(
      new TH1F(Form("ADC%s", extension), "ADC;ADC", 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1F(Form("Samples%s", extension),
                               "N samples;N samples", 1024, -0.5, 1023.5));
  histos.emplace_back(new TH2F(Form("ADCvsSample%s", extension),
                               "ADC vs N samples (all tracks);N samples;ADC",
                               1024, -0.5, 1023.5, 10001, -0.5, 100009.5));

  histos.emplace_back(
      new TH1F(Form("ADC%sMatch", extension), "ADC;ADC", 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1F(Form("Samples%sMatch", extension),
                               "N samples;N samples", 1024, -0.5, 1023.5));
  histos.emplace_back(
      new TH2F(Form("ADCvsSample%sMatch", extension),
               "ADC vs N samples (matched tracks);N samples;ADC", 1024, -0.5,
               1023.5, 10001, -0.5, 100009.5));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillChargeHistos(const std::vector<const mch::Digit *> &digits,
                      gsl::span<TH1 *> histos, int deMin, int deMax) {
  /// fill track charge histograms
  for (const auto digit : digits) {
    if (digit->getDetID() >= deMin && digit->getDetID() <= deMax) {
      histos[0]->Fill(digit->getADC());
      histos[1]->Fill(digit->getNofSamples());
      histos[2]->Fill(digit->getNofSamples(), digit->getADC());
    }
  }
}

//_________________________________________________________________________________________________
void DrawChargeHistos(gsl::span<TH1 *> histos, const char *extension) {
  /// draw track charge histograms

  TCanvas *cHist = new TCanvas(Form("cCharge%s", extension),
                               Form("cCharge%s", extension), 10, 10, 800, 800);
  cHist->Divide(2, 2);
  cHist->cd(1);
  gPad->SetLogy();
  histos[0]->SetStats(false);
  histos[0]->SetLineColor(4);
  histos[0]->Draw();
  histos[3]->SetLineColor(2);
  histos[3]->Draw("same");
  cHist->cd(2);
  gPad->SetLogy();
  histos[1]->SetStats(false);
  histos[1]->SetLineColor(4);
  histos[1]->Draw();
  histos[4]->SetLineColor(2);
  histos[4]->Draw("same");
  cHist->cd(3);
  gPad->SetLogz();
  histos[2]->SetStats(false);
  histos[2]->Draw("colz");
  cHist->cd(4);
  gPad->SetLogz();
  histos[5]->SetStats(false);
  histos[5]->Draw("colz");

  TLegend *lHist = new TLegend(0.5, 0.8, 0.9, 0.9);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0], "all tracks", "l");
  lHist->AddEntry(histos[3], "matched tracks", "l");
  cHist->cd(1);
  lHist->Draw("same");

  static TF1 *fSignal = new TF1("fSignal", signalCut, 0, 1023, 4);
  fSignal->SetParameters(signalParam);
  fSignal->SetLineColor(2);
  static TF1 *fBackground = new TF1("fBackground", backgroundCut, 0, 1023, 4);
  fBackground->SetParameters(backgroundParam);
  fBackground->SetLineColor(4);
  cHist->cd(3);
  fSignal->Draw("same");
  fBackground->Draw("same");
  cHist->cd(4);
  fSignal->Draw("same");
  fBackground->Draw("same");
}

//_________________________________________________________________________________________________
void CreateCorrelationHistos(std::vector<TH1 *> &histos) {
  /// create correlation histograms between number of digits and total
  /// charge

  histos.emplace_back(new TH2F("ChargevsNDigits",
                               "Charge vs N digits;N digits;ADC", 100, 0, 100,
                               10000, 0, 100000));
  histos.emplace_back(new TH2F("ChargevsNDigitsSt1",
                               "Charge vs N digits (St1);N digits;ADC", 100, 0,
                               100, 10000, 0, 100000));
  histos.emplace_back(new TH2F("ChargevsNDigitsSt2",
                               "Charge vs N digits (St2);N digits;ADC", 100, 0,
                               100, 10000, 0, 100000));
  histos.emplace_back(new TH2F("ChargevsNDigitsSt345",
                               "Charge vs N digits (St345);N digits;ADC", 100,
                               0, 100, 10000, 0, 100000));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillCorrelationHistos(const std::vector<const mch::Digit *> &digits,
                           TH1 *hist, double timeRef, int deMin, int deMax) {
  /// fill correlation histograms between number of digits and total
  /// charge

  uint32_t charge(0);
  int nDigits(0);

  for (const auto digit : digits) {
    if (digit->getDetID() < deMin || digit->getDetID() > deMax) {
      continue;
    }
    if (digit->getTime() < timeRef - bcIntegrationRange ||
        digit->getTime() > timeRef + bcIntegrationRange) {
      continue;
    }
    charge += digit->getADC();
    ++nDigits;
  }

  hist->Fill(nDigits, charge);
}

//_________________________________________________________________________________________________
void DrawCorrelationHistos(std::vector<TH1 *> &histos) {
  /// draw correlation histograms between number of digits and total
  /// charge

  TCanvas *cCorr = new TCanvas("cCorr", "cCorr", 10, 10, 800, 800);
  cCorr->Divide(2, 2);
  for (int i = 0; i < 4; ++i) {
    cCorr->cd(i + 1);
    gPad->SetLogz();
    histos[i]->Draw("boxcolz");
  }
}

//_________________________________________________________________________________________________
double signalCut(double *x, double *p) {
  /// function used to select the signal
  double x0 = pow(p[0] / p[2], 1. / p[3]) + p[1];
  if (x[0] < x0) {
    return p[0];
  } else {
    return p[2] * pow(x[0] - p[1], p[3]);
  }
}

//_________________________________________________________________________________________________
double backgroundCut(double *x, double *p) {
  /// function used to select the signal
  double x0 = (p[3] * p[2] - p[1] * p[0]) / (p[3] - p[1]);
  if (x[0] < x0) {
    return p[1] * (x[0] - p[0]);
  } else {
    return p[3] * (x[0] - p[2]);
  }
}

//_________________________________________________________________________________________________
void WriteHistos(TFile *f, const char *dirName,
                 const std::vector<TH1 *> &histos) {
  /// write histograms in the subdirectory dirName

  f->mkdir(dirName, dirName, true);
  f->cd(dirName);

  for (auto h : histos) {
    h->Write();
  }
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
