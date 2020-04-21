// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Math/interface/deltaPhi.h"
// track data formats
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/GeometrySurface/interface/PlaneBuilder.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FastSimulation/CaloGeometryTools/interface/Transform3DPJ.h"
#include "FastSimulation/Event/interface/FSimEvent.h"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include "FastSimulation/Event/interface/FSimVertex.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeGeometry/interface/MagVolumeOutsideValidity.h"
#include "TrackPropagation/RungeKutta/interface/defaultRKPropagator.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TH1F.h"
#include "TTree.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

namespace HGCal_helpers {

class coordinates {
 public:
  coordinates() : x(0), y(0), z(0), eta(100), phi(0) {}
  float x, y, z, eta, phi;
  inline math::XYZTLorentzVectorD toVector() {
    return math::XYZTLorentzVectorD(x, y, z, 0);
  }
};

class simpleTrackPropagator {
 public:
  simpleTrackPropagator(MagneticField const *f)
      : field_(f), prod_(field_, alongMomentum, 5.e-5), absz_target_(0) {
    ROOT::Math::SMatrixIdentity id;
    AlgebraicSymMatrix55 C(id);
    C *= 0.001;
    err_ = CurvilinearTrajectoryError(C);
  }
  void setPropagationTargetZ(const float &z);

  bool propagate(const double px, const double py, const double pz, const double x, const double y,
                 const double z, const float charge, coordinates &coords) const;

  bool propagate(const math::XYZTLorentzVectorD &momentum, const math::XYZTLorentzVectorD &position,
                 const float charge, coordinates &coords) const;

 private:
  simpleTrackPropagator() : field_(0), prod_(field_, alongMomentum, 5.e-5), absz_target_(0) {}
  const RKPropagatorInS &RKProp() const { return prod_.propagator; }
  Plane::PlanePointer targetPlaneForward_, targetPlaneBackward_;
  MagneticField const *field_;
  CurvilinearTrajectoryError err_;
  defaultRKPropagator::Product prod_;
  float absz_target_;
};

void simpleTrackPropagator::setPropagationTargetZ(const float &z) {
  targetPlaneForward_ = Plane::build(Plane::PositionType(0, 0, std::abs(z)), Plane::RotationType());
  targetPlaneBackward_ =
      Plane::build(Plane::PositionType(0, 0, -std::abs(z)), Plane::RotationType());
  absz_target_ = std::abs(z);
}
bool simpleTrackPropagator::propagate(const double px, const double py, const double pz,
                                      const double x, const double y, const double z,
                                      const float charge, coordinates &output) const {
  output = coordinates();

  typedef TrajectoryStateOnSurface TSOS;
  GlobalPoint startingPosition(x, y, z);
  GlobalVector startingMomentum(px, py, pz);
  Plane::PlanePointer startingPlane =
      Plane::build(Plane::PositionType(x, y, z), Plane::RotationType());
  TSOS startingStateP(
      GlobalTrajectoryParameters(startingPosition, startingMomentum, charge, field_), err_,
      *startingPlane);

  TSOS trackStateP;
  if (pz > 0) {
    trackStateP = RKProp().propagate(startingStateP, *targetPlaneForward_);
  } else {
    trackStateP = RKProp().propagate(startingStateP, *targetPlaneBackward_);
  }
  if (trackStateP.isValid()) {
    output.x = trackStateP.globalPosition().x();
    output.y = trackStateP.globalPosition().y();
    output.z = trackStateP.globalPosition().z();
    output.phi = trackStateP.globalPosition().phi();
    output.eta = trackStateP.globalPosition().eta();
    return true;
  }
  return false;
}

bool simpleTrackPropagator::propagate(const math::XYZTLorentzVectorD &momentum,
                                      const math::XYZTLorentzVectorD &position, const float charge,
                                      coordinates &output) const {
  return propagate(momentum.px(), momentum.py(), momentum.pz(), position.x(), position.y(),
                   position.z(), charge, output);
}

}  // HGCal_helpers

class HGCalAnalysis : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
 public:
  //
  // constructors and destructor
  //
  typedef ROOT::Math::Transform3DPJ::Point Point;

  // approximative geometrical values
  static constexpr float hgcalOuterRadius_ = 160.;
  static constexpr float hgcalInnerRadius_ = 25.;

  HGCalAnalysis();
  explicit HGCalAnalysis(const edm::ParameterSet &);
  ~HGCalAnalysis();

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  virtual void beginRun(edm::Run const &iEvent, edm::EventSetup const &) override;
  virtual void endRun(edm::Run const &iEvent, edm::EventSetup const &) override;

 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void endJob() override;
  virtual void fillSimHit(const DetId &detid, const float &fraction, const unsigned int &layer);

  virtual void fillRecHit(const DetId &detid, const float &fraction, const unsigned int &layer,
                          const int &cluster_index_ = -1);

  virtual void fillRecHitHF(const DetId &detid, const float &fraction, const unsigned int &layer,
                          const int &cluster_index_ = -1);

  void clearVariables();

  void retrieveLayerPositions(const edm::EventSetup &, unsigned layers);

  // ---------parameters ----------------------------
  bool readGen_;
  bool storeMoreGenInfo_;
  bool storeGenParticleExtrapolation_;
  bool storeGunParticles_;
  double layerClusterPtThreshold_;
  double propagationPtThreshold_;
  std::string detector_;
  std::string inputTag_HGCalMultiCluster_;
  bool rawRecHits_;
  bool rawSimHits_;
  bool verbose_;

  // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> recHitsEE_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitsFH_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitsBH_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitsNose_;
  edm::EDGetTokenT<HFRecHitCollection> recHitsHF_;

  edm::EDGetTokenT<std::vector<PCaloHit>> simHitsEE_;
  edm::EDGetTokenT<std::vector<PCaloHit>> simHitsFH_;
  edm::EDGetTokenT<std::vector<PCaloHit>> simHitsBH_;

  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticles_;
  edm::EDGetTokenT<edm::HepMCProduct> hev_;
  edm::EDGetTokenT<std::vector<reco::Track>> tracks_;

  TTree *t_;

  ////////////////////
  // event
  //
  edm::RunNumber_t ev_run_;
  edm::LuminosityBlockNumber_t ev_lumi_;
  edm::EventNumber_t ev_event_;
  float vtx_x_;
  float vtx_y_;
  float vtx_z_;

  ////////////////////
  // GenParticles
  //
  std::vector<float> genpart_eta_;
  std::vector<float> genpart_phi_;
  std::vector<float> genpart_pt_;
  std::vector<float> genpart_energy_;
  std::vector<float> genpart_dvx_;
  std::vector<float> genpart_dvy_;
  std::vector<float> genpart_dvz_;
  std::vector<float> genpart_ovx_;
  std::vector<float> genpart_ovy_;
  std::vector<float> genpart_ovz_;
  std::vector<float> genpart_exx_;
  std::vector<float> genpart_exy_;
  std::vector<int> genpart_mother_;
  std::vector<float> genpart_exphi_;
  std::vector<float> genpart_exeta_;
  std::vector<float> genpart_fbrem_;
  std::vector<int> genpart_pid_;
  std::vector<int> genpart_gen_;
  std::vector<int> genpart_reachedEE_;
  std::vector<bool> genpart_fromBeamPipe_;
  std::vector<std::vector<float>> genpart_posx_;
  std::vector<std::vector<float>> genpart_posy_;
  std::vector<std::vector<float>> genpart_posz_;

  ////////////////////
  // reco::GenParticles
  //
  std::vector<float> gen_eta_;
  std::vector<float> gen_phi_;
  std::vector<float> gen_pt_;
  std::vector<float> gen_energy_;
  std::vector<int> gen_charge_;
  std::vector<int> gen_pdgid_;
  std::vector<int> gen_status_;
  std::vector<std::vector<int>> gen_daughters_;

  ////////////////////
  // RecHits
  // associated to layer clusters
  std::vector<float> rechit_eta_;
  std::vector<float> rechit_phi_;
  std::vector<float> rechit_pt_;
  std::vector<float> rechit_energy_;
  std::vector<float> rechit_x_;
  std::vector<float> rechit_y_;
  std::vector<float> rechit_z_;
  std::vector<float> rechit_time_;
  std::vector<float> rechit_thickness_;
  std::vector<int> rechit_layer_;
  std::vector<int> rechit_wafer_u_;
  std::vector<int> rechit_wafer_v_;
  std::vector<int> rechit_cell_u_;
  std::vector<int> rechit_cell_v_;
  std::vector<unsigned int> rechit_detid_;
  std::vector<bool> rechit_isHalf_;
  std::vector<int> rechit_flags_;
  std::vector<int> rechit_cluster2d_;
  std::vector<float> rechit_radius_;

  ////////////////////
  // SimHits
  // associated to sim clusters
  std::vector<float> simhit_eta_;
  std::vector<float> simhit_phi_;
  std::vector<float> simhit_energy_;
  std::vector<float> simhit_x_;
  std::vector<float> simhit_y_;
  std::vector<float> simhit_z_;
  std::vector<int> simhit_layer_;
  std::vector<int> simhit_wafer_u_;
  std::vector<int> simhit_wafer_v_;
  std::vector<int> simhit_cell_u_;
  std::vector<int> simhit_cell_v_;
  std::vector<unsigned int> simhit_detid_;
  std::vector<bool> simhit_isHalf_;
  std::vector<int> simhit_flags_;

  ////////////////////
  // high purity tracks
  //
  std::vector<float> track_eta_;
  std::vector<float> track_phi_;
  std::vector<float> track_pt_;
  std::vector<float> track_energy_;
  std::vector<int> track_charge_;
  std::vector<std::vector<float>> track_posx_;
  std::vector<std::vector<float>> track_posy_;
  std::vector<std::vector<float>> track_posz_;

  // gun particles per vertex
  //
  std::vector<std::vector<int> > gunparticle_id_;
  std::vector<std::vector<float> > gunparticle_energy_;
  std::vector<std::vector<float> > gunparticle_pt_;
  std::vector<std::vector<float> > gunparticle_eta_;
  std::vector<std::vector<float> > gunparticle_phi_;


  ////////////////////
  // helper classes
  //
  unsigned int rechit_index_;
  unsigned int simhit_index_;
  std::map<DetId, const HGCRecHit *> hitmap_;
  std::map<DetId, const HFRecHit *> hfhitmap_;
  std::map<DetId, const PCaloHit *> simhitmap_;

  std::map<DetId, unsigned int> detIdToRecHitIndexMap_;
  std::map<DetId, unsigned int> detIdToSimHitIndexMap_;
  float vz_;  // primary vertex z position
  // to keep track of the 2d clusters stored within the loop on multiclusters
  std::set<edm::Ptr<reco::BasicCluster>> storedLayerClusters_;
  // to keep track of the RecHits stored within the cluster loops
  std::set<DetId> storedRecHits_;
  std::set<DetId> storedSimHits_;
  int algo_;
  HGCalDepthPreClusterer pre_;
  hgcal::RecHitTools recHitTools_;

  // -------convenient tool to deal with simulated tracks
  FSimEvent *mySimEvent_;
  edm::ParameterSet particleFilter_;
  std::vector<float> layerPositions_;


  // and also the magnetic field
  MagneticField const *aField_;

};

HGCalAnalysis::HGCalAnalysis() { ; }

HGCalAnalysis::HGCalAnalysis(const edm::ParameterSet &iConfig)
    : readGen_(iConfig.getParameter<bool>("readGenParticles")),
      storeMoreGenInfo_(iConfig.getParameter<bool>("storeGenParticleOrigin")),
      storeGenParticleExtrapolation_(iConfig.getParameter<bool>("storeGenParticleExtrapolation")),
      storeGunParticles_(iConfig.getParameter<bool>("storeGunParticles")),
      propagationPtThreshold_(iConfig.getUntrackedParameter<double>("propagationPtThreshold", 3.0)),
      detector_(iConfig.getParameter<std::string>("detector")),
      inputTag_HGCalMultiCluster_(iConfig.getParameter<std::string>("inputTag_HGCalMultiCluster")),
      rawRecHits_(iConfig.getParameter<bool>("rawRecHits")),
      rawSimHits_(iConfig.getParameter<bool>("rawSimHits")),
      verbose_(iConfig.getParameter<bool>("verbose")),
      particleFilter_(iConfig.getParameter<edm::ParameterSet>("TestParticleFilter")) {
  // now do what ever initialization is needed
  mySimEvent_ = new FSimEvent(particleFilter_);

  if (detector_ == "all") {
    recHitsEE_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCEERecHits"));
    recHitsFH_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
    recHitsBH_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
    simHitsEE_ = consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsEE"));
    simHitsFH_ = consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsHEback"));
    simHitsBH_ = consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsHEfront"));
    algo_ = 1;
  } else if (detector_ == "EM") {
    recHitsEE_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCEERecHits"));
    simHitsEE_ = consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsEE"));
    algo_ = 2;
  } else if (detector_ == "HAD") {
    recHitsFH_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
    recHitsBH_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCHEBRecHits"));
    simHitsFH_ = consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsHEback"));
    simHitsBH_ = consumes<std::vector<PCaloHit>>(edm::InputTag("g4SimHits", "HGCHitsHEfront"));
    algo_ = 3;
  } else if (detector_ == "HFNose") {
    recHitsNose_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCHFNoseRecHits"));
    recHitsHF_   = consumes<HFRecHitCollection>(iConfig.getUntrackedParameter<string>("HFRecHits","hfreco"));
    algo_ = 4;
  } else if (detector_ == "HF") {
    recHitsHF_   = consumes<HFRecHitCollection>(iConfig.getUntrackedParameter<string>("HFRecHits","hfreco"));
    algo_ = 5;
  }

  hev_ = consumes<edm::HepMCProduct>(edm::InputTag("generatorSmeared"));

  if (readGen_) {
    genParticles_ = consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles"));
  }

  tracks_ = consumes<std::vector<reco::Track>>(edm::InputTag("generalTracks"));

  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> fs;
  fs->make<TH1F>("total", "total", 100, 0, 5.);

  t_ = fs->make<TTree>("hgc", "hgc");

  // event info
  t_->Branch("event", &ev_event_);
  t_->Branch("lumi", &ev_lumi_);
  t_->Branch("run", &ev_run_);
  t_->Branch("vtx_x", &vtx_x_);
  t_->Branch("vtx_y", &vtx_y_);
  t_->Branch("vtx_z", &vtx_z_);

  t_->Branch("genpart_eta", &genpart_eta_);
  t_->Branch("genpart_phi", &genpart_phi_);
  t_->Branch("genpart_pt", &genpart_pt_);
  t_->Branch("genpart_energy", &genpart_energy_);
  t_->Branch("genpart_dvx", &genpart_dvx_);
  t_->Branch("genpart_dvy", &genpart_dvy_);
  t_->Branch("genpart_dvz", &genpart_dvz_);
  if (storeMoreGenInfo_) {
    t_->Branch("genpart_ovx", &genpart_ovx_);
    t_->Branch("genpart_ovy", &genpart_ovy_);
    t_->Branch("genpart_ovz", &genpart_ovz_);
    t_->Branch("genpart_mother", &genpart_mother_);
  }
  if (storeGenParticleExtrapolation_) {
    t_->Branch("genpart_exphi", &genpart_exphi_);
    t_->Branch("genpart_exeta", &genpart_exeta_);
    t_->Branch("genpart_exx", &genpart_exx_);
    t_->Branch("genpart_exy", &genpart_exy_);
  }
  t_->Branch("genpart_fbrem", &genpart_fbrem_);
  t_->Branch("genpart_pid", &genpart_pid_);
  t_->Branch("genpart_gen", &genpart_gen_);
  t_->Branch("genpart_reachedEE", &genpart_reachedEE_);
  t_->Branch("genpart_fromBeamPipe", &genpart_fromBeamPipe_);
  t_->Branch("genpart_posx", &genpart_posx_);
  t_->Branch("genpart_posy", &genpart_posy_);
  t_->Branch("genpart_posz", &genpart_posz_);


  if (readGen_) {
    t_->Branch("gen_eta", &gen_eta_);
    t_->Branch("gen_phi", &gen_phi_);
    t_->Branch("gen_pt", &gen_pt_);
    t_->Branch("gen_energy", &gen_energy_);
    t_->Branch("gen_charge", &gen_charge_);
    t_->Branch("gen_pdgid", &gen_pdgid_);
    t_->Branch("gen_status", &gen_status_);
    t_->Branch("gen_daughters", &gen_daughters_);
  }

  ////////////////////
  // RecHits
  // associated to layer clusters
  t_->Branch("rechit_eta", &rechit_eta_);
  t_->Branch("rechit_phi", &rechit_phi_);
  t_->Branch("rechit_pt", &rechit_pt_);
  t_->Branch("rechit_energy", &rechit_energy_);
  t_->Branch("rechit_x", &rechit_x_);
  t_->Branch("rechit_y", &rechit_y_);
  t_->Branch("rechit_z", &rechit_z_);
  t_->Branch("rechit_time", &rechit_time_);
  t_->Branch("rechit_thickness", &rechit_thickness_);
  t_->Branch("rechit_layer", &rechit_layer_);
  t_->Branch("rechit_wafer_u", &rechit_wafer_u_);
  t_->Branch("rechit_wafer_v", &rechit_wafer_v_);
  t_->Branch("rechit_cell_u", &rechit_cell_u_);
  t_->Branch("rechit_cell_v", &rechit_cell_v_);
  t_->Branch("rechit_detid", &rechit_detid_);
  t_->Branch("rechit_isHalf", &rechit_isHalf_);
  t_->Branch("rechit_flags", &rechit_flags_);
  t_->Branch("rechit_cluster2d", &rechit_cluster2d_);
  t_->Branch("rechit_radius", &rechit_radius_);

  ////////////////////
  // simHits
  // associated to layer clusters
  t_->Branch("simhit_eta", &simhit_eta_);
  t_->Branch("simhit_phi", &simhit_phi_);
  t_->Branch("simhit_energy", &simhit_energy_);
  t_->Branch("simhit_x", &simhit_x_);
  t_->Branch("simhit_y", &simhit_y_);
  t_->Branch("simhit_z", &simhit_z_);
  t_->Branch("simhit_layer", &simhit_layer_);
  t_->Branch("simhit_wafer_u", &simhit_wafer_u_);
  t_->Branch("simhit_wafer_v", &simhit_wafer_v_);
  t_->Branch("simhit_cell_u", &simhit_cell_u_);
  t_->Branch("simhit_cell_v", &simhit_cell_v_);
  t_->Branch("simhit_detid", &simhit_detid_);
  t_->Branch("simhit_isHalf", &simhit_isHalf_);
  t_->Branch("simhit_flags", &simhit_flags_);

  ////////////////////
  // high purity trackstatep
  //
  t_->Branch("track_eta", &track_eta_);
  t_->Branch("track_phi", &track_phi_);
  t_->Branch("track_pt", &track_pt_);
  t_->Branch("track_energy", &track_energy_);
  t_->Branch("track_charge", &track_charge_);
  t_->Branch("track_posx", &track_posx_);
  t_->Branch("track_posy", &track_posy_);
  t_->Branch("track_posz", &track_posz_);

  ////////////////////
  // gun particles
  //
  if (storeGunParticles_) {
    t_->Branch("gunparticle_id", &gunparticle_id_);
    t_->Branch("gunparticle_energy", &gunparticle_energy_);
    t_->Branch("gunparticle_pt", &gunparticle_pt_);
    t_->Branch("gunparticle_eta", &gunparticle_eta_);
    t_->Branch("gunparticle_phi", &gunparticle_phi_);
  }

}
HGCalAnalysis::~HGCalAnalysis() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//
void HGCalAnalysis::clearVariables() {
  ev_run_ = 0;
  ev_lumi_ = 0;
  ev_event_ = 0;
  vtx_x_ = 0;
  vtx_y_ = 0;
  vtx_z_ = 0;

  ////////////////////
  // GenParticles
  //
  genpart_eta_.clear();
  genpart_phi_.clear();
  genpart_pt_.clear();
  genpart_energy_.clear();
  genpart_dvx_.clear();
  genpart_dvy_.clear();
  genpart_dvz_.clear();
  genpart_ovx_.clear();
  genpart_ovy_.clear();
  genpart_ovz_.clear();
  genpart_exx_.clear();
  genpart_exy_.clear();
  genpart_mother_.clear();
  genpart_exphi_.clear();
  genpart_exeta_.clear();
  genpart_fbrem_.clear();
  genpart_pid_.clear();
  genpart_gen_.clear();
  genpart_reachedEE_.clear();
  genpart_fromBeamPipe_.clear();
  genpart_posx_.clear();
  genpart_posy_.clear();
  genpart_posz_.clear();

  ////////////////////
  // reco::GenParticles
  //
  gen_eta_.clear();
  gen_phi_.clear();
  gen_pt_.clear();
  gen_energy_.clear();
  gen_charge_.clear();
  gen_pdgid_.clear();
  gen_status_.clear();
  gen_daughters_.clear();

  ////////////////////
  // SimHits
  // associated to layer clusters
  simhit_eta_.clear();
  simhit_phi_.clear();
  simhit_energy_.clear();
  simhit_x_.clear();
  simhit_y_.clear();
  simhit_z_.clear();
  simhit_layer_.clear();
  simhit_wafer_u_.clear();
  simhit_wafer_v_.clear();
  simhit_cell_u_.clear();
  simhit_cell_v_.clear();
  simhit_detid_.clear();
  simhit_isHalf_.clear();
  simhit_flags_.clear();

  ////////////////////
  // RecHits
  // associated to layer clusters
  rechit_eta_.clear();
  rechit_phi_.clear();
  rechit_pt_.clear();
  rechit_energy_.clear();
  rechit_x_.clear();
  rechit_y_.clear();
  rechit_z_.clear();
  rechit_time_.clear();
  rechit_thickness_.clear();
  rechit_layer_.clear();
  rechit_wafer_u_.clear();
  rechit_wafer_v_.clear();
  rechit_cell_u_.clear();
  rechit_cell_v_.clear();
  rechit_detid_.clear();
  rechit_isHalf_.clear();
  rechit_flags_.clear();
  rechit_cluster2d_.clear();
  rechit_radius_.clear();
  detIdToRecHitIndexMap_.clear();

  ////////////////////
  // high purity tracks
  //
  track_eta_.clear();
  track_phi_.clear();
  track_pt_.clear();
  track_energy_.clear();
  track_charge_.clear();
  track_posx_.clear();
  track_posy_.clear();
  track_posz_.clear();

  ////////////////////
  // gun particles
  //
  gunparticle_id_.clear();
  gunparticle_energy_.clear();
  gunparticle_pt_.clear();
  gunparticle_eta_.clear();
  gunparticle_phi_.clear();
}

void HGCalAnalysis::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  using namespace edm;

  clearVariables();

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleFH;
  Handle<HGCRecHitCollection> recHitHandleBH;
  Handle<HGCRecHitCollection> recHitHandleNose;
  Handle<HFRecHitCollection>  recHitHandleHF;

  Handle<std::vector<PCaloHit>> simHitHandleEE;
  Handle<std::vector<PCaloHit>> simHitHandleFH;
  Handle<std::vector<PCaloHit>> simHitHandleBH;

  Handle<edm::HepMCProduct> hevH;

  iEvent.getByToken(hev_, hevH);

  Handle<std::vector<reco::Track>> trackHandle;

  iEvent.getByToken(tracks_, trackHandle);
  const std::vector<reco::Track> &tracks = *trackHandle;

  HepMC::GenVertex *primaryVertex = *(hevH)->GetEvent()->vertices_begin();
  float vx_ = primaryVertex->position().x() / 10.;  // to put in official units
  float vy_ = primaryVertex->position().y() / 10.;
  vz_ = primaryVertex->position().z() / 10.;
  Point sim_pv(vx_, vy_, vz_);

  // fill the gunparticles per vertex
  if (storeGunParticles_) {
    HepMC::GenEvent::vertex_const_iterator vertex_it;
    for (vertex_it = hevH->GetEvent()->vertices_begin();
         vertex_it != hevH->GetEvent()->vertices_end(); vertex_it++) {
      std::vector<int> gunparticle_id;
      std::vector<float> gunparticle_energy;
      std::vector<float> gunparticle_pt;
      std::vector<float> gunparticle_eta;
      std::vector<float> gunparticle_phi;

      HepMC::GenVertex::particles_out_const_iterator particle_it;
      for (particle_it = (*vertex_it)->particles_out_const_begin();
           particle_it != (*vertex_it)->particles_out_const_end(); particle_it++) {
        gunparticle_id.push_back((*particle_it)->pdg_id());
        gunparticle_energy.push_back((*particle_it)->momentum().e());
        gunparticle_pt.push_back((*particle_it)->momentum().perp());
        gunparticle_eta.push_back((*particle_it)->momentum().eta());
        gunparticle_phi.push_back((*particle_it)->momentum().phi());
      }

      gunparticle_id_.push_back(gunparticle_id);
      gunparticle_energy_.push_back(gunparticle_energy);
      gunparticle_pt_.push_back(gunparticle_pt);
      gunparticle_eta_.push_back(gunparticle_eta);
      gunparticle_phi_.push_back(gunparticle_phi);
    }
  }

  HGCal_helpers::simpleTrackPropagator toHGCalPropagator(aField_);
  toHGCalPropagator.setPropagationTargetZ(layerPositions_[0]);
  std::vector<FSimTrack *> allselectedgentracks;
  unsigned int npart = mySimEvent_->nTracks();
  for (unsigned int i = 0; i < npart; ++i) {
    std::vector<float> xp, yp, zp;
    FSimTrack &myTrack(mySimEvent_->track(i));
    math::XYZTLorentzVectorD vtx(0, 0, 0, 0);

    int reachedEE = 0;  // compute the extrapolations for the particles reaching EE
                        // and for the gen particles
    double fbrem = -1;

    if (std::abs(myTrack.vertex().position().z()) >= layerPositions_[0]) continue;

    unsigned nlayers =   recHitTools_.lastLayerFH();
    if (myTrack.noEndVertex())  // || myTrack.genpartIndex()>=0)
    {
      HGCal_helpers::coordinates propcoords;
      bool reachesHGCal = toHGCalPropagator.propagate(
          myTrack.momentum(), myTrack.vertex().position(), myTrack.charge(), propcoords);
      vtx = propcoords.toVector();

      if (reachesHGCal && vtx.Rho() < hgcalOuterRadius_ && vtx.Rho() > hgcalInnerRadius_) {
        reachedEE = 2;
        double dpt = 0;

        for (int i = 0; i < myTrack.nDaughters(); ++i) dpt += myTrack.daughter(i).momentum().pt();
        if (abs(myTrack.type()) == 11) fbrem = dpt / myTrack.momentum().pt();
      } else if (reachesHGCal && vtx.Rho() > hgcalOuterRadius_)
        reachedEE = 1;

      HGCal_helpers::simpleTrackPropagator indiv_particleProp(aField_);
      for (unsigned il = 0; il < nlayers; ++il) {
        const float charge = myTrack.charge();
        indiv_particleProp.setPropagationTargetZ(layerPositions_[il]);
        HGCal_helpers::coordinates propCoords;
        indiv_particleProp.propagate(myTrack.momentum(), myTrack.vertex().position(), charge,
                                     propCoords);

        xp.push_back(propCoords.x);
        yp.push_back(propCoords.y);
        zp.push_back(propCoords.z);
      }
    } else {
      vtx = myTrack.endVertex().position();
    }
    auto orig_vtx = myTrack.vertex().position();

    allselectedgentracks.push_back(&mySimEvent_->track(i));
    // fill branches
    genpart_eta_.push_back(myTrack.momentum().eta());
    genpart_phi_.push_back(myTrack.momentum().phi());
    genpart_pt_.push_back(myTrack.momentum().pt());
    genpart_energy_.push_back(myTrack.momentum().energy());
    genpart_dvx_.push_back(vtx.x());
    genpart_dvy_.push_back(vtx.y());
    genpart_dvz_.push_back(vtx.z());

    genpart_ovx_.push_back(orig_vtx.x());
    genpart_ovy_.push_back(orig_vtx.y());
    genpart_ovz_.push_back(orig_vtx.z());

    HGCal_helpers::coordinates hitsHGCal;
    toHGCalPropagator.propagate(myTrack.momentum(), orig_vtx, myTrack.charge(), hitsHGCal);

    genpart_exphi_.push_back(hitsHGCal.phi);
    genpart_exeta_.push_back(hitsHGCal.eta);
    genpart_exx_.push_back(hitsHGCal.x);
    genpart_exy_.push_back(hitsHGCal.y);

    genpart_fbrem_.push_back(fbrem);
    genpart_pid_.push_back(myTrack.type());
    genpart_gen_.push_back(myTrack.genpartIndex());
    genpart_reachedEE_.push_back(reachedEE);
    genpart_fromBeamPipe_.push_back(true);

    genpart_posx_.push_back(xp);
    genpart_posy_.push_back(yp);
    genpart_posz_.push_back(zp);
  }

  // associate gen particles to mothers
  genpart_mother_.resize(genpart_posz_.size(), -1);
  for (size_t i = 0; i < allselectedgentracks.size(); i++) {
    const auto tracki = allselectedgentracks.at(i);

    for (size_t j = i + 1; j < allselectedgentracks.size(); j++) {
      const auto trackj = allselectedgentracks.at(j);

      if (!tracki->noMother()) {
        if (&tracki->mother() == trackj) genpart_mother_.at(i) = j;
      }
      if (!trackj->noMother()) {
        if (&trackj->mother() == tracki) genpart_mother_.at(j) = i;
      }
    }
  }

  // make a map detid-rechit
  simhitmap_.clear();
  hitmap_.clear();
  hfhitmap_.clear();

  switch (algo_) {
    case 1: {
      iEvent.getByToken(recHitsEE_, recHitHandleEE);
      iEvent.getByToken(recHitsFH_, recHitHandleFH);
      iEvent.getByToken(recHitsBH_, recHitHandleBH);
      const auto &rechitsEE = *recHitHandleEE;
      const auto &rechitsFH = *recHitHandleFH;
      const auto &rechitsBH = *recHitHandleBH;
      for (unsigned int i = 0; i < rechitsEE.size(); ++i) {
        hitmap_[rechitsEE[i].detid()] = &rechitsEE[i];
      }
      for (unsigned int i = 0; i < rechitsFH.size(); ++i) {
        hitmap_[rechitsFH[i].detid()] = &rechitsFH[i];
      }
      for (unsigned int i = 0; i < rechitsBH.size(); ++i) {
        hitmap_[rechitsBH[i].detid()] = &rechitsBH[i];
      }

      iEvent.getByToken(simHitsEE_, simHitHandleEE);
      iEvent.getByToken(simHitsFH_, simHitHandleFH);
      iEvent.getByToken(simHitsBH_, simHitHandleBH);
      const auto &simhitsEE = *simHitHandleEE;
      const auto &simhitsFH = *simHitHandleFH;
      const auto &simhitsBH = *simHitHandleBH;
      for (unsigned int i = 0; i < simhitsEE.size(); ++i) {
        simhitmap_[simhitsEE[i].id()] = &simhitsEE[i];
      }
      for (unsigned int i = 0; i < simhitsFH.size(); ++i) {
        simhitmap_[simhitsFH[i].id()] = &simhitsFH[i];
      }
      for (unsigned int i = 0; i < simhitsBH.size(); ++i) {
        simhitmap_[simhitsBH[i].id()] = &simhitsBH[i];
      }

      break;
    }
    case 2: {
      iEvent.getByToken(recHitsEE_, recHitHandleEE);
      const HGCRecHitCollection &rechitsEE = *recHitHandleEE;
      for (unsigned int i = 0; i < rechitsEE.size(); i++) {
        hitmap_[rechitsEE[i].detid()] = &rechitsEE[i];
      }

      iEvent.getByToken(simHitsEE_, simHitHandleEE);
      const auto &simhitsEE = *simHitHandleEE;
      for (unsigned int i = 0; i < simhitsEE.size(); ++i) {
        simhitmap_[simhitsEE[i].id()] = &simhitsEE[i];
      }
      break;
    }
    case 3: {
      iEvent.getByToken(recHitsFH_, recHitHandleFH);
      iEvent.getByToken(recHitsBH_, recHitHandleBH);
      const auto &rechitsFH = *recHitHandleFH;
      const auto &rechitsBH = *recHitHandleBH;
      for (unsigned int i = 0; i < rechitsFH.size(); i++) {
        hitmap_[rechitsFH[i].detid()] = &rechitsFH[i];
      }
      for (unsigned int i = 0; i < rechitsBH.size(); i++) {
        hitmap_[rechitsBH[i].detid()] = &rechitsBH[i];
      }

      iEvent.getByToken(simHitsFH_, simHitHandleFH);
      iEvent.getByToken(simHitsBH_, simHitHandleBH);
      const auto &simhitsFH = *simHitHandleFH;
      const auto &simhitsBH = *simHitHandleBH;
      for (unsigned int i = 0; i < simhitsFH.size(); ++i) {
        simhitmap_[simhitsFH[i].id()] = &simhitsFH[i];
      }
      for (unsigned int i = 0; i < simhitsBH.size(); ++i) {
        simhitmap_[simhitsBH[i].id()] = &simhitsBH[i];
      }
      break;
    }
    case 4: {
      iEvent.getByToken(recHitsNose_, recHitHandleNose);
      const auto &rechitsNose = *recHitHandleNose;
      for (unsigned int i = 0; i < rechitsNose.size(); ++i) {
        hitmap_[rechitsNose[i].detid()] = &rechitsNose[i];
      }
      iEvent.getByToken(recHitsHF_, recHitHandleHF);
      const auto &rechitsHF = *recHitHandleHF;
      for (unsigned int i = 0; i < rechitsHF.size(); ++i) {
        hfhitmap_[rechitsHF[i].detid()] = &rechitsHF[i];
      }
      break;
    }
    case 5: {
      iEvent.getByToken(recHitsHF_, recHitHandleHF);
      const auto &rechitsHF = *recHitHandleHF;
      for (unsigned int i = 0; i < rechitsHF.size(); ++i) {
        hfhitmap_[rechitsHF[i].detid()] = &rechitsHF[i];
      }
      break;
    }
    default:
      break;
  }

  rechit_index_ = 0;
  storedRecHits_.clear();
  storedSimHits_.clear();

  // Fill remaining RecHits
  if (rawRecHits_) {
    if (algo_ < 3) {
      const HGCRecHitCollection &rechitsEE = *recHitHandleEE;
      // loop over EE RecHits
      for (HGCRecHitCollection::const_iterator it_hit = rechitsEE.begin(); it_hit < rechitsEE.end();
           ++it_hit) {
        const HGCalDetId detid = it_hit->detid();
        unsigned int layer = recHitTools_.getLayerWithOffset(detid);

        if (storedRecHits_.find(detid) == storedRecHits_.end()) {
          fillRecHit(detid, -1, layer);
        }
      }
    }
    if (algo_ == 1 || algo_ == 3) {
      const HGCRecHitCollection &rechitsFH = *recHitHandleFH;
      // loop over FH RecHits
      for (HGCRecHitCollection::const_iterator it_hit = rechitsFH.begin(); it_hit < rechitsFH.end();
           ++it_hit) {
        const HGCalDetId detid = it_hit->detid();
        unsigned int layer = recHitTools_.getLayerWithOffset(detid);

        if (storedRecHits_.find(detid) == storedRecHits_.end()) {
          fillRecHit(detid, -1, layer);
        }
      }
      const HGCRecHitCollection &rechitsBH = *recHitHandleBH;
      // loop over BH RecHits
      for (HGCRecHitCollection::const_iterator it_hit = rechitsBH.begin(); it_hit < rechitsBH.end();
           ++it_hit) {
        const HGCalDetId detid = it_hit->detid();
        unsigned int layer = recHitTools_.getLayerWithOffset(detid);

        if (storedRecHits_.find(detid) == storedRecHits_.end()) {
          fillRecHit(detid, -1, layer);
        }
      }
    }
    if (algo_ == 4) {
      const HGCRecHitCollection &rechitsNose = *recHitHandleNose;

      // loop over Nose RecHits
      for (HGCRecHitCollection::const_iterator it_hit = rechitsNose.begin(); it_hit < rechitsNose.end();
           ++it_hit) {
        const HGCalDetId detid = it_hit->detid();
        unsigned int layer = recHitTools_.getLayerWithOffset(detid);
        if (storedRecHits_.find(detid) == storedRecHits_.end()) {
          fillRecHit(detid, -1, layer);
        }
      }

      const HFRecHitCollection &rechitsHF = *recHitHandleHF;
       // loop over HF RecHits
      for (HFRecHitCollection::const_iterator it_hit = rechitsHF.begin(); it_hit < rechitsHF.end();
           ++it_hit) {
        const HGCalDetId detid = it_hit->detid();
        unsigned int layer = recHitTools_.getLayerWithOffset(detid);
        if (storedRecHits_.find(detid) == storedRecHits_.end()) {
          fillRecHitHF(detid, -1, layer);
        }
      }
    }

    if (algo_ == 5) {
      const HFRecHitCollection &rechitsHF = *recHitHandleHF;
       // loop over HF RecHits
      for (HFRecHitCollection::const_iterator it_hit = rechitsHF.begin(); it_hit < rechitsHF.end();
           ++it_hit) {
        const HGCalDetId detid = it_hit->detid();
        unsigned int layer = recHitTools_.getLayerWithOffset(detid);
        if (storedRecHits_.find(detid) == storedRecHits_.end()) {
          fillRecHitHF(detid, -1, layer);
        }
      }
     }
  }

  // Fill remaining SimHits
  if (rawSimHits_) {
    if (algo_ < 3) {
      const std::vector<PCaloHit> &simhitsEE = *simHitHandleEE;
      // loop over EE SimHits
      for (std::vector<PCaloHit>::const_iterator it_hit = simhitsEE.begin(); it_hit < simhitsEE.end();
           ++it_hit) {
        const HGCalDetId detid = it_hit->id();
        unsigned int layer = recHitTools_.getLayerWithOffset(detid);

        if (storedSimHits_.find(detid) == storedSimHits_.end()) {
          fillSimHit(detid, -1, layer);
        }
      }
    }
    if (algo_ == 1 || algo_ == 3) {
      const std::vector<PCaloHit> &simhitsFH = *simHitHandleFH;
      // loop over FH SimHits
      for (std::vector<PCaloHit>::const_iterator it_hit = simhitsFH.begin(); it_hit < simhitsFH.end();
           ++it_hit) {
        const HGCalDetId detid = it_hit->id();
        unsigned int layer = recHitTools_.getLayerWithOffset(detid);

        if (storedSimHits_.find(detid) == storedSimHits_.end()) {
          fillSimHit(detid, -1, layer);
        }
      }
      const std::vector<PCaloHit> &simhitsBH = *simHitHandleBH;
      // loop over BH SimHits
      for (std::vector<PCaloHit>::const_iterator it_hit = simhitsBH.begin(); it_hit < simhitsBH.end();
           ++it_hit) {
        const HGCalDetId detid = it_hit->id();
        unsigned int layer = recHitTools_.getLayerWithOffset(detid);

        if (storedSimHits_.find(detid) == storedSimHits_.end()) {
          fillSimHit(detid, -1, layer);
        }
      }
    }
  }

  if (readGen_) {
    Handle<std::vector<reco::GenParticle>> genParticlesHandle;
    iEvent.getByToken(genParticles_, genParticlesHandle);
    for (std::vector<reco::GenParticle>::const_iterator it_p = genParticlesHandle->begin();
         it_p != genParticlesHandle->end(); ++it_p) {
      gen_eta_.push_back(it_p->eta());
      gen_phi_.push_back(it_p->phi());
      gen_pt_.push_back(it_p->pt());
      gen_energy_.push_back(it_p->energy());
      gen_charge_.push_back(it_p->charge());
      gen_pdgid_.push_back(it_p->pdgId());
      gen_status_.push_back(it_p->status());
      std::vector<int> daughters(it_p->daughterRefVector().size(), 0);
      for (unsigned j = 0; j < it_p->daughterRefVector().size(); ++j) {
        daughters[j] = static_cast<int>(it_p->daughterRefVector().at(j).key());
      }
      gen_daughters_.push_back(daughters);
    }
  }

  // loop over tracks
  // prepare for RK propagation
  defaultRKPropagator::Product prod(aField_, alongMomentum, 5.e-5);
  auto &RKProp = prod.propagator;

  for (std::vector<reco::Track>::const_iterator it_track = tracks.begin(); it_track != tracks.end();
       ++it_track) {
    if (!it_track->quality(reco::Track::highPurity)) continue;

    double energy = it_track->pt() * cosh(it_track->eta());

    // save info about reconstructed tracks propoagation to hgcal layers (ony
    // for pt>propagationPtThreshold_ tracks)
    std::vector<float> xp, yp, zp;

    if (it_track->pt() >= propagationPtThreshold_) {
      // Define error matrix
      ROOT::Math::SMatrixIdentity id;
      AlgebraicSymMatrix55 C(id);
      C *= 0.01;
      CurvilinearTrajectoryError err(C);
      typedef TrajectoryStateOnSurface TSOS;

      GlobalPoint startingPosition(it_track->vx(), it_track->vy(), it_track->vz());
      GlobalVector startingMomentum(it_track->px(), it_track->py(), it_track->pz());

      Plane::PlanePointer startingPlane =
          Plane::build(Plane::PositionType(it_track->vx(), it_track->vy(), it_track->vz()),
                       Plane::RotationType());

      TSOS startingStateP(GlobalTrajectoryParameters(startingPosition, startingMomentum,
                                                     it_track->charge(), aField_),
                          err, *startingPlane);

      for (unsigned il = 0; il < layerPositions_.size(); ++il) {
        float xp_curr = 0;
        float yp_curr = 0;
        float zp_curr = 0;

        for (int zside = -1; zside <= 1; zside += 2) {
          // clearly try both sides
          Plane::PlanePointer endPlane = Plane::build(
              Plane::PositionType(0, 0, zside * layerPositions_[il]), Plane::RotationType());
          try {
            /*
            std::cout << "Trying from " <<
            " layer " << il <<
            " starting point " << startingStateP.globalPosition() <<
            std::endl;
            */
            TSOS trackStateP = RKProp.propagate(startingStateP, *endPlane);
            if (trackStateP.isValid()) {
              xp_curr = trackStateP.globalPosition().x();
              yp_curr = trackStateP.globalPosition().y();
              zp_curr = trackStateP.globalPosition().z();

              // std::cout << "Succesfully finished Positive track propagation
              // -------------- with RK: " << trackStateP.globalPosition() <<
              // std::endl;
            }
          } catch (...) {
            std::cout << "MagVolumeOutsideValidity not properly caught!! Lost "
                         "this track "
                      << std::endl;
          }
        }
        xp.push_back(xp_curr);
        yp.push_back(yp_curr);
        zp.push_back(zp_curr);
      }  // closes loop on layers
    }    // closes conditions pt>3

    // save info in tree
    track_pt_.push_back(it_track->pt());
    track_eta_.push_back(it_track->eta());
    track_phi_.push_back(it_track->phi());
    track_energy_.push_back(energy);
    track_charge_.push_back(it_track->charge());
    track_posx_.push_back(xp);
    track_posy_.push_back(yp);
    track_posz_.push_back(zp);

  }  // end loop over tracks

  ev_event_ = iEvent.id().event();
  ev_lumi_ = iEvent.id().luminosityBlock();
  ev_run_ = iEvent.id().run();

  vtx_x_ = vx_;
  vtx_y_ = vy_;
  vtx_z_ = vz_;

  t_->Fill();
}

void HGCalAnalysis::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {
  edm::ESHandle<HepPDT::ParticleDataTable> pdt;
  es.getData(pdt);
  mySimEvent_->initializePdt(&(*pdt));

  recHitTools_.getEventSetup(es);
  retrieveLayerPositions(es, recHitTools_.lastLayerBH());

  edm::ESHandle<MagneticField> magfield;
  es.get<IdealMagneticFieldRecord>().get(magfield);

  aField_ = &(*magfield);
}

void HGCalAnalysis::endRun(edm::Run const &iEvent, edm::EventSetup const &) {}

void HGCalAnalysis::beginJob() { ; }

// ------------ method called once each job just after ending the event loop
// ------------
void HGCalAnalysis::endJob() {}

// ------------ method to be called once
// --------------------------------------------------

void HGCalAnalysis::retrieveLayerPositions(const edm::EventSetup &es, unsigned layers) {
  recHitTools_.getEventSetup(es);

  for (unsigned ilayer = 1; ilayer <= layers; ++ilayer) {
    const GlobalPoint pos = recHitTools_.getPositionLayer(ilayer);
    layerPositions_.push_back(pos.z());
  }
}

void HGCalAnalysis::fillRecHit(const DetId &detid, const float &fraction, const unsigned int &layer,
                               const int &cluster_index_) {
  //std::cout << "in fillRecHit" << std::endl;
  int flags = 0x0;
  if (fraction > 0. && fraction < 1.)
    flags = 0x1;
  else if (fraction < 0.)
    flags = 0x3;
  else if (fraction == 0.)
    flags = 0x2;
  const HGCRecHit *hit = hitmap_[detid];
  const GlobalPoint position = recHitTools_.getPosition(detid);
  std::pair<int, int> wafer;
  std::pair<int, int> cell;

  if (detid.det() == DetId::Forward || detid.det() == DetId::HGCalEE || detid.det() == DetId::HGCalHSi) {
    wafer = recHitTools_.getWafer(detid);
    cell = recHitTools_.getCell(detid);
  }
  else {
    wafer = std::pair<int, int>(std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max());
    cell = std::pair<int, int>(std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max());
  }
  const double cellThickness =
      ((detid.det() == DetId::Forward || detid.det() == DetId::HGCalEE || detid.det() == DetId::HGCalHSi) ? recHitTools_.getSiThickness(detid)
                                            : std::numeric_limits<std::float_t>::max());
  const bool isHalfCell = recHitTools_.isHalfCell(detid);
  const double eta = recHitTools_.getEta(position);
  const double phi = recHitTools_.getPhi(position);
  const double pt = recHitTools_.getPt(position, hit->energy());
  const double radius =
      ((detid.det() == DetId::Forward || detid.det() == DetId::HGCalEE || detid.det() == DetId::HGCalHSi) ? recHitTools_.getRadiusToSide(detid) : -1.);

  // fill the vectors
  rechit_eta_.push_back(eta);
  rechit_phi_.push_back(phi);
  rechit_pt_.push_back(pt);
  rechit_energy_.push_back(hit->energy());
  rechit_layer_.push_back(layer);
  rechit_wafer_u_.push_back(wafer.first);
  rechit_wafer_v_.push_back(wafer.second);
  rechit_cell_u_.push_back(cell.first);
  rechit_cell_v_.push_back(cell.second);
  rechit_detid_.push_back(detid);
  rechit_x_.push_back(position.x());
  rechit_y_.push_back(position.y());
  rechit_z_.push_back(position.z());
  rechit_time_.push_back(hit->time());
  rechit_thickness_.push_back(cellThickness);
  rechit_isHalf_.push_back(isHalfCell);
  rechit_flags_.push_back(flags);
  rechit_cluster2d_.push_back(cluster_index_);
  rechit_radius_.push_back(radius);

  storedRecHits_.insert(detid);
  detIdToRecHitIndexMap_[detid] = rechit_index_;
  ++rechit_index_;
}

void HGCalAnalysis::fillRecHitHF(const DetId &detid, const float &fraction, const unsigned int &layer,
                               const int &cluster_index_) {
  //std::cout << "in fillRecHit" << std::endl;

  // set flag to 0 for HF hits for now
  int flags = 0x0;

  const HFRecHit *hit = hfhitmap_[detid];
  const GlobalPoint position = recHitTools_.getPosition(detid);
  std::pair<int, int> wafer;
  std::pair<int, int> cell;

  wafer = std::pair<int, int>(std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max());
  cell  = std::pair<int, int>(std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max());

  const double cellThickness = std::numeric_limits<std::float_t>::max();

  const bool isHalfCell = recHitTools_.isHalfCell(detid);

  const double eta = recHitTools_.getEta(position);
  const double phi = recHitTools_.getPhi(position);
  const double pt = recHitTools_.getPt(position, hit->energy());
  const double radius = sqrt(position.x()*position.x() + position.y()*position.y());

  // fill the vectors
  rechit_eta_.push_back(eta);
  rechit_phi_.push_back(phi);
  rechit_pt_.push_back(pt);
  rechit_energy_.push_back(hit->energy());
  rechit_layer_.push_back(layer);
  rechit_wafer_u_.push_back(wafer.first);
  rechit_wafer_v_.push_back(wafer.second);
  rechit_cell_u_.push_back(cell.first);
  rechit_cell_v_.push_back(cell.second);
  rechit_detid_.push_back(detid);
  rechit_x_.push_back(position.x());
  rechit_y_.push_back(position.y());
  rechit_z_.push_back(position.z());
  rechit_time_.push_back(hit->time());
  rechit_thickness_.push_back(cellThickness);
  rechit_isHalf_.push_back(isHalfCell);
  rechit_flags_.push_back(flags);
  rechit_cluster2d_.push_back(cluster_index_);
  rechit_radius_.push_back(radius);

  storedRecHits_.insert(detid);
  detIdToRecHitIndexMap_[detid] = rechit_index_;
  ++rechit_index_;
}

void HGCalAnalysis::fillSimHit(const DetId &detid, const float &fraction, const unsigned int &layer) {
  int flags = 0x0;
  if (fraction > 0. && fraction < 1.)
    flags = 0x1;
  else if (fraction < 0.)
    flags = 0x3;
  else if (fraction == 0.)
    flags = 0x2;
  const PCaloHit *hit = simhitmap_[detid];
  const GlobalPoint position = recHitTools_.getPosition(detid);
  std::pair<int, int> wafer;
  std::pair<int, int> cell;

  if (detid.det() == DetId::Forward || detid.det() == DetId::HGCalEE || detid.det() == DetId::HGCalHSi) {
    wafer = recHitTools_.getWafer(detid);
    cell = recHitTools_.getCell(detid);
  }
  else {
    wafer = std::pair<int, int>(std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max());
    cell = std::pair<int, int>(std::numeric_limits<unsigned int>::max(), std::numeric_limits<unsigned int>::max());
  }
  const bool isHalfCell = recHitTools_.isHalfCell(detid);
  const double eta = recHitTools_.getEta(position);
  const double phi = recHitTools_.getPhi(position);

  // fill the vectors
  simhit_eta_.push_back(eta);
  simhit_phi_.push_back(phi);
  simhit_energy_.push_back(hit->energy());
  simhit_layer_.push_back(layer);
  simhit_wafer_u_.push_back(wafer.first);
  simhit_wafer_v_.push_back(wafer.second);
  simhit_cell_u_.push_back(cell.first);
  simhit_cell_v_.push_back(cell.second);
  simhit_detid_.push_back(detid);
  simhit_x_.push_back(position.x());
  simhit_y_.push_back(position.y());
  simhit_z_.push_back(position.z());
  simhit_isHalf_.push_back(isHalfCell);
  simhit_flags_.push_back(flags);

  storedRecHits_.insert(detid);
  detIdToRecHitIndexMap_[detid] = rechit_index_;
  ++rechit_index_;
}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------

void HGCalAnalysis::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

/*
   Surface::RotationType HGCalAnalysis::rotation( const GlobalVector& zDir)
   const
   {
   GlobalVector zAxis = zDir.unit();
   GlobalVector yAxis( zAxis.y(), -zAxis.x(), 0);
   GlobalVector xAxis = yAxis.cross( zAxis);
   return Surface::RotationType( xAxis, yAxis, zAxis);
   }
 */

// define this as a plug-in
DEFINE_FWK_MODULE(HGCalAnalysis);
