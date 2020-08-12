//////////////////////////////////////////////
//
//  Module to dump raw data for ROI training
//
//  mwang@fnal.gov
//  tjyang@fnal.gov
//
//////////////////////////////////////////////

#include <random>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// DUNETPC specific includes
//#include "dune/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"
//#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"
//#include "dune/DuneInterface/AdcTypes.h"
//#include "dune/DuneInterface/SimChannelExtractService.h"

#include "c2numpy.h"

using std::cout;
using std::endl;
using std::ofstream;
using std::string;

struct WireSigInfo {
  int pdgcode;
  std::string genlab;
  std::string procid;
  unsigned int tdcmin;
  unsigned int tdcmax;
  int numel;
  double edep;
};

namespace nnet {
  class RawWaveformDump;
}

class nnet::RawWaveformDump : public art::EDAnalyzer {

public:
  explicit RawWaveformDump(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  RawWaveformDump(RawWaveformDump const&) = delete;
  RawWaveformDump(RawWaveformDump&&) = delete;
  RawWaveformDump& operator=(RawWaveformDump const&) = delete;
  RawWaveformDump& operator=(RawWaveformDump&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  //void reconfigure(fhicl::ParameterSet const & p);

  void beginJob() override;
  void endJob() override;

private:
  std::string fDumpWaveformsFileName;

  std::string
    fSimulationProducerLabel; ///< The name of the producer that tracked simulated particles through the detector
  std::string fDigitModuleLabel; ///< module that made digits

  std::string fSelectGenLabel;
  std::string fSelectProcID;
  int fSelectPDGCode;
  std::string fPlaneToDump;
  double fMinParticleEnergyGeV;
  double fMinEnergyDepositedMeV;
  int fMinNumberOfElectrons;
  int fMaxNumberOfElectrons;
  bool fSaveSignal;
  art::ServiceHandle<geo::Geometry> fgeom;
  art::ServiceHandle<cheat::ParticleInventoryService> PIS;
  //art::ServiceHandle<SimChannelExtractService> m_pscx;
  detinfo::DetectorClocks const* fClks;

  std::default_random_engine rndm_engine;

  c2numpy_writer npywriter;
};

//-----------------------------------------------------------------------
struct genFinder {
private:
  typedef std::pair<int, std::string> track_id_to_string;
  std::vector<track_id_to_string> track_id_map;
  std::set<std::string> generator_names;
  bool isSorted = false;

public:
  void
  sort_now()
  {
    std::sort(this->track_id_map.begin(),
              this->track_id_map.end(),
              [](const auto& a, const auto& b) { return (a.first < b.first); });
    isSorted = true;
  }
  void
  add(const int& track_id, const std::string& gname)
  {
    this->track_id_map.push_back(std::make_pair(track_id, gname));
    generator_names.emplace(gname);
    isSorted = false;
  }
  bool
  has_gen(std::string gname)
  {
    return static_cast<bool>(generator_names.count(gname));
  };
  std::string
  get_gen(int tid)
  {
    if (!isSorted) { this->sort_now(); }
    return std::lower_bound(track_id_map.begin(),
                            track_id_map.end(),
                            tid,
                            [](const auto& a, const auto& b) { return (a.first < b); })
      ->second;
  };
};
std::unique_ptr<genFinder> gf(new genFinder());

//-----------------------------------------------------------------------
nnet::RawWaveformDump::RawWaveformDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fDumpWaveformsFileName(p.get<std::string>("DumpWaveformsFileName", "dumpwaveforms"))
  , fSimulationProducerLabel(p.get<std::string>("SimulationProducerLabel", "largeant"))
  , fDigitModuleLabel(p.get<std::string>("DigitModuleLabel", "daq"))
  , fSelectGenLabel(p.get<std::string>("SelectGenLabel", "ANY"))
  , fSelectProcID(p.get<std::string>("SelectProcID", "ANY"))
  , fSelectPDGCode(p.get<int>("SelectPDGCode", 0))
  , fPlaneToDump(p.get<std::string>("PlaneToDump", "U"))
  , fMinParticleEnergyGeV(p.get<double>("MinParticleEnergyGeV", 0.))
  , fMinEnergyDepositedMeV(p.get<double>("MinEnergyDepositedMeV", 0.))
  , fMinNumberOfElectrons(p.get<int>("MinNumberOfElectrons", 1000))
  , fMaxNumberOfElectrons(p.get<int>("MaxNumberOfElectrons", 100000))
  , fSaveSignal(p.get<bool>("SaveSignal", true))
  , fClks(lar::providerFrom<detinfo::DetectorClocksService>())
{
  if (std::getenv("PROCESS")) { fDumpWaveformsFileName += string(std::getenv("PROCESS")) + "-"; }
  //this->reconfigure(p);
}

//-----------------------------------------------------------------------
void
nnet::RawWaveformDump::beginJob()
{

  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  std::random_device rndm_device; // this will give us our seed
  rndm_engine.seed(rndm_device());

  c2numpy_init(&npywriter, fDumpWaveformsFileName, 50000);
  c2numpy_addcolumn(&npywriter, "evt", C2NUMPY_UINT32);
  c2numpy_addcolumn(&npywriter, "chan", C2NUMPY_UINT32);
  c2numpy_addcolumn(&npywriter, "view", (c2numpy_type)((int)C2NUMPY_STRING + 1));
  c2numpy_addcolumn(&npywriter, "ntrk", C2NUMPY_UINT16);

  for (unsigned int i = 0; i < 5; i++) {
    std::ostringstream name;

    name.str("");
    name << "tid" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_INT32);

    name.str("");
    name << "pdg" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_INT32);

    name.str("");
    name << "gen" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), (c2numpy_type)((int)C2NUMPY_STRING + 6));

    name.str("");
    name << "pid" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), (c2numpy_type)((int)C2NUMPY_STRING + 7));

    name.str("");
    name << "edp" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_FLOAT32);

    name.str("");
    name << "nel" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_UINT32);

    name.str("");
    name << "sti" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_UINT16);

    name.str("");
    name << "stf" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_UINT16);
  }

  for (unsigned int i = 0; i < detprop->ReadOutWindowSize(); i++) {
    std::ostringstream name;
    name << "tck_" << i;
    c2numpy_addcolumn(&npywriter, name.str().c_str(), C2NUMPY_INT16);
  }
}

//-----------------------------------------------------------------------
void
nnet::RawWaveformDump::endJob()
{
  c2numpy_close(&npywriter);
}

//-----------------------------------------------------------------------
void
nnet::RawWaveformDump::analyze(art::Event const& evt)
{
  cout << "Event "
       << " " << evt.id().run() << " " << evt.id().subRun() << " " << evt.id().event() << endl;

  // ... Read in the digit List object(s).
  art::Handle<std::vector<raw::RawDigit>> digitVecHandle;
  evt.getByLabel(fDigitModuleLabel, digitVecHandle);

  if (!digitVecHandle->size()) return;

  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  // ... Use the handle to get a particular (0th) element of collection.
  art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
  unsigned int dataSize = digitVec0->Samples(); //size of raw data vectors
  if (dataSize != detprop->ReadOutWindowSize()) {
    std::cout << "!!!!! Bad dataSize: " << dataSize << std::endl;
    return;
  }

  std::vector<short> rawadc(dataSize); // vector to hold uncompressed adc values later

  // ... Build a map from channel number -> rawdigitVec
  std::map<raw::ChannelID_t, art::Ptr<raw::RawDigit>> rawdigitMap;
  raw::ChannelID_t chnum = raw::InvalidChannelID; // channel number
  for (size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter) {
    art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
    chnum = digitVec->Channel();
    if (chnum == raw::InvalidChannelID) continue;
    rawdigitMap[chnum] = digitVec;
  }

  // ... Read in MC particle list
  art::Handle<std::vector<simb::MCParticle>> particleHandle;
  if (!evt.getByLabel(fSimulationProducerLabel, particleHandle)) {
    throw cet::exception("AnalysisExample")
      << " No simb::MCParticle objects in this event - "
      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  // ... Read in sim channel list
  auto simChannelHandle =
    evt.getValidHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);

  if (!simChannelHandle->size()) return;

  // ... Create a map of track IDs to generator labels
  //Get a list of generator names.
  std::vector<art::Handle<std::vector<simb::MCTruth>>> mcHandles;
  evt.getManyByType(mcHandles);
  std::vector<std::pair<int, std::string>> track_id_to_label;

  for (auto const& mcHandle : mcHandles) {
    const std::string& sModuleLabel = mcHandle.provenance()->moduleLabel();
    art::FindManyP<simb::MCParticle> findMCParts(mcHandle, evt, "largeant");
    std::vector<art::Ptr<simb::MCParticle>> mcParts = findMCParts.at(0);
    for (const art::Ptr<simb::MCParticle> ptr : mcParts) {
      int track_id = ptr->TrackId();
      gf->add(track_id, sModuleLabel);
    }
  }

  std::string dummystr6 = "none  ";
  std::string dummystr7 = "none   ";

  if (fSaveSignal) {
    // .. create a channel number to trackid-wire signal info map
    std::map<raw::ChannelID_t, std::map<int, WireSigInfo>> Ch2TrkWSInfoMap;

    // .. create a track ID to vector of channel numbers (in w/c this track deposited energy) map
    std::map<int, std::vector<raw::ChannelID_t>> Trk2ChVecMap;

    // ... Loop over simChannels
    for (auto const& channel : (*simChannelHandle)) {

      // .. get simChannel channel number
      const raw::ChannelID_t ch1 = channel.Channel();
      if (ch1 == raw::InvalidChannelID) continue;
      if (geo::PlaneGeo::ViewName(fgeom->View(ch1)) != fPlaneToDump[0]) continue;

      bool selectThisChannel = false;

      // .. create a track ID to wire signal info map
      std::map<int, WireSigInfo> Trk2WSInfoMap;

      // ... Loop over all ticks with ionization energy deposited
      auto const& timeSlices = channel.TDCIDEMap();
      for (auto const& timeSlice : timeSlices) {

        auto const& energyDeposits = timeSlice.second;
        auto const tpctime = timeSlice.first;
        unsigned int tdctick = static_cast<unsigned int>(fClks->TPCTDC2Tick(double(tpctime)));
        if (tdctick != tpctime)
          std::cout << "tpctime: " << tpctime << ", tdctick: " << tdctick << std::endl;
        if (tdctick < 0 || tdctick > (dataSize - 1)) continue;

        // ... Loop over all energy depositions in this tick
        for (auto const& energyDeposit : energyDeposits) {

          if (!energyDeposit.trackID) continue;
          int trkid = energyDeposit.trackID;
          simb::MCParticle particle = PIS->TrackIdToMotherParticle(trkid);
          //std::cout << energyDeposit.trackID << " " << trkid << " " << particle.TrackId() << std::endl;

          // .. ignore this energy deposition if incident particle energy below some threshold
          if (particle.E() < fMinParticleEnergyGeV) continue;

          int eve_id = PIS->TrackIdToEveTrackId(trkid);
          if (!eve_id) continue;
          std::string genlab = gf->get_gen(eve_id);

          if (Trk2WSInfoMap.find(trkid) == Trk2WSInfoMap.end()) {
            WireSigInfo wsinf;
            wsinf.pdgcode = particle.PdgCode();
            wsinf.genlab = genlab;
            wsinf.procid = particle.Process();
            wsinf.tdcmin = dataSize - 1;
            wsinf.tdcmax = 0;
            wsinf.edep = 0.;
            wsinf.numel = 0;
            Trk2WSInfoMap.insert(std::pair<int, WireSigInfo>(trkid, wsinf));
          }
          if (tdctick < Trk2WSInfoMap.at(trkid).tdcmin) Trk2WSInfoMap.at(trkid).tdcmin = tdctick;
          if (tdctick > Trk2WSInfoMap.at(trkid).tdcmax) Trk2WSInfoMap.at(trkid).tdcmax = tdctick;
          Trk2WSInfoMap.at(trkid).edep += energyDeposit.energy;
          Trk2WSInfoMap.at(trkid).numel += energyDeposit.numElectrons;
        }
      }

      if (!Trk2WSInfoMap.empty()) {
        for (std::pair<int, WireSigInfo> itmap : Trk2WSInfoMap) {
          if (fSelectGenLabel != "ANY") {
            if (itmap.second.genlab != fSelectGenLabel) continue;
          }
          if (fSelectProcID != "ANY") {
            if (itmap.second.procid != fSelectProcID) continue;
          }
          if (fSelectPDGCode != 0) {
            if (itmap.second.pdgcode != fSelectPDGCode) continue;
          }
          itmap.second.genlab.resize(6, ' ');
          itmap.second.procid.resize(7, ' ');
          if (itmap.second.numel >= fMinNumberOfElectrons &&
              itmap.second.edep >= fMinEnergyDepositedMeV) {
            if (fMaxNumberOfElectrons >= 0 && itmap.second.numel >= fMaxNumberOfElectrons) {
              continue;
            }
            else {
              int trkid = itmap.first;
              if (Trk2ChVecMap.find(trkid) == Trk2ChVecMap.end()) {
                std::vector<raw::ChannelID_t> chvec;
                Trk2ChVecMap.insert(std::pair<int, std::vector<raw::ChannelID_t>>(trkid, chvec));
              }
              Trk2ChVecMap.at(trkid).push_back(ch1);
              selectThisChannel = true;
            }
          }
        } // loop over Trk2WSinfoMap
        if (selectThisChannel) {
          Ch2TrkWSInfoMap.insert(
            std::pair<raw::ChannelID_t, std::map<int, WireSigInfo>>(ch1, Trk2WSInfoMap));
        }
      } // if Trk2WSInfoMap not empty

    } // loop over SimChannels

    // ... Now write out the signal waveforms for each track
    if (!Trk2ChVecMap.empty()) {
      for (auto const& ittrk : Trk2ChVecMap) {
        std::uniform_int_distribution<int> rndm_dist(0, ittrk.second.size() - 1);
        int i =
          rndm_dist(rndm_engine); // randomly select one channel with a signal from this particle
        chnum = ittrk.second[i];

        std::map<raw::ChannelID_t, std::map<int, WireSigInfo>>::iterator itchn;
        itchn = Ch2TrkWSInfoMap.find(chnum);
        if (itchn != Ch2TrkWSInfoMap.end()) {
          auto search = rawdigitMap.find(chnum);
          if (search == rawdigitMap.end()) continue;
          art::Ptr<raw::RawDigit> rawdig = (*search).second;
          raw::Uncompress(rawdig->ADCs(), rawadc, rawdig->GetPedestal(), rawdig->Compression());

          c2numpy_uint32(&npywriter, evt.id().event());
          c2numpy_uint32(&npywriter, chnum);
          c2numpy_string(&npywriter, geo::PlaneGeo::ViewName(fgeom->View(chnum)).c_str());
          c2numpy_uint16(&npywriter,
                         itchn->second.size()); // size of Trk2WSInfoMap, or number of peaks

          // .. write out info for each peak
          unsigned int icnt = 0;
          for (auto const& it : itchn->second) {
            c2numpy_int32(&npywriter, it.first);                  // trackid
            c2numpy_int32(&npywriter, it.second.pdgcode);         // pdgcode
            c2numpy_string(&npywriter, it.second.genlab.c_str()); // genlab
            c2numpy_string(&npywriter, it.second.procid.c_str()); // procid
            c2numpy_float32(&npywriter, it.second.edep);          // edepo
            c2numpy_uint32(&npywriter, it.second.numel);          // numelec
            c2numpy_uint16(&npywriter, it.second.tdcmin);         // stck1
            c2numpy_uint16(&npywriter, it.second.tdcmax);         // stc2
            icnt++;
            if (icnt == 5) break;
          }

          // .. pad with 0's if number of peaks less than 5
          for (unsigned int i = icnt; i < 5; ++i) {
            c2numpy_int32(&npywriter, 0);
            c2numpy_int32(&npywriter, 0);
            c2numpy_string(&npywriter, dummystr6.c_str());
            c2numpy_string(&npywriter, dummystr7.c_str());
            c2numpy_float32(&npywriter, 0.);
            c2numpy_uint32(&npywriter, 0);
            c2numpy_uint16(&npywriter, 0);
            c2numpy_uint16(&npywriter, 0);
          }

          // .. now write out the ADC values
          for (unsigned int itck = 0; itck < dataSize; ++itck) {
            rawadc[itck] -= rawdig->GetPedestal();
            c2numpy_int16(&npywriter, rawadc[itck]);
          }
        }
      }
    }
  }
  else {
    //save noise
    int noisechancount = 0;
    std::map<raw::ChannelID_t, bool> signalMap;
    for (auto const& channel : (*simChannelHandle)) {
      signalMap[channel.Channel()] = true;
    }

    std::vector<short> rawadc(dataSize); // vector to hold uncompressed adc values later

    for (size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter) {
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      if (signalMap[digitVec->Channel()]) continue;
      if (geo::PlaneGeo::ViewName(fgeom->View(digitVec->Channel())) != fPlaneToDump[0]) continue;
      raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->GetPedestal(), digitVec->Compression());
      c2numpy_uint32(&npywriter, evt.id().event());
      c2numpy_uint32(&npywriter, digitVec->Channel());
      c2numpy_string(&npywriter, geo::PlaneGeo::ViewName(fgeom->View(digitVec->Channel())).c_str());
      c2numpy_uint16(&npywriter, 0); //number of peaks
      for (unsigned int i = 0; i < 5; ++i) {
        c2numpy_int32(&npywriter, 0);
        c2numpy_int32(&npywriter, 0);
        c2numpy_string(&npywriter, dummystr6.c_str());
        c2numpy_string(&npywriter, dummystr7.c_str());
        c2numpy_float32(&npywriter, 0.);
        c2numpy_uint32(&npywriter, 0);
        c2numpy_uint16(&npywriter, 0);
        c2numpy_uint16(&npywriter, 0);
      }
      for (unsigned int itck = 0; itck < dataSize; ++itck) {
        rawadc[itck] -= digitVec->GetPedestal();
        c2numpy_int16(&npywriter, rawadc[itck]);
      }
      ++noisechancount;
    }
    std::cout << "Total number of noise channels " << noisechancount << std::endl;
  }
}
DEFINE_ART_MODULE(nnet::RawWaveformDump)
