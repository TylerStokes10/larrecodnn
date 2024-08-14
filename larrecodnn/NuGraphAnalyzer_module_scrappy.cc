////////////////////////////////////////////////////////////////////////
// Class:       NuGraphAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        NuGraphAnalyzer_module.cc
//
// Generated at Mon Nov 20 13:42:17 2023 by Giuseppe Cerati using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// saving output
#include "TTree.h"
#include "art_root_io/TFileService.h"

#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/RecoBase/Hit.h"

class NuGraphAnalyzer;

using std::vector;

class NuGraphAnalyzer : public art::EDAnalyzer {
public:
  explicit NuGraphAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuGraphAnalyzer(NuGraphAnalyzer const&) = delete;
  NuGraphAnalyzer(NuGraphAnalyzer&&) = delete;
  NuGraphAnalyzer& operator=(NuGraphAnalyzer const&) = delete;
  NuGraphAnalyzer& operator=(NuGraphAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  // Declare member data here.
  TTree* _tree;
  int _run, _subrun, _event, _id;
  //float _x_filter;
  float  _pion, _muon, _kaon, _hadron, _shower, _michel, _diffuse;
 // float  _nu, _pdk;
};

NuGraphAnalyzer::NuGraphAnalyzer(fhicl::ParameterSet const& p) : EDAnalyzer{p} // ,
// More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("NuGraphOutput", "NuGraphOutput");
  _tree->Branch("run", &_run, "run/I");
  _tree->Branch("subrun", &_subrun, "subrun/I");
  _tree->Branch("event", &_event, "event/I");
  _tree->Branch("id", &_id, "id/I");
 // _tree->Branch("x_filter", &_x_filter, "x_filter/F");
  _tree->Branch("pion", &_pion, "pion/F");
  _tree->Branch("muon", &_muon, "muon/F");
  _tree->Branch("kaon", &_kaon, "kaon/F");
  _tree->Branch("hadron", &_hadron, "hadron/F");
  _tree->Branch("shower", &_shower, "shower/F");
  _tree->Branch("michel", &_michel, "michel/F");
  _tree->Branch("diffuse", &_diffuse, "diffuse/F");


  _tree->Branch("nu", &_nu, "nu/F");
  _tree->Branch("pdk", &_pdk, "pdk/F");
 }
/*void NuGraphAnalyzer::analyze(art::Event const& e)
{
  // Retrieve the description for accessing indices
  art::Handle<anab::MVADescription<2>> GNNDescription;
  e.getByLabel(art::InputTag("NuGraph", "event"), GNNDescription);

  if (!GNNDescription.isValid()) {
    mf::LogError("NuGraphAnalyzer") << "Failed to retrieve GNNDescription.";
    return;
  }

  // Retrieve hits with associated event-level data
  auto const& hitsWithScores = proxy::getCollection<std::vector<recob::Hit>>(
    e,
    GNNDescription->dataTag(), // Assume this correctly specifies the hit collection
    //proxy::withParallelData<anab::FeatureVector<1>>(art::InputTag("NuGraph", "filter")),
    proxy::withParallelData<anab::FeatureVector<2>>(art::InputTag("NuGraph", "event"))
  );

  std::cout << "Number of hits: " << hitsWithScores.size() << std::endl;

  for (auto& h : hitsWithScores) {
    //const auto& assocFilter = h.get<anab::FeatureVector<1>>();
    const auto& assocEvt = h.get<anab::FeatureVector<2>>();

    _event = e.id().event();
    _subrun = e.id().subRun();
    _run = e.run();
    _id = h.index();
//    _x_filter = assocFilter.at(0);
    _nu = assocEvt.at(GNNDescription->getIndex("nu"));
    _pdk = assocEvt.at(GNNDescription->getIndex("pdk"));

    _tree->Fill();
  }
}
*/
#include <vector>  // Ensure this is included
/*
void NuGraphAnalyzer::analyze(art::Event const& e) {
    art::Handle<std::vector<anab::FeatureVector<2>>> eventDataHandle;
    e.getByLabel(art::InputTag("NuGraph", "event"), eventDataHandle);

    if (!eventDataHandle.isValid() || eventDataHandle->empty()) {
        mf::LogError("NuGraphAnalyzer") << "Failed to retrieve or empty event data.";
        return;
    }

    // Assuming only one entry per event
    const auto& eventData = eventDataHandle->at(0);

    // Directly access indices if MVADescription is not providing labels (assumed fixed indices)
    _nu = eventData.at(0);  // Assuming index 0 for 'nu'
    _pdk = eventData.at(1);  // Assuming index 1 for 'pdk'

    _event = e.id().event();
    _subrun = e.id().subRun();
    _run = e.run();

    std::cout << "Filling Tree with nu: " << _nu << ", pdk: " << _pdk << std::endl;
    _tree->Fill();
}
*/

void NuGraphAnalyzer::analyze(art::Event const& e)
{

  art::Handle<anab::MVADescription<7>> GNNDescription;
  e.getByLabel(art::InputTag("NuGraph", "semantic"), GNNDescription);

  auto const& hitsWithScores = proxy::getCollection<std::vector<recob::Hit>>(
    e,
    GNNDescription->dataTag(), //tag of the hit collection we ran the GNN on
   // proxy::withParallelData<anab::FeatureVector<1>>(art::InputTag("NuGraph", "filter")),
    proxy::withParallelData<anab::FeatureVector<7>>(art::InputTag("NuGraph", "semantic")),
    proxy::withParallelData<anab::FeatureVector<2>>(art::InputTag("NuGraph", "event")));

  std::cout << hitsWithScores.size() << std::endl;
  for (auto& h : hitsWithScores) {
    int event_info_added = 0;
    if(event_info_added == 0)
    {
     //  const auto& assocFilter = h.get<anab::FeatureVector<1>>();
       const auto& assocSemantic = h.get<anab::FeatureVector<7>>();
    //   const auto& assocEvt = h.get<anab::FeatureVector<2>>();
       _event = e.event();
       _subrun = e.subRun();
       _run = e.run();
       _id = h.index();
     //  _x_filter = assocFilter.at(0);
       _pion = assocSemantic.at(GNNDescription->getIndex("pion"));
       _muon = assocSemantic.at(GNNDescription->getIndex("muon"));
       _kaon = assocSemantic.at(GNNDescription->getIndex("kaon"));
       _hadron = assocSemantic.at(GNNDescription->getIndex("hadron"));
       _shower = assocSemantic.at(GNNDescription->getIndex("shower"));
       _michel = assocSemantic.at(GNNDescription->getIndex("michel"));
       _diffuse = assocSemantic.at(GNNDescription->getIndex("diffuse"));
     //  _nu = assocEvt.at(GNNDescription->getIndex("nu"));
      // _pdk = assocEvt.at(GNNDescription->getIndex("pdk"));
       _tree->Fill();
       event_info_added += 1;
    }

    else
    {
     //  const auto& assocFilter = h.get<anab::FeatureVector<1>>();
       const auto& assocSemantic = h.get<anab::FeatureVector<7>>();
       _event = e.event();
       _subrun = e.subRun();
       _run = e.run();
       _id = h.index();
     //  _x_filter = assocFilter.at(0);
       _pion = assocSemantic.at(GNNDescription->getIndex("pion"));
       _muon = assocSemantic.at(GNNDescription->getIndex("muon"));
       _kaon = assocSemantic.at(GNNDescription->getIndex("kaon"));
       _hadron = assocSemantic.at(GNNDescription->getIndex("hadron"));
       _shower = assocSemantic.at(GNNDescription->getIndex("shower"));
       _michel = assocSemantic.at(GNNDescription->getIndex("michel"));
       _diffuse = assocSemantic.at(GNNDescription->getIndex("diffuse"));
 
       _tree->Fill();
  }
 }

}

DEFINE_ART_MODULE(NuGraphAnalyzer)
