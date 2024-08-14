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
  float _x_filter, _MIP, _HIP, _shower, _michel, _diffuse, _ccnue, _ccnumu, _ccnutau, _nc;
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
  _tree->Branch("x_filter", &_x_filter, "x_filter/F");
  _tree->Branch("MIP", &_MIP, "MIP/F");
  _tree->Branch("HIP", &_HIP, "HIP/F");
  _tree->Branch("shower", &_shower, "shower/F");
  _tree->Branch("michel", &_michel, "michel/F");
  _tree->Branch("diffuse", &_diffuse, "diffuse/F");
  _tree->Branch("cc_nue", &_ccnue, "cc_nue/F");
  _tree->Branch("cc_numu", &_ccnumu, "cc_numu/F");
  _tree->Branch("cc_nutau", &_ccnutau, "cc_nutau/F");
  _tree->Branch("nc", &_nc, "nc/F");
}

void NuGraphAnalyzer::analyze(art::Event const& e)
{

  art::Handle<anab::MVADescription<5>> GNNDescription;
  e.getByLabel(art::InputTag("NuGraph", "semantic"), GNNDescription);

  auto const& hitsWithScores = proxy::getCollection<std::vector<recob::Hit>>(
    e,
    GNNDescription->dataTag(), //tag of the hit collection we ran the GNN on
    proxy::withParallelData<anab::FeatureVector<1>>(art::InputTag("NuGraph", "filter")),
    proxy::withParallelData<anab::FeatureVector<5>>(art::InputTag("NuGraph", "semantic")),
    proxy::withParallelData<anab::FeatureVector<4>>(art::InputTag("NuGraph", "event")));

  std::cout << hitsWithScores.size() << std::endl;
  for (auto& h : hitsWithScores) {
    int event_info_added = 0;
    if(event_info_added == 0)
    {
       const auto& assocFilter = h.get<anab::FeatureVector<1>>();
       const auto& assocSemantic = h.get<anab::FeatureVector<5>>();
       const auto& assocEvt = h.get<anab::FeatureVector<4>>();
       _event = e.event();
       _subrun = e.subRun();
       _run = e.run();
       _id = h.index();
       _x_filter = assocFilter.at(0);
       _MIP = assocSemantic.at(GNNDescription->getIndex("MIP"));
       _HIP = assocSemantic.at(GNNDescription->getIndex("HIP"));
       _shower = assocSemantic.at(GNNDescription->getIndex("shower"));
       _michel = assocSemantic.at(GNNDescription->getIndex("michel"));
       _diffuse = assocSemantic.at(GNNDescription->getIndex("diffuse"));
       _ccnue = assocEvt.at(GNNDescription->getIndex("cc_nue"));
       _ccnumu = assocEvt.at(GNNDescription->getIndex("cc_numu"));
       _ccnutau = assocEvt.at(GNNDescription->getIndex("cc_nutau"));
       _nc = assocEvt.at(GNNDescription->getIndex("nc"));
       _tree->Fill();
       event_info_added += 1;
    }

    else
    {
       const auto& assocFilter = h.get<anab::FeatureVector<1>>();
       const auto& assocSemantic = h.get<anab::FeatureVector<5>>();
       _event = e.event();
       _subrun = e.subRun();
       _run = e.run();
       _id = h.index();
       _x_filter = assocFilter.at(0);
       _MIP = assocSemantic.at(GNNDescription->getIndex("MIP"));
       _HIP = assocSemantic.at(GNNDescription->getIndex("HIP"));
       _shower = assocSemantic.at(GNNDescription->getIndex("shower"));
       _michel = assocSemantic.at(GNNDescription->getIndex("michel"));
       _diffuse = assocSemantic.at(GNNDescription->getIndex("diffuse"));
       _tree->Fill();
  }
 }

}

DEFINE_ART_MODULE(NuGraphAnalyzer)
