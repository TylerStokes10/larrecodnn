#include <vector>
#include <iostream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// For ROOT output
#include "TTree.h"
#include "art_root_io/TFileService.h"

#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/RecoBase/Hit.h"

using std::vector;

class NuGraphAnalyzer : public art::EDAnalyzer {
public:
    explicit NuGraphAnalyzer(fhicl::ParameterSet const& p);
    void analyze(art::Event const& e) override;

private:
    TTree* _eventTree;
    TTree* _hitTree;
    int _run, _subrun, _event;
    float _nu, _pdk;  // Event-level data
    int _id;
    float _pion, _muon, _kaon, _hadron, _shower, _michel, _diffuse;  // Hit-level data
};

NuGraphAnalyzer::NuGraphAnalyzer(fhicl::ParameterSet const& p) : EDAnalyzer{p} {
    art::ServiceHandle<art::TFileService> tfs;
    _eventTree = tfs->make<TTree>("EventTree", "Event-level data");
    _eventTree->Branch("run", &_run, "run/I");
    _eventTree->Branch("subrun", &_subrun, "subrun/I");
    _eventTree->Branch("event", &_event, "event/I");
    _eventTree->Branch("nu", &_nu, "nu/F");
    _eventTree->Branch("pdk", &_pdk, "pdk/F");

    _hitTree = tfs->make<TTree>("HitTree", "Hit-level data");
    _hitTree->Branch("id", &_id, "id/I");
    _hitTree->Branch("pion", &_pion, "pion/F");
    _hitTree->Branch("muon", &_muon, "muon/F");
    _hitTree->Branch("kaon", &_kaon, "kaon/F");
    _hitTree->Branch("hadron", &_hadron, "hadron/F");
    _hitTree->Branch("shower", &_shower, "shower/F");
    _hitTree->Branch("michel", &_michel, "michel/F");
    _hitTree->Branch("diffuse", &_diffuse, "diffuse/F");
}

void NuGraphAnalyzer::analyze(art::Event const& e) {
    _event = e.id().event();
    _subrun = e.id().subRun();
    _run = e.run();

    art::Handle<std::vector<anab::FeatureVector<2>>> eventDataHandle;
    e.getByLabel(art::InputTag("NuGraph", "event"), eventDataHandle);
    if (eventDataHandle.isValid() && !eventDataHandle->empty()) {
        _nu = eventDataHandle->at(0).at(1);
        _pdk = eventDataHandle->at(0).at(0);
    } else {
        _nu = 0.0f;
        _pdk = 0.0f;
    }

    // Fill event-level data only once per event
    _eventTree->Fill();
/*
    art::Handle<anab::MVADescription<7>> GNNDescription;
    e.getByLabel(art::InputTag("NuGraph", "semantic"), GNNDescription);
    if (!GNNDescription.isValid()) {
        mf::LogError("NuGraphAnalyzer") << "Failed to retrieve GNNDescription.";
        return;
    }

    auto const& hitsWithScores = proxy::getCollection<std::vector<recob::Hit>>(
        e,
        GNNDescription->dataTag(),
        proxy::withParallelData<anab::FeatureVector<7>>(art::InputTag("NuGraph", "semantic"))
    );

    for (auto& hit : hitsWithScores) {
        _id = hit.index();
        const auto& assocSemantic = hit.get<anab::FeatureVector<7>>();
        _pion = assocSemantic.at(GNNDescription->getIndex("pion"));
        _muon = assocSemantic.at(GNNDescription->getIndex("muon"));
        _kaon = assocSemantic.at(GNNDescription->getIndex("kaon"));
        _hadron = assocSemantic.at(GNNDescription->getIndex("hadron"));
        _shower = assocSemantic.at(GNNDescription->getIndex("shower"));
        _michel = assocSemantic.at(GNNDescription->getIndex("michel"));
        _diffuse = assocSemantic.at(GNNDescription->getIndex("diffuse"));

        // Fill hit-level data per hit
        _hitTree->Fill();
    }
*/}

DEFINE_ART_MODULE(NuGraphAnalyzer)

