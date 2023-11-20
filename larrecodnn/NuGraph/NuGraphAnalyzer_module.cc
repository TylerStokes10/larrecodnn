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
#include "art_root_io/TFileService.h"
#include "TTree.h"

#include "lardataobj/AnalysisBase/MVAOutput.h"

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
  TTree *_tree;
  int _run, _subrun, _event, _id;
  float _x_filter, _MIP, _HIP, _shower, _michel, _diffuse;
};


NuGraphAnalyzer::NuGraphAnalyzer(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
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

}

void NuGraphAnalyzer::analyze(art::Event const& e)
{

  // art::Handle< vector< float > > hfilter;
  art::Handle< vector<anab::FeatureVector<1> > > hfilter;
  e.getByLabel("NuGraph","filter",hfilter);

  // art::Handle< vector< vector<float> > > hsemantic;
  art::Handle< vector<anab::FeatureVector<5> > > hsemantic;
  e.getByLabel("NuGraph","semantic",hsemantic);

  std::cout << hfilter->size() << std::endl;
  for (size_t i=0;i<hfilter->size();i++) {
    _event = e.event();
    _subrun = e.subRun();
    _run = e.run();
    _id = i;
    _x_filter = hfilter->at(i).at(0);
    // _x_filter = hfilter->at(i);
    _MIP = hsemantic->at(i).at(0);
    _HIP = hsemantic->at(i).at(1);
    _shower = hsemantic->at(i).at(2);
    _michel = hsemantic->at(i).at(3);
    _diffuse = hsemantic->at(i).at(4);
    _tree->Fill();
  }

}

DEFINE_ART_MODULE(NuGraphAnalyzer)
