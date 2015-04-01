////////////////////////////////////////////////////////////////////////
// Class:       MCBTDemo
// Module Type: analyzer
// File:        MCBTDemo_module.cc
//
// Generated at Thu Jan  8 08:16:24 2015 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_08_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>
#include "Geometry/Geometry.h"
#include "Simulation/SimChannel.h"
#include "MCBase/MCTrack.h"
#include "MCBase/MCShower.h"
#include "RecoBase/Track.h"
#include "RecoBase/Hit.h"
#include "MCBTAlg.h"

class MCBTDemo : public art::EDAnalyzer {
public:
  explicit MCBTDemo(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCBTDemo(MCBTDemo const &) = delete;
  MCBTDemo(MCBTDemo &&) = delete;
  MCBTDemo & operator = (MCBTDemo const &) = delete;
  MCBTDemo & operator = (MCBTDemo &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;


private:

  // Declare member data here.

};


MCBTDemo::MCBTDemo(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void MCBTDemo::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  art::Handle<std::vector<sim::MCTrack> > mctHandle;
  e.getByLabel("mcreco",mctHandle);

  art::Handle<std::vector<sim::SimChannel> > schHandle;
  e.getByLabel("largeant",schHandle);

  art::Handle<std::vector<recob::Track> > trkHandle;
  e.getByLabel("trackkalmanhit",trkHandle);

  if(!mctHandle.isValid() || !schHandle.isValid() || !trkHandle.isValid()) return;

  // Collect G4 track ID from MCTrack whose energy loss > 100 MeV inside the detector
  std::vector<unsigned int> g4_track_id;
  for(auto const& mct : *mctHandle) {

    if(!mct.size()) continue;

    double dE = (*mct.begin()).Momentum().E() - (*mct.rbegin()).Momentum().E();
    if(dE > 100) g4_track_id.push_back(mct.TrackID());

  }

  if(g4_track_id.size()) {

    art::ServiceHandle<geo::Geometry> geo;
    btutil::MCBTAlg alg_mct(g4_track_id,*schHandle);

    auto sum_mcq_v = alg_mct.MCQSum(2);
    std::cout<<"Total charge contents on W plane:"<<std::endl;
    for(size_t i=0; i<sum_mcq_v.size()-1; ++i)
      std::cout<<" MCTrack " << i << " => " << sum_mcq_v[i] <<std::endl;
    std::cout<<" Others => "<<(*sum_mcq_v.rbegin())<<std::endl;

    // Loop over reconstructed tracks and find charge fraction
    art::FindManyP<recob::Hit> hit_coll_v(trkHandle, e, "trackkalmanhit");

    for(size_t i=0; i<trkHandle->size(); ++i) {

      const std::vector<art::Ptr<recob::Hit> > hit_coll = hit_coll_v.at(i);

      std::vector<btutil::WireRange_t> hits;

      for(auto const& h_ptr : hit_coll) {

	if(geo->ChannelToWire(h_ptr->Channel())[0].Plane != ::geo::kW) continue;

	hits.emplace_back(h_ptr->Channel(),h_ptr->StartTick(), h_ptr->EndTick());

      }

      auto mcq_v = alg_mct.MCQ(hits);

      auto mcq_frac_v = alg_mct.MCQFrac(hits);

      std::cout << "Track " << i << " "
		<< "Y plane Charge from first MCTrack: " << mcq_v[0]
		<< " ... Purity: " << mcq_frac_v[0]
		<< " ... Efficiency: " << mcq_v[0] / sum_mcq_v[0] << std::endl;
    }
  }
}

DEFINE_ART_MODULE(MCBTDemo)
