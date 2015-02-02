////////////////////////////////////////////////////////////////////////
// Class:       RFFHitFinder
// Module Type: producer
// File:        RFFHitFinder_module.cc
//
// Generated at Fri Jan 30 17:27:31 2015 by Wesley Ketchum using artmod
// from cetpkgsupport v1_08_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

namespace hit{
  class RFFHitFinder;
}

namespace hit{

  class RFFHitFinder : public art::EDProducer {
  public:
    explicit RFFHitFinder(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    RFFHitFinder(RFFHitFinder const &) = delete;
    RFFHitFinder(RFFHitFinder &&) = delete;
    RFFHitFinder & operator = (RFFHitFinder const &) = delete;
    RFFHitFinder & operator = (RFFHitFinder &&) = delete;
    
    // Required functions.
    void produce(art::Event & e) override;
    
    // Selected optional functions.
    void reconfigure(fhicl::ParameterSet const & p) override;
    
  private:
    
    // Declare member data here.
    
  };
  

  RFFHitFinder::RFFHitFinder(fhicl::ParameterSet const & p)
  // :
  // Initialize member data here.
  {
    // Call appropriate produces<>() functions here.
  }
  
  void RFFHitFinder::produce(art::Event & e)
  {
    // Implementation of required member function here.
  }
  
  void RFFHitFinder::reconfigure(fhicl::ParameterSet const & p)
  {
    // Implementation of optional member function here.
  }

}

DEFINE_ART_MODULE(hit::RFFHitFinder)
