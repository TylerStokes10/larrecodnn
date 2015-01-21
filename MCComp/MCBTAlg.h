/**
 * \file MCBTAlg.h
 *
 * \ingroup MCBTAlg
 * 
 * \brief Class def header for a class MCBTAlg
 *
 * @author littlejohn, kaleko, terao
 */

/** \addtogroup MCBTAlg

    @{*/
#ifndef RECOTOOL_MCSTBACKTRACKER_H
#define RECOTOOL_MCSTBACKTRACKER_H

#include <iostream>
#include <vector>
#include <map>
#include <set>
//#include "DataFormat/simch.h"
//#include "LArUtil/TimeService.h"
#include "Simulation/SimChannel.h"
#include "Utilities/TimeService.h"
#include "Geometry/Geometry.h"
#include "MCBTAlgConstants.h"
/**
   \class MCBTAlg
   MCBTAlg is meant to back-track reco-ed hits/clusters to MCShower/MCTrack
 */

namespace btutil {

  struct WireRange_t {
    unsigned int ch;
    double start, end;
    WireRange_t() {
      ch    = std::numeric_limits<unsigned int>::max();
      start = end = std::numeric_limits<double>::max();
    }
    WireRange_t(unsigned int c,double s, double e)
    { ch = c; start = s; end = e; }		
  };

  typedef std::vector<double> edep_info_t; // vector of energy deposition
  
  typedef std::vector<edep_info_t> ch_info_t; // vector of time (index) for each edep (value)
  
  class MCBTAlg {
    
  public:
    
    MCBTAlg(const std::vector<unsigned int>& g4_trackid_v,
	    const std::vector<sim::SimChannel>& simch_v);

    /**
       Returns MC charge sum per MCX for a specified plane
     */
    const std::vector<double>& MCQSum(const size_t plane_id) const;

    /**
       Relate Hit => MCShower/MCTrack (called MCX). 
       Returns a vector of double w/ length = # of relevant MCX + 1.
       Each entry is # drifted electrons from each relevant MCX.
       The last element contains a sum of drifted electrons that do not belong
       to any of relevant MCX.
    */      
    std::vector<double> MCQ(const WireRange_t& hit) const;

    /**
       Relate Hit => MCX. 
       Returns a vector of double w/ length = # of relevant MCXs + 1.
       Each entry is a fraction of # drifted electrons within the specified time
       range from each relevant MCX. The last element contains a sum of drifted 
       electrons that do not belong to any of relevant MCX.
    */      
    std::vector<double> MCQFrac(const WireRange_t& hit) const;

    /**
       Relate Cluster => MCX. 
       Returns a vector of double w/ length = # of relevant MCXs + 1.
       Each entry is # drifted electrons from each relevant MCX.
       The last element contains a sum of drifted electrons that do not belong
       to any of relevant MCX.
    */      
    std::vector<double> MCQ(const std::vector<WireRange_t>& hit_v) const;

    /**
       Relate Cluster => MCX. 
       Returns a vector of double w/ length = # of relevant MCXs + 1.
       Each entry is a fraction of # drifted electrons within the specified time
       range from each relevant MCX. The last element contains a sum of drifted 
       electrons that do not belong to any of relevant MCX.
    */      
    std::vector<double> MCQFrac(const std::vector<WireRange_t>& hit_v) const;
      
    size_t Index(const unsigned int g4_track_id) const;
    
  protected:
      
    void Register(unsigned int g4_track_id);

    std::vector<ch_info_t> _event_info;
    std::vector<size_t> _trkid_to_index;
    std::vector<std::vector<double> > _sum_mcq;
    size_t _num_parts;
  };
}
#endif
/** @} */ // end of doxygen group 

