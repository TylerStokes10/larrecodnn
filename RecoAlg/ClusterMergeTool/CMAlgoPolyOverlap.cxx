#ifndef CMALGOPOLYOVERLAP_CXX
#define CMALGOPOLYOVERLAP_CXX

#include "CMAlgoPolyOverlap.h"

namespace cluster {

  CMAlgoPolyOverlap::CMAlgoPolyOverlap()
  {
    // Nothing to be done in the base class
    SetDebug(true);
    SetVerbose(true);
    this->reconfigure();
  }


  void CMAlgoPolyOverlap::reconfigure(){

    //not sure what needs to be reset/reconfigured for this algo
    
  }//end reconfigure function

  
  bool CMAlgoPolyOverlap::Bool(const ClusterParamsAlg &cluster1,
			       const ClusterParamsAlg &cluster2)
  {

    //if either has < 3 sides do not merge!
    if ( (cluster1.GetParams().PolyObject.Size() < 2) or
	 (cluster2.GetParams().PolyObject.Size() < 2) ){
      return false;
    }
    if (_debug and cluster1.GetParams().N_Hits > 10 and cluster2.GetParams().N_Hits > 10) {
      std::cout << "Cluster 1:" << std::endl;
      std::cout << "\tN_Hits: " << cluster1.GetParams().N_Hits << std::endl;
      std::cout << "\tN Sides:" << cluster1.GetParams().PolyObject.Size() << std::endl;
      for (unsigned int n=0; n < cluster1.GetParams().PolyObject.Size(); n++)
	std::cout << "\t\t\t(" << cluster1.GetParams().PolyObject.Point(n).first << ", "
		  << cluster1.GetParams().PolyObject.Point(n).second << ")" << std::endl;
      std::cout << "Cluster 2:" << std::endl;
      std::cout << "\tN_Hits: " << cluster2.GetParams().N_Hits << std::endl;
      std::cout << "\tN Sides:" << cluster2.GetParams().PolyObject.Size() << std::endl;
      for (unsigned int n=0; n < cluster2.GetParams().PolyObject.Size(); n++)
	std::cout << "\t\t\t(" << cluster2.GetParams().PolyObject.Point(n).first << ", "
		  << cluster2.GetParams().PolyObject.Point(n).second << ")" << std::endl;
    }

    //if the two polygons overlap even partially
    //then return true! --> MERGE!
    if ( cluster1.GetParams().PolyObject.PolyOverlapSegments(cluster2.GetParams().PolyObject) ){
      if (_verbose) { std::cout << "Overlap...merging!" << std::endl; }
      return true;
    }
    return false;

  }


}

#endif
