/**
 * @file   LineCluster_module.cc
 * @brief  Cluster finder using cluster crawler algorithm
 * @author Bruce Baller (bballer@fnal.gov)
 * 
 * Generated at Fri Jun  7 09:44:09 2013 by Bruce Baller using artmod 
 * from cetpkgsupport v1_02_00.
 */


// C/C++ standard libraries
#include <string>
#include <utility> // std::unique_ptr<>

// Framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Utilities/InputTag.h"

//LArSoft includes
#include "RecoAlg/CCHitFinderAlg.h"
#include "RecoAlg/ClusterCrawlerAlg.h"

// ... more includes in the implementation section

namespace cluster {
  /**
   * @brief Produces clusters by ClusterCrawler algorithm
   * 
   * Configuration parameters
   * -------------------------
   * 
   * - *HitFinderModuleLabel* (InputTag, mandatory): label of the hits to be
   *   used as input (usually the label of the producing module is enough)
   * - *ClusterCrawlerAlg* (parameter set, mandatory): full configuration for
   *   ClusterCrawlerAlg algorithm
   * 
   */
  class LineCluster: public art::EDProducer {
    
    public:
      explicit LineCluster(fhicl::ParameterSet const & pset);
      virtual ~LineCluster() = default;
      
      void reconfigure(fhicl::ParameterSet const & pset) override;
      void produce(art::Event & evt) override;
      
    private:
      std::unique_ptr<ClusterCrawlerAlg> fCCAlg; // define ClusterCrawlerAlg object
      
      art::InputTag fHitFinderLabel; ///< label of module producing input hits
      std::string fFilteredDataModuleLabel; ///< label of module producing filtered wires
    
  }; // class LineCluster
  
} // namespace cluster

//******************************************************************************
//*** implementation
//***

// C/C++ standard libraries
#include <vector>
#include <memory> // std::move()

// Framework libraries
#include "art/Utilities/Exception.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Assns.h"

//LArSoft includes
#include "SimpleTypesAndConstants/geo_types.h"
#include "Utilities/AssociationUtil.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Vertex.h"
#include "RecoBaseArt/HitCreator.h" // recob::HitCollectionAssociator
#include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "RecoAlg/ClusterParamsImportWrapper.h"


namespace cluster {
  
  //----------------------------------------------------------------------------
  LineCluster::LineCluster(fhicl::ParameterSet const& pset) {
    reconfigure(pset);
    
    // let HitCollectionAssociator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionAssociator::declare_products(*this);
    
    produces< std::vector<recob::Cluster> >();
    produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::EndPoint2D> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::Vertex, unsigned short> >();
    //produces< art::Assns<recob::Cluster, recob::EndPoint2D, unsigned short> >();
  } // LineCluster::LineCluster()
  
  
  //----------------------------------------------------------------------------
  void LineCluster::reconfigure(fhicl::ParameterSet const & pset)
  {
    fHitFinderLabel = pset.get<art::InputTag>("HitFinderModuleLabel");
    // TODO: set default fFilteredDataModuleLabel = ""?
    fFilteredDataModuleLabel = pset.get<std::string>("FilteredDataModuleLabel", "NA");
    
    // this trick avoids double configuration on construction
    if (fCCAlg)
      fCCAlg->reconfigure(pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg"));
    else {
      fCCAlg.reset(new ClusterCrawlerAlg
        (pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg")));
    }
  } // LineCluster::reconfigure()
  
  //----------------------------------------------------------------------------
  void LineCluster::produce(art::Event & evt)
  {
    // fetch the wires needed by CCHitFinder

    // make this accessible to ClusterCrawler_module
    art::ValidHandle< std::vector<recob::Hit>> hitVecHandle
     = evt.getValidHandle<std::vector<recob::Hit>>(fHitFinderLabel);

    // decide whether to use filtered wires
    fCCAlg->ClearFilteredWires();
    std::string isNA ("NA");
    if(isNA.compare(fFilteredDataModuleLabel) != 0) {
      art::ValidHandle< std::vector<recob::Wire>> wireVecHandle
      = evt.getValidHandle<std::vector<recob::Wire>>(fFilteredDataModuleLabel);
      // pass (filtered) recob::wires so that ClusterCrawlerAlg can use it to check
      // for noisy, dead or low-noise wires
      fCCAlg->CheckFilteredWires(*wireVecHandle);
    }

    // look for clusters in all planes
    fCCAlg->RunCrawler(*hitVecHandle);
    
    std::unique_ptr<std::vector<recob::Hit>> FinalHits
      (new std::vector<recob::Hit>(std::move(fCCAlg->YieldHits())));
    
    // shcol contains the hit collection
    // and its associations to wires and raw digits;
    // we get the association to raw digits through wire associations
    recob::HitRefinerAssociator shcol(*this, evt, fHitFinderLabel);
    std::vector<recob::Cluster> sccol;
    std::vector<recob::Vertex> sv3col;
    std::vector<recob::EndPoint2D> sv2col;

    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit> > 
        hc_assn(new art::Assns<recob::Cluster, recob::Hit>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Vertex, unsigned short>> 
        cv_assn(new art::Assns<recob::Cluster, recob::Vertex, unsigned short>);

    std::vector<ClusterCrawlerAlg::ClusterStore> const& Clusters
      = fCCAlg->GetClusters();
    
    std::vector<short> const& inClus = fCCAlg->GetinClus();

// Consistency check
      for(unsigned int icl = 0; icl < Clusters.size(); ++icl) {
        ClusterCrawlerAlg::ClusterStore const& clstr = Clusters[icl];
        if(clstr.ID < 0) continue;
        geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(clstr.CTP);
        unsigned short plane = planeID.Plane;
        for(unsigned short ii = 0; ii < clstr.tclhits.size(); ++ii) {
          unsigned int iht = clstr.tclhits[ii];
          recob::Hit const& theHit = FinalHits->at(iht);
          if(theHit.WireID().Plane != plane) {
            mf::LogError("LineCluster")<<"Cluster-hit plane mis-match "<<theHit.WireID().Plane<<" "<<plane
            <<" in cluster "<<clstr.ID<<" WT "<<clstr.BeginWir<<":"<<(int)clstr.BeginTim<<" cluster CTP "<<clstr.CTP;
            return;
          }
          if(inClus[iht] != clstr.ID) {
            mf::LogError("LineCluster") << "InClus mis-match " << inClus[iht]
            << " ID " << clstr.ID << " in cluster ID " << clstr.ID<<" cluster ProcCode "<<clstr.ProcCode;;
            return;
          }
        } // ii
      } // icl
    
    // make EndPoints (aka 2D vertices)
    std::vector<ClusterCrawlerAlg::VtxStore> const& EndPts = fCCAlg->GetEndPoints();
    art::ServiceHandle<geo::Geometry> geom;
    unsigned int vtxID = 0, end, wire;
    for(ClusterCrawlerAlg::VtxStore const& vtx2: EndPts) {
      if(vtx2.NClusters == 0) continue;
      ++vtxID;
      wire = (0.5 + vtx2.Wire);
      geo::PlaneID plID = ClusterCrawlerAlg::DecodeCTP(vtx2.CTP);
      geo::WireID wID = geo::WireID(plID.Cryostat, plID.TPC, plID.Plane, wire);
      geo::View_t view = geom->View(wID);
      sv2col.emplace_back((double)vtx2.Time,    // Time
                          wID,                  // WireID
                          0,                    // strength - not relevant
                          vtxID,                // ID
                          view,                 // View
                          0);                   // total charge - not relevant
    } // Endpoints
    // convert 2D Vertex vector to unique_ptrs
    std::unique_ptr<std::vector<recob::EndPoint2D> > v2col(new std::vector<recob::EndPoint2D>(std::move(sv2col)));

    // make 3D vertices
    std::vector<ClusterCrawlerAlg::Vtx3Store> const& Vertices = fCCAlg->GetVertices();
    double xyz[3] = {0, 0, 0};
    vtxID = 0;
    for(ClusterCrawlerAlg::Vtx3Store const& vtx3: Vertices) {
      // ignore incomplete vertices
      if(vtx3.Ptr2D[0] < 0) continue;
      if(vtx3.Ptr2D[1] < 0) continue;
      if(vtx3.Ptr2D[2] < 0) continue;
      ++vtxID;
      xyz[0] = vtx3.X;
      xyz[1] = vtx3.Y;
      xyz[2] = vtx3.Z;
      sv3col.emplace_back(xyz, vtxID);
    } // 3D vertices
    // convert Vertex vector to unique_ptrs
    std::unique_ptr<std::vector<recob::Vertex> > v3col(new std::vector<recob::Vertex>(std::move(sv3col)));
    
    // make the clusters and associations
    float sumChg, sumADC;
    unsigned int clsID = 0, nclhits;
    for(unsigned int icl = 0; icl < Clusters.size(); ++icl) {
      ClusterCrawlerAlg::ClusterStore const& clstr = Clusters[icl];
      if(clstr.ID < 0) continue;
      ++clsID;
      sumChg = 0;
      sumADC = 0;
      geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(clstr.CTP);
      unsigned short plane = planeID.Plane;
      nclhits = clstr.tclhits.size();
      std::vector<unsigned int> clsHitIndices;
      // correct the hit indices to refer to the valid hits that were just added
      for(unsigned int itt = 0; itt < nclhits; ++itt) {
        unsigned int iht = clstr.tclhits[itt];
        recob::Hit const& hit = FinalHits->at(iht);
        sumChg += hit.Integral();
        sumADC += hit.SummedADC();
      } // itt
      // get the wire, plane from a hit
      unsigned int iht = clstr.tclhits[0];
      
      geo::View_t view = FinalHits->at(iht).View();
      sccol.emplace_back(
          (float)clstr.BeginWir,  // Start wire
          0,                      // sigma start wire
          clstr.BeginTim,         // start tick
          0,                      // sigma start tick
          clstr.BeginChg,         // start charge
          clstr.BeginAng,         // start angle
          0,                      // start opening angle (0 for line-like clusters)
          (float)clstr.EndWir,    // end wire
          0,                      // sigma end wire
          clstr.EndTim,           // end tick
          0,                      // sigma end tick
          clstr.EndChg,           // end charge
          clstr.EndAng,           // end angle
          0,                      // end opening angle (0 for line-like clusters)
          sumChg,                 // integral
          0,                      // sigma integral
          sumADC,                 // summed ADC
          0,                      // sigma summed ADC
          nclhits,                // n hits
          0,                      // wires over hits
          0,                      // width (0 for line-like clusters)
          clsID,                  // ID
          view,                   // view
          planeID,                // plane
          recob::Cluster::Sentry  // sentry
          );
      // make the cluster - hit association
      if(!util::CreateAssn(
        *this, evt, *hc_assn, sccol.size()-1, clstr.tclhits.begin(), clstr.tclhits.end())
        )
      {
        throw art::Exception(art::errors::InsertFailure)
          <<"Failed to associate hit "<<iht<<" with cluster "<<icl;
      } // exception
      // make the cluster - endpoint associations
      if(clstr.BeginVtx >= 0) {
        end = 0;
        // See if this endpoint is associated with a 3D vertex
        unsigned short vtxIndex = 0;
        for(ClusterCrawlerAlg::Vtx3Store const& vtx3: Vertices) {
          // ignore incomplete vertices
          if(vtx3.Ptr2D[0] < 0) continue;
          if(vtx3.Ptr2D[1] < 0) continue;
          if(vtx3.Ptr2D[2] < 0) continue;
          if(vtx3.Ptr2D[plane] == clstr.BeginVtx) {
            if(!util::CreateAssnD(*this, evt, *cv_assn, clsID - 1, vtxIndex, end))
            {
              throw art::Exception(art::errors::InsertFailure)
                <<"Failed to associate cluster "<<icl<<" with vertex";
            } // exception
            break;
          } // vertex match
          ++vtxIndex;
        } // 3D vertices
      } // clstr.BeginVtx >= 0
      if(clstr.EndVtx >= 0) {
        end = 1;
        // See if this endpoint is associated with a 3D vertex
        unsigned short vtxIndex = 0;
        for(ClusterCrawlerAlg::Vtx3Store const& vtx3: Vertices) {
          // ignore incomplete vertices
          if(vtx3.Ptr2D[0] < 0) continue;
          if(vtx3.Ptr2D[1] < 0) continue;
          if(vtx3.Ptr2D[2] < 0) continue;
          if(vtx3.Ptr2D[plane] == clstr.EndVtx) {
            if(!util::CreateAssnD(*this, evt, *cv_assn, clsID - 1, vtxIndex, end))
            {
              throw art::Exception(art::errors::InsertFailure)
                <<"Failed to associate cluster "<<icl<<" with endpoint";
            } // exception
            break;
          } // vertex match
          ++vtxIndex;
        } // 3D vertices
      } // clstr.BeginVtx >= 0
    } // icl
    
    // convert cluster vector to unique_ptrs
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(std::move(sccol)));

    shcol.use_hits(std::move(FinalHits));
    
    // clean up
    fCCAlg->ClearResults();

    // move the hit collection and the associations into the event:
    shcol.put_into(evt);
    evt.put(std::move(ccol));
    evt.put(std::move(hc_assn));
    evt.put(std::move(v2col));
    evt.put(std::move(v3col));
    evt.put(std::move(cv_assn));

  } // LineCluster::produce()
  
  
  
  //----------------------------------------------------------------------------
  DEFINE_ART_MODULE(LineCluster)
  
} // namespace cluster

