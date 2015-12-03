////////////////////////////////////////////////////////////////////////
// ClusterCrawlerAlg.h
//
// ClusterCrawlerAlg class
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef CLUSTERCRAWLERALG_H
#define CLUSTERCRAWLERALG_H


// C/C++ standard libraries
#include <array>
#include <map>
#include <vector>
#include <memory> // std::move()
#include <utility> // std::pair<>

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 

// LArSoft libraries
#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "RecoAlg/CCHitFinderAlg.h"
#include "RecoAlg/LinFitAlg.h"


namespace cluster {
  
  class ClusterCrawlerAlg {
    public:
    
    // some functions to handle the CTP_t type
    typedef unsigned int CTP_t;
    
    static constexpr unsigned int CTPpad = 1000; // alignment for CTP sub-items
    static CTP_t EncodeCTP
      (unsigned int cryo, unsigned int tpc, unsigned int plane)
      { return cryo * CTPpad*CTPpad + tpc * CTPpad + plane; }
    static CTP_t EncodeCTP(const geo::PlaneID& planeID)
      { return EncodeCTP(planeID.Cryostat, planeID.TPC, planeID.Plane); }
    static geo::PlaneID DecodeCTP(CTP_t CTP)
      { return { CTP / (CTPpad*CTPpad), CTP / CTPpad % CTPpad, CTP % CTPpad }; }
    
    /// @{
    /// @name Data structures for the reconstruction results
    
    /// struct of temporary clusters
    struct ClusterStore {
      short ID;         // Cluster ID. ID < 0 = abandoned cluster
      short ProcCode;   // Processor code for debugging
      short StopCode;   // code for the reason for stopping cluster tracking
      CTP_t CTP;        // Cryostat/TPC/Plane code
      float BeginSlp;   // beginning slope (= DS end = high wire number)
      float BeginSlpErr; // error
      float BeginAng;
      unsigned short BeginWir;   // begin wire
      float BeginTim;   // begin time
      float BeginChg;   // beginning average charge
			float BeginChgNear; // charge near the cluster Begin
			short BeginVtx; 	// ID of the Begin vertex
      float EndSlp;     // end slope (= US end = low  wire number)
      float EndAng;
      float EndSlpErr;  // error
      unsigned short EndWir;     // end wire
      float EndTim;     // ending time
      float EndChg;     // ending average charge
			float EndChgNear; // charge near the cluster End
      short EndVtx;     // ID of the end vertex
      std::vector<unsigned int> tclhits; // hits on the cluster
    }; // ClusterStore

    /// struct of temporary 2D vertices (end points)
    struct VtxStore {
      float Wire;
			float WireErr;
      float Time;
			float TimeErr;
      unsigned short NClusters;  // = 0 for abandoned vertices
			float ChiDOF;
      short Topo; 			// 1 = US-US, 2 = US-DS, 3 = DS-US, 4 = DS-DS, 5 = Star,
												// 6 = hammer, 7 = vtx3clustermatch, 8 = vtx3clustersplit
      CTP_t CTP;
      bool Fixed;                 // Vertex position fixed (should not be re-fit)
    };
    
    /// struct of temporary 3D vertices
    struct Vtx3Store {
      std::array<short, 3> Ptr2D; // pointers to 2D vertices in each plane
      float X;                    // x position
      float XErr;                 // x position error
      float Y;                    // y position
      float YErr;                 // y position error
      float Z;                    // z position
      float ZErr;                 // z position error
      short Wire;                 // wire number for an incomplete 3D vertex
      unsigned short CStat;
      unsigned short TPC;
      unsigned short ProcCode;
    };
    
    ClusterCrawlerAlg(fhicl::ParameterSet const& pset);

    virtual void reconfigure(fhicl::ParameterSet const& pset);
    void RunCrawler(std::vector<recob::Hit> const& srchits);
    
    // initializes the vector of filtered wires
    void ClearFilteredWires();
    void CheckFilteredWires(std::vector<recob::Wire> const& Wires);

    /// @{
    /// @name Result retrieval
    
    std::vector<short> const& GetinClus() const {return inClus; }
    
    /// Returns (and loses) the collection of reconstructed hits
    std::vector<recob::Hit>&& YieldHits() { return std::move(fHits); }
    
    /// Returns a constant reference to the clusters found
    std::vector<ClusterStore> const& GetClusters() const { return tcl; }
    
    /// Returns a constant reference to the 2D end points found
    std::vector<VtxStore> const& GetEndPoints() const { return vtx; }
    
    /// Returns a constant reference to the 3D vertices found
    std::vector<Vtx3Store> const& GetVertices() const { return vtx3; }
    
    
    /// Deletes all the results (might saves memory)
    /// @note The current implementation typically does NOT save memory.
    void ClearResults();
    
    /// @}
        
    /// Comparison for sorting hits by wire and hit multiplet
    static bool SortByMultiplet(recob::Hit const& a, recob::Hit const& b);
    
    
////////////////////////////////////

    private:
    
    
    unsigned short fNumPass;                 ///< number of passes over the hit collection
    std::vector<unsigned short> fMaxHitsFit; ///< Max number of hits fitted
    std::vector<unsigned short> fMinHits;    ///< Min number of hits to make a cluster
    std::vector<unsigned short> fNHitsAve;   ///< number of US hits used to compute fAveChg
                                    ///< set to > 2 to do a charge fit using fNHitsAve hits
    std::vector<float> fChiCut;     ///< stop adding hits to clusters if chisq too high
    std::vector<float> fKinkChiRat; ///< Max consecutive chisq increase for the last 
                                    ///< 3 hits on the cluster
    std::vector<float> fKinkAngCut; ///< kink angle cut made after fKinkChiRat
    std::vector<float> fChgCut;     ///< charge difference cut for adding a hit to a cluster
    std::vector<unsigned short> fMaxWirSkip; ///< max number of wires that can be skipped while crawling
    std::vector<unsigned short> fMinWirAfterSkip; ///< minimum number of hits on consecutive wires
                                    ///< after skipping
    std::vector<bool> fDoMerge;     ///< try to merge clusters?
    std::vector<float> fTimeDelta;  ///< max time difference for matching
    std::vector<float> fMergeChgCut;  ///< max charge ratio for matching
    std::vector<bool> fFindVertices;    ///< run vertexing code after clustering?
    std::vector<bool> fLACrawl;    ///< Crawl Large Angle clusters on pass?
		bool fHammerCluster;					 ///< look for hammer type clusters
    
    bool fuBCode;     ///< patch in MicroBooNE-specific code
    bool fUseChannelFilter;
		
		float fMinAmp;									///< expected minimum signal

    bool fChkClusterDS;
    bool fVtxClusterSplit;
    bool fFindStarVertices;

    // global cuts and parameters 
    float fHitErrFac;   ///< hit time error = fHitErrFac * hit RMS used for cluster fit
    float fMinHitFrac;
    float fLAClusAngleCut;  ///< call Large Angle Clustering code if > 0
		unsigned short fLAClusMaxHitsFit; ///< max hits fitted on a Large Angle cluster
    float fLAClusSlopeCut;
    float fHitMergeChiCut; ///< Merge cluster hit-multiplets if the separation chisq
                             ///< is < cut. Set < 0 for no merging
    float fMergeOverlapAngCut;   ///< angle cut for merging overlapping clusters
    unsigned short fAllowNoHitWire;
		float fVertex2DCut; 	///< 2D vtx -> cluster matching cut (chisq/dof)
    float fVertex2DWireErrCut;
    float fVertex3DCut;   ///< 2D vtx -> 3D vtx matching cut (chisq/dof)

    short fDebugPlane;
    short fDebugWire;  ///< set to the Begin Wire and Hit of a cluster to print
    short fDebugHit;   ///< out detailed information while crawling
    
    // Wires that have been determined by some filter (e.g. NoiseFilter) to be good
    std::vector<geo::WireID> fFilteredWires;
 
    // these variables define the cluster used during crawling
    float clpar[3];     ///< cluster parameters for the current fit with
                        ///< origin at the US wire on the cluster (in clpar[2])
    float clparerr[2];  ///< cluster parameter errors
    float clChisq;     ///< chisq of the current fit
    float fAveChg;  ///< average charge at leading edge of cluster
    float fChgSlp;  ///< slope of the  charge vs wire
    float fAveHitWidth; ///< average width (EndTick - StartTick) of hits
    float fAveDTick;    ///< average dTick of hits added to the cluster
    
    bool prt;
    bool vtxprt;
    unsigned short NClusters;
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;
    
    std::vector<recob::Hit> fHits; ///< our version of the hits
    std::vector<short> inClus;    ///< Hit used in cluster (-1 = obsolete, 0 = free)
    std::vector< ClusterStore > tcl; ///< the clusters we are creating
    std::vector< VtxStore > vtx; ///< the endpoints we are reconstructing
    std::vector< Vtx3Store > vtx3; ///< the 3D vertices we are reconstructing
    
    trkf::LinFitAlg fLinFitAlg;

    float clBeginSlp;  ///< begin slope (= DS end = high wire number)
    float clBeginAng;
    float clBeginSlpErr; 
    unsigned short clBeginWir;  ///< begin wire
    float clBeginTim;  ///< begin time
    float clBeginChg;  ///< begin average charge
		float clBeginChgNear; ///< nearby charge
    float clEndSlp;    ///< slope at the end   (= US end = low  wire number)
    float clEndAng;
    float clEndSlpErr; 
    unsigned short clEndWir;    ///< begin wire
    float clEndTim;    ///< begin time
    float clEndChg;    ///< end average charge
		float clEndChgNear;  ///< nearby charge
    short clStopCode;     ///< code for the reason for stopping cluster tracking
                        ///< 0 = no signal on the next wire
                        ///< 1 = skipped too many occupied/dead wires
                        ///< 2 = failed the fMinWirAfterSkip cut
                        ///< 3 = ended on a kink. Fails fKinkChiRat
                        ///< 4 = failed the fChiCut cut
                        ///< 5 = cluster split by VtxClusterSplit
                        ///< 6 = stop at a vertex
												///< 7 = LA crawl stopped due to slope cut
    short clProcCode;     ///< Processor code = pass number
                        ///< +   10 ChkMerge
                        ///< +   20 ChkMerge with overlapping hits
                        ///< +  100 ChkMerge12
                        ///< +  200 ClusterFix
                        ///< +  300 LACrawlUS
                        ///< +  500 MergeOverlap
                        ///< + 1000 VtxClusterSplit
                        ///< + 2000 failed pass N cuts but passes pass N=1 cuts
                        ///< + 3000 Cluster hits merged
                        ///< + 5000 ChkClusterDS
                        ///< +10000 Vtx3ClusterSplit
    CTP_t clCTP;        ///< Cryostat/TPC/Plane code
		bool clLA;					///< using Large Angle crawling code
    
    unsigned short fFirstWire;    ///< the first wire with a hit
    unsigned int fFirstHit;     ///< first hit used
    unsigned short fLastWire;      ///< the last wire with a hit
    unsigned short cstat;         // the current cryostat
    unsigned short tpc;         // the current TPC
    unsigned short plane;         // the current plane
    unsigned short fNumWires;   // number of wires in the current plane
    unsigned short fMaxTime;    // number of time samples in the current plane
    float fScaleF;     ///< scale factor from Tick/Wire to dx/du
    
    unsigned short pass;

    // vector of pairs of first (.first) and last+1 (.second) hit on each wire
    // in the range fFirstWire to fLastWire. A value of -2 indicates that there
    // are no hits on the wire. A value of -1 indicates that the wire is dead
    std::vector< std::pair<int, int> > WireHitRange;
    
    std::vector<unsigned int> fcl2hits;  ///< vector of hits used in the cluster
    std::vector<float> chifits;   ///< fit chisq for monitoring kinks, etc
    std::vector<short> hitNear;   ///< Number of nearby
                                  ///< hits that were merged have hitnear < 0
    
    std::vector<float> chgNear; ///< charge near a cluster on each wire
		float fChgNearWindow; 		///< window (ticks) for finding nearby charge
		float fChgNearCut;				///< cut on ratio of nearby/cluster charge to
															///< to define a shower-like cluster
    std::vector<float> projDTick;

    std::string fhitsModuleLabel;
    
    // hit multiplets that have been saved before merging.
    std::vector<recob::Hit> unMergedHits;
    // RestoreUnMergedHits should be called before abandoning fcl2hits to
    // restore any hits that may have been merged in the fcl2hits vector
    // Set ntrim < 0 to check all hits
    void RestoreUnMergedClusterHits(short ntrim);
    // Restore a single merged hit
    void RestoreUnMergedHit(unsigned int theHit);
    // ClearUnMergedHits should (optionally) be called after a new
    // temporary (tcl) cluster has been created
    void ClearUnMergedHits();

    // ******** crawling routines *****************

    // Loops over wires looking for seed clusters
    void ClusterLoop();
    // Returns true if the hits on a cluster have a consistent width
    bool ClusterHitsOK(short nHitChk);
    // Finds a hit on wire kwire, adds it to the cluster and re-fits it
    void AddHit(unsigned short kwire, bool& HitOK, bool& SigOK);
    // Finds a hit on wire kwire, adds it to a LargeAngle cluster and re-fits it
    void AddLAHit(unsigned short kwire, bool& ChkCharge, bool& HitOK, bool& SigOK);
    // Fits the cluster hits in fcl2hits to a straight line
    void FitCluster();
    // Fits the charge of the cluster hits in fcl2hits
    void FitClusterChg();
     // Fits the middle of a temporary cluster it1 using hits iht to iht + nhit
    void FitClusterMid(unsigned short it1, unsigned int iht, short nhit);
    // Crawls along a trail of hits UpStream
    void CrawlUS();
    // Crawls along a trail of hits UpStream - Large Angle version
    void LACrawlUS();

    // ************** cluster merging routines *******************

    // Compares two cluster combinations to see if they should be merged
    void ChkMerge();
    // Checks merge for cluster cl2 within the bounds of cl1
    void ChkMerge12(unsigned short it1, unsigned short it2, bool& didit);
    // Merges clusters cl1 and cl2
    void DoMerge(unsigned short it1, unsigned short it2, short ProcCode);
      
    // ************** hit merging routines *******************

    // check the number of nearby hits to the ones added to the current cluster.
    // If there are too many, merge the hits and re-fit
    void ChkClusterNearbyHits(bool prt);
    // merge the hits in a multiplet into one hit
    void MergeHits(const unsigned int theHit, bool& didMerge);
    /// Resets the local index and multiplicity of all the hits in [begin;end[
    void FixMultipletLocalIndices
      (size_t begin, size_t end, short int multiplicity = -1);

    // ************** cluster finish routines *******************

    // Check the fraction of wires with hits
    void CheckClusterHitFrac(bool prt);

    // Try to extend clusters downstream
    void ChkClusterDS();

    // Try to merge overlapping clusters
    void MergeOverlap();
    
    /// Marks the cluster as obsolete and frees hits still associated with it
    void MakeClusterObsolete(unsigned short icl);
    /// Restores an obsolete cluster
    void RestoreObsoleteCluster(unsigned short icl);
    
    /// Removes obsolete hits from hits, updating the indices
    void RemoveObsoleteHits();

    // ************** 2D vertex routines *******************

    // Find 2D vertices
    void FindVertices();
    // Find 2D star topology vertices
    void FindStarVertices();
    // check a vertex (vw, fvt) made with clusters it1, and it2 against the
    // vector of existing clusters
    void ChkVertex(float fvw, float fvt, unsigned short it1, unsigned short it2, short topo);
    // try to attach a cluster to an existing vertex
    void ClusterVertex(unsigned short it2);
    // try to attach a cluster to a specified vertex
    void VertexCluster(unsigned short ivx);
    // Refine cluster ends near vertices
    void RefineVertexClusters(unsigned short ivx);
    // Split clusters that cross a vertex
    bool VtxClusterSplit();
    // returns true if a vertex is encountered while crawling
    bool CrawlVtxChk(unsigned short kwire);
    // returns true if this cluster is between a vertex and another
    // cluster that is associated with the vertex
    bool CrawlVtxChk2();
    // use a vertex constraint to start a cluster
    void VtxConstraint(unsigned short iwire, unsigned int ihit,
      unsigned short jwire, unsigned int& useHit, bool& doConstrain);
    // fit the vertex position
    void FitVtx(unsigned short iv);
    // weight and fit all vertices
    void FitAllVtx(CTP_t inCTP);

    // ************** 3D vertex routines *******************

    // match vertices between planes
    void VtxMatch(geo::TPCID const& tpcid);
    // Match clusters to endpoints using 3D vertex information
    void Vtx3ClusterMatch(geo::TPCID const& tpcid);
    // split clusters using 3D vertex information
    void Vtx3ClusterSplit(geo::TPCID const& tpcid);
		// look for a long cluster that stops at a short cluster in two views
		void FindHammerClusters();

    // ************** utility routines *******************

    // inits everything
    void CrawlInit();
    // fills the wirehitrange vector for the supplied Cryostat/TPC/Plane code
    void GetHitRange(CTP_t CTP);
    // Stores cluster information in a temporary vector
    bool TmpStore();
    // Gets a temp cluster and puts it into the working cluster variables
    void TmpGet(unsigned short it1);
    // Shortens the fcl2hits, chifits, etc vectors by the specified amount
    void FclTrimUS(unsigned short nTrim);
    // Splits a cluster into two clusters at position pos. Associates the
    // new clusters with a vertex
    bool SplitCluster(unsigned short icl, unsigned short pos, unsigned short ivx);
    // Prints cluster information to the screen
    void PrintClusters();
    // check for a signal on all wires between two points
    bool ChkSignal(unsigned short wire1, float time1, unsigned short wire2, float time2);
    // returns an angle-dependent scale factor for weighting fits, etc
    float AngleFactor(float slope);
    // calculate the kink angle between hits 0-2 and 3 - 5 on the leading edge of
    // the cluster under construction
    float EndKinkAngle();
    /// Returns true if there are no duplicates in the hit list for next cluster
    bool CheckHitDuplicates
      (std::string location, std::string marker = "") const;
    // Find the distance of closest approach between the end of a cluster and a (wire,tick) position
    float DoCA(short icl, unsigned short end, float vwire, float vtick);
    // Find the Chisq/DOF between the end of a cluster and a (wire,tick) position
    float ClusterVertexChi(short icl, unsigned short end, unsigned short ivx);
    // Find the Chisq/DOF between a point and a vertex
    float PointVertexChi(float wire, float tick, unsigned short ivx);
    
    /// Returns a pair of first and past-the-last index
    /// of all the contiguous hits belonging to the same multiplet
    std::pair<size_t, size_t> FindHitMultiplet(size_t iHit) const;

    void CheckHitClusterAssociations();
    
    /// Returns whether the two hits belong to the same multiplet
    static bool areInSameMultiplet(recob::Hit const& first_hit, recob::Hit const& second_hit);

    
  }; // class ClusterCrawlerAlg


} // namespace cluster

#endif // ifndef CLUSTERCRAWLERALG_H
