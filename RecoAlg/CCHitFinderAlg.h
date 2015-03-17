////////////////////////////////////////////////////////////////////////
// ClusterCrawlerAlg.h
//
// ClusterCrawlerAlg class
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef CCHITFINDERALG_H
#define CCHITFINDERALG_H

// C/C++ standard libraries
#include <string>
#include <vector>

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 

// LArSoft libraries
#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/Geometry.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"


namespace hit {

  /**
   * @brief Hit findre algorithm designed to work with Cluster Crawler
   * 
   * This algorithm used to store hits in a proprietary `CCHit` data structure.
   * It has now been changed to use `recob::Hit` class directly.
   * It is possible to translate the former into the latter, with one exception,
   * as follows:
   *     
   *     // this is the original CCHit definition
   *     struct CCHit {
   *       float Charge;            // recob::Hit::Integral()
   *       float ChargeErr;         // recob::Hit::SigmaIntegral()
   *       float Amplitude;         // recob::Hit::PeakAmplitude()
   *       float AmplitudeErr;      // recob::Hit::SigmaPeakAmplitude()
   *       float Time;              // recob::Hit::PeakTime()
   *       float TimeErr;           // recob::Hit::SigmaPeakTime()
   *       float RMS;               // recob::Hit::RMS()
   *       float RMSErr;            // dropped
   *       float ChiDOF;            // recob::Hit::GoodnessOfFit()
   *       int   DOF;               // recob::Hit::DegreesOfFreedom()
   *       float ADCSum;            // recob::Hit::SummedADC()
   *       unsigned short WireNum;  // recob::Hit::WireID().Wire
   *       unsigned short numHits;  // recob::Hit::Multiplicity()
   *       unsigned int LoHitID;    // see below
   *       float LoTime;            // recob::Hit::StartTick()
   *       float HiTime;            // recob::Hit::EndTick()
   *       short InClus;            // dropped; see below
   *       geo::WireID WirID;       // recob::Hit::WireID()
   *       recob::Wire const* Wire; // dropped; see below
   *     };
   *     
   * The uncertainty on RMS has been dropped for good.
   * 
   * The `LoHitID` member used to mean the index of the first hit in the "hit
   * train" (that is the set of hits extracted from the same region of
   * interest). That is a concept that is not portable. If your hit list is
   * still the original one as produced by this algorithm, or if at least the
   * hits from the same train are stored sorted and contiguously, for a hit with
   * index `iHit`, the equivalent value of `LoHitID` is
   * `iHit - hit.LocalIndex()`.
   * 
   * There is no pointer to the wire any more in `recob::Hit`. The wire can be
   * obtained through associations, that are typically produced by the art
   * module that runs CCHitFinderAlg (e.g. `CCHitFinder`). The channel ID is
   * also directly available as `recob::Hit::Channel()`.
   */
  
  class CCHitFinderAlg {
  
  public:
    
    std::vector<recob::Hit> allhits;
    
    // struct for passing hit fitting cuts to ClusterCrawler
    struct HitCuts {
      float MinSigInd;
      float MinSigCol;
      float MinRMSInd;
      float MinRMSCol;
      float ChiSplit;
      std::vector<float> ChiNorms;
    };
    HitCuts hitcuts;

    CCHitFinderAlg(fhicl::ParameterSet const& pset);
    virtual ~CCHitFinderAlg();

    void reconfigure(fhicl::ParameterSet const& pset);

    void RunCCHitFinder(std::vector<recob::Wire> const& Wires);
    
    /// Returns (and loses) the collection of reconstructed hits
    std::vector<recob::Hit>&& YieldHits() { return std::move(allhits); }
    
  private:
    
    float fMinSigInd;     ///<Induction signal height threshold 
    float fMinSigCol;     ///<Collection signal height threshold 
    float fMinRMSInd;      ///<Initial rms for induction fit
    float fMinRMSCol;      ///<Initial rms for collection fit
    unsigned short fMaxBumps; // make a crude hit if > MaxBumps are found in the RAT
    unsigned short fMaxXtraHits; // max num of hits in Region Above Threshold
    float fChiSplit;      ///<Estimated noise error on the Signal
    float ChgNorm;     // Area norm for the wire we are working on

    std::vector<float> fChiNorms;
    std::vector<float> fTimeOffsets;
    std::vector<float> fChgNorms;

    raw::ChannelID_t theChannel;
    unsigned short theWireNum;
    unsigned short thePlane;
    float minRMS;
    float minSig;
    float chinorm;
    float timeoff;
    static constexpr float Sqrt2Pi = 2.5066;
    static constexpr float SqrtPi  = 1.7725;


//    bool prt;
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;

    // fit n Gaussians possibly with bounds setting (parmin, parmax)
    void FitNG(unsigned short nGaus, unsigned short npt, float *ticks,
       float *signl);
    // parameters, errors, lower limit, upper limits for FitNG
    std::vector<double> par;
    std::vector<double> parerr;
    std::vector<double> parmin;
    std::vector<double> parmax;
    float chidof;
    int dof;
    std::vector<unsigned short> bumps;
    
    /// exchange data about the originating wire
    class HitChannelInfo_t {
        public:
      recob::Wire const* wire;
      geo::WireID wireID;
      geo::SigType_t sigType;
      
      HitChannelInfo_t
        (recob::Wire const* w, geo::WireID wid, geo::Geometry const& geom);
    }; // HitChannelInfo_t
    
    // make a cruddy hit if fitting fails
    void MakeCrudeHit(unsigned short npt, float *ticks, float *signl);
    // store the hits
    void StoreHits(unsigned short TStart, unsigned short npt, 
      HitChannelInfo_t info, float adcsum
      );

    // study hit finding and fitting
    bool fStudyHits;
    std::vector< short > fUWireRange, fUTickRange;
    std::vector< short > fVWireRange, fVTickRange;
    std::vector< short > fWWireRange, fWTickRange;
    void StudyHits(unsigned short flag, unsigned short npt = 0,
      float *ticks = 0, float *signl = 0, unsigned short tstart = 0);
    std::vector<int> bumpCnt;
    std::vector<int> RATCnt;
    std::vector<float> bumpChi;
    std::vector<float> bumpRMS;
    std::vector<int> hitCnt;
    std::vector<float> hitRMS;
    // use to determine the slope of protons
    std::vector<float> loWire;
    std::vector<float> loTime;
    std::vector<float> hiWire;
    std::vector<float> hiTime;
    bool SelRAT; // set true if a Region Above Threshold should be studied


  }; // class CCHitFinderAlg
  
} // namespace hit

#endif // ifndef CCHITFINDERALG_H
