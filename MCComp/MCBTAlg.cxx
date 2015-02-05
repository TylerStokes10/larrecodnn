#ifndef MCBTALG_CXX
#define MCBTALG_CXX

#include "MCBTAlg.h"

namespace btutil {

  MCBTAlg::MCBTAlg(const std::vector<unsigned int>& g4_trackid_v,
		   const std::vector<sim::SimChannel>& simch_v)
  {
    Reset(g4_trackid_v,simch_v);
  }

  void Reset(const std::vector<unsigned int>& g4_trackid_v,
	     const std::vector<sim::SimChannel>& simch_v)
  {
    _num_parts = 0;
    _sum_mcq.clear();
    _trkid_to_index.clear();
    _event_info.clear();
    // 
    for(auto const& id : g4_trackid_v)
      Register(id);

    art::ServiceHandle<geo::Geometry> geo;
    _num_parts = g4_trackid_v.size() + 1;
    _sum_mcq.resize(geo->Nplanes(),std::vector<double>(_num_parts,0));

    for(auto const& sch : simch_v) {
      
      auto const ch = sch.Channel();
      if(_event_info.size() <= ch) _event_info.resize(ch+1);
      
      auto& ch_info = _event_info[ch];

      size_t plane = geo->ChannelToWire(ch)[0].Plane;

      for(auto const& time_ide : sch.TDCIDEMap()) {
	
	auto const& time  = time_ide.first;
	auto const& ide_v = time_ide.second;
	
	if(ch_info.size() <= time) ch_info.resize(time+1);
	
	auto& edep_info = ch_info[time];
	
	edep_info.resize(_num_parts,0);
	
	for(auto const& ide : ide_v) {
	  
	  size_t index = kINVALID_INDEX;
	  if(ide.trackID < (int)(_trkid_to_index.size())){
	    index = _trkid_to_index[ide.trackID];
	  }
	  if(_num_parts <= index) {
	    (*edep_info.rbegin()) += ide.numElectrons;
	    (*(_sum_mcq[plane]).rbegin()) += ide.numElectrons;
	  }
	  else {
	    edep_info[index] += ide.numElectrons;
	    _sum_mcq[plane][index] += ide.numElectrons;
	  }
	}
      }
    }
  }

  const std::vector<double>& MCBTAlg::MCQSum(const size_t plane_id) const
  {
    if(plane_id > _sum_mcq.size())
      throw MCBTException(Form("Invalid plane requested: %zu",plane_id));
    return _sum_mcq[plane_id];
  }

  std::vector<double> MCBTAlg::MCQ(const WireRange_t& hit) const
  {
    std::vector<double> res(_num_parts,0);
    
    if(_event_info.size() <= hit.ch) return res;
    
    auto const& ch_info = _event_info[hit.ch];

    art::ServiceHandle<util::TimeService> ts;

    size_t start = (size_t)(ts->TPCTick2TDC(hit.start));
    size_t end   = (size_t)(ts->TPCTick2TDC(hit.end))+1;
    for(size_t tick=start; tick<end && tick<ch_info.size(); ++tick) {
      
      auto const& edep_info = ch_info[tick];
      
      if(!edep_info.size()) continue;
      
      for(size_t part_index = 0; part_index<_num_parts; ++part_index)
	
	res[part_index] += edep_info[part_index];
      
    }
    return res;
  }


  std::vector<double> MCBTAlg::MCQFrac(const WireRange_t& hit) const
  { 
    auto res = MCQ(hit);
    if(!res.size()) return res;
    
    double sum = 0;
    for(auto const& v : res) sum += v;
    for(size_t i=0; i<(res.size()-1); ++i)
      res[i] /= (sum - (*res.rbegin()));
    (*res.rbegin()) /= sum;
    return res;
  }

  std::vector<double> MCBTAlg::MCQ(const std::vector<WireRange_t>& hit_v) const
  {
    std::vector<double> res(_num_parts,0);
    for(auto const& h : hit_v) {
      auto tmp_res = MCQ(h);
      for(size_t i=0; i<res.size(); ++i) res[i] += tmp_res[i];
    }
    return res;
  }

  std::vector<double> MCBTAlg::MCQFrac(const std::vector<WireRange_t>& hit_v) const
  {
    auto res = MCQ(hit_v);
    if(!res.size()) return res;
    
    double sum = 0;
    for(auto const& v : res) sum += v;
    for(size_t i=0; i<(res.size()-1); ++i)
      res[i] /= (sum - (*res.rbegin()));
    (*res.rbegin()) /= sum;
    return res;
  }

  size_t MCBTAlg::Index(const unsigned int g4_track_id) const
  { 
    if(g4_track_id >= _trkid_to_index.size()) return kINVALID_INDEX;
    return _trkid_to_index[g4_track_id];
  }

  void MCBTAlg::Register(unsigned int g4_track_id)
  { 
    if(_trkid_to_index.size() <= g4_track_id)
      _trkid_to_index.resize(g4_track_id+1,kINVALID_INDEX);
    
    if(_trkid_to_index[g4_track_id] == kINVALID_INDEX){
      _trkid_to_index[g4_track_id] = _num_parts;
      ++_num_parts;
    }
    std::cout<<"Registered: "<<g4_track_id<<" => "<<_trkid_to_index[g4_track_id]<<std::endl;
    return;
  }

}

#endif  
