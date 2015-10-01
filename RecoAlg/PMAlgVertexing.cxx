////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgVertexing
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), August 2015
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "RecoAlg/PMAlgVertexing.h"

#include "RecoAlg/PMAlg/Utilities.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

pma::PMAlgVertexing::PMAlgVertexing(const fhicl::ParameterSet& pset)
{
	this->reconfigure(pset); 
}
// ------------------------------------------------------

pma::PMAlgVertexing::~PMAlgVertexing(void)
{
	cleanTracks();
}
// ------------------------------------------------------

void pma::PMAlgVertexing::reconfigure(const fhicl::ParameterSet& pset)
{
	fMinTrackLength = pset.get< double >("MinTrackLength");

	//fInputVtxDist2D = pset.get< double >("InputVtxDist2D");
	//fInputVtxDistY = pset.get< double >("InputVtxDistY");
}
// ------------------------------------------------------

void pma::PMAlgVertexing::cleanTracks(void)
{
	for (auto t : fOutTracks) delete t;
	fOutTracks.clear();

	for (auto t : fShortTracks) delete t;
	fShortTracks.clear();

	for (auto t : fEmTracks) delete t;
	fEmTracks.clear();
}
// ------------------------------------------------------

void pma::PMAlgVertexing::collectTracks(
	std::vector< pma::Track3D* >& result)
{
	for (auto t : result) delete t;
	result.clear();

	for (auto t : fOutTracks) result.push_back(t);
	fOutTracks.clear();

	for (auto t : fShortTracks) result.push_back(t);
	fShortTracks.clear();

	for (auto t : fEmTracks) result.push_back(t);
	fEmTracks.clear();
}
// ------------------------------------------------------

void pma::PMAlgVertexing::sortTracks(
	const std::vector< pma::Track3D* >& trk_input)
{
	cleanTracks();

	for (auto t : trk_input)
	{
		double l = t->Length();
		if (t->GetTag() == pma::Track3D::kEmLike)
		{
			if (l > 2 * fMinTrackLength) fOutTracks.push_back(new pma::Track3D(*t));
			else fEmTracks.push_back(new pma::Track3D(*t));
		}
		else
		{
			if (l > fMinTrackLength) fOutTracks.push_back(new pma::Track3D(*t));
			else fEmTracks.push_back(new pma::Track3D(*t));
		}
	}
	mf::LogVerbatim("pma::PMAlgVertexing") << "long tracks: " << fOutTracks.size() << std::endl;
	mf::LogVerbatim("pma::PMAlgVertexing") << "em and short tracks: " << fEmTracks.size() << std::endl;
}
// ------------------------------------------------------

std::vector< pma::VtxCandidate > pma::PMAlgVertexing::firstPassCandidates(void)
{
	std::vector< pma::VtxCandidate > candidates;
	for (size_t t = 0; t < fOutTracks.size() - 1; t++)
	{
		for (size_t u = t + 1; u < fOutTracks.size(); u++)
		{
			pma::VtxCandidate candidate;
			if (!candidate.Add(fOutTracks[t])) break; // no segments with length > thr

			// **************************** try Mse2D / or only Mse ************************************
			if (candidate.Add(fOutTracks[u]) && (sqrt(candidate.Mse()) < 1.0))
			//if (candidate.Add(fOutTracks[u]) && (sqrt(candidate.Mse()) < 2.0) && (candidate.Mse2D() < 1.0))
				candidates.push_back(candidate);
		}
	}
	return candidates;
}

std::vector< pma::VtxCandidate > pma::PMAlgVertexing::secondPassCandidates(void)
{
	std::vector< pma::VtxCandidate > candidates;
	for (size_t t = 0; t < fOutTracks.size(); t++)
		if (fOutTracks[t]->Length() > fMinTrackLength)
	{
		for (size_t u = 0; u < fEmTracks.size(); u++)
		{
			pma::VtxCandidate candidate;
			if (!candidate.Add(fOutTracks[t])) break; // no segments with length > thr

			if (candidate.Add(fEmTracks[u]) && (sqrt(candidate.Mse()) < 1.0))
			{
				candidates.push_back(candidate);
			}
		}
	}
	return candidates;
}

bool pma::PMAlgVertexing::findOneVtx(std::vector< pma::VtxCandidate >& candidates)
{
	bool merged = true;
	while (merged && (candidates.size() > 1))
	{
		merged = false;
		double d_thr = 1.0; // 1.0 = max weighted dist. threshold
		double d, dmin = d_thr;

		size_t k_best, l_best, k = 0;
		while (k < candidates.size() - 1)
		{
			size_t l = k + 1;
			while (l < candidates.size())
			{
				if (candidates[l].Has(candidates[k]))
				{
					candidates[k] = candidates[l];
					candidates.erase(candidates.begin() + l);
					mf::LogVerbatim("pma::PMAlgVertexing") << "removed (k)";
				}
				else if (candidates[k].Has(candidates[l]))
				{
					candidates.erase(candidates.begin() + l);
					mf::LogVerbatim("pma::PMAlgVertexing") << "removed (l)";
				}
				else
				{
					d = candidates[k].Test(candidates[l]);
					if (d < dmin) { dmin = d; k_best = k; l_best = l; }
					l++;
				}
			}
			k++;
		}
		if ((dmin < d_thr) && candidates[k_best].MergeWith(candidates[l_best]))
		{
			mf::LogVerbatim("pma::PMAlgVertexing") << "merged";
			candidates.erase(candidates.begin() + l_best);
			merged = true;
		}
	}

	int s, nmax = 0, c_best = -1;
	double a, amax = 0.0;
	//bool pure_endpoints = false;

	mf::LogVerbatim("pma::PMAlgVertexing") << "*** Vtx candidates: " << candidates.size();
	for (size_t v = 0; v < candidates.size(); v++)
	{
		s = (int)candidates[v].Size(2 * fMinTrackLength);
		a = candidates[v].MaxAngle(1.0);

		if ((s > nmax) || ((s == nmax) && (a > amax)))
		{
			nmax = s; amax = a; c_best = v;
		}

		mf::LogVerbatim("pma::PMAlgVertexing")
			<< "center x:" << candidates[v].Center().X()
			<< " y:" << candidates[v].Center().Y()
			<< " z:" << candidates[v].Center().Z();

		for (size_t i = 0; i < candidates[v].Size(); i++)
			mf::LogVerbatim("pma::PMAlgVertexing")
				<< "     trk:" << i << " "
				<< candidates[v].Track(i).first->size();

		mf::LogVerbatim("pma::PMAlgVertexing")
			<< " dist 3D:" << sqrt(candidates[v].Mse())
			<< " 2D:" << sqrt(candidates[v].Mse2D())
			<< " max ang:" << a;
	}

	if (c_best >= 0)
	{
		return candidates[c_best].JoinTracks(fOutTracks, fEmTracks);
	}
	else return false;
}
// ------------------------------------------------------

size_t pma::PMAlgVertexing::run(
	std::vector< pma::Track3D* >& trk_input)
{
	if (trk_input.size() < 2)
	{
		mf::LogWarning("pma::PMAlgVertexing") << "need min two source tracks!";
		return 0;
	}

	sortTracks(trk_input); // copy input and split by tag/size

	size_t nvtx = 0;
	mf::LogVerbatim("pma::PMAlgVertexing") << "Pass #1:";
	if (fOutTracks.size() > 1)
	{
		bool found = true;
		while (found)
		{
			auto candidates = firstPassCandidates();
			if (candidates.size())
			{
				found = findOneVtx(candidates);
				if (found) nvtx++;
			}
			else found = false;
		}
		mf::LogVerbatim("pma::PMAlgVertexing") << "  " << nvtx << " vertices.";
	}
	else mf::LogVerbatim("pma::PMAlgVertexing") << " ...short tracks only.";

	mf::LogVerbatim("pma::PMAlgVertexing") << "Pass #2:";
	if (fOutTracks.size() && fEmTracks.size())
	{
		bool found = true;
		while (found && fEmTracks.size())
		{
			auto candidates = secondPassCandidates();
			if (candidates.size())
			{
				found = findOneVtx(candidates);
				if (found) nvtx++;
			}
			else found = false;
		}

		mf::LogVerbatim("pma::PMAlgVertexing") << "  " << nvtx << " vertices.";
	}
	else mf::LogVerbatim("pma::PMAlgVertexing") << " ...no tracks.";

	collectTracks(trk_input);

	mergeBrokenTracks(trk_input);

	return nvtx;
}
// ------------------------------------------------------

size_t pma::PMAlgVertexing::run(
	std::vector< pma::Track3D* >& trk_input,
	const std::vector< TVector3 >& vtx_input)
{
	sortTracks(trk_input); // copy input and split by tag/size

	// ....

	//collectTracks(trk_input); // return output in place of (deleted) input

	return 0;
}
// ------------------------------------------------------

std::vector< std::pair<double, double> > pma::PMAlgVertexing::getdQdx(
	const pma::Track3D& trk) const
{
	std::vector< std::pair<double, double> > result;

	unsigned int view = geo::kZ;
	unsigned int nhits = trk.NHits(view);
	unsigned int max = nhits;

	nhits = trk.NHits(geo::kV);
	if (nhits > max) { max = nhits; view = geo::kV; }

	nhits = trk.NHits(geo::kU);
	if (nhits > max) { max = nhits; view = geo::kU; }

	if (max >= 16)
	{
		std::map< size_t, std::vector<double> > dqdx;
		trk.GetRawdEdxSequence(dqdx, view);

		for (size_t i = 0; i < trk.size(); i++)
		{
			auto it = dqdx.find(i);
			if (it != dqdx.end())
			{
				if (it->second[6] > 0.0) // dx > 0
				{
					double dvalue = it->second[5] / it->second[6];
					result.emplace_back(std::pair<double, double>(dvalue, it->second[7]));
				}
			}
		}
	}

	return result;
}
// ------------------------------------------------------

bool pma::PMAlgVertexing::isSingleParticle(pma::Track3D* trk1, pma::Track3D* trk2) const
{
	return false;
}

void pma::PMAlgVertexing::mergeBrokenTracks(std::vector< pma::Track3D* >& trk_input) const
{
	if (trk_input.size() < 2) return;

	mf::LogVerbatim("pma::PMAlgVertexing") << "Find and merge tracks broken by vertices.";
	bool merged = true;
	while (merged)
	{
		merged = false;
		for (size_t t = 0; t < trk_input.size(); t++)
		{
			pma::Track3D* trk1 = trk_input[t];

			pma::Node3D* node = trk1->Nodes().front();
			if (node->Prev())
			{
				pma::Segment3D* seg = static_cast< pma::Segment3D* >(node->Prev());
				pma::Track3D* trk2 = seg->Parent();
				if ((trk1 != trk2) && isSingleParticle(trk2, trk1)) // note: reverse order
				{
				}
			}
			for (size_t n = 0; n < node->NextCount(); n++)
			{
				pma::Segment3D* seg = static_cast< pma::Segment3D* >(node->Next(n));
				pma::Track3D* trk2 = seg->Parent();
				if ((trk1 != trk2) && isSingleParticle(trk1, trk2))
				{
				}
			}

			node = trk1->Nodes().back();
			for (size_t n = 0; n < node->NextCount(); n++)
			{
				pma::Segment3D* seg = static_cast< pma::Segment3D* >(node->Next(n));
				pma::Track3D* trk2 = seg->Parent();
				if ((trk1 != trk2) && isSingleParticle(trk1, trk2))
				{
				}
			}
		}
	}
}
// ------------------------------------------------------

void pma::PMAlgVertexing::splitMergedTracks(std::vector< pma::Track3D* >& trk_input) const
{
	if (trk_input.size() < 1) return;

	mf::LogVerbatim("pma::PMAlgVertexing") << "Find missed vtx by dQ/dx steps along merged tracks.";
	size_t t = 0;
	while (t < trk_input.size())
	{
		t++;
	}
}
// ------------------------------------------------------

std::vector< std::pair< TVector3, std::vector< size_t > > > pma::PMAlgVertexing::getVertices(
	const std::vector< pma::Track3D* >& tracks) const
{
	std::vector< std::pair< TVector3, std::vector< size_t > > > vsel;
	std::vector< pma::Node3D const * > bnodes;

	for (size_t t = 0; t < tracks.size(); t++)
	{
		pma::Track3D const * trk = tracks[t];
		for (auto node : trk->Nodes())
			if (node->IsBranching())
		{
			bool found = false;
			for (size_t n = 0; n < bnodes.size(); n++)
				if (node == bnodes[n])
			{
				vsel[n].second.push_back(t);
				found = true; break;
			}
			if (!found)
			{
				std::vector< size_t > tidx;
				tidx.push_back(t);
				vsel.push_back(std::pair< TVector3, std::vector< size_t > >(node->Point3D(), tidx));
				bnodes.push_back(node);
			}
		}
	}

	return vsel;
}
// ------------------------------------------------------

