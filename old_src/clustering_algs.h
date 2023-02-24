#pragma once

#include "basic_types.h"
#include "curve.h"

enum class ClusterAlg {
	SingleLinkage,
	CompleteLinkage,
	Gonzalez,
};
std::string toString(ClusterAlg cluster_alg);

Clustering computeClustering(Curves const& curves, int k, int l, ClusterAlg cluster_alg);

Clustering singleLinkage(Curves const& curves, int k, int l);
Clustering completeLinkage(Curves const& curves, int k, int l);
Clustering runGonzalez(Curves const& curves, int k, int l);

// assign curves to closest clusters
void updateClustering(Curves const& curves, Clustering& clustering);

distance_t calcDiameter(Curves const& curves, CurveIDs const& curve_ids);
