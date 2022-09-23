#pragma once

#include "curve.h"
#include "id.h"

#include <vector>

struct Cluster
{
	CurveIDs curve_ids;
	Curve center_curve;
	distance_t cost = std::numeric_limits<distance_t>::max();
};

using Clustering = std::vector<Cluster>;
using ClusterID = ID<Cluster>;
using ClusterIDs = std::vector<ClusterID>;
