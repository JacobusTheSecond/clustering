#pragma once

//#include "center_algs.h"
#include "clustering_algs.h"
#include "curve.h"
#include "candidate.h"
/*
Clustering computeCenterClustering(
	Curves const& curves, int k, int l, ClusterAlg cluster_alg, CenterAlg center_alg, int max_rounds = 10);

Clustering approxCover(Curves& curves, double delta,int l,int max_rounds = 10);
*/
Curves greedyCover(Curves& curves, double delta, int l, int max_rounds = 10, bool show = false);
std::vector<std::pair<int,Candidate>> greedyCoverUnsanitizedOutput(Curves& curves, double delta, int l, int max_rounds = 10, bool show = false);
