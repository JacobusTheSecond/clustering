#include "clustering_algs.h"

#include "curve_simplification.h"
#include "frechet_light.h"
#include "matrix.h"
#include "random.h"
#include "union_find.h"

#include <limits>
#include <omp.h>

namespace
{

// TODO: Computes all distances, not only one per pair.
template <typename Comp>
Clustering linkage(Curves const& curves, int k, int l, Comp comp)
{
	// compute all pairwise Fréchet distances
	FrechetLight frechet_light;
	Matrix<distance_t> dist_matrix(curves.size(), curves.size());
	for (CurveID curve_id1 = 0; curve_id1 < curves.size(); ++curve_id1) {
		for (CurveID curve_id2 = 0; curve_id2 < curves.size(); ++curve_id2) {
			if (curve_id1 == curve_id2) {
				dist_matrix(curve_id1, curve_id2) = 0;
			}
			else if (curve_id1 < curve_id2) {
				auto dist = frechet_light.calcDistance(curves[curve_id1], curves[curve_id2]);
				dist_matrix(curve_id1,curve_id2) = dist;
				dist_matrix(curve_id2,curve_id1) = dist;
			}
		}
	}

	// create initial clusters
	CurveIDs base_set(curves.size());
	std::iota(base_set.begin(), base_set.end(), 0);
	UnionFind<CurveID> union_find(base_set);

	// merge clusters until there are exactly k
	while ((int)union_find.getRoots().size() > k) {
		// find two clusters to merge
		distance_t min_dist = std::numeric_limits<distance_t>::max();
		CurveID min_id1 = CurveID(), min_id2 = CurveID();
		for (auto curve_id1: union_find.getRoots()) {
			for (auto curve_id2: union_find.getRoots()) {
				if (curve_id1 == curve_id2) { continue; }

				auto new_dist = dist_matrix(curve_id1, curve_id2);
				if (new_dist < min_dist) {
					min_dist = new_dist;
					min_id1 = curve_id1;
					min_id2 = curve_id2;
				}
			}
		}

		// merge clusters and adapt distances
		min_id1 = union_find.findRoot(min_id1);
		min_id2 = union_find.findRoot(min_id2);
		auto new_root_id = union_find.uniteSets(min_id1, min_id2);
		auto new_child_id = (new_root_id == min_id1 ? min_id2 : min_id1);

		for (auto curve_id: union_find.getRoots()) {
			auto root_dist = dist_matrix(new_root_id, curve_id);
			auto child_dist = dist_matrix(new_child_id, curve_id);
			dist_matrix(new_root_id, curve_id) = comp(root_dist, child_dist);
		}
	}

	// construct the result
	Clustering result(union_find.getRoots().size());

	std::unordered_map<CurveID, std::size_t> to_cluster_id;
	ClusterID cluster_id = 0;
	for (auto curve_id: union_find.getRoots()) {
		to_cluster_id[curve_id] = cluster_id;
		++cluster_id;
	}
	for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
		auto cluster_id = to_cluster_id[union_find.findRoot(curve_id)];
		result[cluster_id].curve_ids.push_back(curve_id);
	}

	// We just take the root curves as centers. They don't have any special meaning,
	// but at least we supply some centers.
	for (auto curve_id: union_find.getRoots()) {
		auto cluster_id = to_cluster_id[curve_id];
		result[cluster_id].center_curve = simplify(curves[curve_id], l);
	}

	return result;
}

} // end anonymous namespace

std::string toString(ClusterAlg cluster_alg) {
	switch(cluster_alg) {
	case ClusterAlg::SingleLinkage: return "SingleLinkage";
	case ClusterAlg::CompleteLinkage: return "CompleteLinkage";
	case ClusterAlg::Gonzalez: return "Gonzalez";
	}
	ERROR("Unknown cluster_alg.");
}

Clustering computeClustering(Curves const& curves, int k, int l, ClusterAlg cluster_alg)
{
	switch (cluster_alg) {
	case ClusterAlg::SingleLinkage:
		return singleLinkage(curves, k, l);
	case ClusterAlg::CompleteLinkage:
		return completeLinkage(curves, k, l);
	case ClusterAlg::Gonzalez:
		return runGonzalez(curves, k, l);
	}

	ERROR("No matching cluster_alg enum passed.");
}

Clustering singleLinkage(Curves const& curves, int k, int l)
{
	auto min = [](distance_t a, distance_t b) { return std::min<distance_t>(a,b); };
	return linkage(curves, k, l, min);
}

Clustering completeLinkage(Curves const& curves, int k, int l)
{
	auto max = [](distance_t a, distance_t b) { return std::max<distance_t>(a,b); };
	return linkage(curves, k, l, max);
}

Clustering runGonzalez(Curves const& curves, int k, int l)
{
	Clustering result;

	FrechetLight frechet_light;
	auto max_dist = std::numeric_limits<distance_t>::max();
	std::vector<distance_t> distances_to_center(curves.size(), max_dist);
	ClusterIDs closest_center(curves.size());

	Random random;
	CurveID center_id = random.getUniformInt(0, curves.size()-1);

    std::cout << "Center #1/"<< k << " found\n";

	// add as center and update closest distances to center
	auto center_curve = simplify(curves[center_id], l);
	result.push_back({{}, center_curve});
    int count = 0;
#pragma omp parallel for private(frechet_light) default(shared)
	for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
        count += 1;
        if ((count)%(curves.size()/10)==0)
            std::cout << "    Progress on Frechet computations:" << ((100*(count))/curves.size()) << "% (" <<count << "/" << curves.size() << ") with "<< omp_get_num_threads() << " threads\n";
		auto& current_dist = distances_to_center[curve_id];
		if (frechet_light.lessThanWithFilters(current_dist, center_curve, curves[curve_id])) {
			current_dist = frechet_light.calcDistance(center_curve, curves[curve_id]);
			closest_center[curve_id] = result.size()-1;
		}
	}

	while ((int)result.size() < k) {
		auto center_it = std::max_element(distances_to_center.begin(), distances_to_center.end());
		auto center_id = std::distance(distances_to_center.begin(), center_it);
		auto center_curve = simplify(curves[center_id], l);
		result.push_back({{}, center_curve});
        std::cout << "Center #" << (result.size()+1) <<"/"<< k << " found\n";
        int count = 0;
#pragma omp parallel for private(frechet_light) default(shared)
		for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
            count +=1;
            if ((count)%(curves.size()/10)==0)
                std::cout << "    Progress on Frechet computations:" << ((100*(count))/curves.size()) << "% (" <<count << "/" << curves.size() << ") with "<< omp_get_num_threads() << " threads\n";
			auto& current_dist = distances_to_center[curve_id];
			if (frechet_light.lessThanWithFilters(current_dist, center_curve, curves[curve_id])) {
				current_dist = frechet_light.calcDistance(center_curve, curves[curve_id]);
				closest_center[curve_id] = result.size()-1;
			}
		}
	}

	for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
		auto cluster_id = closest_center[curve_id];
		result[cluster_id].curve_ids.push_back(curve_id);
	}

	return result;
}

void updateClustering(Curves const& curves, Clustering& clustering)
{
	FrechetLight frechet_light;

	// clear clusters
	for (auto& cluster: clustering) {
		cluster.curve_ids.clear();
	}

    std::cout << "Updating...\n";

	// compute new clusters
    int count = 0;
#pragma omp parallel for shared(curves,clustering,count,std::cout) private(frechet_light) default(none)
	for (CurveID curve_id = 0; curve_id < curves.size(); ++curve_id) {
        count += 1;
        if(count%(curves.size()/10)==0)
            std::cout << "Progress on Frechet updates: " << (100*count)/curves.size() << "% ("<< count<< "/"<<curves.size() <<") with " << omp_get_num_threads() << " threads\n";
		distance_t min_dist = std::numeric_limits<distance_t>::max();
		ClusterID min_cluster_id = ClusterID();
		for (ClusterID cluster_id = 0; cluster_id < clustering.size(); ++cluster_id) {
			auto const& center_curve = clustering[cluster_id].center_curve;
			if (frechet_light.lessThanWithFilters(min_dist, curves[curve_id], center_curve)) {
				min_dist = frechet_light.calcDistance(curves[curve_id], center_curve);
				min_cluster_id = cluster_id;
			}
		}

		clustering[min_cluster_id].curve_ids.push_back(curve_id);
	}
}

distance_t calcDiameter(Curves const& curves, CurveIDs const& curve_ids)
{
	FrechetLight frechet_light;
	distance_t max_distance = 0.;
	for (CurveID curve_id1 = 0; curve_id1 < curves.size(); ++curve_id1) {
		for (CurveID curve_id2 = 0; curve_id2 < curves.size(); ++curve_id2) {
			if (curve_id1 == curve_id2) {
				continue;
			}
			else if (curve_id1 < curve_id2) {
				auto dist = frechet_light.calcDistance(curves[curve_id1], curves[curve_id2]);
				max_distance = std::max(max_distance, dist);
			}
		}
	}

	return max_distance;
}
