#include "center_algs.h"

#include "curve_simplification.h"
#include "frechet_light.h"
#include "frechet_matching.h"
#include "geometric_algs.h"

#include <limits>

namespace
{

bool calcKXCenters(Curves const& curves, Clustering& clustering, int l, C2CDist c2c_dist)
{
	bool found_new_center = false;

	// compute cluster costs in case they haven't been computed yet
	for (auto& cluster: clustering) {
		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist);
		}
	}

	for (auto& cluster: clustering) {
		for (CurveID curve_id1: cluster.curve_ids) {
			auto simplified_curve = simplify(curves[curve_id1], l);
			auto dist = calcC2CDist(curves, simplified_curve, cluster.curve_ids, c2c_dist);
			if (dist < cluster.cost) {
				cluster.center_curve = simplified_curve;
				cluster.cost = dist;
				found_new_center = true;
			}
		}
	}

	return found_new_center;
}

} // end anonymous namespace

std::string toString(CenterAlg center_alg) {
	switch(center_alg) {
	case CenterAlg::kMedian: return "kMedian";
	case CenterAlg::kMeans: return "kMeans";
	case CenterAlg::kCenter: return "kCenter";
	case CenterAlg::FSA: return "FSA";
	}
	ERROR("Unknown center_alg.");
}

distance_t calcC2CDist(
	Curves const& curves, Curve const& center_curve, CurveIDs const& curve_ids, C2CDist c2c_dist)
{
	FrechetLight frechet_light;

	distance_t dist = 0;
#pragma omp parallel for private(frechet_light) shared(curve_ids,center_curve,curves,c2c_dist,dist) default(none)
	for (auto curve_id: curve_ids) {
		auto curve_dist = frechet_light.calcDistance(center_curve, curves[curve_id]);
#pragma omp critical
        {
            switch (c2c_dist) {
                case C2CDist::Median:
                    dist += curve_dist;
                    break;
                case C2CDist::Mean:
                    dist += curve_dist * curve_dist;
                    break;
                case C2CDist::Max:
                    dist = std::max(dist, curve_dist);
                    break;
            }
        }
	}

	return dist;
}

bool computerCenters(Curves const& curves, Clustering& clustering, int l, CenterAlg center_alg)
{
	switch (center_alg) {
	case CenterAlg::kMedian:
		return calcKMedianCenters(curves, clustering, l);
	case CenterAlg::kMeans:
		return calcKMeansCenters(curves, clustering, l);
	case CenterAlg::kCenter:
		return calcKCenterCenters(curves, clustering, l);
	case CenterAlg::FSA:
		return calcFSACenters(curves, clustering, l);
	}

	ERROR("No matching center_alg enum passed.");
}

bool calcKMedianCenters(Curves const& curves, Clustering& clustering, int l)
{
	return calcKXCenters(curves, clustering, l, C2CDist::Median);
}

bool calcKMeansCenters(Curves const& curves, Clustering& clustering, int l)
{
	return calcKXCenters(curves, clustering, l, C2CDist::Mean);
}

bool calcKCenterCenters(Curves const& curves, Clustering& clustering, int l)
{
	return calcKXCenters(curves, clustering, l, C2CDist::Max);
}

// TODO: There is some unnecessary pushing around of data here. Fix that to increase performance.
bool calcFSACenters(Curves const& curves, Clustering& clustering, int l, C2CDist c2c_dist)
{
	bool found_new_center = false;
	FrechetLight frechet_light;

	// compute cluster costs in case they haven't been computed yet
    int count = 0;
#pragma omp parallel for shared(clustering,curves,c2c_dist,count,std::cout) default(none)
	for (auto& cluster: clustering) {
        count += 1;
        if(count%(clustering.size()/100)==0)
            std::cout << "FSA Phase 1 Progress: " << (100*count)/clustering.size() << "%\n";
		if (cluster.cost == std::numeric_limits<distance_t>::max()) {
			cluster.cost = calcC2CDist(curves, cluster.center_curve, cluster.curve_ids, c2c_dist);
		}
	}

    count = 0;
#pragma omp parallel for shared(clustering,curves,c2c_dist,found_new_center,count,std::cout) default(none)
	for (auto& cluster: clustering) {

        count += 1;
        if(count%(clustering.size()/100)==0)
            std::cout << "FSA Phase 2 Progress: " << (100*count)/clustering.size() << "%\n";
		std::vector<Points> matchings;
		auto const& center_curve = cluster.center_curve;
		Curve new_center_curve;

		for (auto curve_id: cluster.curve_ids) {
			auto const& curve = curves[curve_id];
			auto matching = calcMatching(cluster.center_curve, curve);
			matchings.push_back(std::move(matching));
		}

		for (PointID point_id = 0; point_id < center_curve.size(); ++point_id) {
			Points matching_points;
			for (auto const& matching: matchings) {
				matching_points.push_back(matching[point_id]);
			}
			auto min_enclosing_circle = calcMinEnclosingCircle(matching_points);
			new_center_curve.push_back(min_enclosing_circle.center);
		}

		if (center_curve != new_center_curve) {
			auto new_dist = calcC2CDist(curves, new_center_curve, cluster.curve_ids, c2c_dist);
			if (new_dist < cluster.cost) {
				cluster.center_curve = std::move(new_center_curve);
				cluster.cost = new_dist;
				found_new_center = true;
			}
		}
	}

	return found_new_center;
}
