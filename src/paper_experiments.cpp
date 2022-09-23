#include "center_clustering_algs.h"
#include "defs.h"
#include "frechet_light.h"
#include "io.h"

#include <algorithm>
#include <chrono>
#include <string>
#include <unordered_map>
#include <vector>

// timing
using hrc = std::chrono::high_resolution_clock;
using time_point = hrc::time_point;
using ns = std::chrono::nanoseconds;


// forward declarations
void testAlgorithms();
void CompleteLinkageVSOurs();
void findBestKL();
void printCurrentKL();
void computeAllClusterings();
void testConvergence();
void siteClustering();
void exportAllWaypoints();
void recognizePigeons();
void characterClassification();
void characterClassificationNaive();
void plotCharacterCenters();

int main()
{
	// testAlgorithms();
	// CompleteLinkageVSOurs();
	// findBestKL();
	// printCurrentKL();
	computeAllClusterings();
	// testConvergence();
	// siteClustering();
	// exportAllWaypoints();
	// recognizePigeons();
	// characterClassification();
	// characterClassificationNaive();
	plotCharacterCenters();
}

// fixed paths
namespace
{

using Strings = std::vector<std::string>;
Strings const site_directories = {
	"../data/Data_for_Mann_et_al_RSBL 2/Bladon & Church route recapping/bladon heath/",
	"../data/Data_for_Mann_et_al_RSBL 2/Bladon & Church route recapping/church hanborough/",
	"../data/Data_for_Mann_et_al_RSBL 2/Horspath/",
	"../data/Data_for_Mann_et_al_RSBL 2/Weston/"
};
std::vector<Strings> const pigeon_directories = {
	{"a55", "brc", "c17", "c35", "p29", "p39", "p94"},
	{"a94", "c22", "c70", "k77", "l29", "liv", "r47", "s93"},
	{"H22", "H27", "H30", "H35", "H38", "H41", "H42", "H71"},
	{"H23", "H31", "H32", "H34", "H36", "H50", "H58", "H62"}
};

std::string const characters_path = "../data/characters/data/";

std::vector<int> site_ks = { 9, 10, 10, 8 };
std::vector<std::vector<int>> pigeon_ks = {
	{4, 3, 3, 4, 4, 5, 4},
	{4, 3, 5, 4, 4, 3, 4, 5},
	{5, 3, 4, 4, 4, 6, 4, 5},
	{6, 3, 4, 4, 3, 5, 5, 6}
};

std::vector<int> site_ls = { 8, 10, 12, 15 };
std::vector<std::vector<int>> pigeon_ls = {
	{ 11, 10, 8, 8, 10, 14, 11},
	{ 7, 11, 10, 9, 11, 9, 10, 10},
	{ 9, 12, 6, 12, 10, 11, 11, 9},
	{ 10, 9, 11, 13, 11, 11, 11, 12}
};

} // end anonymous namespace

void testAlgorithms()
{
	std::size_t site_id = 0;
	auto header_size = 1;

	auto base_path = site_directories[site_id];
	auto curves = io::readCurves(base_path, header_size);

	std::vector<ClusterAlg> cluster_algs = {
		ClusterAlg::SingleLinkage, ClusterAlg::CompleteLinkage, ClusterAlg::Gonzalez };
	std::vector<CenterAlg> center_algs = {
		CenterAlg::kMedian, CenterAlg::kMeans, CenterAlg::kCenter, CenterAlg::FSA };

	int k = site_ks[site_id];
	int l = site_ls[site_id];
	for (auto cluster_alg: cluster_algs) {
		for (auto center_alg: center_algs) {
			auto start = hrc::now();
			auto clustering = computeCenterClustering(curves, k, l, cluster_alg, center_alg, 1);
			auto time = std::chrono::duration_cast<ns>(hrc::now()-start).count();
			distance_t median = 0, mean = 0, center = 0;
			for (auto const& c: clustering) {
				median += calcC2CDist(curves, c.center_curve, c.curve_ids, C2CDist::Median);
				mean += calcC2CDist(curves, c.center_curve, c.curve_ids, C2CDist::Mean);
				center += calcC2CDist(curves, c.center_curve, c.curve_ids, C2CDist::Max);
			}
			center /= clustering.size();
			std::cout << toString(cluster_alg) << " / " << toString(center_alg) << ": "
					  << mean << " / " << median << " / " << center << std::endl;
			std::cout << "Time: " << time/1000000. << " ms\n" << std::endl;
		}
	}
}

void CompleteLinkageVSOurs()
{
	std::size_t site_id = 0;
	auto header_size = 1;

	auto base_path = site_directories[site_id];
	auto curves = io::readCurves(base_path, header_size);

	std::vector<ClusterAlg> cluster_algs = {
		ClusterAlg::SingleLinkage, ClusterAlg::CompleteLinkage, ClusterAlg::Gonzalez };
	std::vector<CenterAlg> center_algs = {
		CenterAlg::kMedian, CenterAlg::kMeans, CenterAlg::kCenter, CenterAlg::FSA };

	int k = site_ks[site_id];
	int l = site_ls[site_id];

	distance_t center_avg1 = 0.;
	distance_t diamter_avg1 = 0.;
	std::size_t time_avg1 = 0;
	distance_t center_avg2 = 0.;
	distance_t diamter_avg2 = 0.;
	std::size_t time_avg2 = 0;
	for (int i = 0; i < 10; ++i) {
		{
			// compute clustering for complete linkage
			auto start = hrc::now();
			auto clustering = computeClustering(curves, k, l, ClusterAlg::CompleteLinkage);
			computerCenters(curves, clustering, l, CenterAlg::kMeans);
			auto time = std::chrono::duration_cast<ns>(hrc::now()-start).count();
			distance_t max_diameter = 0, center = 0;
			for (auto const& c: clustering) {
				auto diameter = calcDiameter(curves, c.curve_ids);
				max_diameter = std::max(max_diameter, diameter);
				center += calcC2CDist(curves, c.center_curve, c.curve_ids, C2CDist::Max);
			}
			center /= clustering.size();

			diamter_avg1 += max_diameter;
			center_avg1 += center;
			time_avg1 += time;
		}

		{
			// compute clustering for Gonzalez + FSA
			auto start = hrc::now();
			auto clustering = computeCenterClustering(curves, k, l, ClusterAlg::Gonzalez, CenterAlg::FSA, 1);
			auto time = std::chrono::duration_cast<ns>(hrc::now()-start).count();
			distance_t max_diameter = 0, center = 0;
			for (auto const& c: clustering) {
				auto diameter = calcDiameter(curves, c.curve_ids);
				max_diameter = std::max(max_diameter, diameter);
				center += calcC2CDist(curves, c.center_curve, c.curve_ids, C2CDist::Max);
			}
			center /= clustering.size();

			diamter_avg2 += max_diameter;
			center_avg2 += center;
			time_avg2 += time;
		}
	}

	diamter_avg1 /= 10.;
	center_avg1 /= 10.;
	time_avg1 /= 10.;
	diamter_avg2 /= 10.;
	center_avg2 /= 10.;
	time_avg2 /= 10.;

	std::cout << "Complete Linkage with mean centers" << ": "
		<< diamter_avg1 << " / " << center_avg1 << std::endl;
	std::cout << "Time: " << time_avg1/1000000. << " ms\n" << std::endl;

	std::cout << "Gonzalez + FSA" << ": "
		<< diamter_avg2 << " / " << center_avg2 << std::endl;
	std::cout << "Time: " << time_avg2/1000000. << " ms\n" << std::endl;
}

namespace
{

void findBestKL(std::string const& path)
{
	std::cout << "Find best k and l for " << path << "\n";

	auto header_size = 1;
	auto curves = io::readCurves(path, header_size);

	ClusterAlg cluster_alg = ClusterAlg::Gonzalez;
	CenterAlg center_alg = CenterAlg::FSA;

	int const k_min = 1, k_max = 10;
	int const l_min = 3, l_max = 15;
	for (int k = k_min; k <= k_max; ++k) {
		// std::cout << "k = " << k << "\n";
		for (int l = l_min; l <= l_max; ++l) {
			auto clustering = computeCenterClustering(curves, k, l, cluster_alg, center_alg, 10);
			distance_t median = 0, mean = 0, center = 0;
			for (auto const& c: clustering) {
				median += calcC2CDist(curves, c.center_curve, c.curve_ids, C2CDist::Median);
				mean += calcC2CDist(curves, c.center_curve, c.curve_ids, C2CDist::Mean);
				center += calcC2CDist(curves, c.center_curve, c.curve_ids, C2CDist::Max);
			}
			center /= clustering.size();
			// std::cout << "l = " << l << ": " << mean << " / " << median << " / " << center << std::endl;
			std::cout << center << std::endl;
		}
		if (l_min != l_max) { std::cout << "\n"; }
	}

	std::cout << "\n";
}

}

void findBestKL()
{
	for (std::size_t site_id = 0; site_id < site_directories.size(); ++site_id) {
		auto const& site_dir = site_directories[site_id];
		// for the single pigeons
		for (std::size_t pigeon_id = 0; pigeon_id < pigeon_directories[site_id].size(); ++pigeon_id) {
			auto const& pigeon_dir = pigeon_directories[site_id][pigeon_id];
			auto const path = site_dir + pigeon_dir;
			findBestKL(path);
		}
		// for the whole site
		findBestKL(site_dir);
	}
}

void printCurrentKL()
{
	for (std::size_t site_id = 0; site_id < site_directories.size(); ++site_id) {
		// for the single pigeons
		for (std::size_t pigeon_id = 0; pigeon_id < pigeon_directories[site_id].size(); ++pigeon_id) {
			auto const k = pigeon_ks[site_id][pigeon_id];
			auto const l = pigeon_ls[site_id][pigeon_id];
			std::cout << pigeon_directories[site_id][pigeon_id] << ": " << "k = " << k << ", l = " << l << std::endl;
		}
		// for the whole site
		auto const k = site_ks[site_id];
		auto const l = site_ls[site_id];
		std::cout << "site " << site_id << ": " << "k = " << k << ", l = " << l << std::endl;
	}
}

void computeAllClusterings()
{
	auto header_size = 1;

	ClusterAlg cluster_alg = ClusterAlg::Gonzalez;
	CenterAlg center_alg = CenterAlg::FSA;

	for (std::size_t site_id = 0; site_id < site_directories.size(); ++site_id) {
		auto const& site_dir = site_directories[site_id];
		// for the single pigeons
		for (std::size_t pigeon_id = 0; pigeon_id < pigeon_directories[site_id].size(); ++pigeon_id) {
			auto const& pigeon_dir = pigeon_directories[site_id][pigeon_id];
			auto const path = site_dir + pigeon_dir;
			auto curves = io::readCurves(path, header_size);
			auto k = pigeon_ks[site_id][pigeon_id];
			auto l = pigeon_ls[site_id][pigeon_id];
			auto clustering = computeCenterClustering(curves, k, l, cluster_alg, center_alg, 3);
			auto filename = std::to_string(site_id) + "_" + pigeon_dir + ".clustering";
			io::exportClustering(filename, path, clustering, curves);
		}
		// for the whole site
		auto curves = io::readCurves(site_dir, header_size);
		auto k = site_ks[site_id];
		auto l = site_ls[site_id];
		auto clustering = computeCenterClustering(curves, k, l, cluster_alg, center_alg, 3);
		auto filename = std::to_string(site_id) + ".clustering";
		io::exportClustering(filename, site_dir, clustering, curves);
	}
}

namespace
{

void testConvergence(std::string const& base_path, int k, int l)
{
	auto header_size = 1;
	auto curves = io::readCurves(base_path, header_size);

	ClusterAlg cluster_alg = ClusterAlg::Gonzalez;
	CenterAlg center_alg = CenterAlg::FSA;

	auto clustering = computeCenterClustering(curves, k, l, cluster_alg, center_alg, 10);

	// calculate max_filename_id
	std::vector<std::pair<std::string, ClusterID>> filename_id_pairs;
	for (auto cluster_id = 0; cluster_id < (int)clustering.size(); ++cluster_id) {
		auto const& cluster = clustering[cluster_id];
		for (auto curve_id: cluster.curve_ids) {
			filename_id_pairs.emplace_back(curves[curve_id].filename, cluster_id);
		}
	}
	std::sort(filename_id_pairs.begin(), filename_id_pairs.end());
	std::string max_filename = filename_id_pairs.back().first;
	ClusterID max_filename_id = filename_id_pairs.back().second;

	// calculate max_size_id
	auto size_comp = [&](Cluster const& c1, Cluster const& c2) {
		return c1.curve_ids.size() < c2.curve_ids.size();
	};
	auto it = std::max_element(clustering.begin(), clustering.end(), size_comp);
	ClusterID max_size_id = std::distance(clustering.begin(), it);

	if (max_filename_id == max_size_id) {
		// Additionally calculate suffix length
		int i = filename_id_pairs.size()-1;
		while (i >= 0) {
			if (filename_id_pairs[i].second != max_size_id) { break; }
			--i;
		}
		int suffix_length = filename_id_pairs.size()-1-i;
		std::cout << "Yes! Suffix length: " << suffix_length << "/" << curves.size() << "\n";
	}
	else {
		std::cout << "No!" << "\n";
	}
}

}

void testConvergence()
{
	for (std::size_t site_id = 0; site_id < site_directories.size(); ++site_id) {
		auto const& site_dir = site_directories[site_id];
		for (std::size_t pigeon_id = 0; pigeon_id < pigeon_directories[site_id].size(); ++pigeon_id) {
			auto const path = site_dir + pigeon_directories[site_id][pigeon_id];
			auto const k = pigeon_ks[site_id][pigeon_id];
			auto const l = pigeon_ls[site_id][pigeon_id];
			std::cout << pigeon_directories[site_id][pigeon_id] << ": ";
			testConvergence(path, k, l);
		}
	}
}

void siteClustering()
{
	int header_size = 1;
	std::size_t const site_id = 0;

	auto const path = site_directories[site_id];
	auto const curves = io::readCurves(path, header_size);

	ClusterAlg cluster_alg = ClusterAlg::Gonzalez;
	CenterAlg center_alg = CenterAlg::FSA;
	int const k = 9;
	int const l = 12;

	auto clustering = computeCenterClustering(curves, k, l, cluster_alg, center_alg, 1);
	io::exportClustering("out.clustering", path, clustering, curves);
	io::exportCentersGPX("centers.gpx", clustering);
}

namespace
{

void exportWaypoints(std::size_t const site_id)
{
	int header_size = 1;

	ClusterAlg cluster_alg = ClusterAlg::Gonzalez;
	CenterAlg center_alg = CenterAlg::FSA;

	auto const& site_dir = site_directories[site_id];
	Curves majority_center_curves;
	for (std::size_t pigeon_id = 0; pigeon_id < pigeon_directories[site_id].size(); ++pigeon_id) {
		auto const path = site_dir + pigeon_directories[site_id][pigeon_id];
		auto const curves = io::readCurves(path, header_size);
		int const k = pigeon_ks[site_id][pigeon_id];
		int const l = pigeon_ls[site_id][pigeon_id];
		auto clustering = computeCenterClustering(curves, k, l, cluster_alg, center_alg, 10);

		auto size_comp = [&](Cluster const& c1, Cluster const& c2) {
			return c1.curve_ids.size() < c2.curve_ids.size();
		};
		auto it = std::max_element(clustering.begin(), clustering.end(), size_comp);
		majority_center_curves.push_back(std::move(it->center_curve));
	}

	io::exportCurvesGPX("waypoints_" + std::to_string(site_id) + ".gpx", majority_center_curves);
}

} // end anonymous namespace

void exportAllWaypoints()
{
	for (std::size_t site_id = 0; site_id < site_directories.size(); ++site_id) {
		exportWaypoints(site_id);
	}
}

void recognizePigeons()
{
	int header_size = 1;
	ClusterAlg cluster_alg = ClusterAlg::Gonzalez;
	CenterAlg center_alg = CenterAlg::FSA;

	FrechetLight frechet_light;

	for (std::size_t site_id = 0; site_id < site_directories.size(); ++site_id) {
		auto const& site_path = site_directories[site_id];
		std::size_t number_of_pigeons = pigeon_directories[site_id].size();

		// read curves
		std::vector<Curves> pigeon_curves;
		for (std::size_t pigeon_id = 0; pigeon_id < number_of_pigeons; ++pigeon_id) {
			auto const path = site_path + pigeon_directories[site_id][pigeon_id];
			auto new_curves = io::readCurves(path, header_size);
			pigeon_curves.push_back(std::move(new_curves));
		}
		
		// learning
		Curves centers;
		for (std::size_t pigeon_id = 0; pigeon_id < number_of_pigeons; ++pigeon_id) {
			auto const k = pigeon_ks[site_id][pigeon_id];
			auto const l = pigeon_ls[site_id][pigeon_id];
			auto clustering = computeCenterClustering(
				pigeon_curves[pigeon_id], k, l, cluster_alg, center_alg, 10);
			for (auto& cluster: clustering) {
				centers.push_back(std::move(cluster.center_curve));
			}

			// for simplicity, fill up centers using the first such that there are k in total.
			if (clustering.empty()) { ERROR("The clustering has to contain at least one center."); }
			auto representative = centers.back();
			while (centers.size()%k != 0) { centers.push_back(representative); }
		}

		// testing on all trajectories
		std::size_t correct = 0;
		std::size_t wrong = 0;
		for (std::size_t pigeon_id = 0; pigeon_id < number_of_pigeons; ++pigeon_id) {
			for (auto const& curve: pigeon_curves[pigeon_id]) {
				// classify curve
				distance_t min_dist = std::numeric_limits<distance_t>::max();
				std::size_t min_id = 0;
				for (CurveID center_id = 0; center_id < centers.size(); ++center_id) {
					auto const& center = centers[center_id];
					auto dist = frechet_light.calcDistance(curve, center);
					if (dist < min_dist) {
						min_dist = dist;
						min_id = center_id;
					}
				}

				// check if classification is correct
				auto const k = pigeon_ks[site_id][pigeon_id];
				if (pigeon_id == min_id/k) ++correct; else ++wrong;
			}
		}

		std::cout << "Site " << site_id << std::endl;
		std::cout << "Correct classifications: " << correct << "\n";
		std::cout << "Wrong classifications: " << wrong << "\n";
		std::cout << "Success rate: " << float(correct)/(correct+wrong) << "\n\n";
	}
}

void characterClassification()
{
	int header_size = 0;
	ClusterAlg cluster_alg = ClusterAlg::Gonzalez;
	CenterAlg center_alg = CenterAlg::FSA;
	int const k = 6;
	int const l = 12;
	int const k_xv = 5; // the k of the cross-validation

	auto curves = io::readCurves(characters_path, header_size);
	std::random_shuffle(curves.begin(), curves.end());
	FrechetLight frechet_light;

	// put curves into own vectors
	std::unordered_map<char, Curves> char_to_curves;
	for (auto const& curve: curves) {
		char_to_curves[curve.filename.front()].push_back(curve);
	}
	std::vector<std::vector<Curves>> curves_by_char;
	for (auto const& pair: char_to_curves) {
		curves_by_char.emplace_back(k_xv);
		auto& curves = pair.second;
		for (CurveID i = 0; i < curves.size(); ++i) {
			curves_by_char.back()[i%k_xv].push_back(std::move(curves[i]));
		}
	}

	int correct = 0;
	int wrong = 0;
	std::size_t time = 0;
	for (int i = 0; i < k_xv; ++i) {
		Curves centers;
		for (auto const& curves_slots: curves_by_char) {
			Curves training_curves;
			for (int j = 0; j < k_xv; ++j) {
				if (i == j) { continue; }
				training_curves.insert(
					training_curves.end(), curves_slots[j].begin(), curves_slots[j].end());
			}
			auto clustering = computeCenterClustering(
				training_curves, k, l, cluster_alg, center_alg, 1);
			for (auto const& cluster: clustering) {
				centers.push_back(std::move(cluster.center_curve));
			}

			// for simplicity, fill up centers using the first such that there are k in total.
			if (clustering.empty()) { ERROR("The clustering has to contain at least one center."); }
			auto representative = centers.back();
			while (centers.size()%k != 0) { centers.push_back(representative); }
		}

		// cross-validate on remaining data
		auto start = hrc::now();
		for (std::size_t char_id = 0; char_id < curves_by_char.size(); ++char_id) {
			auto const& curves = curves_by_char[char_id][i];
			for (auto const& curve: curves) {
				// classify curve
				distance_t min_dist = std::numeric_limits<distance_t>::max();
				std::size_t min_id = 0;
				for (CurveID center_id = 0; center_id < centers.size(); ++center_id) {
					auto const& center = centers[center_id];
					if (frechet_light.lessThanWithFilters(min_dist, curve, center)) {
						auto dist = frechet_light.calcDistance(curve, center);
						min_dist = dist;
						min_id = center_id;
					}
				}

				// check if classification is correct
				if (char_id == min_id/k) ++correct; else ++wrong;
			}
		}
		time += std::chrono::duration_cast<ns>(hrc::now()-start).count();
	}

	std::cout << "Time: " << time/1000000. << " ms\n";
	std::cout << "Correct classifications: " << correct << "\n";
	std::cout << "Wrong classifications: " << wrong << "\n";
	std::cout << "Success rate: " << float(correct)/(correct+wrong) << "\n";
}

void characterClassificationNaive()
{	int header_size = 0;
	int const k_xv = 10; // the k of the cross-validation

	auto curves = io::readCurves(characters_path, header_size);
	std::random_shuffle(curves.begin(), curves.end());
	FrechetLight frechet_light;

	// put curves into own vectors
	std::unordered_map<char, Curves> char_to_curves;
	for (auto const& curve: curves) {
		char_to_curves[curve.filename.front()].push_back(curve);
	}
	std::vector<std::vector<Curves>> curves_by_char;
	for (auto const& pair: char_to_curves) {
		curves_by_char.emplace_back(k_xv);
		auto& curves = pair.second;
		for (CurveID i = 0; i < curves.size(); ++i) {
			curves_by_char.back()[i%k_xv].push_back(std::move(curves[i]));
		}
	}

	int correct = 0;
	int wrong = 0;
	auto start = hrc::now();
	for (int i = 0; i < k_xv; ++i) {
		// cross-validate on remaining data
		for (std::size_t char_id = 0; char_id < curves_by_char.size(); ++char_id) {
			auto const& curves = curves_by_char[char_id][i];
			for (auto const& curve: curves) {
				// classify curve
				distance_t min_dist = std::numeric_limits<distance_t>::max();
				std::size_t min_char_id = 0;
				for (std::size_t c = 0; c < curves_by_char.size(); ++c) {
					for (int j = 0; j < k_xv; ++j) {
						if (i == j) { continue; }
						for (auto const& test_curve: curves_by_char[c][j]) {
							if (frechet_light.lessThanWithFilters(min_dist, curve, test_curve)) {
								auto dist = frechet_light.calcDistance(curve, test_curve);
								min_dist = dist;
								min_char_id = c;
							}
						}
					}
				}

				// check if classification is correct
				if (min_char_id == char_id) ++correct; else ++wrong;
			}
		}
	}
	std::size_t time = std::chrono::duration_cast<ns>(hrc::now()-start).count();

	std::cout << "Time: " << time/1000000. << " ms\n";
	std::cout << "Correct classifications: " << correct << "\n";
	std::cout << "Wrong classifications: " << wrong << "\n";
	std::cout << "Success rate: " << float(correct)/(correct+wrong) << "\n";
}

void plotCharacterCenters()
{
	int header_size = 0;
	ClusterAlg cluster_alg = ClusterAlg::Gonzalez;
	CenterAlg center_alg = CenterAlg::FSA;
	int const k = 6;
	int const l = 12;
	int const k_xv = 5; // the k of the cross-validation

	auto curves = io::readCurves(characters_path, header_size);
	std::random_shuffle(curves.begin(), curves.end());
	FrechetLight frechet_light;

	// put curves into own vectors
	std::unordered_map<char, Curves> char_to_curves;
	for (auto const& curve: curves) {
		char_to_curves[curve.filename.front()].push_back(curve);
	}
	std::vector<std::vector<Curves>> curves_by_char;
	for (auto const& pair: char_to_curves) {
		curves_by_char.emplace_back(k_xv);
		auto& curves = pair.second;
		for (CurveID i = 0; i < curves.size(); ++i) {
			curves_by_char.back()[i%k_xv].push_back(std::move(curves[i]));
		}
	}

	for (int i = 0; i < k_xv; ++i) {
		for (auto const& curves_slots: curves_by_char) {
			Curves training_curves;
			for (int j = 0; j < k_xv; ++j) {
				if (i == j) { continue; }
				training_curves.insert(
					training_curves.end(), curves_slots[j].begin(), curves_slots[j].end());
			}
			auto clustering = computeCenterClustering(
				training_curves, k, l, cluster_alg, center_alg, 1);
			auto const filename = curves_slots[i].front().filename[0] + std::to_string(i) + ".clustering";
			// io::exportCenters(filename, clustering);
			io::exportClustering(filename, characters_path, clustering, training_curves);
		}
	}
}
