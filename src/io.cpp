#include "io.h"

#include <fstream>
#include <sstream>

namespace
{

void readCurve(std::ifstream& curve_file, Curve& curve, int header_size = 0)
{
	// Read everything into a stringstream.
	std::stringstream ss;
	ss << curve_file.rdbuf();

	auto ignore_count = std::numeric_limits<std::streamsize>::max();

	// ignore header lines
	for (int i = 0; i < header_size; ++i) {
		ss.ignore(ignore_count, '\n');
	}

	std::string x_str, y_str;
	while (ss >> x_str >> y_str) {
		distance_t x, y;
		x = std::stod(x_str);
		y = std::stod(y_str);

		ss.ignore(ignore_count, '\n');
		// ignore duplicate coordinates
		if (curve.size() && curve.back().x() == x && curve.back().y() == y) {
			continue;
		}
		curve.push_back({x, y});
	}
}

} // end anonymous namespace

namespace io
{

Curves readCurves(std::string const& base_path, int header_size)
{
	// read filenames of curve files
	std::string curve_data_file = base_path + "/dataset.txt";
	std::ifstream file(curve_data_file);
	std::vector<std::string> curve_filenames;
	if (file.is_open()) {
		std::string line;
		while (std::getline(file, line)) {
			curve_filenames.push_back(line);
		}
	}
	else {
		ERROR("The curve data file could not be opened: " << curve_data_file);
	}

	// read curves
	Curves curve_data;
	curve_data.reserve(curve_filenames.size());

	for (auto const& curve_filename: curve_filenames) {
		std::ifstream curve_file(base_path + "/" + curve_filename);
		if (curve_file.is_open()) {
			curve_data.emplace_back();
			readCurve(curve_file, curve_data.back(), header_size);
			curve_data.back().filename = curve_filename;

			if (curve_data.back().empty()) { curve_data.pop_back(); }
		}
		else {
			ERROR("A curve file could not be opened: " << curve_filename);
		}
	}

	return curve_data;
}


void exportSubcurve(std::string const& filename, Curve curve, ParamPoint s, ParamPoint t){

    std::ofstream file(filename);
    if (!file.is_open()) {
        ERROR("Couldn't open clustering file: " << filename);
    }
    for(auto coord:curve.eval(s)){
        file << coord << " ";
    }
    file << "\n";

    int step = (s>t)?-1:1;

    for(int i = s.id+std::max(step,0);(s>t)?i>t.id:i<=t.id;i+=step){
        for(auto coord:curve[i]){
            file << coord << " ";
        }
        file << "\n";
    }
    for(auto coord:curve.eval(t)){
        file << coord << " ";
    }
    file << "\n";
}

void exportSubcurve(std::string const& filename, Curve curve, ParamPoint s, ParamPoint t, int interpol){

    std::ofstream file(filename);
    if (!file.is_open()) {
        ERROR("Couldn't open clustering file: " << filename);
    }
    for(int i=0;i<=interpol;++i){
        for(auto coord:curve.eval(s.interpol(t,i,interpol))){
            file << coord << " ";
        }
        file << "\n";
    }
}

void exportCurve(std::string const& curve_data_file,Curve c){

    std::ofstream file(curve_data_file);
    if (!file.is_open()) {
        ERROR("Couldn't open clustering file: " << curve_data_file);
    }
    for(auto point:c){
        for(auto coord:point){
            file << coord << " ";
        }
        file << "\n";
    }
}

void exportClustering(std::string const& filename, std::string const& base_path, Clustering const& clustering, Curves const& curves)
{
	std::ofstream file(filename);
	if (!file.is_open()) {
		ERROR("Couldn't open clustering file: " << filename);
	}

	std::size_t clustering_size = 0;
	for (auto const& cluster: clustering) {
		if (!cluster.curve_ids.empty()) { ++clustering_size; }
	}

	// header
	file << clustering_size << "\n";
	file << base_path << "\n";

	// clusters
	file.precision(10);
	for (auto const& cluster: clustering) {
		if (cluster.curve_ids.empty()) { continue; }

		// write center
		for (auto const& point: cluster.center_curve) {
            for(auto const& coord : point)
			    file << coord << " ";
            file << "\n";
		}
		file << "-------------------\n";

		// write cluster curves
		/*for (auto curve_id: cluster.curve_ids) {
			file << curves[curve_id].filename << " ";
		}*/
		file << "\n";
	}
}

void exportCenters(std::string const& filename, Clustering const& clustering)
{
	std::ofstream file(filename);
	if (!file.is_open()) {
		ERROR("Couldn't open centers file: " << filename);
	}

	// header
	file << clustering.size() << "\n";
	file << "\n"; // write empty base path

	// clusters
	file.precision(10);
	for (auto const& cluster: clustering) {
		// write center
		for (auto const& point: cluster.center_curve) {
            for(auto const& coord : point)
                file << coord << " ";
		}
		file << "\n";
		file << "\n"; // write emtpy curve ids string
	}

}

void exportCentersGPX(std::string const& filename, Clustering const& clustering)
{
	std::ofstream file(filename);
	if (!file.is_open()) {
		ERROR("Couldn't open clustering file: " << filename);
	}

	std::size_t clustering_size = 0;
	for (auto const& cluster: clustering) {
		if (!cluster.curve_ids.empty()) { ++clustering_size; }
	}

	// header
	file << "<?xml version=\"1.0\" encoding=\"utf-8\" standalone=\"yes\"?>\n";
	file << "<gpx version=\"1.1\">\n";

	// clusters
	file.precision(10);
	for (auto const& cluster: clustering) {
		if (cluster.curve_ids.empty()) { continue; }

		// write center
		file << "<trk>\n";
		file << "<trkseg>\n";
		for (auto const& point: cluster.center_curve) {
			file << "<trkpt lat=\"" << point.y() << "\" lon=\"" << point.x() << "\"></trkpt>\n";
		}
		file << "</trkseg>\n";
		file << "</trk>\n";
	}

	// footer
	file << "</gpx>\n";
}

void exportCurvesGPX(std::string const& filename, Curves const& curves)
{
	std::ofstream file(filename);
	if (!file.is_open()) {
		ERROR("Couldn't open file: " << filename);
	}

	// header
	file << "<?xml version=\"1.0\" encoding=\"utf-8\" standalone=\"yes\"?>\n";
	file << "<gpx version=\"1.1\">\n";

	// clusters
	file.precision(10);
	for (auto const& curve: curves) {
		file << "<rte>\n";
		for (auto const& point: curve) {
			file << "<rtept lat=\"" << point.y() << "\" lon=\"" << point.x() << "\"></rtept>\n";
		}
		file << "</rte>\n";
	}

	// footer
	file << "</gpx>\n";
}

} // end namespace io
