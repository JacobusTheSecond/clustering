#pragma once

#include "basic_types.h"
#include "curve.h"


//labelset for 86_1,86_2 and 86_4
enum Label {_labelstart,transition,walk,jump,punch,leg_kick,squat,run,stand,arm_up,drink,stretch,slap,turn,_labelend};

namespace io
{

Curves readCurves(std::string const& curve_data_file, int header_size = 0);
void exportClustering(std::string const& filename, std::string const& base_path, Clustering const& clustering, Curves const& curves);
void exportCurve(std::string const& curve_data_file,Curve);
void exportCenters(std::string const& filename, Clustering const& clustering);
void exportCentersGPX(std::string const& filename, Clustering const& clustering);
void exportCurvesGPX(std::string const& filename, Curves const& curves);
void exportSubcurve(std::string const& filename, Curve curve, ParamPoint s, ParamPoint t);
void exportSubcurve(std::string const& filename, Curve curve, ParamPoint s, ParamPoint t, int interpol);
std::vector<std::pair<Label, int>> readHM36GroundTruth(const std::string &ground_truth_data_file, int labels);

} // end namespace io
