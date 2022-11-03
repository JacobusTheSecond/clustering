#pragma once

#include "basic_types.h"
#include "curve.h"

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

} // end namespace io
