//
// Created by Jacobus Conradi on 17.04.23.
//

#include <fstream>
#include "Curve.h"
#include "geometry_types.h"

Curve::Curve(const Points& points)
        : points(points), prefix_length(points.size()), weights(points.size(),0.5f)
{
    if (points.empty()) { return; }

    auto const& front = points.front();
    extreme_points = {front,front};
    prefix_length[0] = 0;

    for (PointID i = 1; i < points.size(); ++i)
    {
        auto segment_distance = points[i - 1].dist(points[i]);
        prefix_length[i] = prefix_length[i - 1] + segment_distance;

        extreme_points.min = extreme_points.min.min(points[i]);
        extreme_points.max = extreme_points.max.max(points[i]);
    }
}

void Curve::push_back(Point const& point, double w)
{
    if(size() == 0)
        extreme_points = {point,point};
    if (prefix_length.size()) {
        auto segment_distance = points.back().dist(point);
        prefix_length.push_back(prefix_length.back() + segment_distance);
    }
    else {
        prefix_length.push_back(0);
    }

    extreme_points.min = extreme_points.min.min(point);
    extreme_points.max = extreme_points.max.max(point);

    points.push_back(point*1.0);
    weights.push_back(w);
}

auto Curve::getExtremePoints() const -> ExtremePoints const&
{
    return extreme_points;
}

distance_t Curve::getUpperBoundDistance(Curve const& other) const
{
    auto const& extreme1 = this->getExtremePoints();
    auto const& extreme2 = other.getExtremePoints();

    Point min_point = extreme1.min.min(extreme2.min);
    Point max_point = extreme1.max.max(extreme2.max);

    return min_point.dist(max_point);
}

Curve::Curve(std::string _filename,dimensions_t d) {
    name = _filename;
    std::ifstream file(name);
    std::string line;
    while(file){
        double number;
        Point p(0);
        for (int i=0;i<d;i++){
            file >> number;
            p.push_back(number);
        }
        //somehow the last point is added 96 times without this... very weird
        if(file) {
            push_back(p);
            weights.push_back(0.5);
        }
    }

    prefix_length = std::vector<distance_t>(points.size());
    prefix_length[0] = 0;
    for(int i=1;i<points.size();++i){
        prefix_length[i] = prefix_length[i-1] + (points[i-1].dist(points[i]));
    }
}

void Curve::clear() {
    points.clear();
    weights.clear();
    prefix_length.clear();
}

Curve::Curve(const Curve& c, std::vector<int> times){
    for(int startindex = 0; startindex<c.size()-1;++startindex){
        Point s = c[startindex];
        Point t = c[startindex+1];
        int stime = times[startindex];
        int ttime = times[startindex+1];
        for(int i=stime;i<ttime;++i){
            push_back(s + (t-s)*((double)(i-stime)/(double)(ttime-stime)));
            weights.push_back(0.5);
        }
    }
    push_back(c.back());
    weights.push_back(0.5);
}

Point Curve::eval(CPoint t) {
    Point ps = operator[](t.getPoint());
    Point pt = operator[](t.getPoint()+1);
    return ps + (pt - ps)*t.getFraction();
}

Curve Curve::constructSubcurve(CPoint s, CPoint t) {
    Curve result;

    result.push_back(eval(s), evalWeight(s));

    int step = (s>t)?-1:1;

    for(PointID i = s.getPoint()+std::max(step,0);(s>t)?i>t.getPoint():i<=t.getPoint();i+=step){
        result.push_back(operator[](i),weights[i]);
    }
    result.push_back(eval(t), evalWeight(t));
    return result;
}

Curve::Curve(const py::array_t<distance_t> &in, const std::string &name)
    : points(in.request().shape[0], in.request().ndim > 1 ? in.request().shape[1] : 1),
    name{name},prefix_length(in.request().shape[0]),weights(in.request().shape[0],0.5f) {
    const auto array_dim = in.ndim();

    if (array_dim > 2){
        std::cerr << "A Curve requires a 1- or 2-dimensional numpy array of type " << typeid(distance_t).name() << "." << std::endl;
        std::cerr << "Current dimensions: " << array_dim << std::endl;
        std::cerr << "WARNING: constructed empty curve" << std::endl;
        return;
    }

    if (array_dim == 2) {
#pragma omp parallel for simd default(none) shared(in)
        for (auto i = 0; i < size(); ++i)
            for (dimensions_t j = 0; j < dimensions(); ++j) {
                points[i][j] = *in.data(i, j);
            }

    } else {
#pragma omp parallel for simd default(none) shared(in)
        for (auto i = 0; i < size(); ++i) {
            points[i][0] = *in.data(i);
        }
    }

    if (empty()) {
        std::cerr << "WARNING: constructed empty curve" << std::endl;
        return;
    }

    auto const& front = points.front();
    extreme_points = {front,front};
    prefix_length[0] = 0;

    for (PointID i = 1; i < points.size(); ++i)
    {
        auto segment_distance = points[i - 1].dist(points[i]);
        prefix_length[i] = prefix_length[i - 1] + segment_distance;

        extreme_points.min = extreme_points.min.min(points[i]);
        extreme_points.max = extreme_points.max.max(points[i]);
    }
}


std::string Curve::repr() const {
    std::stringstream ss;
    ss << "Curve '" << name << "' of complexity " << size() << std::flush;
    return ss.str();
}

//std::string Curves::repr() const {
//    std::stringstream ss;
//    ss << "Curves collection with " << size() << " curves" << std::flush;
//    return ss.str();
//}

std::string Curve::str() const {
    std::stringstream ss;
    ss << name << std::endl;
    for(auto p : points){
        ss << p.str() << std::endl;
    }
    return ss.str();
}

//std::string Curves::str() const {
//    std::stringstream ss;
//    ss << *this;
//    return ss.str();
//}
