#ifndef CLUSTERING_CURVE_H
#define CLUSTERING_CURVE_H

#include "geometry_types.h"
#include <pybind11/pybind11.h>


class Curve{
public:
    Curve() = default;
    Curve(const Points& points);
    Curve(std::string,dimensions_t);
    Curve(const Curve& c, std::vector<int> times);
    Curve(const py::array_t<distance_t> &in, const std::string &name = "unnamed curve");


    Point eval(CPoint);

    dimensions_t dimensions(){
        assert(size() >0);
        return operator[](0).dimensions();
    }

    std::size_t size() const { return points.size(); }
    bool empty() const { return points.empty(); }
    Point const& operator[](PointID i) const { return points[i]; }
    bool operator==(Curve const& other) const {
        return std::equal(points.cbegin(), points.cend(), other.points.cbegin(), other.points.cend());
    }
    bool operator!=(Curve const& other) const {
        return !(*this == other);
    }
    Point interpolate_at(CPoint const& pt) const  {
        assert(pt.getFraction() >= 0. && pt.getFraction() <= 1.);
        assert((pt.getPoint() < points.size()-1 || (pt.getPoint() == points.size()-1 && pt.getFraction() == 0.)));
        return pt.getFraction() == 0. ? points[pt.getPoint()] : points[pt.getPoint()]*(1.-pt.getFraction()) + points[pt.getPoint()+1]*pt.getFraction();
    }
    distance_t curve_length(PointID i, PointID j) const
    { return prefix_length[j] - prefix_length[i]; }

    distance_t subcurve_length(CPoint s, CPoint t) const{
        if(s.getPoint() == t.getPoint()){
            distance_t length = ((t.getFraction() - s.getFraction())*curve_length(s.getPoint(),s.getPoint()+1));
            return length;
        }else{
            distance_t length = ( (1-s.getFraction())*curve_length(s.getPoint(), s.getPoint()+1)+
                                  curve_length(s.getPoint() + 1, t.getPoint()) +
                                  (t.getFraction()-0)*curve_length(t.getPoint(),t.getPoint()+1));
            distance_t l = length;
            return length;
        }
    }

    Point front() const { return points.front(); }
    Point back() const { return points.back(); }

    void push_back(Point const& point, double w = 0.5);

    void clear();

    Points::const_iterator begin() { return points.begin(); }
    Points::const_iterator end() { return points.end(); }
    Points::const_iterator begin() const { return points.cbegin(); }
    Points::const_iterator end() const { return points.cend(); }

    //std::string filename;

    struct ExtremePoints { Point min, max; };
    ExtremePoints const& getExtremePoints() const;
    distance_t getUpperBoundDistance(Curve const& other) const;
    Curve constructSubcurve(CPoint s, CPoint t);


    std::vector<distance_t> prefix_length;

    void assignWeights(std::vector<double> ws){
        assert(ws.size() == size());
        weights = ws;
        maxW = 0;
        for(auto w : weights){
            if(w>maxW){
                maxW = w;
            }
        }
    };

    template <typename func> void assignWeights(func f){
        weights.clear();
        for(auto p : points){
            weights.push_back(f(p));
        }
        maxW = 0;
        for(auto w : weights){
            if(w>maxW){
                maxW = w;
            }
        }
    }

    double weight(int i){
        return weights[i];
    }

    double maxWeight(){
        return maxW;
    }

    inline auto as_ndarray() const {
        py::list l;
        for (const Point &elem : *this) {
            l.append(elem.as_ndarray());
        }
        return py::array_t<distance_t>(l);
    }

    std::string str() const;

    std::string repr() const;

    std::string get_name() const{return name;};

    void set_name(std::string n) {name = n;};

private:
    double maxW;
    Points points;
    ExtremePoints extreme_points;
    std::vector<double> weights;
    std::string name;
    double evalWeight(CPoint t) {
        double ws = weights[t.getPoint()];
        double wt = weights[t.getPoint() + 1];
        return ws + (wt - ws) * t.getFraction();
    };




};
/*
class Curves {
public:
    std::vector<Curve> res;

    Curve& get(int i){
        return res[i];
    }

    Curve& operator[](int i){
        return res[i];
    }

    Curve& front(){
        return res.front();
    }

    Curve& back(){
        return res.back();
    }

    void clear(){
        res.clear();
    }

    void push_back(Curve c){
        res.push_back(c);
    }

    bool empty(){
        return res.empty();
    }

    void add(Curve& val){
        res.push_back(val);
    }

    int size(){
        return res.size();
    }

    inline auto as_ndarray() const {
        py::list l;
        for (const Curve &elem : res) {
            l.append(elem.as_ndarray());
        }
        return py::array_t<distance_t>(l);
    }
};
*/
 using Curves = std::vector<Curve>;
/*
class Curves : private std::vector<Curve>{
public:
    using Parent = std::vector<Curve>;

    using Parent::operator[];
    using Parent::Parent;
    using Parent::size;
    using Parent::push_back;
    using Parent::emplace_back;

    void add(Curve c){
        push_back(c);
    }

    inline auto as_ndarray() const {
        py::list l;
        for (const Curve &elem : *this) {
            l.append(elem.as_ndarray());
        }
        return py::array_t<distance_t>(l);
    }
};
*/

#endif //CLUSTERING_CURVE_H
