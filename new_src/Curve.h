//
// Created by Jacobus Conradi on 17.04.23.
//

#ifndef CLUSTERING_CURVE_H
#define CLUSTERING_CURVE_H

#include "geometry_types.h"


class Curve
{
public:
    Curve() = default;
    Curve(const Points& points);
    Curve(std::string,dimensions_t);
    Curve(const Curve& c, std::vector<int> times);

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

    std::string filename;

    struct ExtremePoints { Point min, max; };
    ExtremePoints const& getExtremePoints() const;
    distance_t getUpperBoundDistance(Curve const& other) const;
    Curve constructSubcurve(CPoint s, CPoint t);


    std::vector<distance_t> prefix_length;

    void assignWeights(std::vector<double> ws){
        assert(ws.size() == size());
        weights = ws;
    };

    template <typename func> void assignWeights(func f){
        weights.clear();
        for(auto p : points){
            weights.push_back(f(p));
        }
    }

    double weight(int i){
        return weights[i];
    }

private:
    Points points;
    ExtremePoints extreme_points;
    std::vector<double> weights;
    double evalWeight(CPoint t) {
        double ws = weights[t.getPoint()];
        double wt = weights[t.getPoint() + 1];
        return ws + (wt - ws) * t.getFraction();
    };




};
using Curves = std::vector<Curve>;


#endif //CLUSTERING_CURVE_H
