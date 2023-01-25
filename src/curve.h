#pragma once

#include "defs.h"
#include "geometry_basics.h"
#include "id.h"

typedef struct BreakpointCandidate{
    int startIdx;
    Point start;
    int endIdx;
    Point end;
} BreakpointCandidate;

// Represents a trajectory. Additionally to the points given in the input file,
// we also store the length of any prefix of the trajectory.
class Curve
{
public:
    Curve() = default;
    Curve(const Points& points);
    Curve(std::string,dimensions_t);
    Curve(const Curve& c, std::vector<int> times);

    Point eval(ParamPoint);

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

    distance_t subcurve_length(ParamPoint s, ParamPoint t) const{
        if(s.id == t.id){
            distance_t length = ((t.t - s.t)*curve_length(s.id,s.id+1));
            return length;
        }else{
            distance_t length = ( (1-s.t)*curve_length(s.id, s.id+1)+
                                  curve_length(s.id + 1, t.id) +
                                  (t.t-0)*curve_length(t.id,t.id+1));
            distance_t l = length;
            return length;
        }
    }

    Point front() const { return points.front(); }
    Point back() const { return points.back(); }

    void push_back(Point const& point);

    void clear();

	Points::const_iterator begin() { return points.begin(); }
	Points::const_iterator end() { return points.end(); }
	Points::const_iterator begin() const { return points.cbegin(); }
	Points::const_iterator end() const { return points.cend(); }
	
	std::string filename;

	struct ExtremePoints { Point min, max; };
	ExtremePoints const& getExtremePoints() const;
	distance_t getUpperBoundDistance(Curve const& other) const;
    Curve constructSubcurve(ParamPoint s, ParamPoint t);


    std::vector<distance_t> prefix_length;

private:
    Points points;
	ExtremePoints extreme_points;



};
using Curves = std::vector<Curve>;

//std::ostream& operator<<(std::ostream& out, const Curve& curve);
