//
// Created by Jacobus Conradi on 17.04.23.
//

#ifndef CLUSTERING_GEOMETRY_TYPES_H
#define CLUSTERING_GEOMETRY_TYPES_H

#include <utility>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "id.h"
#include "defs.h"
#include "basic_types.h"

using distance_t = double;
using dimensions_t = int;

//
// Point
//

class Point : public std::vector<distance_t>{
public:
    inline Point() : std::vector<distance_t>(2){};
    inline Point(const dimensions_t d) : std::vector<distance_t>(d,0) {}
    inline Point(distance_t x, distance_t y){
        std::vector<distance_t>::push_back(x);
        std::vector<distance_t>::push_back(y);
    }

    inline Point(distance_t x, distance_t y,distance_t z){
        std::vector<distance_t>::push_back(x);
        std::vector<distance_t>::push_back(y);
        std::vector<distance_t>::push_back(z);
    }

    inline dimensions_t dimensions() const {
        return size();
    }

    inline distance_t get(const dimensions_t i) const {
        return std::vector<distance_t>::operator[](i);
    }

    inline void set(const dimensions_t i, distance_t val) {
        std::vector<distance_t>::operator[](i) = val;
    }

    inline distance_t x() const{
        assert(std::vector<distance_t>::size()>1);
        return std::vector<distance_t>::operator[](0);
    }

    inline distance_t y() const{
        assert(std::vector<distance_t>::size()>1);
        return std::vector<distance_t>::operator[](1);
    }

#pragma omp declare simd
    inline const distance_t& operator[](const dimensions_t i) const {
        return std::vector<distance_t>::operator[](i);
    }

#pragma omp declare simd
    inline distance_t& operator[](const dimensions_t i) {
        return std::vector<distance_t>::operator[](i);
    }

    inline Point& operator+=(const Point &point) {
#pragma omp simd
        for (dimensions_t i = 0; i < dimensions(); ++i){
            operator[](i) += point[i];
        }
        return *this;
    }

    inline Point& operator-=(const Point &point) {
#pragma omp simd
        for (dimensions_t i = 0; i < dimensions(); ++i){
            operator[](i) -= point[i];
        }
        return *this;
    }

    inline Point& operator/=(const distance_t distance) {
#pragma omp simd
        for (dimensions_t i = 0; i < dimensions(); ++i){
            operator[](i) /= distance;
        }
        return *this;
    }

    inline Point operator+(const Point &point) const {
        Point result = *this;
#pragma omp simd
        for (dimensions_t i = 0; i < dimensions(); ++i){
            result[i] += point[i];
        }
        return result;
    }

    inline Point operator-(const Point &point) const {
        Point result = *this;
#pragma omp simd
        for (dimensions_t i = 0; i < dimensions(); ++i){
            result[i] -= point[i];
        }
        return result;
    }

    template<typename T>
    inline Point operator*(const T mult) const {
        Point result = *this;
#pragma omp simd
        for (dimensions_t i = 0; i < dimensions(); ++i){
            result[i] *= mult;
        }
        return result;
    }

    inline distance_t operator*(const Point &p) const {
        assert(dimensions() == p.dimensions());
        distance_t result = 0;
//#pragma omp simd reduction(+: result)
        for (dimensions_t i = 0; i < dimensions(); ++i) {
            result += operator[](i) * p[i];
        }
        return result;
    }

    inline Point operator/(const distance_t dist) const {
        Point result = *this;
#pragma omp simd
        for (dimensions_t i = 0; i < dimensions(); ++i){
            result[i] /= dist;
        }
        return result;
    }

    inline distance_t dist_sqr(const Point &point) const {
        distance_t result = 0, temp;
#pragma omp simd private(temp) reduction(+: result)
        for (dimensions_t i = 0; i < dimensions(); ++i){
            temp = operator[](i) - point[i];
            result += temp * temp;
        }
        return result;
    }

    inline distance_t dist(const Point &point) const {
        return std::sqrt(dist_sqr(point));
    }

    inline distance_t dot(const Point& point) const{
        distance_t result = 0;
#pragma omp simd reduction(+: result)
        for(dimensions_t i=0;i<dimensions();++i){
            result += operator[](i) * point[i];
        }
        return result;
    }

    inline distance_t length_sqr() const {
        distance_t result = 0;
#pragma omp simd reduction(+: result)
        for (dimensions_t i = 0; i < dimensions(); ++i){
            result += operator[](i) * operator[](i);
        }
        return result;
    }

    inline distance_t length() const {
        return std::sqrt(length_sqr());
    }

    inline Point max(Point const& other) const {
        assert(size() == other.size());
        Point result(size());
#pragma omp simd
        for(dimensions_t i = 0; i< dimensions();++i){
            result[i] = std::max(operator[](i),other[i]);
        }
        return result;
    }

    inline Point min(const Point &other) const  {
        assert(size() == other.size());
        Point result(size());
#pragma omp simd
        for(dimensions_t i = 0; i< dimensions();++i){
            result[i] = std::min(operator[](i),other[i]);
        }
        return result;
    }

    inline bool operator==(Point const& other) const{
        int result = 0;
#pragma omp simd reduction(+: result)
        for (dimensions_t i = 0; i < dimensions(); ++i){
            result += operator[](i) != other[i];
        }
        return result==0;
    }
    inline bool operator!=(Point const& other) const{
        return !(operator==(other));
    }
};
using Points = std::vector<Point>;
using PointID = ID<Point>;

//
// Interval
//

struct Interval
{
    distance_t begin;
    distance_t end;

    Interval()
            : begin(1.),
              end(0.) {}

    Interval(distance_t begin, distance_t end)
            : begin(begin),
              end(end) {}

    bool is_empty() const { return begin > end; }
    bool intersects(Interval const& other) const
    {
        if (is_empty() || other.is_empty()) { return false; }

        return (other.begin >= begin && other.begin <= end) ||
               (other.end >= begin && other.end <= end) ||
               (other.begin <= begin && other.end >= end);
    }

    double clip(double t) const{
        assert(!is_empty());
        if (t<begin)
            return begin;
        if (t>end)
            return end;
        return t;
    }

    //intersects two intervals
    Interval clip(Interval t) const{
        if(is_empty()){
            return {begin,end};
        }
        return {clip(t.begin),clip(t.end)};
    }

    bool contains(distance_t t) const{
        return begin <= t && t <= end;
    }
};

//
// CPoint
//

class CPoint {
private:
    PointID point;
    distance_t fraction;

    void normalize() {
        assert(fraction >= 0. && fraction <= 1.);
        if (fraction == 1.) {
            fraction = 0.;
            ++point;
        }
    }
public:
    CPoint(PointID point, distance_t fraction)
            : point(point), fraction(fraction)
    {
        //normalize();
    }
    CPoint() : point(std::numeric_limits<PointID::IDType>::max()), fraction(0.) {}

    PointID getPoint() const { return point; }
    distance_t getFraction() const { return fraction; }
    distance_t convert() const { return (distance_t) point + fraction; }
    void setPoint(PointID point) { this->point = point; }
    void setFraction(distance_t frac) { fraction = frac; /*normalize();*/ }

    bool operator<(CPoint const& other) const {
        return point < other.point || (point == other.point && fraction < other.fraction);
    }
    bool operator<=(CPoint const& other) const {
        return point < other.point || (point == other.point && fraction <= other.fraction);
    }
    bool operator>(CPoint const& other) const {
        return point > other.point || (point == other.point && fraction > other.fraction);
    }
    bool operator>=(CPoint const& other) const {
        return point > other.point || (point == other.point && fraction >= other.fraction);
    }
    bool operator==(CPoint const& other) const {
        return point == other.point && fraction == other.fraction;
    }
    bool operator!=(CPoint const& other) const {
        return point != other.point or fraction != other.fraction;
    }
    bool operator<(PointID other) const {
        return point < other;
    }
    bool operator>(PointID other) const {
        return point > other || (point == other && fraction > 0.);
    }
    bool operator<=(PointID other) const {
        return point < other || (point == other && fraction == 0.);
    }
    bool operator>=(PointID other) const {
        return point >= other;
    }
    bool operator==(PointID other) const {
        return point == other && fraction == 0.;
    }
    bool operator!=(size_t other) const {
        return !(point == other);
    }
    /*
    CPoint operator+(distance_t other) const {
        assert(other <= 1.);
        PointID p = point; distance_t f = fraction + other;
        if (f > 1.) {
            ++p; f -= 1.;
        }
        return CPoint(p, f);
    }
    CPoint operator-(distance_t other) const {
        assert(other <= 1.);
        PointID p = point; distance_t f = fraction - other;
        if (f < 0.) {
            --p; f += 1.;
        }
        return CPoint(p, f);
    }
     */
    CPoint ceil() const {
        return fraction > 0 ? CPoint(point + 1, 0.) : CPoint(point, 0.);
    }
    CPoint floor() const {
        return CPoint(point, 0.);
    }
    std::string to_string() const {
        //return std::to_string( (double) point + fraction);
        std::stringstream stream;
        stream << std::fixed << std::setprecision(10) << (double) point + fraction;
        return stream.str();
    }

    CPoint interpol(CPoint _end, int numerator, int denominator) {
        double start = (double)(getPoint()) + getFraction();
        double end = (double)(_end.getPoint()) + _end.getFraction();
        double length = end - start;
        double interpoler = (length * numerator)/denominator;

        PointID t = getPoint();
        double f = getFraction() + interpoler;
        while(f <= -EPSILON){
            t -= 1;
            f += 1.0;
        }
        while(f >= 1.0+EPSILON){
            t += 1;
            f -= 1.0;
        }
        return {t,f};
    }

    friend std::ostream& operator<<(std::ostream& out, const CPoint& p);
};

//TODO: CPair

typedef std::vector<std::pair<Label,CPoint>> ParamLabeling;

struct CInterval;
using CIntervals = std::vector<CInterval>;
using CIntervalsID = ID<CIntervals>;
using CIntervalID = std::size_t;
using CIntervalIDs = std::vector<CIntervalID>;

using CPoints = std::vector<CPoint>;

using CPosition = std::array<CPoint, 2>;
using CPositions = std::vector<CPosition>;

using CurveID = std::size_t;
using CurveIDs = std::vector<CurveID>;

struct CInterval
{
    CPoint begin;
    CPoint end;

    CurveID fixed_curve = -1;

    CInterval(CPoint begin, CPoint end, CurveID fixed_curve)
            : begin(std::move(begin)), end(std::move(end)), fixed_curve(fixed_curve) {}

    CInterval()
            : begin(std::numeric_limits<PointID::IDType>::max(), 0.),
              end(std::numeric_limits<PointID::IDType>::lowest(), 0.) {}

    CInterval(CInterval const& other) = default;

    CInterval(CPoint  begin, CPoint  end)
            : begin(std::move(begin)), end(std::move(end)) {}

    CInterval(PointID point1, distance_t fraction1, PointID point2, distance_t fraction2)
            : begin(point1, fraction1), end(point2, fraction2) {}

    CInterval(PointID begin, PointID end)
            : begin(begin, 0.), end(end, 0.) {}


    bool operator<(CInterval const& other) const {
        return begin < other.begin || (begin == other.begin && end < other.end);
    }

    bool is_empty() const { return end < begin; }

    void make_empty() {
        begin = {std::numeric_limits<PointID::IDType>::max(), 0.};
        end = {std::numeric_limits<PointID::IDType>::lowest(), 0.};
    }
    void clamp(CPoint const& min, CPoint const& max) {
        begin = std::max(min, begin);
        end = std::min(max, end);
    }

    bool contains(CPoint const& t) const{
        return begin <= t && t <= end;
    }


    CPoint getBegin() const {
        return begin;
    }

    CurveID getCurveIndex() const {
        return fixed_curve;
    }

    CPoint getEnd() const {
        return end;
    }
};

class IntersectionAlgorithm
{
public:
    static constexpr distance_t eps = 1e-8;

    /*
     * Returns which section of the line segment from line_start to line_end is inside the circle given by circle_center and radius.
     * If the circle and line segment do not intersect, the result is the empty Interval (and outer is the empty Interval, too).
     * Otherwise the result is an interval [x,y] such that the distance at x and at y is at most the radius, i.e., [x,y] is a subset of the free interval.
     * The optional output "outer" is an interval strictly containing the free interval.
     * In other words, "outer" is an interval [x',y'] containing [x,y] such that x-x', y'-y <= eps and:
     * If x = 0 then x' = -eps, while if x > 0 then the distance at x' is more than the radius.
     * If y = 1 then y' = 1+eps, while if y < 1 then the distance at y' is more than the radius.
     */
    static Interval intersection_interval(Point circle_center, distance_t radius, Point line_start, Point line_end, Interval * outer = nullptr);
    //intersection on uv with pull around st
    static Interval intersection_interval(const Point& s, const Point& t, distance_t r, const Point& u, const Point& v, Interval * outer = nullptr);
    //fast underlying primitive
    static Interval intersection_interval(double a1, double b1, double b2, double c1, double c2, double c3, double r, Interval *outer = nullptr);
private:
    IntersectionAlgorithm() {} // Make class static-only
    static inline bool smallDistanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center, distance_t radius_sqr);
    static inline distance_t distanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center);

    static constexpr distance_t save_eps = 0.5 * eps;
    static constexpr distance_t save_eps_half = 0.25 * eps;

    static inline bool smallDistanceAt(distance_t hatt, double d1, double d2, double d3, double d4, double d5, double d6, distance_t radius_sqr);
};

#endif //CLUSTERING_GEOMETRY_TYPES_H
