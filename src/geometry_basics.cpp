#pragma once

#include <omp.h>
#include "geometry_basics.h"

namespace
{

template<typename T>
T pow2(T d) { return std::pow(d, 2); }

} // end anonymous namespace

//
// intersection_interval
//

inline bool IntersectionAlgorithm::smallDistanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center, distance_t radius_sqr) {
    double val = circle_center.dist_sqr(line_start * (1. - interpolate) + line_end * interpolate);
    return val <= radius_sqr;
}
// assume circle_center = p_1 + (p_2 - p_1) * t
// d1 = t*t*(p_2 - p_1).len_sqr = c3
// d2 = t*(line_start - p_1).dot( (p_2 - p_1) ) = c2
// d3 = (line_start - p_1).len_sqr = c1
// d4 = (p_1 - line_start).dot( (line_end - line_start) ) = -b1
// d5 = t*(p_2 - p_1).dot( (line_end - line_start) ) = b2
// d6 = (line_end - line_start).len_sqr = a
bool IntersectionAlgorithm::smallDistanceAt(distance_t hatt, double d1, double d2,double d3, double d4, double d5, double d6, distance_t radius_sqr) {
    //return a -2*interpolate*b + interpolate*interpolate*c <= radius_sqr;
    double val = d1 - 2*d2 + d3 -2*hatt*(d4 + d5) + hatt*hatt*d6;
    return val <= radius_sqr;
}

inline distance_t IntersectionAlgorithm::distanceAt(distance_t interpolate, Point line_start, Point line_end, Point circle_center) {
	return circle_center.dist_sqr(line_start * (1. - interpolate) + line_end * interpolate);
}

Interval IntersectionAlgorithm::intersection_interval(Point circle_center, distance_t radius, Point line_start, Point line_end, Interval * outer /* = nullptr*/)
{
    // The line can be represented as line_start + lambda * v
    //const Point v = line_end - line_start;
	const distance_t rad_sqr = radius * radius;

    const Point u = line_end - line_start, v = line_start - circle_center;
    const distance_t ulen_sqr = u.length_sqr(), vlen_sqr = v.length_sqr();

    // Find points p = line_start + lambda * v with
    //     dist(p, circle_center) = radius
    // <=> sqrt(p.x^2 + p.y^2) = radius
    // <=> p.x^2 + p.y^2 = radius^2
    // <=> (line_start.x + lambda * v.x)^2 + (line_start.y + lambda * v.y)^2 = radius^2
    // <=> (line_start.x^2 + 2 * line_start.x * lambda * v.x + lambda^2 * v.x^2) + (line_start.y^2 + 2 * line_start.y * lambda * v.y + lambda^2 * v.y^2) = radius^2
    // <=> lambda^2 * (v.x^2 + v.y^2) + lambda * (2 * line_start.x * v.x + 2 * line_start.y * v.y) + line_start.x^2 + line_start.y^2) - radius^2 = 0
    // let a := v.x^2 + v.y^2, 
	// let b := line_start.x * v.x + line_start.y * v.y, 
	// let c := line_start.x^2 + line_start.y^2 - radius^2
    // <=> lambda^2 * a + lambda * 2 b + c = 0
    // <=> lambda^2 + (2 b / a) * lambda + (c / a) = 0
    // <=> lambda1/2 = - (b / a) +/- sqrt((b / a)^2 - c / a)
	
    const distance_t a = ulen_sqr;
    const distance_t b = v*u;
    const distance_t c = vlen_sqr - rad_sqr;

	distance_t mid = - b / a;
    distance_t discriminant = pow2(mid) - c / a;

	const bool smallDistAtZero = smallDistanceAt(0., line_start, line_end, circle_center, rad_sqr);
	const bool smallDistAtOne = smallDistanceAt(1., line_start, line_end, circle_center, rad_sqr);
	bool smallDistAtMid = smallDistanceAt(mid, line_start, line_end, circle_center, rad_sqr);
	
	if (smallDistAtZero && smallDistAtOne) {
		if (outer != nullptr) { *outer = Interval(-eps, 1. + eps); }
		return Interval(0, 1);
	}
	
	if (!smallDistAtMid && smallDistAtZero) {
		mid = 0.;
		smallDistAtMid = true;
	}
	else if (!smallDistAtMid && smallDistAtOne) {
		mid = 1.;
		smallDistAtMid = true;
	}
	
	// Here we need the guarantee that if the free interval has length at least eps
	// then at mid the distance is <=radius
	// This is an assumption about the precision of distance_t computations
	// All remaining rules are free of such assumptions! 
	// (except for trivial ones like this: x + y and x - y have distance at most 2y up to negligible error)
    if (!smallDistAtMid) {
		if (outer != nullptr) { *outer = Interval(); }
		return Interval(); // no intersection;
    }
	
	if (mid <= 0. and !smallDistAtZero) {
		if (outer != nullptr) { *outer = Interval(); }
		return Interval();
	}
	if (mid >= 1. and !smallDistAtOne) {
		if (outer != nullptr) { *outer = Interval(); }
		return Interval();
	}
	
	discriminant = std::max<distance_t>(discriminant, 0.);
	distance_t sqrt_discr = 0.;
	bool sqrt_discr_computed = false;
	distance_t begin, end;
	
	if (smallDistAtZero) {
		begin = 0.;
		if (outer != nullptr) { outer->begin = -eps; }
	}
	else {
		sqrt_discr = std::sqrt(discriminant);
		sqrt_discr_computed = true;
		
		const distance_t lambda1 = mid - sqrt_discr;
		const distance_t innershift = std::min<distance_t>(lambda1 + save_eps_half, std::min<distance_t>(1., mid));
		const distance_t outershift = lambda1 - save_eps_half;
		if (innershift >= outershift && smallDistanceAt(innershift, line_start, line_end, circle_center, rad_sqr) && !smallDistanceAt(outershift, line_start, line_end, circle_center, rad_sqr)) {
			begin = innershift;
			if (outer != nullptr) { outer->begin = outershift; }
		}
		else {
			distance_t left = 0., right = std::min<distance_t>(mid, 1.);
			// invariants throughout binary search:
			//  * !smallDistanceAt(left)
			//  * smallDistanceAt(right)
			//  * 0 <= left <= right <= min(mid,1)
			// Clearly this is stays true after an iteration.
			// Why is it true in the beginning?
			// If smallDistanceAt(0.) then begin would already be set (fourth rule).
			// If !smallDistanceAt(right), then either !smallDistanceAt(mid), contradicting the very first rule, 
			//  or mid >= 1. and smallDistanceAt(1.), contradicting the third rule.
			// Finally, since !smallDistanceAt(left) we cannot have mid <= 0 by the second rule. Thus, right = min(mid,1) >= 0. = left
			while (right - left > save_eps) {
				distance_t m = 0.5 * (left + right);
				if (smallDistanceAt(m, line_start, line_end, circle_center, rad_sqr)) { right = m; }
				else { left = m; }
			}
			begin = right;
			if (outer != nullptr) { outer->begin = left; }
		}
	}
	
	if (smallDistAtOne) {
		end = 1.;
		if (outer != nullptr) { outer->end = 1. + eps; }
	}
	else {
		if (!sqrt_discr_computed) {
			sqrt_discr = std::sqrt(discriminant);
		}
		
		const distance_t lambda2 = mid + sqrt_discr;
		const distance_t innershift = std::max<distance_t>(lambda2 - save_eps_half, std::max<distance_t>(0., mid));
		const distance_t outershift = lambda2 + save_eps_half;
		if (innershift <= outershift && smallDistanceAt(innershift, line_start, line_end, circle_center, rad_sqr) && !smallDistanceAt(outershift, line_start, line_end, circle_center, rad_sqr)) {
			end = innershift;
			if (outer != nullptr) { outer->end = outershift; }
		}
		else {
			distance_t left = std::max<distance_t>(mid, 0.), right = 1.;
			// invariants throughout binary search:
			//  * smallDistanceAt(left)
			//  * !smallDistanceAt(right)
			//  * max(mid,0) <= left <= right <= 1
			while (right - left > save_eps) {
				distance_t m = 0.5 * (left + right);
				if (smallDistanceAt(m, line_start, line_end, circle_center, rad_sqr)) { left = m; }
				else { right = m; }
			}
			end = left;
			if (outer != nullptr) { outer->end = right; }
		}
	}
	
	assert(smallDistanceAt(begin, line_start, line_end, circle_center, rad_sqr));
	assert(smallDistanceAt(end, line_start, line_end, circle_center, rad_sqr));
	assert(0. <= begin && begin <= end && end <= 1.);
	
	assert(outer == nullptr || outer->begin < 0. || !smallDistanceAt(outer->begin, line_start, line_end, circle_center, rad_sqr));
	assert(outer == nullptr || outer->end > 1. || !smallDistanceAt(outer->end, line_start, line_end, circle_center, rad_sqr));
	assert(outer == nullptr || (outer->begin < begin && begin - outer->begin <= eps));
	assert(outer == nullptr || (outer->end > end && outer->end - end <= eps));
	
	assert(outer == nullptr || (outer->begin <= 1.));
	assert(outer == nullptr || (outer->end >= 0.));

    //Point p = line_start + (line_end-line_start)*begin;
	
	// These asssertions fail - use exact arithmetic to make it work??
	//assert(begin - eps < 0. || !smallDistanceAt(begin - eps, line_start, line_end, circle_center, rad_sqr));
	//assert(end + eps > 1. || !smallDistanceAt(end + eps, line_start, line_end, circle_center, rad_sqr));
	
	return Interval{ begin, end };
}

// originally: circle_center, radius, line_start, line_end
//
// assume: circle_center = p_1 + (p_2 - p_1)*t
//
// now:
// a = (line_end - line_start).len_sqr
// b1 = (line_start - p_1).dot( (line_end - line_start) )
// b2 = t*(p_2 - p_1).dot( (line_end - line_start) )
// c1 = (line_start - p_1).len_sqr
// c2 = t*(line_start - p_1).dot( (p_2 - p_1) )
// c3 = t*t*(p_2 - p_1).len_sqr
// d1 = t*t*(p_2 - p_1).len_sqr = c3
// d2 = t*(p_1 - line_start).dot( (p_2 - p_1) ) = c2
// d3 = (line_start - p_1).len_sqr = c1
// d4 = (p_1 - line_start).dot( (line_end - line_start) ) = -b1
// d5 = t*(p_2 - p_1).dot( (line_end - line_start) ) = b2
// d6 = (line_end - line_start).len_sqr = a
Interval IntersectionAlgorithm::intersection_interval(double a1, double b1, double b2, double c1, double c2, double c3, double r, Interval *outer) {
    distance_t rr = r*r;

    double a = a1;
    double b = b1 - b2;
    double c = c1 - 2*c2 + c3 - rr;
    double d1 = c3;
    double d2 = c2;
    double d3 = c1;
    double d4 = -b1;
    double d5 = b2;
    double d6 = a;


    distance_t mid = - b / a;
    distance_t discriminant = pow2(mid) - c / a;

    const bool smallDistAtZero = smallDistanceAt(0.,d1,d2,d3,d4,d5,d6,rr);
    const bool smallDistAtOne = smallDistanceAt(1.,d1,d2,d3,d4,d5,d6,rr);
    bool smallDistAtMid = smallDistanceAt(mid,d1,d2,d3,d4,d5,d6,rr);

    if (smallDistAtZero && smallDistAtOne) {
        if (outer != nullptr) { *outer = Interval(-eps, 1. + eps); }
        return Interval(0, 1);
    }

    if (!smallDistAtMid && smallDistAtZero) {
        mid = 0.;
        smallDistAtMid = true;
    }
    else if (!smallDistAtMid && smallDistAtOne) {
        mid = 1.;
        smallDistAtMid = true;
    }

    // Here we need the guarantee that if the free interval has length at least eps
    // then at mid the distance is <=radius
    // This is an assumption about the precision of distance_t computations
    // All remaining rules are free of such assumptions!
    // (except for trivial ones like this: x + y and x - y have distance at most 2y up to negligible error)
    if (!smallDistAtMid) {
        if (outer != nullptr) { *outer = Interval(); }
        return Interval(); // no intersection;
    }

    if (mid <= 0. and !smallDistAtZero) {
        if (outer != nullptr) { *outer = Interval(); }
        return Interval();
    }
    if (mid >= 1. and !smallDistAtOne) {
        if (outer != nullptr) { *outer = Interval(); }
        return Interval();
    }

    discriminant = std::max<distance_t>(discriminant, 0.);
    distance_t sqrt_discr = 0.;
    bool sqrt_discr_computed = false;
    distance_t begin, end;

    if (smallDistAtZero) {
        begin = 0.;
        if (outer != nullptr) { outer->begin = -eps; }
    }
    else {
        sqrt_discr = std::sqrt(discriminant);
        sqrt_discr_computed = true;

        const distance_t lambda1 = mid - sqrt_discr;
        const distance_t innershift = std::min<distance_t>(lambda1 + save_eps_half, std::min<distance_t>(1., mid));
        const distance_t outershift = lambda1 - save_eps_half;
        if (innershift >= outershift && smallDistanceAt(innershift,d1,d2,d3,d4,d5,d6,rr) && !smallDistanceAt(outershift,d1,d2,d3,d4,d5,d6,rr)) {
            begin = innershift;
            if (outer != nullptr) { outer->begin = outershift; }
        }
        else {
            distance_t left = 0., right = std::min<distance_t>(mid, 1.);
            // invariants throughout binary search:
            //  * !smallDistanceAt(left)
            //  * smallDistanceAt(right)
            //  * 0 <= left <= right <= min(mid,1)
            // Clearly this is stays true after an iteration.
            // Why is it true in the beginning?
            // If smallDistanceAt(0.) then begin would already be set (fourth rule).
            // If !smallDistanceAt(right), then either !smallDistanceAt(mid), contradicting the very first rule,
            //  or mid >= 1. and smallDistanceAt(1.), contradicting the third rule.
            // Finally, since !smallDistanceAt(left) we cannot have mid <= 0 by the second rule. Thus, right = min(mid,1) >= 0. = left
            while (right - left > save_eps) {
                distance_t m = 0.5 * (left + right);
                if (smallDistanceAt(m,d1,d2,d3,d4,d5,d6,rr)) { right = m; }
                else { left = m; }
            }
            begin = right;
            if (outer != nullptr) { outer->begin = left; }
        }
    }

    if (smallDistAtOne) {
        end = 1.;
        if (outer != nullptr) { outer->end = 1. + eps; }
    }
    else {
        if (!sqrt_discr_computed) {
            sqrt_discr = std::sqrt(discriminant);
        }

        const distance_t lambda2 = mid + sqrt_discr;
        const distance_t innershift = std::max<distance_t>(lambda2 - save_eps_half, std::max<distance_t>(0., mid));
        const distance_t outershift = lambda2 + save_eps_half;
        if (innershift <= outershift && smallDistanceAt(innershift,d1,d2,d3,d4,d5,d6,rr) && !smallDistanceAt(outershift,d1,d2,d3,d4,d5,d6,rr)) {
            end = innershift;
            if (outer != nullptr) { outer->end = outershift; }
        }
        else {
            distance_t left = std::max<distance_t>(mid, 0.), right = 1.;
            // invariants throughout binary search:
            //  * smallDistanceAt(left)
            //  * !smallDistanceAt(right)
            //  * max(mid,0) <= left <= right <= 1
            while (right - left > save_eps) {
                distance_t m = 0.5 * (left + right);
                if (smallDistanceAt(m,d1,d2,d3,d4,d5,d6,rr)) { left = m; }
                else { right = m; }
            }
            end = left;
            if (outer != nullptr) { outer->end = right; }
        }
    }

    assert(smallDistanceAt(begin,d1,d2,d3,d4,d5,d6,rr));
    assert(smallDistanceAt(end,d1,d2,d3,d4,d5,d6,rr));
    assert(0. <= begin && begin <= end && end <= 1.);

    assert(outer == nullptr || outer->begin < 0. || !smallDistanceAt(outer->begin,d1,d2,d3,d4,d5,d6,rr));
    assert(outer == nullptr || outer->end > 1. || !smallDistanceAt(outer->end,d1,d2,d3,d4,d5,d6,rr));
    assert(outer == nullptr || (outer->begin < begin && begin - outer->begin <= eps));
    assert(outer == nullptr || (outer->end > end && outer->end - end <= eps));

    assert(outer == nullptr || (outer->begin <= 1.));
    assert(outer == nullptr || (outer->end >= 0.));

    //Point p = line_start + (line_end-line_start)*begin;

    // These asssertions fail - use exact arithmetic to make it work??
    //assert(begin - eps < 0. || !smallDistanceAt(begin - eps, line_start, line_end, circle_center, rad_sqr));
    //assert(end + eps > 1. || !smallDistanceAt(end + eps, line_start, line_end, circle_center, rad_sqr));

    return Interval{ begin, end };
}

Interval IntersectionAlgorithm::intersection_interval(const Point& s, const Point& t, distance_t r, const Point& u, const Point& v, Interval *outer) {
    assert(s!=t);
    Point n = t-s;
    Point vs = v-s;
    Point us = u-s;
    double pvs = (vs*n)/(n*n),pus = (us*n)/(n*n);
    Interval outerCylinderInterval;
    Interval cylinderInterval = intersection_interval(s-s,r,us - (n*pus),vs - (n*pvs),&outerCylinderInterval);
    Interval outerStartInterval;
    Interval startInterval = intersection_interval(s,r,u,v,&outerStartInterval);
    Interval outerEndInterval;
    Interval endInterval = intersection_interval(t,r,u,v,&outerEndInterval);

    Interval result;
    //verify start
    Point pCylinder = u + (v-u)*cylinderInterval.begin;
    if(intersection_interval(pCylinder,r,s,t).is_empty()) {
        result.begin = std::min(startInterval.begin, endInterval.begin);
        if(outer != nullptr)
            outer->begin = std::min(outerStartInterval.begin,outerEndInterval.begin);
    }else{
        result.begin = cylinderInterval.begin;
        if(outer != nullptr)
            outer->begin = outerCylinderInterval.begin;
    }


    Point qCylinder = u + (v-u)*cylinderInterval.end;
    if(intersection_interval(qCylinder,r,s,t).is_empty()) {
        result.end = std::max(startInterval.end, endInterval.end);
        if(outer != nullptr)
            outer->end = std::max(outerStartInterval.end,outerEndInterval.end);
    }else{
        result.end = cylinderInterval.end;
        if(outer != nullptr)
            outer->end = outerCylinderInterval.end;
    }

    return result;

}

Ellipse segmentsToEllipse(Point const& a1, Point const& b1, Point const& a2, Point const& b2, distance_t distance)
{
	Ellipse e;

	auto const pi = std::atan2(0, -1);

	// Check if segments are parallel
	auto dir1 = b1 - a1;
	auto dir2 = b2 - a2;
	dir1 /= (dir1.length());
	dir2 /= (dir2.length_sqr());
	if (std::abs(dir1*dir2) >= 0.999) {
		e.invalidate();
		return e;
	}

	// First assign coefficients of x^2, y^2, ...
	// Those are calculated by hand
    //   A = a1*a1 - 2 a1*b1 + b1*b1
    //   B = - 2 a1*a2 + 2 a1*b2 + 2 a2*b1 - 2 b1*b2
    //   C = a2*a2 - 2 a2*b2 + b2*b2
    //   D = 2 a1*b1 - 2 a1*b2 - 2 b1*b1 + 2 b1*b2
    //   E = - 2 a2*b1 + 2 a2*b2 + 2 b1*b2 - 2 b2*b2
    //   F = b1*b1 - 2 b1*b2 + 2 b2*b2 - dd

    //
	// OLD STUFF
    //

    auto a1a1 = a1*a1,a1a2 = a1*a2,a1b1 = a1*b1, a1b2 = a1*b2;
    auto a2a2 = a2*a2,a2b1 =a2*b1, a2b2 = a2*b2;
    auto b1b1 = b1*b1, b1b2 = b1*b2;
    auto b2b2 = b2*b2;

    auto A = a1a1 - 2*a1b1 + b1b1;
    auto B = -2*a1a2 + 2*a1b2 + 2*a2b1 - 2*b1b2;
    auto C = a2a2 - 2*a2b2 + b2b2;
    auto D=2*a1b1 - 2*a1b2 - 2*b1b1 + 2*b1b2;
    auto E = -2*a2b1 + 2*a2b2 + 2*b1b2 - 2*b2b2;
    auto F = b1b1 - 2*b1b2 + 2*b2b2 - pow2(distance);

    //auto A = pow2(a1.x) - 2*a1.x*b1.x + pow2(a1.y) - 2*a1.y*b1.y + pow2(b1.x) + pow2(b1.y);
    //auto B = -2*a1.x*a2.x + 2*a1.x*b2.x - 2*a1.y*a2.y + 2*a1.y*b2.y + 2*a2.x*b1.x + 2*a2.y*b1.y - 2*b1.x*b2.x - 2*b1.y*b2.y;
    //auto C = pow2(a2.x) - 2*a2.x*b2.x + pow2(a2.y) - 2*a2.y*b2.y + pow2(b2.x) + pow2(b2.y);
    //auto D = 2*a1.x*b1.x - 2*a1.x*b2.x + 2*a1.y*b1.y - 2*a1.y*b2.y - 2*pow2(b1.x) + 2*b1.x*b2.x - 2*pow2(b1.y) + 2*b1.y*b2.y;
    //auto E = -2*a2.x*b1.x + 2*a2.x*b2.x - 2*a2.y*b1.y + 2*a2.y*b2.y + 2*b1.x*b2.x + 2*b1.y*b2.y - 2*pow2(b2.x) - 2*pow2(b2.y);
    //auto F = pow2(b1.x) - 2*b1.x*b2.x + pow2(b1.y) - 2*b1.y*b2.y + pow2(b2.x) + pow2(b2.y) - pow2(distance);

	// This should not fail if they are not parallel
	assert(pow2(B) - 4*A*C <= 0.0);

	// Now convert those to the form of the Ellipse type
	// See: https://en.wikipedia.org/wiki/Ellipse#Canonical_form
	e.center.set(0,(2*C*D - B*E)/(pow2(B) - 4*A*C));
	e.center.set(1,(2*A*E - B*D)/(pow2(B) - 4*A*C));
	e.width = -std::sqrt(2*(A*pow2(E) + C*pow2(D) - B*D*E + (pow2(B) - 4*A*C)*F)*(A+C+std::sqrt(pow2(A-C) + pow2(B))))/(pow2(B) - 4*A*C);
	e.height = -std::sqrt(2*(A*pow2(E) + C*pow2(D) - B*D*E + (pow2(B) - 4*A*C)*F)*(A+C-std::sqrt(pow2(A-C) + pow2(B))))/(pow2(B) - 4*A*C);

	if (B == 0) {
		e.alpha = A < C ? 0 : pi/2;
	}
	else {
		e.alpha = std::atan((C-A-std::sqrt(pow2(A-C) + pow2(B)))/B);
	}
	e.alpha = e.alpha*360/(2*pi);

	return e;
}


Point Point::min(const Point &other) const  {
    if(size() != other.size())
        std::cout << size() << " " << other.size() << "\n";
    Point result(size());
#pragma omp simd
    for(dimensions_t i = 0; i< dimensions();++i){
        result[i] = std::min(operator[](i),other[i]);
    }
    return result;
}

