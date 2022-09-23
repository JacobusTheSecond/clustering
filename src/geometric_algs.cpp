#include "geometric_algs.h"

#include <algorithm>
#include <random>

//
// All the functions in this class work with squared radii to increase performance.
// Keep this in mind when reusing this code!
//

namespace
{

bool inCircle(Point const& p, Circle const& circle) {
	const distance_t eps = 0.00000000000001;
	return p.dist_sqr(circle.center) <= circle.radius + eps;
}

Circle circleFromTwo(Point p1, Point p2)
{
	auto center  = (p1 + p2)/2.;
	auto radius = p1.dist_sqr(center);
	return Circle(center, radius);
}

Circle calcCircle(Point p1, Point p2, Point p3)
{
	static constexpr distance_t epsilon = 1e-10;

	// compute normal vectors and their offset vectors
	auto q1 = (p1 + p2)/2.;
	auto q2 = (p2 + p3)/2.;
	auto n1 = Point((p2 - p1).y(), -1*(p2 - p1).x());
	auto n2 = Point((p3 - p2).y(), -1*(p3 - p2).x());

	// Avoid division by zero. Due to non-colinearity n1.x is non-zero after the swap.
	if (n1.x() == 0) {
		std::swap(n1, n2);
		std::swap(q1, q2);
	}
	assert(n1.x() != 0);

	// main calculation
	auto a = q1.y() - q2.y() + (q2.x() - q1.x())*n1.y()/n1.x();
	auto b = n2.y() - (n2.x()*n1.y())/n1.x();

	// if the points p1, p2, p3 are colinear
	if (std::abs(b) < epsilon) {
		// rather naive, but well...
		if (inCircle(p1, circleFromTwo(p2,p3))) { return circleFromTwo(p2,p3); }
		if (inCircle(p2, circleFromTwo(p1,p3))) { return circleFromTwo(p1,p3); }
		return circleFromTwo(p1,p2);
	}

	auto i2 = a/b;

	auto center = q2 + n2*i2;
	auto radius = center.dist_sqr(p1);
	return Circle(center, radius);
}

} // end anonymous namespace

// XXX: Theoretically we should random shuffle after all except the most inner
// current_circle update. However, I think practically this is slower and an
// initial random shuffle suffices. I didn't test this though.
Circle calcMinEnclosingCircle(Points points)
{
	if (points.size() == 0) {
		return Circle();
	}
	if (points.size() == 1) {
		return Circle(points.front(), 0.);
	}

	std::shuffle(points.begin(), points.end(), std::mt19937(std::random_device()()));

	auto current_circle = circleFromTwo(points[0], points[1]);
	for (PointID i = 2; i < points.size(); ++i) {
		auto p_i = points[i];
		if (inCircle(p_i, current_circle)) { continue; }

		// from here: p_i not in current_circle
		current_circle = circleFromTwo(p_i, points[0]);
		for (PointID j = 1; j < i; ++j) {
			auto p_j = points[j];
			if (inCircle(p_j, current_circle)) { continue; }

			// from here: p_j not in current_circle
			current_circle = circleFromTwo(p_i, p_j);
			for (PointID k = 0; k < j; ++k) {
				auto p_k = points[k];
				if (inCircle(p_k, current_circle)) { continue; }

				// from here: p_k not in current_circle
				current_circle = calcCircle(p_i, p_j, p_k);
			}
		}
	}

	current_circle.radius = std::sqrt(current_circle.radius);
	return current_circle;
}
