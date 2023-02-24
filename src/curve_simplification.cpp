#include "curve_simplification.h"

#include "defs.h"
#include "frechet_light.h"

namespace
{

Curve simplify(Curve const& curve, distance_t distance)
{
	Curve simplified_curve({curve.front()});
	Curve prefix_curve({curve.front()});

	FrechetLight frechet_light;
	for (PointID id = 1; id < curve.size()-1; ++id) {
		auto const& point = curve[id];
		prefix_curve.push_back(point);
		auto line_segment = Curve({prefix_curve.front(), prefix_curve.back()});
		bool less_than = frechet_light.lessThanWithFilters(distance, line_segment, prefix_curve);
		if (!less_than) {
			simplified_curve.push_back(point);
			prefix_curve = Curve({point});
		}
	}
	simplified_curve.push_back(curve.back());

	return simplified_curve;
}

} // end anonymous namespace

// triest to find an l-simplification with a small distance
Curve simplify(Curve const& curve, int l)
{
	static constexpr distance_t epsilon = 1e-8;

	assert(l >= 2);
	if ((int)curve.size() < l) { return curve; }

	distance_t min = 0.;
	distance_t max = curve.getUpperBoundDistance(curve);

	Curve simplified_curve = curve;
	while (max-min > epsilon) {
		auto split = (max + min)/2.;
		simplified_curve = simplify(curve, split);
		if ((int)simplified_curve.size() <= l) {
			max = split;
		}
		else {
			min = split;
		}
	}

	return simplify(curve, max);
}

bool _less_than_or_equal_edge(const distance_t dist, const Point &p1, const Point &p2, const Curve &c,int start,int end) {
    assert(start <= end && end < c.size() && 0<=start);
    if(p1.dist_sqr(c[start])>dist*dist || p2.dist_sqr(c[end])>dist*dist)
        return false;
    distance_t min=0;
    for (int i = start+1;i<end;++i){
        Interval freeInterval = IntersectionAlgorithm::intersection_interval(c[i],dist,p1,p2);
        if (freeInterval.is_empty() || freeInterval.begin > freeInterval.end)
            return false;
        min = std::max(min,freeInterval.begin);
    }
    return true;
}

Curve good_simplification(const Curve &c, const distance_t delta, std::vector<int>* indices) {
    std::vector<int> I;
    if(indices != nullptr)
        indices->push_back(0);
    I.push_back(0);
    int length = 0;
    for (int i=1;i<c.size();++i){
        while(length > 0 && _less_than_or_equal_edge(2*delta,c[I[length-1]],c[i],c,I[length-1],i)){
            I.pop_back();
            if(indices != nullptr)
                indices->pop_back();
            length -= 1;
        }
        if(c[I[length]].dist_sqr(c[i])>=delta*delta/9.0){
            I.push_back(i);
            if(indices != nullptr)
                indices->push_back(i);
            length += 1;
        }
    }
    Curve simplification;
    for (int i=0;i<=length;++i){
        simplification.push_back(c[I[i]]);
        simplification.prefix_length.pop_back();
        //simplification.prefix_length.push_back(i);
        //simplification.prefix_length.push_back(I[i]);
        simplification.prefix_length.push_back(c.prefix_length[I[i]]);
    }

    std::cout << "Complexity: " << simplification.size() << "\n";
    return simplification;
}

Curve good_simplification(const Curve &c, const distance_t delta, FrameLabeling & groundthruth, ParamLabeling * simplifiedLabels) {
    assert(groundthruth[groundthruth.size()-1].second == c.size());
    std::vector<int> I;
    I.push_back(0);
    int length = 0;
    for (int i=1;i<c.size();++i){
        while(length > 0 && _less_than_or_equal_edge(2*delta,c[I[length-1]],c[i],c,I[length-1],i)){
            I.pop_back();
            length -= 1;
        }
        if(c[I[length]].dist_sqr(c[i])>=delta*delta/9.0){
            I.push_back(i);
            length += 1;
        }
    }
    Curve simplification;
    int gt_counter = 0;
    std::pair<Label,int>& cur = groundthruth[gt_counter];
    for (int i=0;i<=length;++i){
        simplification.push_back(c[I[i]]);
        simplification.prefix_length.pop_back();
        //simplification.prefix_length.push_back(i);
        //simplification.prefix_length.push_back(I[i]);
        simplification.prefix_length.push_back(c.prefix_length[I[i]]);
        while(i > 0 && cur.second-1 > I[i-1] && cur.second-1 <= I[i]){
            simplifiedLabels->push_back({cur.first,{i-1,(((double)(cur.second-1-I[i-1]))/(I[i]-I[i-1]))}});
            gt_counter++;
            if(gt_counter < groundthruth.size())
                cur = groundthruth[gt_counter];
            else {
                assert(i == length);
                break;
            }
        }
    }

    std::cout << "Complexity: " << simplification.size() << "\n";
    return simplification;
}

Curve good_simplification(const Curve &c, const distance_t delta, FrameLabeling &groundthruth,
                          ParamLabeling *simplifiedLabels, int maxSkip) {

    assert(groundthruth[groundthruth.size()-1].second == c.size());
    std::vector<int> I;
    I.push_back(0);
    int length = 0;
    for (int i=1;i<c.size();++i){
        bool cantSkip = false;
        while(length > 0 && i - I[length-1] < maxSkip && _less_than_or_equal_edge(2*delta,c[I[length-1]],c[i],c,I[length-1],i)){
            I.pop_back();
            length -= 1;
        }
        if((length == 0 && i-I[length] >= maxSkip) || (length > 0 && i - I[length-1] >= maxSkip)){
            cantSkip = true;
        }
        if(c[I[length]].dist_sqr(c[i])>=delta*delta/9.0 || cantSkip){
            I.push_back(i);
            length += 1;
        }
    }
    Curve simplification;
    int gt_counter = 0;
    std::pair<Label,int>& cur = groundthruth[gt_counter];
    for (int i=0;i<=length;++i){
        simplification.push_back(c[I[i]]);
        simplification.prefix_length.pop_back();
        //simplification.prefix_length.push_back(i);
        //simplification.prefix_length.push_back(I[i]);
        simplification.prefix_length.push_back(c.prefix_length[I[i]]);
        while(i > 0 && cur.second-1 > I[i-1] && cur.second-1 <= I[i]){
            simplifiedLabels->push_back({cur.first,{i-1,(((double)(cur.second-1-I[i-1]))/(I[i]-I[i-1]))}});
            gt_counter++;
            if(gt_counter < groundthruth.size())
                cur = groundthruth[gt_counter];
            else {
                assert(i == length);
                break;
            }
        }
    }

    std::cout << "Complexity: " << simplification.size() << "\n";
    return simplification;
}

Curve good_simplification(const Curve &c, const distance_t delta, int maxSkip) {

    std::vector<int> I;
    I.push_back(0);
    int length = 0;
    for (int i=1;i<c.size();++i){
        bool cantSkip = false;
        while(length > 0 && i - I[length-1] < maxSkip && _less_than_or_equal_edge(2*delta,c[I[length-1]],c[i],c,I[length-1],i)){
            I.pop_back();
            length -= 1;
        }
        if((length == 0 && i-I[length] >= maxSkip) || (length > 0 && i - I[length-1] >= maxSkip)){
            cantSkip = true;
        }
        if(c[I[length]].dist_sqr(c[i])>=delta*delta/9.0 || cantSkip){
            I.push_back(i);
            length += 1;
        }
    }
    Curve simplification;
    for (int i=0;i<=length;++i){
        simplification.push_back(c[I[i]]);
        simplification.prefix_length.pop_back();
        //simplification.prefix_length.push_back(i);
        //simplification.prefix_length.push_back(I[i]);
        simplification.prefix_length.push_back(c.prefix_length[I[i]]);
    }

    std::cout << "Complexity: " << simplification.size() << "\n";
    return simplification;
}
