
#ifndef CLUSTERING_CURVESIMPLIFICATION_H
#define CLUSTERING_CURVESIMPLIFICATION_H


#include "Curve.h"

class CurveSimplifier {
private:
    distance_t delta;
    int samplingRate = std::numeric_limits<int>::max();
    int curveLengthRule = 0;
    std::vector<int> times;num_threads
    ParamLabeling outputLabelings;
    static bool

    _less_than_or_equal_edge(const distance_t dist, const Point &p1, const Point &p2, const Curve &c,int start,int end) {
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
public:
    explicit CurveSimplifier(distance_t _delta){delta = _delta;};
    CurveSimplifier(CurveSimplifier const &simplifier) =default;
    void setSubsampling(int rate){samplingRate = rate;};
    void unsetSubsampling(){samplingRate = std::numeric_limits<int>::max();};
    void setDelta(distance_t _delta){delta = _delta;};
    std::vector<int> getTimes(){
        return times;
    }

    Curve simplify(Curve& c, FrameLabeling* inputLabelings=nullptr){
        assert(inputLabelings == nullptr || (inputLabelings->back().second==c.size()));
        times.clear();
        times.push_back(0);
        int length = 0;
        for (int i=1;i<c.size();++i){
            bool cantSkip = false;
            while(length > 0 && i - times[length-1] < samplingRate && _less_than_or_equal_edge(2*delta,c[times[length-1]],c[i],c,times[length-1],i)){
                times.pop_back();
                length -= 1;
            }
            if((length == 0 && i-times[length] >= samplingRate) || (length > 0 && i - times[length-1] >= samplingRate)){
                cantSkip = true;
            }
            if(c[times[length]].dist_sqr(c[i])>=delta*delta/9.0 || cantSkip){
                times.push_back(i);
                length += 1;
            }
        }
        if(length == 0){
            times.push_back(c.size()-1);
            length +=1;
        }
        Curve simplification;
        simplification.set_name(c.get_name());

        std::vector<double> weights;
        //TODO incorporate rule
        for (int i=0;i<=length;++i) {
            simplification.push_back(c[times[i]]);
            simplification.prefix_length.pop_back();
            //simplification.prefix_length.push_back(i);
            //simplification.prefix_length.push_back(times[i]);
            simplification.prefix_length.push_back(c.prefix_length[times[i]]);
            weights.push_back(c.weight(times[i]));
        }
        simplification.assignWeights(weights);
        if(inputLabelings != nullptr) {
            int gt_counter = 0;
            std::pair<Label, int> &cur = (*inputLabelings)[gt_counter];
            for (int i = 0; i <= length; ++i) {
                while (i > 0 && cur.second - 1 > times[i - 1] && cur.second - 1 <= times[i]) {
                    outputLabelings.push_back(
                            {cur.first, {i - 1, (((double) (cur.second - 1 - times[i - 1])) / (times[i] - times[i - 1]))}});
                    gt_counter++;
                    if (gt_counter < inputLabelings->size())
                        cur = (*inputLabelings)[gt_counter];
                    else {
                        assert(i == length);
                        break;
                    }
                }
            }
        }
        return simplification;
    }
    //TODO setepsilon?

    std::vector<std::pair<Label, CPoint>> &getSimplifiedGTs() {
        return outputLabelings;
    }
};


#endif //CLUSTERING_CURVESIMPLIFICATION_H
