
#include "center_clustering_algs.h"
bool cmpLeftLower(CInterval l, CInterval r){
    return (l.getCurveIndex() < r.getCurveIndex()) || (l.getCurveIndex() == r.getCurveIndex() && (l.getBegin() < r.getBegin()));
}

bool cmpLeftUpper(CInterval l, CInterval r){
    return (l.getCurveIndex() < r.getCurveIndex()) || (l.getCurveIndex() == r.getCurveIndex() && (l.getBegin() < r.getEnd()));
}

bool cmpRightUpper(CInterval l, CInterval r){
    return (l.getCurveIndex() < r.getCurveIndex()) || (l.getCurveIndex() == r.getCurveIndex() && (l.getEnd() < r.getEnd()));
}
//TODO: something here is fucked.
void updateCandidate(Curves& curves, Candidate& c, std::vector<CInterval> covering, const std::vector<double>& suffixLengths, int roundID){
    double newLength = 0;

    for(auto matching : c.matching){

        auto l = std::upper_bound(covering.begin(),covering.end(),matching,cmpLeftUpper);
        int lIdx = l - covering.begin();
        //lIdx points to the first interval that can lie right of s

        auto r = std::upper_bound(covering.begin(),covering.end(),matching, cmpRightUpper);
        int rIdx = r - covering.begin();
        //rIdx points to the first interval that can lie right of t

        bool leftContains = (l!=covering.end()) && (l->getCurveIndex() == matching.getCurveIndex()) && (l->contains(matching.getBegin()));
        bool rightContains = (r!=covering.end()) && (r->getCurveIndex() == matching.getCurveIndex()) && (r->contains(matching.getEnd()));

        if(leftContains && rightContains && lIdx == rIdx)
            continue;

        newLength += curves[matching.getCurveIndex()].subcurve_length(matching.getBegin(),matching.getEnd());

        if(leftContains){
            newLength -= curves[matching.getCurveIndex()].subcurve_length(matching.getBegin(),l->getEnd());
            lIdx ++;
        }
        //now lIdx points to the first interval that lies strictly to the right of s

        if(rightContains){
            newLength -= curves[matching.getCurveIndex()].subcurve_length(r->getBegin(),matching.getEnd());
            //rIdx --;
        }
        //now rIdx points to the last inerval that lies strictly to the left of t
        if(rIdx > lIdx)
            newLength -= (suffixLengths[rIdx] - suffixLengths[lIdx]);
    }
    if(c.semiUpdatedCoverLength + EPSILON < newLength){
        std::cout << "  ???" << c.semiUpdatedCoverLength << " -> " << newLength <<" with error " << newLength - c.semiUpdatedCoverLength <<"???"<<std::endl;
        //updateCandidate(curves,c,covering,suffixLengths,roundID);
    }
    c.semiUpdatedCoverLength = newLength;
    c.roundOfUpdate = roundID;
}

void printFirst50(CandidateSetPQ& cs){
    std::vector<Candidate> temp;
    for(int i=0;i<5;++i){
        Candidate c = cs.top();
        cs.pop();
        temp.push_back(c);
        std::cout << "(" << c.roundOfUpdate << ";" << c.semiUpdatedCoverLength << ") ";
    }
    std::cout << "\n";
    for(Candidate c : temp)
        cs.push(c);
}

double lengthOfUncovered(Curves curves, std::vector<Candidate> candidateSet){
    //step 1: merge
    std::vector<CInterval> presorted;
    for(auto c:candidateSet){
        presorted.insert(presorted.end(),c.matching.begin(),c.matching.end());
    }
    std::sort(presorted.begin(), presorted.end(), cmpLeftLower);


    double uncovered = 0;
    std::pair<int, CPoint> pcur = {0, {0, 0}};

    std::vector<CInterval> covering;

    if(!presorted.empty()) {
        CInterval cur = presorted[0];
        for (int covI = 1; covI < presorted.size(); covI++) {
            CInterval next = presorted[covI];
            if ((next.getCurveIndex() == cur.getCurveIndex()) && (next.getBegin() < cur.getEnd())) {
                cur.end = std::max(next.getEnd(), cur.getEnd());
            } else {
                covering.push_back(cur);
                cur = next;
            }
        }
        covering.push_back(cur);

        //now covering is a disjoint set of intervals
        for (auto cov: covering) {
            while (pcur.first < cov.getCurveIndex()) {
                uncovered += curves[pcur.first].subcurve_length(pcur.second,
                                                                {(int) (curves[pcur.first].size()) - 2, 1.0});
                pcur = {pcur.first + 1, {0, 0}};
            }
            if (cov.getCurveIndex() == pcur.first) {
                uncovered += std::max(0.0,curves[cov.getCurveIndex()].subcurve_length(pcur.second, cov.getBegin()));
                pcur = {pcur.first, cov.end};
                if((CPoint){(int)(curves[pcur.first].size()-2),1.0} <= pcur.second ){
                    pcur = {pcur.first + 1, {0, 0}};
                }
            }
        }
    }
    while(pcur.first < curves.size()){
        uncovered += curves[pcur.first].subcurve_length(pcur.second,{(int)(curves[pcur.first].size())-2,1.0});
        pcur = {pcur.first+1,{0,0}};
    }
    return uncovered;
}
/*
Curves greedyCoverAlreadySimplified(Curves& curves, double delta, int l, int max_rounds, bool show){
    std::vector<Candidate> bestResultVisualizer = greedyCoverUnsanitizedOutput(curves,delta,l,max_rounds,show,[=](const Candidate& a){return a.getEnd().getPoint() >a.getBegin().getPoint() + l/4;});
    Curves bestresult;

    for (auto &sub: bestResultVisualizer) {
        bestresult.push_back(
                curves[sub.getCurveIndex()].constructSubcurve(sub.getBegin(), sub.getEnd()));
    }

    std::cout << "\nImportances: ";
    for(auto c : bestResultVisualizer){
        std::cout << c.importance << " ";
    }
    std::cout << "\n";
    for(int count = 0;count < ((20<bestResultVisualizer.size())?20:bestResultVisualizer.size());count++) {
        Candidate bestCandidate = bestResultVisualizer[count];
        Candidate bC = bestCandidate;
        for (int i = 0; i < bC.visualMatching.size(); ++i) {
            auto matching = bC.visualMatching[i];
            io::exportSubcurve(
                    "/Users/styx/data/curveclustering/results/cluster/matching" + std::to_string(count) +"/interval" + std::to_string(i) + ".txt",
                    curves[matching.getCurveIndex()], matching.getBegin(), matching.end, 100);
        }
    }
    return bestresult;
}
 */