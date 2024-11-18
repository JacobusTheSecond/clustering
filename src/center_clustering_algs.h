#pragma once

//#include "center_algs.h"
//#include "../src/clustering_algs.h"
#include "Curve.h"
#include "Candidate.h"
#include "io.h"
#include "CurveSimplification.h"
#include "SWatch.h"
#include "FreespaceVisualizer.h"
#include <algorithm>

/*
Clustering computeCenterClustering(
	Curves const& curves, int k, int l, ClusterAlg cluster_alg, CenterAlg center_alg, int max_rounds = 10);

Clustering approxCover(Curves& curves, double delta,int l,int max_rounds = 10);
*/

class Cluster{
public:
    Cluster(CInterval c, CIntervals& m){
        center = c;
        matching = m;
    }
    CInterval getCenter() const{return center;};
    CIntervals& getMatching(){return matching;};
    int size() const{return matching.size();};
    inline auto as_ndarray() const {
        py::list l;
        for (auto &elem : matching) {
            l.append(elem.as_ndarray());
        }
        return py::array_t<distance_t>(l);
    }

    CInterval const& operator[](PointID i) const { return matching[i]; }

    void merge(Cluster& other) {
        // merge cluster centers
        center = CInterval::Union(center, other.center);
        // merge matchings
        for (int i = 0; i < other.matching.size(); i++) {
            // try to insert other.matching[i]
            bool inserted = false;
            for (int j = 0; j < matching.size(); j++) {
                if (matching[j].intersects(other.matching[i])) {
                    matching[j] = CInterval::Union(matching[j], other.matching[i]);
                    inserted = true;
                    break;
                }
            }
            if (!inserted) {
                matching.push_back(other.matching[i]);
            }
        }
    }

private:
    CInterval center;
    CIntervals matching;
};

class ClusteringResult : public std::vector<Cluster>{
public:
    ClusteringResult(std::vector<Candidate>& cs){
        for(auto c:cs){
            CInterval center{c.getBegin(),c.getEnd(),c.getCurveIndex()};
            emplace_back(center,c.matching);
        }
    }
    Cluster const& get(PointID i) const { return operator[](i); }
    int len(){return (*this).size();}
};

Curves greedyCoverAlreadySimplified(Curves &curves, double delta, int l, int max_rounds = 10, bool show = false);

void
updateCandidate(Curves &curves, Candidate &c, std::vector<CInterval> covering, const std::vector<double> &suffixLengths,
                int roundID);

double lengthOfUncovered(Curves curves, std::vector<Candidate> candidateSet);

bool cmpLeftLower(CInterval l, CInterval r);

//template<typename func> Curves greedyCover(Curves& curves, double delta, int l, );

class CurveClusterer {
private:
    bool withSubsamplingRule;
    int subsamplingFrequency;

    bool withSort = false;

    //visual stuff
    bool showFreespaces;


    double simplificationDelta;
    double freespaceDelta;


    //TODO: make freespaces persistent and updatable
    //TODO: unsimplify candidate infos

    bool hasGTs;
    //std::vector<std::pair<Label,int>> groundTruth;

    std::vector<std::vector<int>> times;

    //correct Map to Base stuff
    std::vector<std::vector<CPoint>> vertexMaps;


public:

    double getSimpDelta() const{return simplificationDelta;};
    double getFreeDelta() const{return freespaceDelta;};

    //TODO: make these accessible in a better way, for filters...
    Curves simplifiedCurves;
    Curves unsimplifiedCurves;
    std::vector<ParamLabeling> simplifiedGTs;
    std::vector<int> simpIDtoOriginID;

    CurveClusterer(int samplingRate, bool showFlag) {
        updateFlags(samplingRate, showFlag);
    };

    CurveClusterer() {
        updateFlags(-1, false);
    };

    //TODO: this seems wrong
    CPoint mapSimplificationToBase(int curveIdx, CPoint q) {
        auto x = std::lower_bound(vertexMaps[curveIdx].begin(),vertexMaps[curveIdx].end(),q);
        int i = x - vertexMaps[curveIdx].begin();

        if(i>=vertexMaps[curveIdx].size()){
            //this should mean that we simply map to the very last vertex
            return {vertexMaps[curveIdx].size()-2,1.0};
        }

        if(i==0){
            return {0,0.0};
        }else{
            auto s = vertexMaps[curveIdx][i-1];
            auto t = vertexMaps[curveIdx][i];

            // assert(t>s);
            // std::cout << (t.getFraction() - s.getFraction()) << std::endl;
            if(t.getPoint() == s.getPoint()){
                //the typical case
                if(t>s) {
                    return {i - 1, (q.getFraction() - s.getFraction()) / (t.getFraction() - s.getFraction())};
                }else{
                    return {i-1,0.0};
                }
            }else if(t.getPoint() == s.getPoint()+1){
                return {i - 1, (q.getFraction() - s.getFraction()) / (1.0 - s.getFraction())};
            }else{
                assert(false);
            }

        }
    }

    int mergeOverlappingClusters(ClusteringResult& clusters, double overlapRatio = 0.7) {
        int sizeBefore = clusters.size();
        bool changes;
        do {
            changes = false;
            for (int i = 0; i < clusters.size() && !changes; i++) {
                for (int j = i + 1; j < clusters.size() && !changes; j++) {
                    if (clusters[i].getCenter().getCurveIndex() != clusters[j].getCenter().getCurveIndex())
                        continue;
                    double totalOverlap = 0;
                    double totalLengthI = 0;
                    double totalLengthJ = 0;
                    // calculate total matching overlap
                    for (const CInterval& matching0 : clusters[i].getMatching()) {
                        for (const CInterval& matching1 : clusters[j].getMatching()) {
                            if (matching0.getCurveIndex() != matching1.getCurveIndex())
                                continue;

                            double start0 = mapSimplificationToBase(matching0.getCurveIndex(), matching0.begin).convert();
                            double end0 = mapSimplificationToBase(matching0.getCurveIndex(), matching0.end).convert();
                            double start1 = mapSimplificationToBase(matching1.getCurveIndex(), matching1.begin).convert();
                            double end1 = mapSimplificationToBase(matching1.getCurveIndex(), matching1.end).convert();
                            
                            totalOverlap += std::max(0.0, std::min(end0, end1) - std::max(start0, start1));
                            totalLengthI += end0 - start0;
                            totalLengthJ += end1 - start1;
                        }
                    }
                    
                    // std::cout << totalOverlap << " " << totalLengthI << std::endl;
                    if (totalOverlap / totalLengthI >= overlapRatio || totalOverlap / totalLengthJ >= overlapRatio) {
                        clusters[i].merge(clusters[j]);
                        clusters.erase(clusters.begin() + j);
                        changes = true;
                    }
                }
            }
        } while (changes);

        return sizeBefore - clusters.size();
    }

    void updateFlags(int samplingRate, bool showFlag) {
        if (samplingRate < 0) {
            subsamplingFrequency = -1;
            withSubsamplingRule = false;
        } else {
            subsamplingFrequency = samplingRate;
            withSubsamplingRule = true;
        }
        showFreespaces = showFlag;
    }

    void setWithSort(bool sort) {
        withSort = sort;
    }

    int test(Curves &curves){
        return curves.size();
    }

    void initCurves(Curves &curves, double _simpDelta, double _freeDelta=-1) {
        CurveSimplifier cs(_simpDelta);
        if (subsamplingFrequency >= 0) {
            cs.setSubsampling(subsamplingFrequency);
        } else {
            cs.unsetSubsampling();
        }
        hasGTs = false;
        simplificationDelta = _simpDelta;
        double const guarantee = 1.0 + 1.0 + 2 * (7.0 / 3.0);
        freespaceDelta = (_freeDelta<0?guarantee*_simpDelta:_freeDelta);
        simplifiedCurves.clear();
        times.clear();
        simpIDtoOriginID.clear();
        vertexMaps.clear();
        int total = 0;

        Curves localSimps;

        std::vector<std::vector<int>> localTimes;
        std::vector<int> localIDMap;
        std::vector<std::vector<CPoint>> localVertexMap;

        std::vector<int> cur(omp_get_max_threads(),0);
        std::vector<int> len(omp_get_max_threads(),0);

        for(int i=0;i<omp_get_max_threads();i++){
            std::cout << "["<<i<<"/"<<omp_get_max_threads()<<"] simplified " << cur[i] << " with length " << len[i] <<".     "<<std::endl;
        }
        std::cout << "Progress: " << total << "/" << curves.size()<<std::endl;

#pragma omp declare reduction (mergeSimps : Curves : omp_out.insert(omp_out.end(), std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))
#pragma omp declare reduction (mergeTimes : std::vector<std::vector<int>> : omp_out.insert(omp_out.end(), std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))
#pragma omp declare reduction (mergeIDMap : std::vector<int> : omp_out.insert(omp_out.end(), std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))
#pragma omp declare reduction (mergeVMap : std::vector<std::vector<CPoint>> : omp_out.insert(omp_out.end(), std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))

//#pragma omp parallel for default(none) shared(l,filter,upAggregated) reduction(merge: localSet)
#pragma omp parallel for default(none) shared(curves, std::cout, total,cur,len) firstprivate(cs) reduction(mergeSimps: localSimps) reduction(mergeTimes: localTimes) reduction(mergeIDMap: localIDMap) reduction(mergeVMap: localVertexMap) schedule(dynamic)
        for (int i = 0; i < curves.size(); ++i) {
            Curve curve = curves[i];
            auto simp = cs.simplify(curve);
//#pragma omp critical
            {
                total++;
                localSimps.push_back(simp);
                localIDMap.push_back(i);
                localTimes.push_back(cs.getTimes());
                cur[omp_get_thread_num()]=i;
                len[omp_get_thread_num()]=localSimps.back().size();
#pragma omp critical
                {
                    if (((total + 1) % 10 == 0) || total == curves.size() ) {
                        printf("\033[0G\033[%dA",omp_get_max_threads()+1);
                        for(int j=0;j<omp_get_max_threads();j++){
                            std::cout << "["<<j<<"/"<<omp_get_max_threads()<<"] simplified " << cur[j] << " with length " << len[j] <<".     " << std::endl;
                        }
                        std::cout << "Progress: " << total << "/" << curves.size()<<std::endl;
                    }
                }
                //std::cout << "Simplified " << i << " with length " << localSimps.back().size() << " (" << total
                //          << "/" << curves.size() << ")" << std::endl;


                //build vertexMapping

                localVertexMap.emplace_back();
                int sj = 0;
                for (int ci=0;ci<curve.size();ci++){
                    if(sj + 2 < localTimes.back().size() && localTimes.back()[sj+1]<ci){
                        sj ++;
                    }
                    auto interval = IntersectionAlgorithm::intersection_interval(curve[ci],7*simplificationDelta/3,simp[sj],simp[sj+1]);
                    assert(!interval.is_empty());
                    if(localVertexMap.back().empty()){
                        localVertexMap.back().push_back({sj,interval.begin});
                    }else{
                        localVertexMap.back().push_back(std::max(localVertexMap.back().back(),{sj,interval.begin}));
                    }
                }
            };
        }
        simplifiedCurves = localSimps;
        unsimplifiedCurves = localSimps;
        times = localTimes;
        simpIDtoOriginID = localIDMap;
        vertexMaps = localVertexMap;
        Curves curves2(simplifiedCurves.size());
        for(int i = 0; i < curves.size(); i++) {
            unsimplifiedCurves[i] = curves[simpIDtoOriginID[i]];
            //std::cout << simpIDtoOriginID[i] << std::endl;
        }
        //simplifiedCurves = curves2;

        }

    void initCurvesWithGT(Curves &curves, double _simpDelta, std::vector<FrameLabeling> GTs, double _freeDelta=-1) {
        assert(curves.size() == GTs.size());

        assert(curves.size() == GTs.size());
        CurveSimplifier cs(_simpDelta);
        if (subsamplingFrequency >= 0) {
            cs.setSubsampling(subsamplingFrequency);
        } else {
            cs.unsetSubsampling();
        }
        hasGTs = true;
        simplificationDelta = _simpDelta;
        double const guarantee = 1.0 + 1.0 + 2 * (7.0 / 3.0);
        freespaceDelta = (_freeDelta<0?guarantee*_simpDelta:_freeDelta);
        simplifiedCurves.clear();
        simplifiedGTs.clear();
        simpIDtoOriginID.clear();
        vertexMaps.clear();
        times.clear();
        int total = 0;


        Curves localSimps;
        std::vector<std::vector<int>> localTimes;
        std::vector<int> localIDMap;
        std::vector<std::vector<CPoint>> localVertexMap;
        std::vector<std::vector<std::pair<Label,CPoint>>> localsimpGTs;

        std::vector<int> cur(omp_get_max_threads(),0);
        std::vector<int> len(omp_get_max_threads(),0);

        for(int i=0;i<omp_get_max_threads();i++){
            std::cout << "["<<i<<"/"<<omp_get_max_threads()<<"] simplified " << cur[i] << " with length " << len[i] <<".     " << std::endl;
        }
        std::cout << "Progress: " << total << "/" << curves.size()<<std::endl;

//#pragma omp parallel for default(none) shared(curves, std::cout, GTs) firstprivate(cs) schedule(dynamic)
#pragma omp declare reduction (mergeSimps : Curves : omp_out.insert(omp_out.end(), std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))
#pragma omp declare reduction (mergeTimes : std::vector<std::vector<int>> : omp_out.insert(omp_out.end(), std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))
#pragma omp declare reduction (mergeIDMap : std::vector<int> : omp_out.insert(omp_out.end(), std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))
#pragma omp declare reduction (mergeVMap : std::vector<std::vector<CPoint>> : omp_out.insert(omp_out.end(), std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))
#pragma omp declare reduction (mergesimGT : std::vector<std::vector<std::pair<Label,CPoint>>> : omp_out.insert(omp_out.end(), std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))

//#pragma omp parallel for default(none) shared(l,filter,upAggregated) reduction(merge: localSet)
#pragma omp parallel for default(none) shared(curves, std::cout,GTs,cur,len,total) firstprivate(cs) reduction(mergeSimps: localSimps) reduction(mergeTimes: localTimes) reduction(mergeIDMap: localIDMap) reduction(mergeVMap: localVertexMap) reduction(mergesimGT:localsimpGTs) schedule(dynamic)

        for (int i = 0; i < curves.size(); ++i) {

            total++;
            Curve curve = curves[i];
            auto simp = cs.simplify(curve,&(GTs[i]));
            localSimps.push_back(simp);
            localIDMap.push_back(i);
            localTimes.push_back(cs.getTimes());
            localsimpGTs.push_back(cs.getSimplifiedGTs());

            cur[omp_get_thread_num()]=i;
            len[omp_get_thread_num()]=localSimps.back().size();
#pragma omp critical
            {
                if (((total + 1) % 10 == 0) || total == curves.size() ) {
                    printf("\033[0G\033[%dA",omp_get_max_threads()+1);
                    for(int j=0;j<omp_get_max_threads();j++){
                        std::cout << "["<<j<<"/"<<omp_get_max_threads()<<"] simplified " << cur[j] << " with length " << len[j] <<".     " << std::endl;
                    }
                    std::cout << "Progress: " << total << "/" << curves.size()<<std::endl;
                }
            }


            //build vertexMapping

            localVertexMap.emplace_back();
            int sj = 0;
            for (int ci=0;ci<curve.size();ci++){
                if(sj + 2 < localTimes.back().size() && localTimes.back()[sj+1]<ci){
                    sj ++;
                }
                auto interval = IntersectionAlgorithm::intersection_interval(curve[ci],7*simplificationDelta/3,simp[sj],simp[sj+1]);
                assert(!interval.is_empty());
                if(localVertexMap.back().empty()){
                    localVertexMap.back().push_back({sj,interval.begin});
                }else{
                    localVertexMap.back().push_back(std::max(localVertexMap.back().back(),{sj,interval.begin}));
                }
            }
        }
        simplifiedCurves = localSimps;
        times = localTimes;
        simpIDtoOriginID = localIDMap;
        vertexMaps = localVertexMap;
        simplifiedGTs = localsimpGTs;

        
    }

    template<typename func>
    ClusteringResult greedyCover(int l, int rounds, func filter,bool withShow=false, long long* size = nullptr) {

        assert(not simplifiedCurves.empty());

        //Curves curves;
        //Curves bestresult;
        std::vector<Candidate> bestResultVisualizer;
        CandidateSetPQ cs = CandidateSetPQ(simplifiedCurves, freespaceDelta);
        if(withShow){
            SparseFreeSpacesVisualizer sfsv(cs.sparsefreespaces);
            sfsv.show();
        }
        //cs.computeCandidates(l);
        cs.ultrafastComputeSmall(l, filter);

        if(size != nullptr){
            *size = 0;
            while(not cs.empty()){
                *size += cs.top().matching.size();
                cs.pop();
            }
            //*size = (int)(cs.size());
            return bestResultVisualizer;
        }


        //io::exportSubcurve("/Users/styx/data/curveclustering/results/bestcandidate.txt",curves[cs.top().first],cs.top().second.getStart(),cs.top().second.getEnd());

        for (int r = 0; r < rounds; ++r) {

            bool lastRound = r==rounds-1;

            std::vector<Candidate> result;
            std::vector<CInterval> covering;
            std::vector<double> suffixLengths;

            std::cout << "Round " << r<<std::endl;

            int internal_max = 500001;
            for (int i = 0; i < internal_max; ++i) {
                bool lastInternalRound = i==internal_max-1;
                auto ID = 0;
                auto roundID = i;

                //first verify that we need to find another center, otherwise output solution

                if (lastInternalRound ||
                    lengthOfUncovered(simplifiedCurves, result) <= EPSILON) {
                    //std::cout << "\nTrying to refine... ";
                    int deletecount = 0;
                    for (int igni = result.size() - 1; igni >= 0; --igni) {
                        //i is the index to be ignored
                        std::vector<Candidate> temp;
                        for (int j = 0; j < result.size(); ++j) {
                            if (igni != j) {
                                temp.push_back(result[j]);
                            }
                        }
                        double importance = lengthOfUncovered(simplifiedCurves, temp);
                        if (importance <= EPSILON) {
                            //std::cout << " ( " << igni << " , " << importance << " )";
                            result.erase(result.begin() + igni);
                            deletecount += 1;
                        } else {
                            result[igni].importance = importance;
                        }
                    }
                    if (deletecount == 0)
                        //std::cout << "nothing";
                        //std::cout << " greedely deleted\n";
                        for (int igni = result.size() - 1; igni >= 0; --igni) {
                            //i is the index to be ignored
                            std::vector<Candidate> temp;
                            for (int j = 0; j < result.size(); ++j) {
                                if (igni != j) {
                                    temp.push_back(result[j]);
                                }
                            }
                            double importance = lengthOfUncovered(simplifiedCurves, temp);
                            result[igni].importance = importance;
                        }
                    //std::sort(result.begin(), result.end(), [](auto& a, auto& b){return a.second.importance > b.second.importance;});
                    std::cout << "\nSolution of size " << result.size() << " found. ";
                    if (bestResultVisualizer.empty() || bestResultVisualizer.size() > result.size()) {
                        if (!bestResultVisualizer.empty()) {
                            std::cout << "This improves on the current best by "
                                      << bestResultVisualizer.size() - result.size() << "! ";
                        }
                        bestResultVisualizer.clear();
                        for (const auto &sub: result) {
                            bestResultVisualizer.push_back(sub);
                        }
                    }
                    break;

                }

                //update top, until the roundID matches
                while (cs.top().roundOfUpdate != roundID) {
                    ID++;
                    Candidate c = cs.top();
                    cs.pop();
                    updateCandidate(simplifiedCurves, c, covering, suffixLengths, roundID);
                    cs.push(c);
                }

                //if top index is too light, we also stop
                if(cs.top().semiUpdatedCoverLength <= EPSILON){
                    //std::cout << "\nTrying to refine... ";
                    int deletecount = 0;
                    for (int igni = result.size() - 1; igni >= 0; --igni) {
                        //i is the index to be ignored
                        std::vector<Candidate> temp;
                        for (int j = 0; j < result.size(); ++j) {
                            if (igni != j) {
                                temp.push_back(result[j]);
                            }
                        }
                        double importance = lengthOfUncovered(simplifiedCurves, temp);
                        if (importance <= EPSILON) {
                            //std::cout << " ( " << igni << " , " << importance << " )";
                            result.erase(result.begin() + igni);
                            deletecount += 1;
                        } else {
                            result[igni].importance = importance;
                        }
                    }
                    if (deletecount == 0)
                        //std::cout << "nothing";
                        //std::cout << " greedely deleted\n";
                        for (int igni = result.size() - 1; igni >= 0; --igni) {
                            //i is the index to be ignored
                            std::vector<Candidate> temp;
                            for (int j = 0; j < result.size(); ++j) {
                                if (igni != j) {
                                    temp.push_back(result[j]);
                                }
                            }
                            double importance = lengthOfUncovered(simplifiedCurves, temp);
                            result[igni].importance = importance;
                        }
                    //std::sort(result.begin(), result.end(), [](auto& a, auto& b){return a.second.importance > b.second.importance;});
                    std::cout << "\nSolution of size " << result.size() << " found. ";
                    if (bestResultVisualizer.empty() || bestResultVisualizer.size() > result.size()) {
                        if (!bestResultVisualizer.empty()) {
                            std::cout << "This improves on the current best by "
                                      << bestResultVisualizer.size() - result.size() << "! ";
                        }
                        bestResultVisualizer.clear();
                        for (const auto &sub: result) {
                            bestResultVisualizer.push_back(sub);
                        }
                    }
                    break;

                }

                //otherwise add top element to cs


                //printFirst50(cs);

                double addedweight;
                //cs.showCovering(result);
                //std::cout << cs.top().second.semiUpdatedCoverLength << " ";
                //std::pair<int,Candidate> ctemp = cs.top();
                //std::cout << cs.top().second.matchings.size();
                //updateCandidate(curves, ctemp.second, covering, suffixLengths, roundID);


                //cs.top() will be added to c
                if (i == 0) {
                    std::vector<Candidate> temp;
                    for (int j = 0; j < r; ++j) {
                        temp.push_back(cs.top());
                        cs.pop();
                    }
                    Candidate c = cs.top();
                    result.push_back(c);

                    cs.pop();
                    addedweight = c.semiUpdatedCoverLength;
                    c.semiUpdatedCoverLength = 0;
                    cs.push(c);
                    for (const auto &ctemp: temp)
                        cs.push(ctemp);
                } else {
                    Candidate c = cs.top();
                    result.push_back(c);

                    //io::exportSubcurve("/Users/styx/data/curveclustering/results/bestcandidate.txt",curves[0],cs.top().getStart(),cs.top().second.getEnd());
                    //cs.showCovering(result);

                    cs.pop();
                    addedweight = c.semiUpdatedCoverLength;
                    c.semiUpdatedCoverLength = 0;
                    cs.push(c);
                }

                std::cout << "\033[0G"<<"Identified center #"<<i<<std::flush;

                //std::cout << " added weight: " << addedweight << ". Updated " << ID << " lengths\n";


                Candidate c = result.back();

                //update covering
                std::vector<CInterval> temp(covering);
                temp.insert(temp.end(), c.matching.begin(), c.matching.end());
                std::sort(temp.begin(), temp.end(), cmpLeftLower);
                covering.clear();
                suffixLengths.clear();
                suffixLengths.push_back(0);
                CInterval cur = temp[0];
                for (int covI = 1; covI < temp.size(); covI++) {
                    CInterval next = temp[covI];
                    if ((next.getCurveIndex() == cur.getCurveIndex()) && (next.getBegin() < cur.end)) {
                        cur.end = std::max(next.end, cur.end);
                    } else {
                        covering.push_back(cur);
                        suffixLengths.push_back(
                                suffixLengths.back() +
                                simplifiedCurves[cur.getCurveIndex()].subcurve_length(cur.getBegin(), cur.end));
                        cur = next;
                    }
                }
                covering.push_back(cur);
                suffixLengths.push_back(
                        suffixLengths.back() +
                        simplifiedCurves[cur.getCurveIndex()].subcurve_length(cur.getBegin(), cur.end));
                //if((i+1)%100==0)
                //cs.showCovering(result);
            }
            std::cout << std::endl;
            //std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(swatchinsert.elapsed()).count() << std::endl;
            //std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(swatchupdate.elapsed()).count() << std::endl;
            //std::cout << "Cleaning up for next round ... ";
            if(!lastRound)
                cs.resetWeights();
            //std::cout << "Done";
        }
        if (withSort) {
            std::sort(bestResultVisualizer.begin(), bestResultVisualizer.end(),
                      [](const Candidate &a, const Candidate &b) {
                          return (a.matching[0].getCurveIndex() < b.matching[0].getCurveIndex()) ||
                                 ((a.matching[0].getCurveIndex() == b.matching[0].getCurveIndex()) &&
                                  (a.matching[0].getBegin() < b.matching[0].getBegin()));
                      });
        }
        if (showFreespaces) {
#ifdef HASVISUAL
            cs.showCovering(bestResultVisualizer);
#endif
        }
        ClusteringResult cr(bestResultVisualizer);
        return cr;
    }
};

template<typename func>
std::vector<Candidate>
greedyCoverUnsanitizedOutput(Curves &curves, double delta, int l, int max_rounds, bool show, func filter) {
    Curves bestresult;
    stdc::SWatch swatchinsert;
    stdc::SWatch swatchupdate;
    std::vector<Candidate> bestResultVisualizer;
    double const guarantee = 1.0 + 1.0 + 2 * (7.0 / 3.0);
    CandidateSetPQ cs = CandidateSetPQ(curves, guarantee * delta);
    //cs.computeCandidates(l);
    cs.ultrafastComputeSmall(l, filter);


    //io::exportSubcurve("/Users/styx/data/curveclustering/results/bestcandidate.txt",curves[cs.top().first],cs.top().second.getStart(),cs.top().second.getEnd());

    for (int r = 0; r < max_rounds; ++r) {

        //swatchinsert.reset();
        //swatchupdate.reset();

        std::vector<Candidate> result;
        std::vector<CInterval> covering;
        std::vector<double> suffixLengths;

        std::cout << "Round " << r<<std::endl;

        for (int i = 0; i < 500001; ++i) {
            bool lastInternalRound = i==500000;
            auto ID = 0;
            auto roundID = i;

            if (lastInternalRound || lengthOfUncovered(curves, result) <= EPSILON) {
                //std::cout << "\nTrying to refine... ";
                int deletecount = 0;
                for (int igni = result.size() - 1; igni >= 0; --igni) {
                    //i is the index to be ignored
                    std::vector<Candidate> temp;
                    for (int j = 0; j < result.size(); ++j) {
                        if (igni != j) {
                            temp.push_back(result[j]);
                        }
                    }
                    double importance = lengthOfUncovered(curves, temp);
                    if (importance <= EPSILON) {
                        //std::cout << " ( " << igni << " , " << importance << " )";
                        result.erase(result.begin() + igni);
                        deletecount += 1;
                    } else {
                        result[igni].importance = importance;
                    }
                }
                if (deletecount == 0)
                    //std::cout << "nothing";
                    //std::cout << " greedely deleted\n";
                    for (int igni = result.size() - 1; igni >= 0; --igni) {
                        //i is the index to be ignored
                        std::vector<Candidate> temp;
                        for (int j = 0; j < result.size(); ++j) {
                            if (igni != j) {
                                temp.push_back(result[j]);
                            }
                        }
                        double importance = lengthOfUncovered(curves, temp);
                        result[igni].importance = importance;
                    }
                //std::sort(result.begin(), result.end(), [](auto& a, auto& b){return a.second.importance > b.second.importance;});
                std::cout << " Solution of size " << result.size() << " found. ";
                if (bestresult.empty() || bestresult.size() > result.size()) {
                    if (!bestresult.empty()) {
                        std::cout << "This improves on the current best by " << bestresult.size() - result.size()
                                  << "! ";
                    }
                    bestresult.clear();
                    bestResultVisualizer.clear();
                    for (const auto &sub: result) {
                        bestresult.push_back(
                                curves[sub.getCurveIndex()].constructSubcurve(sub.getBegin(), sub.getEnd()));
                        bestResultVisualizer.push_back(sub);
                    }
                }
                break;

            }

            //update top, until the roundID matches
            while (cs.top().roundOfUpdate != roundID) {
                ID++;
                Candidate c = cs.top();
                cs.pop();
                swatchupdate.start();
                updateCandidate(curves, c, covering, suffixLengths, roundID);
                swatchupdate.stop();
                swatchinsert.start();
                cs.push(c);
                swatchinsert.stop();
            }

            if (cs.top().semiUpdatedCoverLength <= EPSILON) {
                //std::cout << "\nTrying to refine... ";
                int deletecount = 0;
                for (int igni = result.size() - 1; igni >= 0; --igni) {
                    //i is the index to be ignored
                    std::vector<Candidate> temp;
                    for (int j = 0; j < result.size(); ++j) {
                        if (igni != j) {
                            temp.push_back(result[j]);
                        }
                    }
                    double importance = lengthOfUncovered(curves, temp);
                    if (importance <= EPSILON) {
                        //std::cout << " ( " << igni << " , " << importance << " )";
                        result.erase(result.begin() + igni);
                        deletecount += 1;
                    } else {
                        result[igni].importance = importance;
                    }
                }
                if (deletecount == 0)
                    //std::cout << "nothing";
                    //std::cout << " greedely deleted\n";
                    for (int igni = result.size() - 1; igni >= 0; --igni) {
                        //i is the index to be ignored
                        std::vector<Candidate> temp;
                        for (int j = 0; j < result.size(); ++j) {
                            if (igni != j) {
                                temp.push_back(result[j]);
                            }
                        }
                        double importance = lengthOfUncovered(curves, temp);
                        result[igni].importance = importance;
                    }
                //std::sort(result.begin(), result.end(), [](auto& a, auto& b){return a.second.importance > b.second.importance;});
                std::cout << " Solution of size " << result.size() << " found. ";
                if (bestresult.empty() || bestresult.size() > result.size()) {
                    if (!bestresult.empty()) {
                        std::cout << "This improves on the current best by " << bestresult.size() - result.size()
                                  << "! ";
                    }
                    bestresult.clear();
                    bestResultVisualizer.clear();
                    for (const auto &sub: result) {
                        bestresult.push_back(
                                curves[sub.getCurveIndex()].constructSubcurve(sub.getBegin(), sub.getEnd()));
                        bestResultVisualizer.push_back(sub);
                    }
                }
                break;

            }

            //printFirst50(cs);

            double addedweight;
            //cs.showCovering(result);
            //std::cout << cs.top().second.semiUpdatedCoverLength << " ";
            //std::pair<int,Candidate> ctemp = cs.top();
            //std::cout << cs.top().second.matchings.size();
            //updateCandidate(curves, ctemp.second, covering, suffixLengths, roundID);


            //cs.top() will be added to c
            if (i == 0) {
                std::vector<Candidate> temp;
                for (int j = 0; j < r; ++j) {
                    temp.push_back(cs.top());
                    cs.pop();
                }
                Candidate c = cs.top();
                result.push_back(c);

                cs.pop();
                addedweight = c.semiUpdatedCoverLength;
                c.semiUpdatedCoverLength = 0;
                cs.push(c);
                for (const auto &ctemp: temp)
                    cs.push(ctemp);
            } else {
                Candidate c = cs.top();
                result.push_back(c);

                //io::exportSubcurve("/Users/styx/data/curveclustering/results/bestcandidate.txt",curves[0],cs.top().getStart(),cs.top().second.getEnd());
                //cs.showCovering(result);

                cs.pop();
                addedweight = c.semiUpdatedCoverLength;
                c.semiUpdatedCoverLength = 0;
                cs.push(c);
            }

            std::cout << "\033[0G"<<"Identified center #"<<i<<std::flush;

            //std::cout << " added weight: " << addedweight << ". Updated " << ID << " lengths\n";


            Candidate c = result.back();

            //update covering
            std::vector<CInterval> temp(covering);
            temp.insert(temp.end(), c.matching.begin(), c.matching.end());
            std::sort(temp.begin(), temp.end(), cmpLeftLower);
            covering.clear();
            suffixLengths.clear();
            suffixLengths.push_back(0);
            CInterval cur = temp[0];
            for (int covI = 1; covI < temp.size(); covI++) {
                CInterval next = temp[covI];
                if ((next.getCurveIndex() == cur.getCurveIndex()) && (next.getBegin() < cur.end)) {
                    cur.end = std::max(next.end, cur.end);
                } else {
                    covering.push_back(cur);
                    suffixLengths.push_back(
                            suffixLengths.back() +
                            curves[cur.getCurveIndex()].subcurve_length(cur.getBegin(), cur.end));
                    cur = next;
                }
            }
            covering.push_back(cur);
            suffixLengths.push_back(
                    suffixLengths.back() + curves[cur.getCurveIndex()].subcurve_length(cur.getBegin(), cur.end));
            //if((i+1)%100==0)
            //cs.showCovering(result);
        }
        //std::cout << std::endl;
        //std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(swatchinsert.elapsed()).count() << std::endl;
        //std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(swatchupdate.elapsed()).count() << std::endl;
        //std::cout << "Cleaning up for next round ... ";
        if(r + 1 < max_rounds)
            cs.resetWeights();
        //std::cout << "Done";
    }
    std::sort(bestResultVisualizer.begin(), bestResultVisualizer.end(),
              [](const Candidate &a, const Candidate &b) {
                  return (a.matching[0].getCurveIndex() < b.matching[0].getCurveIndex()) ||
                         ((a.matching[0].getCurveIndex() == b.matching[0].getCurveIndex()) &&
                          (a.matching[0].getBegin() < b.matching[0].getBegin()));
              });
    if (show) {
#ifdef HASVISUAL
        cs.showCovering(bestResultVisualizer);
#endif
    }
    return bestResultVisualizer;
}
