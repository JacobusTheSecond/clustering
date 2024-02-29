#pragma once

//#include "center_algs.h"
//#include "../src/clustering_algs.h"
#include "Curve.h"
#include "Candidate.h"
#include "io.h"
#include "CurveSimplification.h"
#include "SWatch.h"
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

    Curves unsimplifiedCurves;

    double delta;


    //TODO: make freespaces persistent and updatable
    //TODO: unsimplify candidate infos

    bool hasGTs;
    //std::vector<std::pair<Label,int>> groundTruth;

    std::vector<std::vector<int>> times;


public:

    double getDelta() const{return delta;};

    //TODO: make these accessible in a better way, for filters...
    Curves simplifiedCurves;
    std::vector<ParamLabeling> simplifiedGTs;
    std::vector<int> simpIDtoOriginID;

    CurveClusterer(int samplingRate, bool showFlag) {
        updateFlags(samplingRate, showFlag);
    };

    CurveClusterer() {
        updateFlags(-1, false);
    };

    //TODO: this seems wrong
    CPoint mapSimplificationToBase(int curveIdx, CPoint t) {
        CurveID p = times[curveIdx][t.getPoint()];
        distance_t frac = t.getFraction();
        while (frac > 1.0 + EPSILON) {
            p += 1;
            frac -= 1.0;
        }
        return {p, frac};
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

    void initCurves(Curves &curves, double _delta) {
        CurveSimplifier cs(_delta);
        if (subsamplingFrequency >= 0) {
            cs.setSubsampling(subsamplingFrequency);
        } else {
            cs.unsetSubsampling();
        }
        hasGTs = false;
        delta = _delta;
        simplifiedCurves.clear();
        times.clear();
        simpIDtoOriginID.clear();
        int total = 0;
//#pragma omp parallel for default(none) shared(curves, std::cout, total) firstprivate(cs) schedule(dynamic)
        for (int i = 0; i < curves.size(); ++i) {
            Curve curve = curves[i];
            auto simp = cs.simplify(curve);
//#pragma omp critical
            {
                total++;
                times.emplace_back();
                simplifiedCurves.push_back(simp);
                simpIDtoOriginID.push_back(i);
                std::cout << "Simplified " << i << " with length " << simplifiedCurves.back().size() << " (" << total
                          << "/" << curves.size() << ")" << std::endl;
                times.back() = cs.getTimes();
                assert(times.size() == simplifiedCurves.size());
            };
        }
    }

    void initCurvesWithGT(Curves &curves, double _delta, std::vector<FrameLabeling> GTs) {
        assert(curves.size() == GTs.size());

        assert(curves.size() == GTs.size());
        CurveSimplifier cs(_delta);
        if (subsamplingFrequency >= 0) {
            cs.setSubsampling(subsamplingFrequency);
        } else {
            cs.unsetSubsampling();
        }
        hasGTs = true;
        delta = _delta;
        simplifiedCurves.clear();
        simplifiedGTs.clear();
        simpIDtoOriginID.clear();
        times.clear();
        for (int i = 0; i < curves.size(); ++i) {
            times.emplace_back();
            Curve curve = curves[i];
            simplifiedGTs.emplace_back();
            simplifiedCurves.push_back(cs.simplify(curve, &(GTs[i])));
            simpIDtoOriginID.push_back(i);
            times[i] = cs.getTimes();
            simplifiedGTs[i] = cs.getSimplifiedGTs();
        }
    }

    template<typename func>
    ClusteringResult greedyCover(int l, int rounds, func filter,long long* size = nullptr) {

        assert(not simplifiedCurves.empty());

        //Curves curves;
        //Curves bestresult;
        std::vector<Candidate> bestResultVisualizer;
        double const guarantee = 1.0 + 1.0 + 2 * (7.0 / 3.0);
        CandidateSetPQ cs = CandidateSetPQ(simplifiedCurves, guarantee * delta);
        //cs.computeCandidates(l);
        cs.ultrafastCompute(l, filter);

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

            std::cout << "Round " << r;

            int internal_max = 500001;
            for (int i = 0; i < internal_max; ++i) {
                bool lastInternalRound = i==internal_max-1;
                auto ID = 0;
                auto roundID = i;
                std::cout << "."<<std::flush;
                //update top, until the roundID matches
                while (cs.top().roundOfUpdate != roundID) {
                    ID++;
                    Candidate c = cs.top();
                    cs.pop();
                    updateCandidate(simplifiedCurves, c, covering, suffixLengths, roundID);
                    cs.push(c);
                }
                if (lastInternalRound ||
                    lengthOfUncovered(simplifiedCurves, result) <= EPSILON ||
                    cs.top().semiUpdatedCoverLength <= EPSILON) {
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
    cs.ultrafastCompute(l, filter);


    //io::exportSubcurve("/Users/styx/data/curveclustering/results/bestcandidate.txt",curves[cs.top().first],cs.top().second.getStart(),cs.top().second.getEnd());

    for (int r = 0; r < max_rounds; ++r) {

        //swatchinsert.reset();
        //swatchupdate.reset();

        std::vector<Candidate> result;
        std::vector<CInterval> covering;
        std::vector<double> suffixLengths;

        std::cout << "Round " << r;

        for (int i = 0; i < 501; ++i) {
            bool lastInternalRound = i==500;
            auto ID = 0;
            auto roundID = i;
            std::cout << ".";
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
            if (lastInternalRound || lengthOfUncovered(curves, result) <= EPSILON || cs.top().semiUpdatedCoverLength <= EPSILON) {
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
