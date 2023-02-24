#pragma once

//#include "center_algs.h"
//#include "../src/clustering_algs.h"
#include "curve.h"
#include "candidate.h"
#include "SWatch.h"
#include "io.h"
#include "curve_simplification.h"
#include <algorithm>

/*
Clustering computeCenterClustering(
	Curves const& curves, int k, int l, ClusterAlg cluster_alg, CenterAlg center_alg, int max_rounds = 10);

Clustering approxCover(Curves& curves, double delta,int l,int max_rounds = 10);
*/
Curves greedyCoverAlreadySimplified(Curves& curves, double delta, int l, int max_rounds = 10, bool show = false);
void updateCandidate(Curves& curves, Candidate& c, std::vector<Subcurve> covering, std::vector<double> suffixLengths, int roundID);
double lengthOfUncovered(Curves curves, std::vector<Candidate> candidateSet);
bool cmpLeftLower(Subcurve l, Subcurve r);

//template<typename func> Curves greedyCover(Curves& curves, double delta, int l, );

class CurveClusterer{
private:
    bool withSubsamplingRule;
    int subsamplingFrequency;

    //visual stuff
    bool showFreespaces;

    Curves unsimplifiedCurves;

    double delta;


    //TODO: make freespaces persistent and updatable
    //TODO: unsimplify candidate infos

    bool hasGTs;
    //std::vector<std::pair<Label,int>> groundTruth;



public:

    //TODO: make these accessible in a better way, for filters...
    Curves simplifiedCurves;
    std::vector<ParamLabeling> simplifiedGTs;

    CurveClusterer(int samplingRate, bool showFlag){
        updateFlags(samplingRate,showFlag);
    };

    void updateFlags(int samplingRate, bool showFlag){
        if (samplingRate<0){
            subsamplingFrequency = -1;
            withSubsamplingRule = false;
        }else{
            subsamplingFrequency = samplingRate;
            withSubsamplingRule = true;
        }
        showFreespaces = showFlag;
    }

    void initCurves(Curves& curves, double _delta){
        hasGTs = false;
        delta = _delta;
        simplifiedCurves.clear();
        for(auto curve : curves){
            if(subsamplingFrequency >= 0){
                simplifiedCurves.push_back(good_simplification(curve,delta,subsamplingFrequency));
            }else{
                simplifiedCurves.push_back(good_simplification(curve,delta));
            }
        }
    }

    void initCurves(Curves& curves, double _delta, std::vector<FrameLabeling> GTs){
        assert(curves.size() == GTs.size());
        hasGTs = true;
        delta = _delta;
        simplifiedCurves.clear();
        simplifiedGTs.clear();
        for(int i=0;i<curves.size();++i){
            Curve curve = curves[i];
            simplifiedGTs.emplace_back();
            if(subsamplingFrequency >= 0){
                simplifiedCurves.push_back(good_simplification(curve,delta,GTs[i],&(simplifiedGTs[i]),subsamplingFrequency));
            }else{
                simplifiedCurves.push_back(good_simplification(curve,delta,GTs[i],&(simplifiedGTs[i])));
            }
        }
    }

    template<typename func> std::vector<Candidate> greedyCover(int l, int rounds, func filter){

        assert(not simplifiedCurves.empty());

        //Curves curves;
        //Curves bestresult;
        std::vector<Candidate> bestResultVisualizer;
        double const guarantee = 1.0+1.0+2*(7.0/3.0);
        CandidateSetPQ cs = CandidateSetPQ(simplifiedCurves,guarantee*delta);
        //cs.computeCandidates(l);
        cs.ultrafastComputeCandidates(l,filter);


        //io::exportSubcurve("/Users/styx/data/curveclustering/results/bestcandidate.txt",curves[cs.top().first],cs.top().second.getStart(),cs.top().second.getEnd());

        for(int r = 0;r<rounds;++r) {

            std::vector<Candidate> result;
            std::vector<Subcurve> covering;
            std::vector<double> suffixLengths;

            std::cout << "\nRound " << r;

            for (int i = 0; i < 501; ++i) {
                auto ID = 0;
                auto roundID = i;
                std::cout << ".";
                //update top, until the roundID matches
                while (cs.top().roundOfUpdate != roundID) {
                    ID++;
                    Candidate c = cs.top();
                    cs.pop();
                    updateCandidate(simplifiedCurves, c, covering, suffixLengths, roundID);
                    cs.push(c);
                }
                if (lengthOfUncovered(simplifiedCurves,result) <= EPSILON || cs.top().semiUpdatedCoverLength <= EPSILON) {
                    std::cout << "\nTrying to refine... ";
                    int deletecount = 0;
                    for(int igni=result.size()-1;igni>=0;--igni){
                        //i is the index to be ignored
                        std::vector<Candidate> temp;
                        for(int j=0;j<result.size();++j){
                            if(igni != j){
                                temp.push_back(result[j]);
                            }
                        }
                        double importance = lengthOfUncovered(simplifiedCurves,temp);
                        if(importance <= EPSILON) {
                            std::cout << " ( " << igni << " , " << importance << " )";
                            result.erase(result.begin() + igni);
                            deletecount += 1;
                        }else{
                            result[igni].importance = importance;
                        }
                    }
                    if (deletecount == 0)
                        std::cout << "nothing";
                    std::cout << " greedely deleted\n";
                    for(int igni=result.size()-1;igni>=0;--igni){
                        //i is the index to be ignored
                        std::vector<Candidate> temp;
                        for(int j=0;j<result.size();++j){
                            if(igni != j){
                                temp.push_back(result[j]);
                            }
                        }
                        double importance = lengthOfUncovered(simplifiedCurves,temp);
                        result[igni].importance = importance;
                    }
                    //std::sort(result.begin(), result.end(), [](auto& a, auto& b){return a.second.importance > b.second.importance;});
                    std::cout << " Solution of size " << result.size() << " found. ";
                    if (bestResultVisualizer.empty() || bestResultVisualizer.size() > result.size()) {
                        if(!bestResultVisualizer.empty()){
                            std::cout << "This improves on the current best by " << bestResultVisualizer.size() - result.size() << "! ";
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
                std::vector<Subcurve> temp(covering);
                temp.insert(temp.end(), c.matchings.begin(), c.matchings.end());
                std::sort(temp.begin(), temp.end(), cmpLeftLower);
                covering.clear();
                suffixLengths.clear();
                suffixLengths.push_back(0);
                Subcurve cur = temp[0];
                for (int covI = 1; covI < temp.size(); covI++) {
                    Subcurve next = temp[covI];
                    if ((next.curveIdx == cur.curveIdx) && (next.start < cur.end)) {
                        cur.end = std::max(next.end, cur.end);
                    } else {
                        covering.push_back(cur);
                        suffixLengths.push_back(
                                suffixLengths.back() + simplifiedCurves[cur.curveIdx].subcurve_length(cur.start, cur.end));
                        cur = next;
                    }
                }
                covering.push_back(cur);
                suffixLengths.push_back(
                        suffixLengths.back() + simplifiedCurves[cur.curveIdx].subcurve_length(cur.start, cur.end));
                //if((i+1)%100==0)
                //cs.showCovering(result);
            }
            std::cout << std::endl;
            //std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(swatchinsert.elapsed()).count() << std::endl;
            //std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(swatchupdate.elapsed()).count() << std::endl;
            std::cout << "Cleaning up for next round ... ";
            cs.resetWeights();
            std::cout << "Done";
        }
        std::sort(bestResultVisualizer.begin(),bestResultVisualizer.end(),
                  [](const Candidate& a,const Candidate& b){
                      return (a.matchings[0].curveIdx < b.matchings[0].curveIdx) ||
                             ( (a.matchings[0].curveIdx == b.matchings[0].curveIdx) && (a.matchings[0].start < b.matchings[0].start));
                  });
        if(showFreespaces) {
#ifdef HASVISUAL
            cs.showCovering(bestResultVisualizer);
#endif
        }
        return bestResultVisualizer;
    }


};

template<typename func> std::vector<Candidate> greedyCoverUnsanitizedOutput(Curves& curves, double delta, int l, int max_rounds, bool show, func filter){
    Curves bestresult;
    stdc::SWatch swatchinsert;
    stdc::SWatch swatchupdate;
    std::vector<Candidate> bestResultVisualizer;
    double const guarantee = 1.0+1.0+2*(7.0/3.0);
    CandidateSetPQ cs = CandidateSetPQ(curves,guarantee*delta);
    //cs.computeCandidates(l);
    cs.ultrafastComputeCandidates(l,filter);


    //io::exportSubcurve("/Users/styx/data/curveclustering/results/bestcandidate.txt",curves[cs.top().first],cs.top().second.getStart(),cs.top().second.getEnd());

    for(int r = 0;r<max_rounds;++r) {

        //swatchinsert.reset();
        //swatchupdate.reset();

        std::vector<Candidate> result;
        std::vector<Subcurve> covering;
        std::vector<double> suffixLengths;

        std::cout << "\nRound " << r;

        for (int i = 0; i < 501; ++i) {
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
            if (lengthOfUncovered(curves,result) <= EPSILON || cs.top().semiUpdatedCoverLength <= EPSILON) {
                std::cout << "\nTrying to refine... ";
                int deletecount = 0;
                for(int igni=result.size()-1;igni>=0;--igni){
                    //i is the index to be ignored
                    std::vector<Candidate> temp;
                    for(int j=0;j<result.size();++j){
                        if(igni != j){
                            temp.push_back(result[j]);
                        }
                    }
                    double importance = lengthOfUncovered(curves,temp);
                    if(importance <= EPSILON) {
                        std::cout << " ( " << igni << " , " << importance << " )";
                        result.erase(result.begin() + igni);
                        deletecount += 1;
                    }else{
                        result[igni].importance = importance;
                    }
                }
                if (deletecount == 0)
                    std::cout << "nothing";
                std::cout << " greedely deleted\n";
                for(int igni=result.size()-1;igni>=0;--igni){
                    //i is the index to be ignored
                    std::vector<Candidate> temp;
                    for(int j=0;j<result.size();++j){
                        if(igni != j){
                            temp.push_back(result[j]);
                        }
                    }
                    double importance = lengthOfUncovered(curves,temp);
                    result[igni].importance = importance;
                }
                //std::sort(result.begin(), result.end(), [](auto& a, auto& b){return a.second.importance > b.second.importance;});
                std::cout << " Solution of size " << result.size() << " found. ";
                if (bestresult.empty() || bestresult.size() > result.size()) {
                    if(!bestresult.empty()){
                        std::cout << "This improves on the current best by " << bestresult.size() - result.size() << "! ";
                    }
                    bestresult.clear();
                    bestResultVisualizer.clear();
                    for (const auto &sub: result) {
                        bestresult.push_back(
                                curves[sub.getIndex()].constructSubcurve(sub.getStart(), sub.getEnd()));
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
            std::vector<Subcurve> temp(covering);
            temp.insert(temp.end(), c.matchings.begin(), c.matchings.end());
            std::sort(temp.begin(), temp.end(), cmpLeftLower);
            covering.clear();
            suffixLengths.clear();
            suffixLengths.push_back(0);
            Subcurve cur = temp[0];
            for (int covI = 1; covI < temp.size(); covI++) {
                Subcurve next = temp[covI];
                if ((next.curveIdx == cur.curveIdx) && (next.start < cur.end)) {
                    cur.end = std::max(next.end, cur.end);
                } else {
                    covering.push_back(cur);
                    suffixLengths.push_back(
                            suffixLengths.back() + curves[cur.curveIdx].subcurve_length(cur.start, cur.end));
                    cur = next;
                }
            }
            covering.push_back(cur);
            suffixLengths.push_back(
                    suffixLengths.back() + curves[cur.curveIdx].subcurve_length(cur.start, cur.end));
            //if((i+1)%100==0)
            //cs.showCovering(result);
        }
        std::cout << std::endl;
        std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(swatchinsert.elapsed()).count() << std::endl;
        std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(swatchupdate.elapsed()).count() << std::endl;
        std::cout << "Cleaning up for next round ... ";
        cs.resetWeights();
        std::cout << "Done";
    }
    std::sort(bestResultVisualizer.begin(),bestResultVisualizer.end(),
              [](const Candidate& a,const Candidate& b){
                  return (a.matchings[0].curveIdx < b.matchings[0].curveIdx) ||
                         ( (a.matchings[0].curveIdx == b.matchings[0].curveIdx) && (a.matchings[0].start < b.matchings[0].start));
              });
    if(show) {
#ifdef HASVISUAL
        cs.showCovering(bestResultVisualizer);
#endif
    }
    return bestResultVisualizer;
}
