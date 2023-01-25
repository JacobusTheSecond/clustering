#pragma once

//#include "center_algs.h"
#include "clustering_algs.h"
#include "curve.h"
#include "candidate.h"
#include "SWatch.h"
#include "io.h"
#include <algorithm>

/*
Clustering computeCenterClustering(
	Curves const& curves, int k, int l, ClusterAlg cluster_alg, CenterAlg center_alg, int max_rounds = 10);

Clustering approxCover(Curves& curves, double delta,int l,int max_rounds = 10);
*/
Curves greedyCover(Curves& curves, double delta, int l, int max_rounds = 10, bool show = false);
void updateCandidate(Curves& curves, Candidate& c, std::vector<std::pair<int,Subcurve>> covering, std::vector<double> suffixLengths, int roundID);
double lengthOfUncovered(Curves curves, std::vector<std::pair<int,Candidate>> candidateSet);
bool cmpLeftLower(std::pair<int,Subcurve> l, std::pair<int,Subcurve> r);

template<typename func> std::vector<std::pair<int,Candidate>> greedyCoverUnsanitizedOutput(Curves& curves, double delta, int l, int max_rounds, bool show, func filter){
    Curves bestresult;
    stdc::SWatch swatchinsert;
    stdc::SWatch swatchupdate;
    std::vector<std::pair<int,Candidate>> bestResultVisualizer;
    double const guarantee = 1.0+1.0+2*(7.0/3.0);
    CandidateSetPQ cs = CandidateSetPQ(curves,guarantee*delta);
    //cs.computeCandidates(l);
    cs.ultrafastComputeCandidates(l,filter);


    //io::exportSubcurve("/Users/styx/data/curveclustering/results/bestcandidate.txt",curves[cs.top().first],cs.top().second.getStart(),cs.top().second.getEnd());

    for(int r = 0;r<max_rounds;++r) {

        //swatchinsert.reset();
        //swatchupdate.reset();

        std::vector<std::pair<int,Candidate>> result;
        std::vector<std::pair<int,Subcurve>> covering;
        std::vector<double> suffixLengths;

        std::cout << "\nRound " << r;

        for (int i = 0; i < 501; ++i) {
            auto ID = 0;
            auto roundID = i;
            std::cout << ".";
            //update top, until the roundID matches
            while (cs.top().second.roundOfUpdate != roundID) {
                ID++;
                std::pair<int, Candidate> c = cs.top();
                cs.pop();
                swatchupdate.start();
                updateCandidate(curves, c.second, covering, suffixLengths, roundID);
                swatchupdate.stop();
                swatchinsert.start();
                cs.push(c);
                swatchinsert.stop();
            }
            if (lengthOfUncovered(curves,result) <= EPSILON || cs.top().second.semiUpdatedCoverLength <= EPSILON) {
                std::cout << "\nTrying to refine... ";
                int deletecount = 0;
                for(int igni=result.size()-1;igni>=0;--igni){
                    //i is the index to be ignored
                    std::vector<std::pair<int,Candidate>> temp;
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
                        result[igni].second.importance = importance;
                    }
                }
                if (deletecount == 0)
                    std::cout << "nothing";
                std::cout << " greedely deleted\n";
                for(int igni=result.size()-1;igni>=0;--igni){
                    //i is the index to be ignored
                    std::vector<std::pair<int,Candidate>> temp;
                    for(int j=0;j<result.size();++j){
                        if(igni != j){
                            temp.push_back(result[j]);
                        }
                    }
                    double importance = lengthOfUncovered(curves,temp);
                    result[igni].second.importance = importance;
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
                                curves[sub.first].constructSubcurve(sub.second.getStart(), sub.second.getEnd()));
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
                std::vector<std::pair<int, Candidate>> temp;
                for (int j = 0; j < r; ++j) {
                    temp.push_back(cs.top());
                    cs.pop();
                }
                std::pair<int, Candidate> c = cs.top();
                result.push_back(c);

                cs.pop();
                addedweight = c.second.semiUpdatedCoverLength;
                c.second.semiUpdatedCoverLength = 0;
                cs.push(c);
                for (const auto &ctemp: temp)
                    cs.push(ctemp);
            } else {
                std::pair<int, Candidate> c = cs.top();
                result.push_back(c);

                //io::exportSubcurve("/Users/styx/data/curveclustering/results/bestcandidate.txt",curves[0],cs.top().second.getStart(),cs.top().second.getEnd());
                //cs.showCovering(result);

                cs.pop();
                addedweight = c.second.semiUpdatedCoverLength;
                c.second.semiUpdatedCoverLength = 0;
                cs.push(c);
            }

            //std::cout << " added weight: " << addedweight << ". Updated " << ID << " lengths\n";


            std::pair<int, Candidate> c = result.back();

            //update covering
            std::vector<std::pair<int, Subcurve>> temp(covering);
            temp.insert(temp.end(), c.second.matchings.begin(), c.second.matchings.end());
            std::sort(temp.begin(), temp.end(), cmpLeftLower);
            covering.clear();
            suffixLengths.clear();
            suffixLengths.push_back(0);
            std::pair<int, Subcurve> cur = temp[0];
            for (int covI = 1; covI < temp.size(); covI++) {
                std::pair<int, Subcurve> next = temp[covI];
                if ((next.first == cur.first) && (next.second.start < cur.second.end)) {
                    cur.second.end = std::max(next.second.end, cur.second.end);
                } else {
                    covering.push_back(cur);
                    suffixLengths.push_back(
                            suffixLengths.back() + curves[cur.first].subcurve_length(cur.second.start, cur.second.end));
                    cur = next;
                }
            }
            covering.push_back(cur);
            suffixLengths.push_back(
                    suffixLengths.back() + curves[cur.first].subcurve_length(cur.second.start, cur.second.end));
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
              [](const std::pair<int,Candidate>& a,const std::pair<int,Candidate>& b){
                  return (a.second.matchings[0].first < b.second.matchings[0].first) ||
                         ( (a.second.matchings[0].first == b.second.matchings[0].first) && (a.second.matchings[0].second.start < b.second.matchings[0].second.start));
              });
    if(show) {
#ifdef HASVISUAL
        cs.showCovering(bestResultVisualizer);
#endif
    }
    return bestResultVisualizer;
}
