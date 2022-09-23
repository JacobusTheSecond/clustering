//
// Created by Jacobus Conradi on 04.02.22.
//

#ifndef CLUSTERING_CANDIDATE_H
#define CLUSTERING_CANDIDATE_H

#include "freespace.h"
//#include <gmp.h>
//#include <gmpxx.h>
#include <queue>
struct Subcurve{

    Subcurve(ParamPoint s, ParamPoint t) : start(s), end(t){};

    bool isTrueSubcurve() const{return start < end;};

    bool contains(ParamPoint t) const{if(isTrueSubcurve())return start <= t && t <= end;return start >= t && t >= end;};;;

    ParamPoint start;
    ParamPoint end;

    ParamPoint getStart() const{return start;}

    ParamPoint getEnd() const{return end;}
};

class Candidate : Subcurve{
public:
    int ID;
    //mpz_t wback;
    //mpz_t w;
    //mpz_t cW;
    Candidate(ParamPoint s, ParamPoint t) : Subcurve(s, t){};// {mpz_init_set_ui(w,1);mpz_init_set_ui(wback,1);mpz_init(cW);};
    int computeCovering(FreeSpace& fs, int curveindex, Curves& curves);
    void resetCoverLength();
    std::vector<std::pair<int,Subcurve>> matchings;
    std::vector<std::pair<int,Subcurve>> visualMatchings;

    bool isDominatedBy(Candidate& other);

    //priority queue stuff
    distance_t optimisticCoverLength;
    distance_t semiUpdatedCoverLength;
    int roundOfUpdate = 0;
    int debugID = -1;

    //others
    ParamPoint getStart() const {return start; };
    ParamPoint getEnd() const {return end; };
    bool inline operator==(Candidate& other){
        return other.start == start && other.end == end;
    };
    bool covers(int curveId, ParamPoint t);
    static std::vector<Candidate> uncompressAndComputeCovering(std::vector<FreeSpace>& freespaces,ParamPoint start, std::vector<ParamPoint> ends, Curves& curves, int threadID=0);
};

class CandidateSet : std::vector<std::vector<Candidate>>{
private:
    Curves curves;
    double delta;
    std::vector<std::vector<FreeSpace>> freespaces;
public:
    //TODO: multithread
    //int totalWeight;
    //mpz_t tW;
    int totalCount = 0;
    CandidateSet(Curves& c,double d);
    //void computeCandidates();
    //void computeCandidates(int l);
    void showCandidates();
    void showCandidate(std::pair<int,int> t);
    Candidate & getCandidate(std::pair<int,int> t){return this->operator[](t.first)[t.second];};
    std::pair<int,ParamPoint> findNonCovered(std::vector<std::pair<int,int>> indices);
    void updateConcurrentWeights();
    void updateConcurrentWeightsForRedraw();
    void freshenWeights();
    void reweight(std::pair<int,ParamPoint> t);
    //std::pair<int,int> draw(mpz_t rand);
    //std::pair<int,int> drawWithoutRedraw(mpz_t rand);
    void showCovering(std::vector<std::pair<int,int>> indices);
    void resetWeights();
};

auto cmpPQ = [](const std::pair<int,Candidate>& left, const std::pair<int,Candidate>& right) { return (left.second.semiUpdatedCoverLength) < (right.second.semiUpdatedCoverLength); };
class CandidateSetPQ : std::priority_queue<std::pair<int,Candidate>,std::vector<std::pair<int,Candidate>>,decltype(cmpPQ)>{
private:
    Curves curves;
    double delta;
    std::vector<std::vector<FreeSpace>> freespaces;
public:
    using std::priority_queue<std::pair<int,Candidate>,std::vector<std::pair<int,Candidate>>,decltype(cmpPQ)>::top;
    using std::priority_queue<std::pair<int,Candidate>,std::vector<std::pair<int,Candidate>>,decltype(cmpPQ)>::pop;
    using std::priority_queue<std::pair<int,Candidate>,std::vector<std::pair<int,Candidate>>,decltype(cmpPQ)>::push;
    CandidateSetPQ(Curves& c,double d);
    void computeCandidates(int l);
    void ultrafastComputeCandidates(int l, int minL = 0);
    void showCovering(std::vector<std::pair<int,Candidate>> candidates);
    void resetWeights();

};

#endif //CLUSTERING_CANDIDATE_H
