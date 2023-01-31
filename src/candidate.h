//
// Created by Jacobus Conradi on 04.02.22.
//

#ifndef CLUSTERING_CANDIDATE_H
#define CLUSTERING_CANDIDATE_H

#include "freespace.h"
//#include <gmp.h>
//#include <gmpxx.h>
#include <queue>
#include <omp.h>

struct Subcurve{

    Subcurve(ParamPoint s, ParamPoint t, int idx) : start(s), end(t), curveIdx(idx){};

    bool isTrueSubcurve() const{return start < end;};

    bool contains(ParamPoint t) const{if(isTrueSubcurve())return start <= t && t <= end;return start >= t && t >= end;};

    ParamPoint start;
    ParamPoint end;

    int curveIdx;

    ParamPoint getStart() const{return start;}

    ParamPoint getEnd() const{return end;}

    int getIndex() const{return curveIdx;}
};

class Candidate : Subcurve{
public:
    Candidate(ParamPoint s, ParamPoint t, int idx) : Subcurve(s, t, idx){};
    int computeCovering(FreeSpace& fs, int curveindex, Curves& curves);
    void resetCoverLength();
    std::vector<Subcurve> matchings;
    std::vector<Subcurve> visualMatchings;

    bool isDominatedBy(Candidate& other);


    ParamPoint getStart() const{return start;}

    ParamPoint getEnd() const{return end;}

    int getIndex() const{return curveIdx;}

    //priority queue stuff
    distance_t optimisticCoverLength;
    distance_t semiUpdatedCoverLength;
    distance_t importance;
    int roundOfUpdate = 0;
    int debugID = -1;

    //others
    bool inline operator==(Candidate& other){
        return other.start == start && other.end == end;
    };
    static std::vector<Candidate> uncompressAndComputeCovering(std::vector<FreeSpace>& freespaces,ParamPoint start, std::vector<ParamPoint> ends, Curves& curves, int curveIdx, int threadID);
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

auto cmpPQ = [](const Candidate& left, const Candidate& right) { return (left.semiUpdatedCoverLength) < (right.semiUpdatedCoverLength); };
class CandidateSetPQ : std::priority_queue<Candidate,std::vector<Candidate>,decltype(cmpPQ)>{
private:
    using Parent = std::priority_queue<Candidate,std::vector<Candidate>,decltype(cmpPQ)>;
    Curves curves;
    double delta;
    std::vector<std::vector<FreeSpace>> freespaces;
public:
    using Parent::top;
    using Parent::pop;
    using Parent::push;
    CandidateSetPQ(Curves& c,double d);
    void computeCandidates(int l);
    void ultrafastComputeCandidates(int l, int minL = 0);
    void showCovering(std::vector<Candidate> candidates);
    void resetWeights();
    template<typename func> void ultrafastComputeCandidates(int l, func filter) {
        std::cout << "Identifying Important Ys";

//#pragma parallel for
        for(auto & col : freespaces){
            std::cout << ".";
            for(auto& fs : col){
                //with better identify function
                fs.ultrafastIdentify(l);
            }
        }
        std::cout << " Done" << std::endl;


        //FreeSpacesVisualizer fsv(freespaces);
        //fsv.show();

        //std::vector<std::pair<int,Candidate>> baseCandidates;
        std::vector<std::tuple<int,ParamPoint,std::vector<ParamPoint>>> compressedCandidates;

        for (int baseCurveIdx = 0;baseCurveIdx<curves.size();baseCurveIdx++){
            for(int idx1=0;idx1<curves.size();idx1++){
                FreeSpace& fs1 = freespaces[baseCurveIdx][idx1];

                /**
                 * First up
                 * */

                std::cout << "Processing upstarts... ";
                int upcount = 0;
                //identify starts
                std::vector<std::pair<ParamPoint,ParamPoint>> startsAndLimits;
                for(auto ups : fs1.upStarts){
                    if((upcount+1)%100==0 || upcount == fs1.upStarts.size()-1){
                        std::cout << "\33[2K\rProcessing upstarts... "<< upcount+1 << "/" << fs1.upStarts.size();
                    }
                    upcount++;
                    auto limits = fs1.findExtremePointsFromPointWrapper(Direction::Up,Direction::Right,ups,{ups.y.id+l-1,1.0},false);
                    bool trivialStart = ups.y.t == 0.0;
                    if(!trivialStart) {
                        startsAndLimits.push_back({{ups.y},
                                                   {limits[0].y}});
                    }
                    else{
                        if(limits[0].y.id == ups.y.id + l - 1)
                            startsAndLimits.push_back({{ups.y},
                                                       {limits[0].y}});
                    }
                }

                sort(startsAndLimits.begin(),startsAndLimits.end());
                startsAndLimits.erase(std::unique(startsAndLimits.begin(),startsAndLimits.end()),startsAndLimits.end());

                std::cout << "\33[2K\rProcessing upstarts... "<< upcount << "/" << fs1.upStarts.size() << " ... Done";
                std::cout <<std::endl<< "Generated " << startsAndLimits.size() << " Start-Limit-Pairs" << std::endl;
                std::cout << "Generating candidates...";
                int upignores = 0;
                int upcandidatecount = 0;
#pragma omp parallel for shared(startsAndLimits,idx1,baseCurveIdx,upcandidatecount,std::cout,upignores,compressedCandidates) default(none)
                //for(auto sal : startsAndLimits){
                for(int sali=0;sali<startsAndLimits.size();++sali){
                    auto sal = startsAndLimits[sali];
                    if((upcandidatecount+1)%100==0 || upcandidatecount == startsAndLimits.size()-1){
#pragma omp critical (block1)
                        {
                            std::cout << "\33[2K\rGenerating candidates... " << upcandidatecount + 1 << "/" << startsAndLimits.size()
                                      << " pairs processed";
                        }
                    }
                    upcandidatecount++;
                    ParamPoint start = sal.first;
                    ParamPoint limit = sal.second;

                    std::vector<ParamPoint> localEnds;

                    for(int idx2=idx1;idx2<curves.size();idx2++) {
                        FreeSpace &fs2 = freespaces[baseCurveIdx][idx2];
                        //iterate over all cells, where we could potentially start at
                        int startOnT2 = 0;

                        bool firstIteration = true;

                        while (startOnT2 < fs2.xSize()) {

                            //first zoom to first nonempty cell
                            Cell *cell = fs2.getCell(startOnT2, start.id);

                            if (!firstIteration) {
                                //then we have handled the cell, and we need to advance
                                startOnT2++;
                                if (startOnT2 >= fs2.xSize())
                                    break;
                                cell = fs2.getCell(startOnT2, start.id);
                                while (cell->left.begin <= start.t && cell->left.end >= start.t) {
                                    startOnT2++;
                                    if (startOnT2 >= fs2.xSize())
                                        break;
                                    cell = fs2.getCell(startOnT2, start.id);
                                }
                            }

                            if (startOnT2 >= fs2.xSize())
                                break;

                            firstIteration = false;

                            //these breaks can be done way better...
                            while (cell->is_empty() || !(cell->topPair.second.y >= start.t &&
                                                         cell->bottomPair.first.y <= start.t)) {
                                startOnT2++;
                                if (startOnT2 >= fs2.xSize())
                                    break;
                                cell = fs2.getCell(startOnT2, start.id);
                            }
                            if (startOnT2 >= fs2.xSize())
                                break;

                            //cell now points to our startcell

                            CellPoint leftMostInT2Cell = cell->getExtremePointAt<Direction::Left>(start.t);

                            FreeSpacePoint s = {{startOnT2,leftMostInT2Cell.x},{start.id,leftMostInT2Cell.y}};

                            std::vector<FreeSpacePoint> localMaxima;

                            fs2.findExtremePointsFromPointWrapper(Direction::Right,Direction::Up,s,limit,true,omp_get_thread_num(),&localMaxima);
                            {
                                for (auto lm: localMaxima) {
                                    if(start.t==0){
                                        if(lm.y.id == limit.id){
                                            localEnds.push_back(lm.y);
                                            //baseCandidates.push_back({baseCurveIdx, {start, lm.y}});
                                        }else{
                                            upignores++;
                                        }
                                    }else {
                                        localEnds.push_back(lm.y);
                                        //baseCandidates.push_back({baseCurveIdx, {start, lm.y}});
                                    }
                                }
                            }
                        }
                    }
                    //assert(!localEnds.empty());
                    std::sort(localEnds.begin(),localEnds.end());
                    localEnds.erase(std::unique(localEnds.begin(),localEnds.end()),localEnds.end());
                    //assert(!localEnds.empty());
                    if(!localEnds.empty() && localEnds[0] == start){
                        //ignore perfectly straight candidates, as these should be covered by other direction
                        localEnds.erase(localEnds.begin());
                    }
#pragma omp critical (block2)
                    {
                        if(!localEnds.empty())
                            compressedCandidates.emplace_back(baseCurveIdx,start,localEnds);
                    }
                }
                std::cout << "\33[2K\rGenerating candidates... " << startsAndLimits.size() << "/" << startsAndLimits.size()
                          << " pairs processed\n";
                /**
                 * Then down
                 * */

                //identify starts
                std::cout << "Processing downstarts... ";
                int downcount = 0;
                startsAndLimits.clear();
                for(auto downs : fs1.downStarts){
                    if((downcount+1)%100==0 || downcount == fs1.downStarts.size()-1){
                        std::cout << "\33[2K\rProcessing downstarts... "<< downcount+1 << "/" << fs1.downStarts.size();
                    }
                    downcount++;
                    auto limits = fs1.findExtremePointsFromPointWrapper(Direction::Down,Direction::Right,downs,{downs.y.id-l+1,0.0},false);
                    bool trivialStart = downs.y.t == 1.0;
                    if(!trivialStart) {
                        startsAndLimits.push_back({{downs.y},
                                                   {limits[0].y}});
                    }
                    else{
                        if(limits[0].y.id == downs.y.id - l + 1)
                            startsAndLimits.push_back({{downs.y},
                                                       {limits[0].y}});
                    }
                }

                sort(startsAndLimits.begin(),startsAndLimits.end());
                startsAndLimits.erase(std::unique(startsAndLimits.begin(),startsAndLimits.end()),startsAndLimits.end());
                std::vector<std::pair<ParamPoint,ParamPoint>> temp;
                std::cout << "\33[2K\rProcessing downstarts... "<< downcount << "/" << fs1.downStarts.size() << " ... Done";
                std::cout <<std::endl<< "Generated " << startsAndLimits.size() << " Start-Limit-Pairs" << std::endl;
                if(!startsAndLimits.empty()) {
                    ParamPoint cur = startsAndLimits[0].first;
                    ParamPoint best = startsAndLimits[0].second;
                    for (auto sal: startsAndLimits) {
                        if (sal.first != cur) {
                            temp.emplace_back(cur, best);
                            cur = sal.first;
                            best = sal.second;
                        } else {
                            if (best > sal.second) {
                                best = sal.second;
                            }
                        }
                    }
                    temp.emplace_back(cur, best);
                    startsAndLimits = temp;
                }
                std::cout << "Purging dominated limits resulting in " << startsAndLimits.size() << " Start-Limit-Pairs" << std::endl;
                std::cout << "Generating candidates...";
                int downignores;
                int downcandidatecount = 0;
#pragma omp parallel for shared(startsAndLimits,idx1,baseCurveIdx,downcandidatecount,std::cout,downignores,compressedCandidates) default(none)
                //for(auto sal : startsAndLimits){
                for(int sali=0;sali<startsAndLimits.size();++sali){
                    auto sal = startsAndLimits[sali];
                    if((downcandidatecount+1)%100==0 || downcandidatecount == startsAndLimits.size()-1){
#pragma omp critical (block3)
                        {
                            std::cout << "\33[2K\rGenerating candidates... " << downcandidatecount + 1 << "/" << startsAndLimits.size()
                                      << " pairs processed";
                        }
                    }
                    downcandidatecount++;
                    ParamPoint start = sal.first;
                    ParamPoint limit = sal.second;

                    std::vector<ParamPoint> localEnds;

                    for (int idx2 = idx1; idx2 < curves.size(); idx2++) {
                        FreeSpace &fs2 = freespaces[baseCurveIdx][idx2];
                        //iterate over all cells, where we could potentially start at
                        int startOnT2 = 0;

                        bool firstIteration = true;

                        while (startOnT2 < fs2.xSize()) {

                            //first zoom to first nonempty cell
                            Cell *cell = fs2.getCell(startOnT2, start.id);

                            if (!firstIteration) {
                                //then we have handled the cell, and we need to advance
                                startOnT2++;
                                if (startOnT2 >= fs2.xSize())
                                    break;
                                cell = fs2.getCell(startOnT2, start.id);
                                while (cell->left.begin <= start.t && cell->left.end >= start.t) {
                                    startOnT2++;
                                    if (startOnT2 >= fs2.xSize())
                                        break;
                                    cell = fs2.getCell(startOnT2, start.id);
                                }
                            }

                            if (startOnT2 >= fs2.xSize())
                                break;

                            firstIteration = false;

                            //these breaks can be done way better...
                            while (cell->is_empty() || !(cell->topPair.second.y >= start.t &&
                                                         cell->bottomPair.first.y <= start.t)) {
                                startOnT2++;
                                if (startOnT2 >= fs2.xSize())
                                    break;
                                cell = fs2.getCell(startOnT2, start.id);
                            }
                            if (startOnT2 >= fs2.xSize())
                                break;

                            //cell now points to our startcell

                            CellPoint leftMostInT2Cell = cell->getExtremePointAt<Direction::Left>(start.t);

                            FreeSpacePoint s = {{startOnT2, leftMostInT2Cell.x},
                                                {start.id,  leftMostInT2Cell.y}};

                            std::vector<FreeSpacePoint> localMaxima;

                            fs2.findExtremePointsFromPointWrapper(Direction::Right, Direction::Down, s, limit, true, omp_get_thread_num(),
                                                                  &localMaxima);

                            for (auto lm: localMaxima) {
                                if(start.t==0){
                                    if(lm.y.id == limit.id){
                                        localEnds.push_back(lm.y);
                                        //baseCandidates.push_back({baseCurveIdx, {start, lm.y}});
                                    }else{
                                        downignores++;
                                    }
                                }else {
                                    localEnds.push_back(lm.y);
                                    //baseCandidates.push_back({baseCurveIdx, {start, lm.y}});
                                }
                            }
                        }
                    }
                    //assert(!localEnds.empty());
                    std::sort(localEnds.begin(),localEnds.end());
                    localEnds.erase(std::unique(localEnds.begin(),localEnds.end()),localEnds.end());
                    std::reverse(localEnds.begin(),localEnds.end());
                    //assert(!localEnds.empty());
                    if(!localEnds.empty() && localEnds[0] == start){
                        //ignore perfectly straight candidates, as these should be covered by other direction
                        localEnds.erase(localEnds.begin());
                    }
#pragma omp critical (block4)
                    {
                        if(!localEnds.empty())
                            compressedCandidates.emplace_back(baseCurveIdx,start,localEnds);
                    }
                }

                std::cout << "\33[2K\rGenerating candidates... " << startsAndLimits.size() << "/" << startsAndLimits.size()
                          << " pairs processed\n";
            }
        }
        std::cout << "Generated " << compressedCandidates.size() << " many compressed candidates."<< std::endl << "Augmenting candidates with covering data...";
        int progress = 0;
        std::vector<Candidate> prepruned;
        //alternative:
        int removes = 0;
        //TODO: fix race condicitions as one thread may reset another threads freespace
#pragma omp parallel for shared(compressedCandidates,freespaces,progress,std::cout, prepruned, removes, filter) default(none)
        //for(auto& cC : compressedCandidates){
        for(int cCi=0;cCi<compressedCandidates.size();++cCi){
            auto cC = compressedCandidates[cCi];
            progress++;
            if(progress % 100==0)
                std::cout << "\33[2K\rUncompressing candidates and augmenting with covering data... " << progress << "/" << compressedCandidates.size() << ", resulting in " << size() << " candidates so far. Removed "<< removes << " filtered candidates so far... ";
            if(std::get<2>(cC).empty())
                continue;
            std::vector<Candidate> result = Candidate::uncompressAndComputeCovering(freespaces[std::get<0>(cC)],std::get<1>(cC),std::get<2>(cC),curves, std::get<0>(cC), omp_get_thread_num());
            for(auto& c: result){
                if(filter(c)) {
#pragma omp critical (block6)
                    {
                        push(c);
                    }
                }else{
#pragma omp critical (block7)
                    {
                        removes++;
                    }
                }
            }
        }

        std::cout << " Done.\nCandidate generation done.";
    }


};

#endif //CLUSTERING_CANDIDATE_H
