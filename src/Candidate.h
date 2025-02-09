//
// Created by Jacobus Conradi on 18.04.23.
//

#ifndef CLUSTERING_CANDIDATE_H
#define CLUSTERING_CANDIDATE_H

#include "geometry_types.h"
#include "Curve.h"
#include "SparseFreespace.h"
#include <queue>
#include <omp.h>


class tinyCandidate:CInterval{
public:
    std::vector<CInterval> matching;

    using CInterval::CInterval;
    using CInterval::getCurveIndex;
    using CInterval::getEnd;
    using CInterval::getBegin;

    distance_t optimisticCoverLength=0;
    distance_t semiUpdatedCoverLength=0;
    distance_t importance=0;
    int roundOfUpdate = 0;

    void resetCoverLength() {
        semiUpdatedCoverLength = optimisticCoverLength;
        roundOfUpdate = 0;
    }

};

class Candidate:CInterval{
public:
    //need explicit copy cosntructor for std::move
    //Candidate(Candidate& r) = default;

    std::vector<CInterval> matching;
    std::vector<CInterval> visualMatching;

    using CInterval::CInterval;
    using CInterval::getCurveIndex;
    using CInterval::getEnd;
    using CInterval::getBegin;

    distance_t optimisticCoverLength=0;
    distance_t semiUpdatedCoverLength=0;
    distance_t importance=0;
    int roundOfUpdate = 0;

    void resetCoverLength() {
        semiUpdatedCoverLength = optimisticCoverLength;
        roundOfUpdate = 0;
    }
};

class SparseFreeSpaces:std::vector<std::vector<SparseFreespace>>{
private:
    using Parent = std::vector<std::vector<SparseFreespace>>;//THIS IS POTENTIALLY NOT SORTED
public:
    using Parent::begin;
    using Parent::end;
    using Parent::size;
    using Parent::operator[];
    SparseFreeSpaces(Curves& c, double d, int tc=1);
};

//typedef std::vector<std::vector<SparseFreespace>> SparseSparseFreeSpaces;


static auto cmpPQ = [](const Candidate& left, const Candidate& right) { return (left.semiUpdatedCoverLength) < (right.semiUpdatedCoverLength); };
class CandidateSetPQ : std::priority_queue<Candidate,std::vector<Candidate>,decltype(cmpPQ)>{
private:
    using Parent = std::priority_queue<Candidate,std::vector<Candidate>,decltype(cmpPQ)>;
    Curves curves;
    double delta;

    std::vector<Candidate> uncompressCandidate(CurveID bIndex, CPoint start, const CPoints& ends, int threadID=0);
public:
    using Parent::top;
    using Parent::pop;
    using Parent::push;
    using Parent::empty;
    using Parent::size;
    using Parent::swap;
    CandidateSetPQ(Curves& c,double d);

    void showFreespaces();


    //std::vector<std::vector<SparseFreespace>> freespaces;
    SparseFreeSpaces sparsefreespaces;
    template<typename func>
    void ultrafastComputeSmall(int l, func filter){
        for (auto &fss: sparsefreespaces) {
            for (auto &fs: fss) {
                fs.identifyStartsAndEnds();
            }
        }

        std::vector<CCPoint> upAggregated,downAggregated;
        for (auto &fss: sparsefreespaces) {
            for (auto &fs: fss) {
                for(auto p : fs.upStarts){
                    upAggregated.emplace_back(fs.BID,p.getPoint(),p.getFraction());
                }
                for(auto p : fs.upEnds){
                    upAggregated.emplace_back(fs.BID,p.getPoint(),p.getFraction());
                }
                for(auto p : fs.downStarts){
                    downAggregated.emplace_back(fs.BID,p.getPoint(),p.getFraction());
                }
                for(auto p : fs.downEnds){
                    downAggregated.emplace_back(fs.BID,p.getPoint(),p.getFraction());
                }
                for(int i = 0;i<fs.ySize()-1;i++){
                    upAggregated.emplace_back(fs.BID,i,0);
                    upAggregated.emplace_back(fs.BID,i,1);
                    downAggregated.emplace_back(fs.BID,i,0);
                    downAggregated.emplace_back(fs.BID,i,1);
                }
            }
        }

        std::sort(upAggregated.begin(), upAggregated.end());
        upAggregated.erase(std::unique(upAggregated.begin(), upAggregated.end()), upAggregated.end());

        std::sort(downAggregated.begin(), downAggregated.end());
        downAggregated.erase(std::unique(downAggregated.begin(), downAggregated.end()), downAggregated.end());

        std::cout << upAggregated.size() << " many critical up points identified.\n" << std::flush;
        std::cout << downAggregated.size() << " many critical down points identified.\n" << std::flush;

        //std::vector<Candidate> cans;

        //UP!
//#pragma omp parallel for default(none) shared(upAggregated)

        std::vector<Candidate> localSet;
#pragma omp declare reduction (merge : std::vector<Candidate> : omp_out.insert(omp_out.end(), std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))
#pragma omp parallel for default(none) shared(l,filter,upAggregated) reduction(merge: localSet)
        for(int i=0;i<upAggregated.size();i++){//auto start : upAggregated){
            auto start = upAggregated[i];
            std::vector<CPoint> ends;

            //populate ends
            for(int offset = 1;i+offset < upAggregated.size();offset *= 2){
                auto potentialEnd = upAggregated[i+offset];
                if(potentialEnd.getCurve() == start.getCurve() && potentialEnd.getPoint() <= start.getPoint()+l){
                    ends.push_back(potentialEnd.getCPoint());
                }else{
                    break;
                }
            }

            if(ends.empty()){
                continue;
            }

            //uncompress
            auto uncompressedCandidates = uncompressCandidate(start.getCurve(),start.getCPoint(),ends,omp_get_thread_num());


            //Step 5: for every candidate
            //      filter stage 2?
            //      compute the actual matching
            //      compute the weight of that matching
            for (auto ucC: uncompressedCandidates) {
                if (!ucC.visualMatching.empty()) {
                    CInterval cur;
                    for (auto m: ucC.visualMatching) {
                        if (cur.is_empty()) {
                            cur = m;
                            continue;
                        }
                        if (cur.fixed_curve == m.fixed_curve && cur.contains(m.begin)) {
                            cur.end = m.end;
                        } else {
                            ucC.matching.push_back(cur);
                            ucC.optimisticCoverLength += curves[cur.fixed_curve].subcurve_length(cur.begin, cur.end);
                            cur = m;
                        }
                    }
                    ucC.matching.push_back(cur);
                    ucC.optimisticCoverLength += curves[cur.fixed_curve].subcurve_length(cur.begin, cur.end);
                    ucC.semiUpdatedCoverLength = ucC.optimisticCoverLength;
                    if (filter(ucC)) {
                        //cans.push_back(ucC);
                        localSet.push_back(ucC);
                    }
                } else {
                    assert(false);
                }
            }
        }

        this->c.insert(this->c.end(),std::make_move_iterator(localSet.begin()),std::make_move_iterator(localSet.end()));
        std::cout << "Generated " << size() << " many up candidates. \n"<<std::flush;
        //std::cout << "comp: " << comp<<" , "<<(double)(comp)/(double)(count)<<"\n"<<std::flush;

        //DOWN
        std::vector<Candidate> localSet2;

#pragma omp parallel for default(none) shared(l,filter,downAggregated) reduction(merge: localSet2)
        for(int i=0;i<downAggregated.size();i++){
            auto start = downAggregated[i];
            std::vector<CPoint> ends;

            //populate ends
            for(int offset = 1;i-offset >= 0;offset *= 2){
                auto potentialEnd = downAggregated[i-offset];
                if(potentialEnd.getCurve() == start.getCurve() && potentialEnd.getPoint() == start.getPoint()){
                    ends.push_back(potentialEnd.getCPoint());
                }else{
                    break;
                }
            }

            if(ends.empty()){
                continue;
            }

            //uncompress
            auto uncompressedCandidates = uncompressCandidate(start.getCurve(),start.getCPoint(),ends,omp_get_thread_num());


            //Step 5: for every candidate
            //      filter stage 2?
            //      compute the actual matching
            //      compute the weight of that matching
            for (auto ucC: uncompressedCandidates) {
                if (!ucC.visualMatching.empty()) {
                    CInterval cur;
                    for (auto m: ucC.visualMatching) {
                        if (cur.is_empty()) {
                            cur = m;
                            continue;
                        }
                        if (cur.fixed_curve == m.fixed_curve && cur.contains(m.begin)) {
                            cur.end = m.end;
                        } else {
                            ucC.matching.push_back(cur);
                            ucC.optimisticCoverLength += curves[cur.fixed_curve].subcurve_length(cur.begin, cur.end);
                            cur = m;
                        }
                    }
                    ucC.matching.push_back(cur);
                    ucC.optimisticCoverLength += curves[cur.fixed_curve].subcurve_length(cur.begin, cur.end);
                    ucC.semiUpdatedCoverLength = ucC.optimisticCoverLength;
                    if (filter(ucC)) {
                        //cans.push_back(ucC);
                        localSet2.push_back(ucC);
                    }
                } else {
                    assert(false);
                }
            }
        }

        this->c.insert(this->c.end(),std::make_move_iterator(localSet2.begin()),std::make_move_iterator(localSet2.end()));
        std::make_heap(this->c.begin(), this->c.end(),cmpPQ);

        std::cout << "Generated " << size() << " many up and down candidates.\n"<<std::flush;
        //std::sort(cans.begin(), cans.end(), [](Candidate &l, Candidate &r) { return l.getBegin() < r.getBegin(); });

    }


    template<typename func> void ultrafastCompute(int l, func filter) {
        ultrafastComputeSmall(l,filter);
        //Step 1: prepare freespaces
        for (auto &fss: sparsefreespaces) {
            for (auto &fs: fss) {
                fs.identifyStartsAndEnds();
            }
        }



        //Step 2: aggregate all starts and ends, and remove duplicates
        std::vector<CPoints> upStarts, downStarts, upEnds, downEnds;
        for (auto &fss: sparsefreespaces) {
            upStarts.emplace_back();
            downStarts.emplace_back();
            upEnds.emplace_back();
            downEnds.emplace_back();
            for (auto &fs: fss) {
                upStarts.back().insert(upStarts.back().end(), fs.upStarts.begin(), fs.upStarts.end());
                upEnds.back().insert(upEnds.back().end(), fs.upEnds.begin(), fs.upEnds.end());
                downStarts.back().insert(downStarts.back().end(), fs.downStarts.begin(), fs.downStarts.end());
                downEnds.back().insert(downEnds.back().end(), fs.downEnds.begin(), fs.downEnds.end());
            }
        }
        int uS = 0, uE = 0, dS = 0, dE = 0;
        for (auto &starts: upStarts) {
            uS += starts.size();
            std::sort(starts.begin(), starts.end());
            starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
        }
        for (auto &starts: downStarts) {
            dS += starts.size();
            std::sort(starts.begin(), starts.end());
            starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
        }
        for (auto &ends: upEnds) {
            uE += ends.size();
            std::sort(ends.begin(), ends.end());
            ends.erase(std::unique(ends.begin(), ends.end()), ends.end());
        }
        for (auto &ends: downEnds) {
            dE += ends.size();
            std::sort(ends.begin(), ends.end());
            ends.erase(std::unique(ends.begin(), ends.end()), ends.end());
        }

        std::cout << "Generated " << uS << " many upstarts and " << uE << " many upends\n"<<std::flush;
        std::cout << "Generated " << dS << " many downstarts and " << dE << " many downends\n"<<std::flush;


        //Step 3: for every start
        //      for every applicable end
        //              create compressed candidate
        int upcC = 0, downcC = 0;
        int upavgC = 0, downavgC = 0;
        std::vector<std::tuple<ID<Curve>, CPoint, std::vector<CPoint>::iterator, std::vector<CPoint>::iterator>> compressedCandidates; //<curveID, startpoint, endpoints>
        for (int id = 0; id < curves.size(); id++) {
            //std::cout << "Processing " << id << "/" << curves.size() << " up compressed candidates" << std::endl;
            for (auto &start: upStarts[id]) {
                auto leftBound = std::lower_bound(upEnds[id].begin(), upEnds[id].end(), start);
                CPoint cutoff = {start.getPoint() + l, 0.0};
                auto rightBound = std::upper_bound(upEnds[id].begin(), upEnds[id].end(), cutoff);
                compressedCandidates.emplace_back(id, start, leftBound, rightBound);
                upcC++;
                upavgC += rightBound - leftBound;
            }
        }

        //std::cout << "Done with up!\n";

#pragma omp parallel for default(none) shared(curves,std::cout,upStarts,upEnds,l,upcC,upavgC,compressedCandidates) schedule(dynamic)
        for (int id = 0; id < curves.size(); id++) {
            for (PointID start = 0; start == 0 || start+l < curves[id].size(); start += 1) {
                CPoint lowercutoff = {std::min((unsigned int) (curves[id].size() - 2), start + l - 1), 0.0};
                auto leftBound = std::lower_bound(upEnds[id].begin(), upEnds[id].end(), lowercutoff);
                CPoint uppercutoff = {std::min((unsigned int) (curves[id].size() - 2), start + l), 0.0};
                auto rightBound = std::upper_bound(upEnds[id].begin(), upEnds[id].end(), uppercutoff);
                if (leftBound != rightBound) {
                    CPoint s = {start, 0.0};
#pragma omp critical
//                    {
                    compressedCandidates.emplace_back(id, s, leftBound, rightBound);
  //              };
                }
            }
        }
        //std::cout << "Done with max length up candidates!\n";

        //TODO: these are downends
//#pragma omp parallel for default(none) shared(curves,std::cout,downStarts,downEnds,l,downcC,downavgC,compressedCandidates) schedule(dynamic)
        for (int id = 0; id < curves.size(); id++) {
//            std::cout << "Processing " << id << "/" << curves.size() << " down compressed candidates" << std::endl;
            for (auto &start: downStarts[id]) {
                CPoint cutoff = {start.getPoint(), 0.0};
                auto leftBound = std::lower_bound(downEnds[id].begin(), downEnds[id].end(), cutoff);
                auto rightBound = std::upper_bound(downEnds[id].begin(), downEnds[id].end(), start);
//#pragma omp critical
//                {
                    compressedCandidates.emplace_back(id, start, leftBound, rightBound);
                    downcC++;
                    downavgC += rightBound - leftBound;
//                };
            }
        }

        std::cout << "Generated " << upcC << " many compressed upcandidates with average size "
                  << (double) (upavgC) / upcC << "  and " << downcC
                  << " many compressed down candidates with average size " << (double) (downavgC) / downcC << " \n";


        std::vector<Candidate> cans;
        //Step 4: for every start
        //      uncompress every candidate that has not been filtered
        for (auto cC: compressedCandidates) {
            CPoints ends(get<2>(cC), get<3>(cC));
            if (ends.empty() || ends.back() > get<1>(cC)) {
                //upcandidate
                if (ends.empty() || (get<1>(cC).getPoint() + l - 1) != ends.back().getPoint() ||
                    ends.back().getFraction() != 1.0) {
                    ends.emplace_back(
                            std::min((unsigned int) (curves[get<0>(cC)].size() - 2), get<1>(cC).getPoint() + l - 1),
                            1.0);
                }
            } else {
                //downcandidate
                //canonical end might not be reachable...

                //if((get<1>(cC).getPoint()) != ends.back().getPoint() || ends.back().getFraction() != 0.0){
                //    ends.emplace_back(get<1>(cC).getPoint(),0.0);
                //}
            }
            auto uncompressedCandidates = uncompressCandidate(get<0>(cC), get<1>(cC), ends);

            //Step 5: for every candidate
            //      filter stage 2?
            //      compute the actual matching
            //      compute the weight of that matching
            for (auto ucC: uncompressedCandidates) {
                if (!ucC.visualMatching.empty()) {
                    CInterval cur;
                    for (auto m: ucC.visualMatching) {
                        if (cur.is_empty()) {
                            cur = m;
                            continue;
                        }
                        if (cur.fixed_curve == m.fixed_curve && cur.contains(m.begin)) {
                            cur.end = m.end;
                        } else {
                            ucC.matching.push_back(cur);
                            ucC.optimisticCoverLength += curves[cur.fixed_curve].subcurve_length(cur.begin, cur.end);
                            cur = m;
                        }
                    }
                    ucC.matching.push_back(cur);
                    ucC.optimisticCoverLength += curves[cur.fixed_curve].subcurve_length(cur.begin, cur.end);
                    ucC.semiUpdatedCoverLength = ucC.optimisticCoverLength;
                    if (filter(ucC)) {
                        cans.push_back(ucC);
                        push(ucC);
                    }
                } else {
                    assert(false);
                }
            }
        }
        std::cout << "Generated " << size() << " many candidates out of " << compressedCandidates.size()
                  << " many compressed Candidates.\n";
        std::sort(cans.begin(), cans.end(), [](Candidate &l, Candidate &r) { return l.getBegin() < r.getBegin(); });
    }

    void resetWeights() {
        std::vector<Candidate> temp;
        while(!empty()) {
            temp.push_back(top());
            pop();
        }
        for(auto& c : temp){
            c.resetCoverLength();
            push(c);
        }
    }

    void showCovering(std::vector<Candidate> candidates);
};


#endif //CLUSTERING_CANDIDATE_H
