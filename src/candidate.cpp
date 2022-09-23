//
// Created by Jacobus Conradi on 04.02.22.
//

#include <omp.h>

#include <utility>
#include "candidate.h"
#ifdef HASVISUAL
#include "free_space_visualizer.h"
#endif
#include <algorithm>

//TODO: improve zoop
int Candidate::computeCovering(FreeSpace &fs, int curveIndex, Curves& curves) {

    //FreeSpaceVisualizer fsv = FreeSpaceVisualizer(fs);

    int count = 0;

    int startOnT = 0;

    bool firstIteration = true;

    while (startOnT < fs.xSize()) {

        //first zoom to first nonempty cell
        Cell *cell = fs.getCell(startOnT, start.id);

        if(!firstIteration){
            //then we have handled the cell, and we need to advance
            startOnT++;
            if (startOnT >= fs.xSize())
                break;
            cell = fs.getCell(startOnT, start.id);
            while(cell->left.begin <= start.t && cell->left.end >= start.t){
                startOnT++;
                if (startOnT >= fs.xSize())
                    break;
                cell = fs.getCell(startOnT, start.id);
            }
        }

        firstIteration = false;


        if (startOnT >= fs.xSize())
            break;

        //these breaks can be done way better...
        while (cell->is_empty() || !(cell->topPair.second.y >= start.t &&
                                     cell->bottomPair.first.y <= start.t)) {
            startOnT++;
            if (startOnT >= fs.xSize())
                break;
            cell = fs.getCell(startOnT, start.id);
        }
        if (startOnT >= fs.xSize())
            break;

        //cell now points to our startcell

        CellPoint leftMostInT2Cell = cell->getExtremePointAt<Direction::Left>(start.t);
        FreeSpacePoint leftMostInCellFSP = {{startOnT, leftMostInT2Cell.x},
                                            {start.id,  leftMostInT2Cell.y}};


        std::vector<FreeSpacePoint> localExtremes;
        if(isTrueSubcurve())
            localExtremes = fs.findExtremePointsFromPointWrapper(
                Direction::Up, Direction::Right, leftMostInCellFSP, end,false);
        else
            localExtremes = fs.findExtremePointsFromPointWrapper(Direction::Down, Direction::Right,leftMostInCellFSP,end,false);

        if(localExtremes.back().y != end){
            continue;
        }

        //TODO: fix this, but for visualization, it is better to push back every item
        //I think the line matchings.back().second.end = localExtremes.back().x; does not properly update ...second.end
        /*if(!matchings.empty() && matchings.back().first == curveIndex && matchings.back().second.end >= leftMostInCellFSP.x && matchings.back().second.end <= localExtremes.back().x)
            matchings.back().second.end = localExtremes.back().x;
        else {*/
        if(!matchings.empty() && matchings.back().first == curveIndex && matchings.back().second.end >= localExtremes.back().x){
            //pass
        }else{
            matchings.push_back({curveIndex, {leftMostInCellFSP.x, localExtremes.back().x}});
            count++;
            //fsv.show(leftMostInCellFSP,localExtremes.back(),start,end);
        }

    }

    //compute optimisticCoverLength
    optimisticCoverLength = 0;
    for (auto matching : matchings){
        auto s = matching.second.start;
        auto t = matching.second.end;
        optimisticCoverLength += curves[matching.first].subcurve_length(s,t);
        }
    semiUpdatedCoverLength = optimisticCoverLength;
    return count;
}

//mimic a<b
bool comp(const std::pair<int,Subcurve>& a,const std::pair<int,Subcurve> & b){
    return (a.first < b.first) || (a.first == b.first && a.second.end < b.second.start);
}
bool Candidate::covers(int curveId, ParamPoint t) {
    std::pair<int,Subcurve> val = {curveId,{t,t}};
    return std::binary_search(matchings.begin(), matchings.end(),val,comp);
}

void Candidate::resetCoverLength() {
    semiUpdatedCoverLength = optimisticCoverLength;
    roundOfUpdate = 0;
}

std::vector<Candidate> Candidate::uncompressAndComputeCovering(std::vector<FreeSpace> &freespaces, ParamPoint start,
                                                               std::vector<ParamPoint> ends, Curves& curves, int threadID) {
    std::vector<Candidate> result;
    result.reserve(ends.size());
    for(auto& end : ends){
        result.emplace_back(start,end);
    }
    assert(!ends.empty());

    assert(start != ends.back());
    bool willBeTrueSubcurves = start < ends.back();

    for(int tidx = 0;tidx<freespaces.size();tidx++) {

        FreeSpace& fs = freespaces[tidx];

        int startOnT = 0;

        bool firstIteration = true;

        while (startOnT < fs.xSize()) {

            //first zoom to first nonempty cell
            Cell *cell = fs.getCell(startOnT, start.id);

            if (!firstIteration) {
                //then we have handled the cell, and we need to advance
                startOnT++;
                if (startOnT >= fs.xSize())
                    break;
                cell = fs.getCell(startOnT, start.id);
                while (cell->left.begin <= start.t && cell->left.end >= start.t) {
                    startOnT++;
                    if (startOnT >= fs.xSize())
                        break;
                    cell = fs.getCell(startOnT, start.id);
                }
            }

            firstIteration = false;


            if (startOnT >= fs.xSize())
                break;

            //these breaks can be done way better...
            while (cell->is_empty() || !(cell->topPair.second.y >= start.t &&
                                         cell->bottomPair.first.y <= start.t)) {
                startOnT++;
                if (startOnT >= fs.xSize())
                    break;
                cell = fs.getCell(startOnT, start.id);
            }
            if (startOnT >= fs.xSize())
                break;

            //cell now points to our startcell

            CellPoint leftMostInT2Cell = cell->getExtremePointAt<Direction::Left>(start.t);
            FreeSpacePoint leftMostInCellFSP = {{startOnT, leftMostInT2Cell.x},
                                                {start.id, leftMostInT2Cell.y}};


            std::vector<FreeSpacePoint> localExtremes;

            if (willBeTrueSubcurves)
                localExtremes = fs.findExtremePointsFromPointUntilEndsWrapper(Direction::Up,leftMostInCellFSP,ends, threadID);
            else
                localExtremes = fs.findExtremePointsFromPointUntilEndsWrapper(Direction::Down,leftMostInCellFSP,ends, threadID);

            for(int i=0;i<ends.size();++i){
                Candidate& c = result[i];
                FreeSpacePoint& end = localExtremes[i];
                if(end.x.id == -1){
                    continue;
                }else{
                    //for visualization we push regardless. This possibly infers errors in the weight computation...
                    c.visualMatchings.push_back({tidx,{leftMostInCellFSP.x,end.x}});

                    if(c.matchings.empty() || c.matchings.back().first != tidx || c.matchings.back().second.end < leftMostInCellFSP.x){
                        c.matchings.push_back({tidx,{leftMostInCellFSP.x,end.x}});
                    }else{
                        c.matchings.back().second.end = std::max(c.matchings.back().second.end,end.x);
                    }
                }
            }
        }
    }

    //compute optimisticCoverLength
    for(auto& r : result) {
        r.optimisticCoverLength = 0;
        for (auto matching: r.matchings) {
            auto s = matching.second.start;
            auto t = matching.second.end;
            r.optimisticCoverLength += curves[matching.first].subcurve_length(s, t);
        }
        r.semiUpdatedCoverLength = r.optimisticCoverLength;
    }
    return result;
}

bool Candidate::isDominatedBy(Candidate &other) {
    int oMI = 0;
    for(auto& m : matchings){
        //go to first maching of other, whose end is further right
        while(oMI < other.matchings.size() && (other.matchings[oMI].first < m.first || (other.matchings[oMI].first == m.first && other.matchings[oMI].second.end < m.second.end)))
            oMI++;
        if(oMI >= other.matchings.size())
            return false;
        //now other.end lies right of this.end
        if( !(m.first == other.matchings[oMI].first && other.matchings[oMI].second.start <= m.second.start))
            return false;
    }
    return true;
}

CandidateSet::CandidateSet(Curves &c, double d) : curves(c) {

    delta = d;

    for(auto & curve_i : curves) {
        freespaces.emplace_back();
        freespaces.back().reserve(curves.size());
    }
    std::cout << "Computing Freespaces";
#pragma omp parallel for default(none) shared(freespaces,curves,std::cout)
    for(int i=0;i<freespaces.size();++i){
        Curve& curve_i = curves[i];
        std::vector<FreeSpace>& col = freespaces[i];
        for(auto & curve_j : curves){
            col.emplace_back(curve_i,curve_j,delta);
        }
        //freespaces.push_back(col);
        std::cout << ".";
    }
    std::cout << " Done\n";
}

//mimic a<b at first coordinate
bool comp2(std::pair<int,Subcurve> a, std::pair<int,Subcurve> b){
    return (a.first < b.first) || (a.first == b.first && a.second.start < b.second.start);
}

ParamPoint halfPointInbetween(ParamPoint a, ParamPoint b){
    assert(a<b);
    //we can interpret a Parampoint as a.id being the whole number, and a.t the decinal places
    bool increaseDigit = a.t + b.t >= 1.0;
    int twoTId = a.id + b.id + (increaseDigit?1:0);
    double twoTT = a.t + b.t - (increaseDigit?1.0:0.0);
    bool uneven = twoTId % 2 != 0;
    ParamPoint r{};
    r.id = twoTId / 2;
    r.t = twoTT / 2.0 + (uneven?0.5:0.0);
    assert(a < r && r < b);
    return r;
}

std::pair<int,ParamPoint> CandidateSet::findNonCovered(std::vector<std::pair<int,int>> indices) {
    std::vector<std::pair<int,Subcurve>> matchings;
    int count = 0;
    for(auto index : indices){
        count += getCandidate(index).matchings.size();
    }
    matchings.reserve(count);

    for(auto index : indices){
        matchings.insert(matchings.end(),getCandidate(index).matchings.begin(),getCandidate(index).matchings.end());
    }
    std::sort(matchings.begin(), matchings.end(),comp2);

    //sweep through the flattend array, finding the last interval of the initial overlap
    int currentId = 0;
    ParamPoint currentT = {0,0.0};
    for(auto & matching:matchings){
        //go to next curve
        if(currentId != matching.first){
            if(currentT.id == curves[currentId].size()-2 && currentT.t >= 1.0){
                currentId++;
                currentT = {0,0.0};
                continue;
            }else{
                ParamPoint max = {static_cast<int>(curves[currentId].size()-2),1.0};
                return {currentId, halfPointInbetween(currentT, max)};
            }
        }
        if(currentT < matching.second.start){
            return {currentId, halfPointInbetween(currentT, matching.second.start)};
        }
        currentT = std::max(matching.second.end,currentT);
    }
    if(currentT.t < 1.0){
        return {currentId, halfPointInbetween(currentT,{currentT.id,1.0})};
    }
    if(currentT.id != static_cast<int>(curves[currentId].size()-2)){
        return {currentId, halfPointInbetween(currentT,{static_cast<int>(curves[currentId].size()-2),1.0})};
    }
    if(currentId != curves.size()-1){
        return {currentId+1,{0,0.0}};
    }
    return {-1,{-1,-1.0}};
}
/*
void CandidateSet::reweight(std::pair<int,ParamPoint> t) {
    mpz_t preTotalWeight;
    mpz_init_set(preTotalWeight,tW);
    int count=0;
    //totalWeight = 0;
    mpz_set_ui(tW,0);
    for (auto & col : *this){
        for(auto & candidate : col){
            if(candidate.covers(t.first,t.second)){
                mpz_mul_ui(candidate.w,candidate.w,2);
                mpz_set(candidate.wback,candidate.w);
                //candidate.weight*=2;
                //assert(mpz_cmp_ui(candidate.w,candidate.weight) == 0);
                count++;
            }
            //candidate.concurrentWeight = totalWeight + candidate.weight;
            mpz_add(candidate.cW,tW,candidate.w);
            //assert(mpz_cmp_ui(candidate.cW,candidate.concurrentWeight) == 0);
            //totalWeight += candidate.weight;
            mpz_add(tW,tW,candidate.w);
            //assert(mpz_cmp_ui(tW,totalWeight) == 0);
        }
    }
    //std::cout << "Updated " << count << " candidates, increasing the total weight from " << preTotalWeight << " to " << tW << "\n";
}

void CandidateSet::updateConcurrentWeights() {
    //totalWeight = 0;
    mpz_set_ui(tW,0);
    for(auto & col : *this){
        for(auto & candidate : col){
            //candidate.concurrentWeight = totalWeight + candidate.weight;
            mpz_add(candidate.cW, tW, candidate.w);
            //assert(mpz_cmp_ui(candidate.cW,candidate.concurrentWeight) == 0);
            //totalWeight += candidate.weight;
            mpz_add(tW,tW,candidate.w);
            //assert(mpz_cmp_ui(tW,totalWeight) == 0);
        }
    }
}
*/
/*bool comp3(std::vector<Candidate>& v, long long rand){
    return v.back().concurrentWeight < rand;
}
bool comp4(Candidate& c, long long rand){
    return c.concurrentWeight < rand;
}*/
/*std::pair<int, int> CandidateSet::draw(long long int rand) {
    assert(rand <= totalWeight);
    auto it = std::lower_bound(begin(),end(),rand,comp3);
    assert(it < end());
    int idx = it - begin();
    auto it2 = std::lower_bound(this->operator[](idx).begin(), this->operator[](idx).end(),rand,comp4);
    return {idx,it2-this->operator[](idx).begin()};
}*/

void CandidateSet::showCandidate(std::pair<int, int> t) {
        /*
    int baseCurveId = t.first;
    std::vector<FreeSpaceVisualizer> fsvs;
    for(FreeSpace & fs : freespaces[baseCurveId])
        fsvs.emplace_back(fs);
    Candidate & candidate = getCandidate(t);
    for(int i =0; i<fsvs.size();++i) {
        //i is target curve index
        fsvs[i].showCandidate(candidate, i);
    }
         */
}

bool comp5(std::pair<int,int> a, std::pair<int,int> b){
    return a.first < b.first;
}
void CandidateSet::showCovering(std::vector<std::pair<int, int>> indices) {
#ifdef HASVISUAL
    FreeSpacesVisualizer fsv = FreeSpacesVisualizer(freespaces);
    fsv.showCandidates(*this,indices);
#endif

}

CandidateSetPQ::CandidateSetPQ(Curves &c, double d) : Parent(cmpPQ), curves(c), delta(d) {
    for(auto & curve_i : curves) {
        freespaces.emplace_back();
        freespaces.back().reserve(curves.size());
    }
    std::cout << "Computing Freespaces";
    int num_threads = omp_get_max_threads();
#pragma omp parallel for default(none) shared(freespaces,curves,std::cout,num_threads)
    for(int i=0;i<freespaces.size();++i){
        Curve& curve_i = curves[i];
        std::vector<FreeSpace>& col = freespaces[i];
        for(auto & curve_j : curves){
            col.emplace_back(curve_i,curve_j,delta,num_threads);
        }
        //freespaces.push_back(col);
        std::cout << ".";
    }
    std::cout << " Done\n";
}
/*
void CandidateSetPQ::computeCandidates(int l) {
    //l is length limit
    std::chrono::steady_clock::time_point totalbegin = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point ysbegin = std::chrono::steady_clock::now();
    std::cout << "Identifying Important Ys\n";

//#pragma parallel for
    for(auto & col : freespaces){
        for(auto& fs : col){
            //with better identify function
            //fs.ultrafastIdentify(l);
            fs.identifyImportantUpDownYs(l);
        }
    }

    FreeSpacesVisualizer fsv(freespaces);
    fsv.show();


    std::chrono::steady_clock::time_point ysend = std::chrono::steady_clock::now();
    std::cout << "Done. Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(ysend - ysbegin).count() << "[µs]" << std::endl;

    std::vector<char> progressvector;
    for(int i=0;i<curves.size();++i)
        progressvector.push_back(' ');

    //prepare for multithreading

//#define SHOW_SINGLETHREADED
    int totalcount = 0;
#ifndef SHOW_SINGLETHREADED
#pragma omp parallel for schedule(guided,1) default(none) shared(l,std::cout, totalcount,progressvector)
    for(int baseCurveIdx = 0; baseCurveIdx < freespaces.size();baseCurveIdx++){
#else
        int loopID = 0;
    for(int baseCurveIdx = 0; baseCurveIdx < freespaces.size();baseCurveIdx++){
#endif
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

#pragma omp critical
        {
            progressvector[baseCurveIdx] = '.';
            std::cout << "[";
            for (char c : progressvector)
                std::cout << c << "|";
            std::cout << "\b]" << " Local time difference = " << 0 << "[µs]" << std::endl;
        }
        std::vector<Candidate> basecurveCandidates;
        auto col = freespaces[baseCurveIdx];
        for(int t1Idx = 0; t1Idx < col.size();t1Idx++){
            FreeSpace& ff1 = col[t1Idx];
            for(int t2Idx = 0; t2Idx < col.size();t2Idx++){
                FreeSpace& ff2 = col[t2Idx];
                //generate candidates

                int Bsize = ff1.ySize();
                int T1size = ff1.xSize();
                int T2size = ff2.xSize();

                assert(Bsize == ff2.ySize());

#ifdef SHOW_SINGLETHREADED
                FreeSpaceVisualizer fsv1 = FreeSpaceVisualizer(ff1);
                fsv1.show();
                FreeSpaceVisualizer fsv2 = FreeSpaceVisualizer(ff2);
                fsv2.show();
#endif

                //iterate over edges on B
                for(int edgeIdxB = 0; edgeIdxB < Bsize;edgeIdxB++){
                    //identify important Y (and corresponding cell)
                    for(int edgeIdxT1 = 0; edgeIdxT1 < T1size;edgeIdxT1++){

                        ///                             ///
                        /// FIRST OF GO TOWARDS THE TOP ///
                        ///                             ///

                        if(!ff1.getCell(edgeIdxT1, edgeIdxB)->importantDownYs.empty()){

                            //TODO: endpoints doesnt need to be a vector, since we really only ever need the last item

                            std::vector<std::pair<ParamPoint, ParamPoint>> endpoints;

                            //now FOR EVERY START

                            Cell *startCellInFF1 = ff1.getCell(edgeIdxT1, edgeIdxB);

                            //we assume, that importantYs is sorted

                            for (std::pair<CellPoint, FreeSpacePoint> paramPointPair: startCellInFF1->importantDownYs) {

                                std::vector<ParamPoint> localEndpoints;

                                CellPoint paramPoint = paramPointPair.first;
                                ParamPoint ystart = {edgeIdxB, paramPoint.y};

                                //Identify limits (in y-direction) for B
                                //first upwards:
                                FreeSpacePoint importantStart = {{edgeIdxT1, paramPoint.x}, ystart};
                                std::vector<FreeSpacePoint> up_extremes = ff1.findExtremePointsFromPointWrapper(Direction::Up,
                                                                                                                Direction::Right,
                                                                                                                importantStart);

                                ParamPoint limit = std::min(up_extremes[0].y,{ystart.id+l-1,1.0});

                                //now we step through T2, finding the first cell, that has a point between ystart and limit

                                int startOnT2 = 0;

                                bool firstIteration = true;

                                while (startOnT2 < ff2.xSize()) {

                                    //first zoom to first nonempty cell
                                    Cell *cell = ff2.getCell(startOnT2, edgeIdxB);

                                    if(!firstIteration){
                                        //then we have handled the cell, and we need to advance
                                        startOnT2++;
                                        if (startOnT2 >= ff2.xSize())
                                            break;
                                        cell = ff2.getCell(startOnT2, edgeIdxB);
                                        while(cell->left.begin <= importantStart.y.t && cell->left.end >= importantStart.y.t){
                                            startOnT2++;
                                            if (startOnT2 >= ff2.xSize())
                                                break;
                                            cell = ff2.getCell(startOnT2, edgeIdxB);
                                        }
                                    }

                                    if (startOnT2 >= ff2.xSize())
                                        break;

                                    firstIteration = false;

                                    //these breaks can be done way better...
                                    while (cell->is_empty() || !(cell->topPair.second.y >= importantStart.y.t &&
                                                                 cell->bottomPair.first.y <= importantStart.y.t)) {
                                        startOnT2++;
                                        if (startOnT2 >= ff2.xSize())
                                            break;
                                        cell = ff2.getCell(startOnT2, edgeIdxB);
                                    }
                                    if (startOnT2 >= ff2.xSize())
                                        break;

                                    //cell now points to our startcell

#ifdef SHOW_SINGLETHREADED
                                    loopID++;
#endif


                                    CellPoint leftMostInT2Cell = cell->getExtremePointAt<Direction::Left>(importantStart.y.t);
                                    FreeSpacePoint leftMostInCellFSP = {{startOnT2, leftMostInT2Cell.x},
                                                                        {edgeIdxB,  leftMostInT2Cell.y}};

                                    assert(leftMostInCellFSP.y >= importantStart.y && leftMostInCellFSP.y <= limit);


                                    std::vector<FreeSpacePoint> localExtremes = ff2.findExtremePointsFromPointWrapper(
                                            Direction::Right, Direction::Up, leftMostInCellFSP, limit);

                                    ParamPoint maximumOnT2 = localExtremes.back().x;

                                    //if it is too short and on the same edge, we stop

                                    if(ystart.id == localExtremes.back().y.id && std::abs(ystart.t - localExtremes.back().y.t) < EPSILON){
                                        continue;
                                    }

                                    //check if localEndpoints contains this endpoint
                                    bool contains = false;
                                    for(auto lE : localEndpoints)
                                        if(localExtremes.back().y == lE){
                                            contains = true;
                                            break;
                                        }

                                    //TODO: add other ways of duplicate removal

                                    //add maximum to endpoints
                                    if (!contains && (endpoints.empty() || endpoints.back().second < maximumOnT2)) {
                                        localEndpoints.push_back(localExtremes.back().y);
                                        endpoints.emplace_back(localExtremes.back().y, localExtremes.back().x);
                                        basecurveCandidates.emplace_back(ystart,localExtremes.back().y);
#ifdef SHOW_SINGLETHREADED
                                        basecurveCandidates.back().ID = loopID;
#endif
                                    }
                                }
                            }
                        }///end of topscope

                        ///                             ///
                        /// THEN GO TOWARDS THE BOTTOM  ///
                        ///                             ///

                        if(!ff1.getCell(edgeIdxT1, edgeIdxB)->importantUpYs.empty()){

                            //TODO: endpoints doesnt need to be a vector, since we really only ever need the last item

                            std::vector<std::pair<ParamPoint, ParamPoint>> endpoints;

                            //now FOR EVERY START

                            Cell *startCellInFF1 = ff1.getCell(edgeIdxT1, edgeIdxB);

                            //we assume, that importantYs is sorted

                            for (std::pair<CellPoint, FreeSpacePoint> paramPointPair: startCellInFF1->importantUpYs) {

                                std::vector<ParamPoint> localEndpoints;

                                CellPoint paramPoint = paramPointPair.first;
                                ParamPoint ystart = {edgeIdxB, paramPoint.y};

                                //Identify limits (in y-direction) for B
                                //now downwards:
                                FreeSpacePoint importantStart = {{edgeIdxT1, paramPoint.x}, ystart};
                                std::vector<FreeSpacePoint> down_extremes = ff1.findExtremePointsFromPointWrapper(Direction::Down,
                                                                                                                  Direction::Right,
                                                                                                                  importantStart);

                                ParamPoint limit = std::max(down_extremes[0].y,{ystart.id-l+1,0.0});

                                //now we step through T2, finding the first cell, that has a point between ystart and limit

                                int startOnT2 = 0;

                                bool firstIteration = true;

                                while (startOnT2 < ff2.xSize()) {

                                    //first zoom to first nonempty cell
                                    Cell *cell = ff2.getCell(startOnT2, edgeIdxB);

                                    if(!firstIteration){
                                        //then we have handled the cell, and we need to advance
                                        startOnT2++;
                                        if (startOnT2 >= ff2.xSize())
                                            break;
                                        cell = ff2.getCell(startOnT2, edgeIdxB);
                                        while(cell->left.begin <= importantStart.y.t && cell->left.end >= importantStart.y.t){
                                            startOnT2++;
                                            if (startOnT2 >= ff2.xSize())
                                                break;
                                            cell = ff2.getCell(startOnT2, edgeIdxB);
                                        }
                                    }

                                    if (startOnT2 >= ff2.xSize())
                                        break;

                                    firstIteration = false;

                                    //these breaks can be done way better...
                                    while (cell->is_empty() || !(cell->topPair.second.y >= importantStart.y.t &&
                                                                 cell->bottomPair.first.y <= importantStart.y.t)) {
                                        startOnT2++;
                                        if (startOnT2 >= ff2.xSize())
                                            break;
                                        cell = ff2.getCell(startOnT2, edgeIdxB);
                                    }
                                    if (startOnT2 >= ff2.xSize())
                                        break;

                                    //cell now points to our startcell

#ifdef SHOW_SINGLETHREADED
                                    loopID++;
#endif


                                    CellPoint leftMostInT2Cell = cell->getExtremePointAt<Direction::Left>(importantStart.y.t);
                                    FreeSpacePoint leftMostInCellFSP = {{startOnT2, leftMostInT2Cell.x},
                                                                        {edgeIdxB,  leftMostInT2Cell.y}};

                                    assert(leftMostInCellFSP.y <= importantStart.y && leftMostInCellFSP.y >= limit);


                                    std::vector<FreeSpacePoint> localExtremes = ff2.findExtremePointsFromPointWrapper(
                                            Direction::Right, Direction::Down, leftMostInCellFSP, limit);

                                    ParamPoint maximumOnT2 = localExtremes.back().x;

                                    //if it is too short and on the same edge, we stop

                                    if(ystart.id == localExtremes.back().y.id && std::abs(ystart.t - localExtremes.back().y.t) < EPSILON){
                                        continue;
                                    }

                                    //check if localEndpoints contains this endpoint
                                    bool contains = false;
                                    for(auto lE : localEndpoints)
                                        if(localExtremes.back().y == lE){
                                            contains = true;
                                            break;
                                        }

                                    //TODO: add other ways of duplicate removal

                                    //add maximum to endpoints
                                    if (!contains && (endpoints.empty() || endpoints.back().second < maximumOnT2)) {
                                        localEndpoints.push_back(localExtremes.back().y);
                                        endpoints.emplace_back(localExtremes.back().y, localExtremes.back().x);
                                        basecurveCandidates.emplace_back(ystart,localExtremes.back().y);

#ifdef SHOW_SINGLETHREADED
                                        basecurveCandidates.back().ID = loopID;
#endif
                                    }
                                }
                            }
                        }///end of topscope

                    }
                }
            }
        }

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
#pragma omp critical
        {
            progressvector[baseCurveIdx] = 'o';
            std::cout << "[";
            for (char c : progressvector)
                std::cout << c << "|";
            std::cout << "\b]" << " Local time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[µs]" << std::endl;
        }

        begin = std::chrono::steady_clock::now();

        //Now all candidates on B are computed and are in back()
//#pragma omp parallel for default(none) shared(basecurveCandidates,col,baseCurveIdx)
        for(Candidate& candidate : basecurveCandidates){
            for(int tIdx = 0; tIdx < col.size();tIdx++)
                candidate.computeCovering(col[tIdx],tIdx, curves);
#pragma omp critical
        {
            push({baseCurveIdx, candidate});
            //if(baseCurveIdx == 0 && candidate.getStart().id == 0){
            //    std::cout << candidate.getStart().id << "," << candidate.getStart().t << " -> " <<candidate.getEnd().id << "," << candidate.getEnd().t << " @ " << candidate.optimisticCoverLength <<  "\n";
           // }
        }
        }
#pragma omp critical
        end = std::chrono::steady_clock::now();
#pragma omp critical
        {
            progressvector[baseCurveIdx] = '0';
            std::cout << "[";
            for (char c : progressvector)
                std::cout << c << "|";
            std::cout << "\b]" << " Local time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[µs]" << std::endl;
        }
    }
    std::chrono::steady_clock::time_point totalend = std::chrono::steady_clock::now();
    std::cout << "Total amount of candidates: " << totalcount << ". Time difference = "<< std::chrono::duration_cast<std::chrono::milliseconds>(totalend - totalbegin).count() << "[µs]" << std::endl;
}
*/
void CandidateSetPQ::showCovering(std::vector<std::pair<int, Candidate>> candidates) {
#ifdef HASVISUAL
    FreeSpacesVisualizer fsv = FreeSpacesVisualizer(freespaces);
    fsv.showCandidates(std::move(candidates));
#endif
}

bool cmpLength(const std::pair<int,Candidate>& l, const std::pair<int,Candidate>& r){
    return l.second.semiUpdatedCoverLength < r.second.semiUpdatedCoverLength;
}

void CandidateSetPQ::ultrafastComputeCandidates(int l, int minL) {
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
    std::vector<std::pair<int,Candidate>> prepruned;
    //alternative:
    int removes = 0;
    //TODO: fix race condicitions as one thread may reset another threads freespace
#pragma omp parallel for shared(compressedCandidates,freespaces,progress,std::cout, prepruned, minL, removes) default(none)
    //for(auto& cC : compressedCandidates){
    for(int cCi=0;cCi<compressedCandidates.size();++cCi){
        auto cC = compressedCandidates[cCi];
        progress++;
        if(progress % 100==0)
            std::cout << "\33[2K\rUncompressing candidates and augmenting with covering data... " << progress << "/" << compressedCandidates.size() << ", resulting in " << size() << " candidates so far. Removes "<< removes << " short candidates so far... ";
        //prune candidates
        for(int i=0;i<std::get<2>(cC).size();++i){
            if(std::abs(std::get<2>(cC)[i].id - std::get<1>(cC).id) < minL){
                std::get<2>(cC).erase(std::get<2>(cC).begin() + i);
                i--;
#pragma omp critical (block5)
                {
                    removes++;
                }
            }
        }
        if(std::get<2>(cC).empty())
            continue;
        std::vector<Candidate> result = Candidate::uncompressAndComputeCovering(freespaces[std::get<0>(cC)],std::get<1>(cC),std::get<2>(cC),curves, omp_get_thread_num());
        for(auto& c: result){
#pragma omp critical (block6)
            {
                push({std::get<0>(cC), c});
            }
        }
    }
    /*
    std::cout << "Attempting to prune starting with " << prepruned.size() << " many candidates...\n";
    std::sort(prepruned.begin(),prepruned.end(), cmpLength);
    std::reverse(prepruned.begin(),prepruned.end());
    for(int i=0;i<prepruned.size();++i){
        auto& bigger = prepruned[i];
        for(int j=i+1;j<prepruned.size();++j){
            auto& smaller = prepruned[j];
            if(smaller.second.isDominatedBy(bigger.second)){
                prepruned.erase(prepruned.begin() + j);
                j--;
            }
        }
    }*/
    //for(const auto& pp: prepruned)
    //    push(pp);
    /*
#pragma omp parallel for shared(baseCandidates,freespaces,progress,std::cout) default(none)
    for(auto& basecandidate : baseCandidates) {
        progress++;
        if(progress % 1000==0)
            std::cout << "\33[2K\rAugmenting candidates with covering data... " << progress << "/" << baseCandidates.size();
        for(int tIdx = 0; tIdx < curves.size();tIdx++)
            basecandidate.second.computeCovering(freespaces[basecandidate.first][tIdx],tIdx, curves);
        #pragma omp critical
        {
            push(basecandidate);
        }
    }
     */

    std::cout << " Done.\nCandidate generation done.";

}

void CandidateSetPQ::resetWeights() {
    std::vector<std::pair<int,Candidate>> temp;
    while(!empty()) {
        temp.push_back(top());
        pop();
    }
    for(auto& c : temp){
        c.second.resetCoverLength();
        push(c);
    }
}
