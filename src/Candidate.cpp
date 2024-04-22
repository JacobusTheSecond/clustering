//
// Created by Jacobus Conradi on 18.04.23.
//

#include "Candidate.h"
#include "FreespaceVisualizer.h"

CandidateSetPQ::CandidateSetPQ(Curves &c, double d) : Parent(cmpPQ), curves(c), delta(d), sparsefreespaces(c,d,omp_get_max_threads()) {
    //std::cout << "Computing Freespaces";
    //int num_threads = omp_get_max_threads();
    //sparsefreespaces()
    //std::cout << " Done\n";
    //FreeSpacesVisualizer fsv(freespaces);
    //fsv.show();
}

CPoints propagateUpAndIntersect(SparseFreespace& sfs, int y, int x, distance_t startheight, const CPoints& ends, CurveID tIndex, int threadID=0){

    SparseGridCell<std::unique_ptr<Cell>>* startcell = sfs.cell(y,x);

    std::vector<SparseGridCell<std::unique_ptr<Cell>>*> resetList;
    resetList.push_back(startcell);
    std::deque<SparseGridCell<std::unique_ptr<Cell>>*> nextList;

    int leastUpdateableStartIndex = 0;
    CPoints result(ends.size(),{0,-1});

    bool initialIteration = true;

    while(initialIteration || !nextList.empty()){
        Interval fromLeft;
        Interval fromBelow;
        SparseGridCell<std::unique_ptr<Cell>>* cur;
        if(initialIteration){
            fromLeft = {startheight,1.0};
            fromBelow = {1.0,0.0};
            cur = startcell;
        }else{
            cur = nextList.front();
            nextList.pop_front();
            fromLeft = (cur->leftId!= -1)?sfs.cell(cur->y,cur->leftId)->data->toRight(threadID):Interval(1.0,0.0);
            fromBelow = (cur->downId!= -1)?sfs.cell(cur->y-1,cur->downId)->data->toAbove(threadID):Interval(1.0,0.0);
        }
        initialIteration = false;
        assert(!(fromBelow.is_empty() && fromLeft.is_empty()));
        cur->data->toRight(threadID) = (fromBelow.is_empty())?cur->data->right.clip({fromLeft.begin,1.0}):cur->data->right;
        cur->data->toAbove(threadID) = (fromLeft.is_empty())?cur->data->top.clip({fromBelow.begin,1.0}):cur->data->top;

        if(cur->upId != -1 && !cur->data->toAbove(threadID).is_empty() && cur->y < ends.back().getPoint()){
            nextList.push_back(sfs.cell(cur->y+1,cur->upId));
            resetList.push_back(sfs.cell(cur->y+1,cur->upId));
        }
        if(cur->rightId != -1 && !cur->data->toRight(threadID).is_empty() //has something to the right to propagate to
            && (nextList.empty() || nextList.front()!=sfs.cell(cur->y,cur->rightId)))  //and the next entry in the nextList is not already the upcoming cell
        {
            nextList.push_front(sfs.cell(cur->y,cur->rightId));
            resetList.push_back(sfs.cell(cur->y,cur->rightId));
        }

        //compute interval of values, that can be updated
        Interval fromBelowApplicableHeights = fromBelow.is_empty()?Interval(1.0,0.0):
                (cur->data->topPair.second.x>=fromBelow.begin?
                    Interval(0.0,cur->data->topPair.second.y): //either the topmost point lies to the right, and is thus reacahble
                    Interval(0.0,cur->data->topMostAt(fromBelow.begin).y)); //or the topmost point lies to the left, and thus the left boundary realizes the topmost point
        Interval fromLeftApplicableHeights = fromLeft.is_empty()?Interval(1.0,0.0):Interval(fromLeft.begin,cur->data->topPair.first.y);
        //the bounding interval containing both
        Interval boundingvalues = fromBelowApplicableHeights.is_empty()?fromLeftApplicableHeights:(fromLeftApplicableHeights.is_empty()?Interval(1.0,0.0):Interval(0.0,fromLeftApplicableHeights.end));
        CPoint cutoff = {cur->y,boundingvalues.end};
        for(int i=leastUpdateableStartIndex;i<ends.size() && ends[i]<=cutoff;i++){
            if(ends[i].getPoint() < cur->y){
                leastUpdateableStartIndex+=1;
                continue;
            }
            assert(ends[i].getPoint() == cur->y);
            if(fromLeftApplicableHeights.contains(ends[i].getFraction()) || fromBelowApplicableHeights.contains(ends[i].getFraction())){
                result[i] = std::max(result[i],{cur->x,cur->data->rightMostAt(ends[i].getFraction()).x});
            }
        }
    }
    for(auto cell : resetList){
        cell->data->toRight(threadID) = Interval();
        cell->data->toAbove(threadID) = Interval();
    }
    return result;
}


CPoints propagateDownAndIntersect(SparseFreespace& sfs, int y, int x, distance_t startheight, const CPoints& ends, CurveID tIndex, int threadID=0){
    SparseGridCell<std::unique_ptr<Cell>>* startcell = sfs.cell(y,x);
    std::vector<SparseGridCell<std::unique_ptr<Cell>>*> resetList;
    resetList.push_back(startcell);
    std::deque<SparseGridCell<std::unique_ptr<Cell>>*> nextList;

    assert(std::all_of(ends.begin(),ends.end(),[&startcell](CPoint c){return c.getPoint() == startcell->y;}));

    int leastUpdateableStartIndex = ((int)ends.size())-1;
    CPoints result(ends.size(),{0,-1});

    bool initialIteration = true;

    while(initialIteration || !nextList.empty()){
        Interval fromLeft;
        Interval fromBelow;
        SparseGridCell<std::unique_ptr<Cell>>* cur;
        if(initialIteration){
            fromLeft = {0.0,startheight};
            cur = startcell;
        }else{
            cur = nextList.front();
            nextList.pop_front();
            fromLeft = (cur->leftId!= -1)?sfs.cell(cur->y,cur->leftId)->data->toRight(threadID):Interval(1.0,0.0);
        }
        initialIteration = false;
        assert(!(fromBelow.is_empty() && fromLeft.is_empty()));
        cur->data->toRight(threadID) = cur->data->right.clip({0.0,fromLeft.end});

        if(cur->rightId != -1 && !cur->data->toRight(threadID).is_empty())
        {nextList.push_front(sfs.cell(cur->y,cur->rightId)); resetList.push_back(sfs.cell(cur->y,cur->rightId));}

        //compute interval of values, that can be updated
        Interval fromLeftApplicableHeights = fromLeft.is_empty()?Interval(1.0,0.0):Interval(cur->data->bottomPair.first.y,fromLeft.end);
        for(int i=leastUpdateableStartIndex;i>=0 && ends[i].getFraction()>=fromLeftApplicableHeights.begin;i--){
            //assert(ends[i].getPoint() == cur->y);
            if(ends[i].getFraction() > fromLeft.end){
                leastUpdateableStartIndex-=1;
                continue;
            }
            if(fromLeftApplicableHeights.contains(ends[i].getFraction())){
                result[i] = std::max(result[i],{cur->x,cur->data->rightMostAt(ends[i].getFraction()).x});
            }
        }
    }
    for(auto cell : resetList){
        cell->data->toRight(threadID) = Interval();
    }
    return result;
}

std::vector<Candidate> CandidateSetPQ::uncompressCandidate(CurveID bIndex, CPoint start, const CPoints& ends,int threadID) {
    //for every (leftmost) intersection of that start in the freespace
    //      determine the maximal x-coordinate at every end-height simultaneously
    //      append this pair into the correct candidates visual matching
    std::vector<SparseFreespace>& fss = sparsefreespaces[bIndex];
    std::vector<Candidate> result;
    //initialize empty candidates
    for(auto end:ends){
        result.emplace_back(start,end,bIndex);
    }
    for(auto& fs : fss){
        //SparseFreespace& fs = fss[tIndex];
        //traverse to first cell in this free space, where the start intersects
        CurveID tId = fs.TID;
        assert(fs.BID == bIndex);
        for(auto xIdx = 0;xIdx<fs.row(start.getPoint()).size();xIdx++){
            auto cell = fs.cell(start.getPoint(),xIdx);
            if(cell->data->bottomPair.first.y <= start.getFraction() && start.getFraction() <= cell->data->topPair.first.y){
                //cell intersects, so we can traverse
                auto limits = (ends[0]<start)?propagateDownAndIntersect(fs,start.getPoint(),xIdx,start.getFraction(),ends,tId,threadID):propagateUpAndIntersect(fs,start.getPoint(),xIdx,start.getFraction(),ends,tId,threadID);
                for(int i=0;i< limits.size();i++){
                    auto limit = limits[i];
                    if(limit != CPoint(0,-1.0)){
                        result[i].visualMatching.push_back({{cell->x,cell->data->leftMostAt(start.getFraction()).x},limit,tId});
                    }
                }
                //propagate through freespace and keep track
            }
            while(cell->rightId != -1 && cell->data->right.begin <= start.getFraction() && start.getFraction() <= cell->data->right.end){
                xIdx++;
                cell = fs.cell(start.getPoint(),xIdx);
            }
        }
    }
    for(auto& c : result){
        std::sort(c.visualMatching.begin(),c.visualMatching.end(),[](CInterval& a, CInterval& b){return a.getCurveIndex() < b.getCurveIndex() || ((a.getCurveIndex() == b.getCurveIndex()) && (a.begin < b.begin));});
    }
    return result;
}

void CandidateSetPQ::showCovering(std::vector<Candidate> candidates) {
#ifdef HASVISUAL
    SparseFreeSpacesVisualizer fsv(sparsefreespaces);
    fsv.showCandidates(std::move(candidates));
#endif
}

void CandidateSetPQ::showFreespaces() {
#ifdef HASVISUAL
    SparseFreeSpacesVisualizer fsv(sparsefreespaces);
    fsv.show();
#endif
}

SparseFreeSpaces::SparseFreeSpaces(Curves &curves, double delta, int threadcount):Parent(curves.size()) {
    //Step1: figure out boundingboxes of curves
    std::vector<std::vector<int>> extremalIndices(curves.size());

    dimensions_t dims = curves.front().front().dimensions();

    std::vector<Point> directions;

    if(dims == 2){
        for(int i=0;i<12;i++){
            double s = sin(M_PI*i/6);
            double c = cos(M_PI*i/6);
            directions.emplace_back(c,s);
        }
    }else if(dims == 3){
        for(int i=0;i<12;i++){
            double s = sin(M_PI*i/6);
            double c = cos(M_PI*i/6);
            //z plane
            directions.emplace_back(c,s,0);
            if(i%3!=0) {
                //y plane
                directions.emplace_back(c,0,s);
                //x plane
                directions.emplace_back(0,c,s);
            }
        }
        directions.emplace_back(0,0,1);
        directions.emplace_back(0,0,-1);
        double v = 1/sqrt(3);
        for(int i=0;i<7;i++){
            double x = ((i&1)==0?-1:1)*v;
            double y = ((i&2)==0?-1:1)*v;
            double z = ((i&4)==0?-1:1)*v;
            directions.emplace_back(x,y,z);
        }
    }else{
        Point base(dims);
        for(int i=0;i<dims;i++){
            directions.push_back(base);
            directions.back()[i] = 1;
            directions.push_back(base);
            directions.back()[i] = -1;
        }
    }

    for(int i=0;i<curves.size();i++){
        for(auto dir : directions){
            extremalIndices[i].push_back(0);
        }
        for(int j=1;j<curves[i].size();j++){
            for(int k = 0;k<directions.size();k++){
                Point dir = directions[k];
                if(dir*(curves[i][j]) > dir*(curves[i][extremalIndices[i][k]])){
                    extremalIndices[i][k] = j;
                }
            }
        }
    }
    std::vector<std::pair<int,int>> intersectionpairs;
    //Step2: look for intersections
#pragma omp parallel for default(none) shared(curves,directions,extremalIndices,delta,std::cout,threadcount,intersectionpairs)
    for(int i=0;i<curves.size();i++){
        int size = 0;
        for(int j=0;j<=i;j++){
            //check if one of the directions seperate i from j and mirrored
            bool trulyseparated = false;
            for(int k=0;k<directions.size();k++){
                bool separated = true;
                for(int l=0;l<curves[j].size();l++){
                    if(curves[j][l]*directions[k]<curves[i][extremalIndices[i][k]]*directions[k]+delta){
                        separated = false;
                        break;
                    }
                }
                if(separated){
                    trulyseparated = true;
                    break;
                }
            }
            for(int k=0;k<directions.size();k++){
                bool separated = true;
                for(int l=0;l<curves[i].size();l++){
                    if(curves[i][l]*directions[k]<curves[j][extremalIndices[j][k]]*directions[k]+delta){
                        separated = false;
                        break;
                    }
                }
                if(separated){
                    trulyseparated = true;
                    break;
                }
            }

            if(!trulyseparated){
                size++;
#pragma omp critical
                {
                    intersectionpairs.emplace_back(i, j);
                };
            }
            /*
            //check for collision between i and j, mirror to (j,i)
            bool intersecting = false;
            //this check is highly directional and not quite stable, as it simply checks, if any vertex of the one polygon lives inside the other
            for(int k=0;k<directions.size();k++){
                bool intersectionpredicate = true;
                for(int k2=0;k2<directions.size();k2++){
                    //check if v_ik is INSIDE the polygon spanned by v_j{K2} (but not good)
                    if(directions[k2]*curves[i][extremalIndices[i][k]] - delta > directions[k2]*curves[j][extremalIndices[j][k2]]) {
                        intersectionpredicate = false;
                        break;
                    }
                }
                if(intersectionpredicate){
                    intersecting = true;
                    break;
                }
            }
            if(!intersecting) {
                for (int k = 0; k < directions.size(); k++) {
                    bool intersectionpredicate = true;
                    for (int k2 = 0; k2 < directions.size(); k2++) {
                        //check if v_jk is INSIDE the polygon spanned by v_i{K2} (but not good)
                        if (directions[k2] * curves[j][extremalIndices[j][k]] - delta >
                            directions[k2] * curves[i][extremalIndices[i][k2]]) {
                            intersectionpredicate = false;
                            break;
                        }
                    }
                    if (intersectionpredicate) {
                        intersecting = true;
                        break;
                    }
                }
            }

            if(intersecting){
                size++;
#pragma omp critical
                {
                    intersectionpairs.emplace_back(i, j);
                };
            }
             */
        }
    }
    std::vector<int> sizes(curves.size(),0);
    for(auto pair : intersectionpairs){
        int i = pair.first;
        int j = pair.second;
        sizes[i]++;
        if(i!=j){
            sizes[j]++;
        }
    }

    for(int i=0;i<curves.size();i++){
        this->operator[](i).reserve(sizes[i]);
    }
    int num = 0;
#pragma omp parallel for default(none) shared(intersectionpairs,std::cout,curves,delta,threadcount,num)
    for(auto pair:intersectionpairs){
        num ++;
        int i = pair.first;
        int j = pair.second;
        std::cout <<"Populating " << i << " " << j << " {" << num << "/"<< intersectionpairs.size()<<"} (" << curves[i].get_name() << " and " << curves[j].get_name() << ")\n"<<std::flush;
        //SparseFreespace(Curve& B, Curve& T, distance_t _delta, int threadcount = 1, CurveID BID = -1, CurveID TID = -1);

        SparseFreespace sfs(curves[i], curves[j], delta, threadcount, i, j);
#pragma omp critical
        {
            this->operator[](i).emplace_back(std::move(sfs));
        };
        if (i != j) {
            SparseFreespace sfs2(curves[j], curves[i], delta, threadcount, j, i);
#pragma omp critical
            {
                this->operator[](j).emplace_back(std::move(sfs2));
            };
        }
    }
}
