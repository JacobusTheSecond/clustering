//
// Created by Jacobus Conradi on 17.04.23.
//

#include "SparseFreespace.h"

CellPoint SparseCell::rightMostAt(double y, CellPoint* outer) const {
    assert(0<=y && y<=1);
    Interval out;
    Interval in = IntersectionAlgorithm::intersection_interval(lra1,lrb1,y*lrb2,lrc1,y*lrc2,y*y*lrc3,delta, &out);

    assert(!in.is_empty());
    if(outer != nullptr)
        outer->x = out.end, outer->y = y;
    return {in.end,y};
}

CellPoint SparseCell::leftMostAt(double y, CellPoint* outer) const {
    assert(0<=y && y<=1);
    Interval out;
    Interval in = IntersectionAlgorithm::intersection_interval(lra1,lrb1,y*lrb2,lrc1,y*lrc2,y*y*lrc3,delta, &out);

    assert(!in.is_empty());
    if(outer != nullptr)
        outer->x = out.begin, outer->y = y;
    return {in.begin,y};
}

CellPoint SparseCell::topMostAt(double x, CellPoint* outer) const {
    assert(0<=x && x<=1);
    Interval out;
    Interval in = IntersectionAlgorithm::intersection_interval(tba1,tbb1,x*tbb2,tbc1,x*tbc2,x*x*tbc3,delta,&out);
    assert(!in.is_empty());
    if(outer != nullptr)
        outer->x = x, outer->y = out.end;
    return {x,in.end};
}

CellPoint SparseCell::bottomMostAt(double x, CellPoint* outer) const {
    assert(0<=x && x<=1);
    Interval out;
    Interval in = IntersectionAlgorithm::intersection_interval(tba1,tbb1,x*tbb2,tbc1,x*tbc2,x*x*tbc3,delta,&out);
    assert(!in.is_empty());
    if(outer != nullptr)
        outer->x = x, outer->y = out.begin;
    return {x,in.begin};
}

SparseCell::SparseCell(Point &a, Point &b, Point &c, Point &d, distance_t delta, int tc):/*a(a),b(b),c(c),d(d),*/Cell(tc),delta(delta) {
    Interval intersection = IntersectionAlgorithm::intersection_interval(c,d,delta,a,b);
    if(intersection.is_empty()){
        isEmpty = true;
    }else {
        isEmpty = false;

        Point dc = d - c, ca = c - a, ba = b - a, ac = a - c;

        //topbottom values
        tba1 = (dc).length_sqr();
        tbb1 = (ca).dot(dc);
        tbb2 = (ba).dot(dc);
        tbc1 = (ca).length_sqr();
        tbc2 = (ca).dot(ba);
        tbc3 = (ba).length_sqr();

        //leftright values
        lra1 = (ba).length_sqr();
        lrb1 = (ac).dot(ba);
        lrb2 = (dc).dot(ba);
        lrc1 = (ac).length_sqr();
        lrc2 = (ac).dot(dc);
        lrc3 = (dc).length_sqr();

        //boundary update
        left = IntersectionAlgorithm::intersection_interval(tba1, tbb1, 0.0 * tbb2, tbc1, 0.0 * tbc2, 0.0 * 0.0 * tbc3,
                                                            delta);
        right = IntersectionAlgorithm::intersection_interval(tba1, tbb1, 1.0 * tbb2, tbc1, 1.0 * tbc2, 1.0 * 1.0 * tbc3,
                                                             delta);
        bottom = IntersectionAlgorithm::intersection_interval(lra1, lrb1, 0.0 * lrb2, lrc1, 0.0 * lrc2,
                                                              0.0 * 0.0 * lrc3, delta);
        top = IntersectionAlgorithm::intersection_interval(lra1, lrb1, 1.0 * lrb2, lrc1, 1.0 * lrc2, 1.0 * 1.0 * lrc3,
                                                           delta);


        Interval clipper = Interval(0.0, 1.0);
        //intersection.begin is x-coord of leftmost point, intersection.end is rightmost
        double x1 = (left.is_empty() ? intersection.begin : 0.0);
        double x2 = (right.is_empty() ? intersection.end : 1.0);
        Interval intersection2;
        Interval intersection3;
        IntersectionAlgorithm::intersection_interval(tba1, tbb1, x1 * tbb2, tbc1, x1 * tbc2, x1 * x1 * tbc3, delta,
                                                     &intersection2);
        IntersectionAlgorithm::intersection_interval(tba1, tbb1, x2 * tbb2, tbc1, x2 * tbc2, x2 * x2 * tbc3, delta,
                                                     &intersection3);
        assert(!intersection2.is_empty());
        assert(!intersection3.is_empty());
        leftPair = {CellPoint(x1, clipper.clip(intersection2.begin)), CellPoint(x1, clipper.clip(intersection2.end))};
        rightPair = {CellPoint(x2, clipper.clip(intersection3.begin)), CellPoint(x2, clipper.clip(intersection3.end))};


        //compute Topmost
        Interval intersection4 = IntersectionAlgorithm::intersection_interval(a, b, delta, c, d);
        double y1 = (bottom.is_empty() ? intersection4.begin : 0.0);
        double y2 = (top.is_empty() ? intersection4.end : 1.0);
        Interval intersection5, intersection6;
        IntersectionAlgorithm::intersection_interval(lra1, lrb1, y1 * lrb2, lrc1, y1 * lrc2, y1 * y1 * lrc3, delta,
                                                     &intersection5);
        IntersectionAlgorithm::intersection_interval(lra1, lrb1, y2 * lrb2, lrc1, y2 * lrc2, y2 * y2 * lrc3, delta,
                                                     &intersection6);
        //bottomMost = Point(intersection5.begin,intersection4.begin);
        //topMost = Point(intersection6.end,intersection4.end);
        bottomPair = {CellPoint(clipper.clip(intersection5.begin), y1), CellPoint(clipper.clip(intersection5.end), y1)};
        topPair = {CellPoint(clipper.clip(intersection6.begin), y2), CellPoint(clipper.clip(intersection6.end), y2)};
    }
}

bool close_linesegments(const Point& p1, const Point& p2, const Point& q1, const Point& q2, distance_t delta){
    Interval intersection = IntersectionAlgorithm::intersection_interval(p1,p2,delta,q1,q2);
    return !intersection.is_empty();
}

SparseFreespace::SparseFreespace(Curve &B, Curve &T, distance_t delta, int threadcount, CurveID BID, CurveID TID):nx(T.size()-1),ny(B.size()-1),Parent(B.size()-1),BID(BID),TID(TID), delta(delta){
  //  int n = 0;
#pragma omp parallel for default(none) shared(T,B,delta, threadcount,std::cout)
    for(int y=0;y<ny;y++){
        for(int x=0;x<nx;x++){
            //generate x-y-cell
//#pragma omp critical
//                {
//                n += 1;
//            }
                //SparseCell c = SparseCell((Point &) T[x],(Point &) T[x + 1],(Point &) B[y],(Point &) B[y + 1],delta, threadcount);
            MinkowskiCell c = MinkowskiCell((Point &) T[x],(Point &) T[x + 1],(Point &) B[y],(Point &) B[y + 1],T.weight(x)*delta,T.weight(x+1)*delta,B.weight(y)*delta,B.weight(y+1)*delta, threadcount);
                //std::cout << T.weight(x)*delta << T.weight(x+1)*delta << B.weight(y)*delta << B.weight(y+1)*delta<<std::endl;
                //SparseGridCell<std::unique_ptr<Cell>> sgc(std::make_unique<SparseCell>(c),x,y);
            if(!c.is_empty()){
                row(y).emplace_back(std::make_unique<MinkowskiCell>(c), x, y);
            }
        }
    }

    //link in x-direction
#pragma omp parallel for default(none)
    for(int y=0;y<ny;y++){
        for(int xidx=0;xidx+1<row(y).size();xidx++){
            if(cell(y,xidx)->x+1 == cell(y,xidx+1)->x && !cell(y,xidx)->data->right.is_empty()){
                operator[](y).operator[](xidx).rightId = xidx+1;//&(operator[](y).operator[](xidx+1));
                operator[](y).operator[](xidx+1).leftId = xidx;//&(operator[](y).operator[](xidx));
            }
        }
    }

    //link in y-direction
//#pragma omp parallel for default(none)
#pragma omp parallel for default(none)
    for(int y=0;y<ny-1;y++){
        if(row(y+1).empty())
            continue;
        int nextXidx = 0;
        for(int xidx=0;xidx<row(y).size();xidx++){
            while(nextXidx < row(y+1).size() && cell(y+1,nextXidx)->x < cell(y,xidx)->x){
                nextXidx++;
            }
            if(nextXidx >= row(y+1).size())
            {
                break;
            }

            if( cell(y+1,nextXidx)->x == cell(y,xidx)->x && !cell(y,xidx)->data->top.is_empty()){
                operator[](y).operator[](xidx).upId = nextXidx;//&(operator[](y+1).operator[](nextXidx));
                operator[](y+1).operator[](nextXidx).downId = xidx;// &(operator[](y).operator[](xidx));
            }
        }
    }
    //std::cout << n << std::endl;
}

int SparseFreespace::xSize() {
    return nx;
}

int SparseFreespace::ySize() {
    return ny;
}

void SparseFreespace::identifyStarts() {
    upStarts.clear();
    downStarts.clear();
    Interval numericalclip(EPSILON,1.0-EPSILON);
    for(int y=0;y<ny;y++){
        if(row(y).empty()){
            continue;
        }

        bool requireTrivialUpStart = false;

        std::vector<double> potentialUpStarts,potentialDownStarts;
        //inducer, cur, induced_by_rightmost
        std::vector<std::tuple<double,double,bool>> potentialUpPairs;
        std::vector<std::tuple<double,double,bool>> potentialDownPairs;
        for(int xidx = row(y).size()-1;xidx>=0;xidx--){
            auto cellwrapper = cell(y,xidx);

            //NEW VERSION
            //UP
            //if we are at a rightmost cell, append the rightmost point
            if(cellwrapper->rightId == -1){
                potentialUpPairs.emplace_back(cellwrapper->data->rightPair.second.y,cellwrapper->data->rightPair.second.y,true);
            }

            //if the rightmost point in the current cell cant be reached and we can go further left, we might
            //be in a monotonicity case. append this as well
            if(cellwrapper->leftId != -1 && cellwrapper->data->rightPair.second.y < cellwrapper->data->leftPair.first.y){
                potentialUpPairs.emplace_back(cellwrapper->data->leftPair.first.y,cellwrapper->data->leftPair.first.y,false);
            }

            //if we can escape the cell to above, append this case too
            if(cellwrapper->upId != -1){
                potentialUpPairs.emplace_back(1,1,true);
            }

            std::vector<std::tuple<double,double,bool>> tempUpPairs;
            for(auto& start : potentialUpPairs){
                double inducer = std::get<0>(start);
                double cur = std::get<1>(start);
                bool induced_by_rightmost = std::get<2>(start);

                if(!induced_by_rightmost){
                    if(cellwrapper->data->leftPair.second.y < cur){
                        cur = cellwrapper->data->leftPair.second.y;
                        induced_by_rightmost = true;
                    }else if(cellwrapper->data->leftPair.first.y > cur){
                        tempUpPairs.emplace_back(cellwrapper->data->leftPair.first.y,cellwrapper->data->leftPair.first.y,false);
                    }
                }
                if(induced_by_rightmost){

                    //if we can go further left
                    if(cellwrapper->leftId != -1){
                        //push coordinate down if possible, or otherwise push the coordinate
                        if(cellwrapper->data->leftPair.second.y <= cur){
                            tempUpPairs.emplace_back(inducer,cellwrapper->data->leftPair.second.y,induced_by_rightmost);
                        }else if(cellwrapper->data->leftPair.first.y > cur){
                            //push the value unless inducer = cur
                            if(cur < inducer){
                                if(cur > EPSILON){
                                    upStarts.emplace_back(y,cur);
                                }else{
                                    upStarts.emplace_back(y,0);
                                }
                                if(inducer <= 1-EPSILON){
                                    upEnds.emplace_back(y,inducer);
                                }
                            }
                        }else{
                            tempUpPairs.emplace_back(start);
                        }
                    }else{
                        //then we are in a leftmost cell, and we need to push everything that is below the second value
                        if(cellwrapper->data->leftPair.second.y < cur){
                            upStarts.emplace_back(y,cellwrapper->data->leftPair.first.y);
                            if(inducer <= 1-EPSILON){
                                upEnds.emplace_back(y,inducer);
                            }
                        }else if(cellwrapper->data->leftPair.first.y >= cur){
                            //push the value unless inducer = cur
                            if(cur < inducer){
                                if(cur > EPSILON){
                                    upStarts.emplace_back(y,cur);
                                }else{
                                    upStarts.emplace_back(y,0);
                                }
                                if(inducer <= 1-EPSILON){
                                    upEnds.emplace_back(y,inducer);
                                }
                            }
                        }
                    }
                }
            }
            potentialUpPairs = tempUpPairs;


            //DOWN
            //if we are at a rightmost cell, append the rightmost point
            if(cellwrapper->rightId == -1){
                potentialDownPairs.emplace_back(cellwrapper->data->rightPair.first.y,cellwrapper->data->rightPair.first.y,true);
            }

            //if the rightmost point in the current cell cant be reached and we can go further left, we might
            //be in a monotonicity case. append this as well
            if(cellwrapper->leftId != -1 && cellwrapper->data->rightPair.first.y > cellwrapper->data->leftPair.second.y){
                potentialDownPairs.emplace_back(cellwrapper->data->leftPair.second.y,cellwrapper->data->leftPair.second.y,false);
            }

            std::vector<std::tuple<double,double,bool>> tempDownPairs;
            for(auto& start : potentialDownPairs){
                double inducer = std::get<0>(start);
                double cur = std::get<1>(start);
                bool induced_by_rightmost = std::get<2>(start);

                if(!induced_by_rightmost){
                    if(cellwrapper->data->leftPair.first.y > cur){
                        cur = cellwrapper->data->leftPair.first.y;
                        induced_by_rightmost = true;
                    }else if(cellwrapper->data->leftPair.second.y <= cur){
                        tempDownPairs.emplace_back(cellwrapper->data->leftPair.second.y,cellwrapper->data->leftPair.second.y,false);
                    }
                }
                if(induced_by_rightmost){

                    //if we can go further left
                    if(cellwrapper->leftId != -1){
                        //push coordinate down if possible, or otherwise push the coordinate
                        if(cellwrapper->data->leftPair.first.y >= cur){
                            tempDownPairs.emplace_back(inducer,cellwrapper->data->leftPair.first.y,induced_by_rightmost);
                        }else if(cellwrapper->data->leftPair.second.y < cur){
                            //push the value unless inducer = cur
                            if(cur > inducer){
                                if(cur < 1-EPSILON){
                                    downStarts.emplace_back(y,cur);
                                }else{
                                    downStarts.emplace_back(y,1);
                                }
                                if(inducer >= EPSILON){
                                    downEnds.emplace_back(y,inducer);
                                }
                            }
                        }else{
                            tempDownPairs.emplace_back(start);
                        }
                    }else{
                        //then we are in a leftmost cell, and we need to push everything that is below the second value
                        if(cellwrapper->data->leftPair.first.y > cur){
                            downStarts.emplace_back(y,cellwrapper->data->leftPair.second.y);
                            if(inducer >= EPSILON){
                                downEnds.emplace_back(y,inducer);
                            }
                        }else if(cellwrapper->data->leftPair.second.y <= cur){
                            //push the value unless inducer = cur
                            if(cur > inducer){
                                if(cur < 1-EPSILON){
                                    downStarts.emplace_back(y,cur);
                                }else{
                                    downStarts.emplace_back(y,1);
                                }
                                if(inducer >= EPSILON){
                                    downEnds.emplace_back(y,inducer);
                                }
                            }
                        }
                    }
                }
            }
            potentialDownPairs = tempDownPairs;

/*
            if(cellwrapper->left == nullptr){

                //push leftmost point if it isnt near the boundary
                if(numericalclip.contains(cellwrapper->data.leftPair.first.y)){
                    upStarts.emplace_back(y,cellwrapper->data.leftPair.first.y);
                }else if((cellwrapper->data.leftPair.first.y < EPSILON) && (not (cellwrapper->down != nullptr && cellwrapper->down->data.leftPair.second.y < 1.0-EPSILON))){
                    //if the leftmost point is below EPSILON, we check, if we can go further down. If this is the case, then we dont need to handle this trivial start
                    requireTrivialUpStart = true;
                }

                //push all potential up starts, that lie BELOW the leftmost point
                for(auto up : potentialUpStarts){
                    if(numericalclip.contains(up) && up < cellwrapper->data.leftPair.first.y)
                        upStarts.emplace_back(y,up);
                }

                //handle downstarts similarly
                if(numericalclip.contains(cellwrapper->data.leftPair.second.y)){
                    downStarts.emplace_back(y,cellwrapper->data.leftPair.second.y);
                }

                for(auto down : potentialDownStarts){
                    if(numericalclip.contains(down) && down > cellwrapper->data.leftPair.second.y)
                        downStarts.emplace_back(y,down);
                }
                potentialUpStarts.clear();
                potentialDownStarts.clear();
            }else{
                std::vector<double> temporaryUpStarts,temporaryDownStarts;
                temporaryUpStarts.push_back(cellwrapper->data.left.end);
                double upcutoff = cellwrapper->data.left.end;
                for(auto up:potentialUpStarts){
                    if(up < upcutoff){
                        temporaryUpStarts.push_back(up);
                    }
                }
                potentialUpStarts = temporaryUpStarts;


                temporaryDownStarts.push_back(cellwrapper->data.left.begin);
                double downcutoff = cellwrapper->data.left.begin;
                for(auto down:potentialDownStarts){
                    if(down > downcutoff){
                        temporaryDownStarts.push_back(down);
                    }
                }
                potentialDownStarts = temporaryDownStarts;

            }
*/
        }

        //push trivial starts
        //if(requireTrivialUpStart)
        //    upStarts.emplace_back(y,0);
        //downStarts.emplace_back(y,1);
    }
    std::sort(upStarts.begin(),upStarts.end());
    upStarts.erase(std::unique(upStarts.begin(),upStarts.end()),upStarts.end());

    std::sort(upEnds.begin(),upEnds.end());
    upEnds.erase(std::unique(upEnds.begin(),upEnds.end()),upEnds.end());

    std::sort(downStarts.begin(),downStarts.end());
    downStarts.erase(std::unique(downStarts.begin(),downStarts.end()),downStarts.end());

    std::sort(downEnds.begin(),downEnds.end());
    downEnds.erase(std::unique(downEnds.begin(),downEnds.end()),downEnds.end());
}

void SparseFreespace::identifyEnds() {
    upEnds.clear();
    downEnds.clear();
    Interval numericalclip(EPSILON,1.0-EPSILON);
    for (int y = 0; y < ny; y++) {
        if (row(y).empty()) {
            continue;
        }

        bool requireTrivialUpEnd = false;

        std::vector<double> potentialUpEnds, potentialDownEnds;
        for (int xidx = 0; xidx < row(y).size(); xidx++) {
            auto cellwrapper = cell(y, xidx);

            if (cellwrapper->rightId == -1) {

                //if we are in a trivial end-situation, dont push the endcoordinate yet, otherwise push it
                if (numericalclip.contains(cellwrapper->data->rightPair.second.y)) {
                    upEnds.emplace_back(y, cellwrapper->data->rightPair.second.y);
                }else if((cellwrapper->data->rightPair.second.y > 1.0 - EPSILON) && (not (cellwrapper->upId != -1 && cell(y+1,cellwrapper->upId)->data->rightPair.first.y < EPSILON))){
                    requireTrivialUpEnd = true;
                }

                //push all potential up end, that lie above the rightmost point
                for (auto up: potentialUpEnds) {
                    if (numericalclip.contains(up)&& up > cellwrapper->data->rightPair.second.y)
                        upEnds.emplace_back(y, up);
                }

                //handle downstarts similarly
                if (numericalclip.contains(cellwrapper->data->rightPair.first.y)) {
                    downEnds.emplace_back(y, cellwrapper->data->rightPair.first.y);
                }

                for (auto down: potentialDownEnds) {
                    if (numericalclip.contains(down) && down < cellwrapper->data->rightPair.first.y)
                        downEnds.emplace_back(y, down);
                }
                potentialUpEnds.clear();
                potentialDownEnds.clear();
            } else {
                std::vector<double> temporaryUpEnds, temporaryDownEnds;
                temporaryUpEnds.push_back(cellwrapper->data->right.begin);
                double upcutoff = cellwrapper->data->right.begin;
                for (auto up: potentialUpEnds) {
                    if (up > upcutoff) {
                        temporaryUpEnds.push_back(up);
                    }
                }
                potentialUpEnds = temporaryUpEnds;


                temporaryDownEnds.push_back(cellwrapper->data->right.end);
                double downcutoff = cellwrapper->data->right.end;
                for (auto down: potentialDownEnds) {
                    if (down < downcutoff) {
                        temporaryDownEnds.push_back(down);
                    }
                }
                potentialDownEnds = temporaryDownEnds;

            }

        }

        if(requireTrivialUpEnd)
            upEnds.emplace_back(y, 1);
        downEnds.emplace_back(y, 0);
    }
    std::sort(upEnds.begin(),upEnds.end());
    std::sort(downEnds.begin(),downEnds.end());
}

MinkowskiCell::MinkowskiCell(Point &a, Point &b, Point &c, Point &d, distance_t ra, distance_t rb, distance_t rc,
                             distance_t rd, int tc): Cell(tc),ra(ra),rb(rb),rc(rc),rd(rd){
    //acac,acab,accd,abab,abcd;
    //caca,cacd,caab,cdcd,cdab;
    Interval intersection = MinkowskiIntersectionAlgorithm::edgeEdgeIntersection(c,d,rc,rd,a,b,ra,rb,nullptr);
    if(intersection.is_empty()){
        isEmpty = true;
    }else{
        isEmpty = false;
        Point ac = c-a;
        Point ab = b-a;
        Point cd = d-c;
        Point ca = a-c;

        acac = ac.dot(ac);
        acab = ac.dot(ab);
        accd = ac.dot(cd);
        abab = ab.dot(ab);
        abcd = ab.dot(cd);
        caca = ca.dot(ca);
        cacd = ca.dot(cd);
        caab = ca.dot(ab);
        cdcd = cd.dot(cd);
        cdab = cd.dot(ab);

        left = MinkowskiIntersectionAlgorithm::pointEdgeIntersection(a,ra,c,d,rc,rd,nullptr);
        right = MinkowskiIntersectionAlgorithm::pointEdgeIntersection(b,rb,c,d,rc,rd,nullptr);
        bottom = MinkowskiIntersectionAlgorithm::pointEdgeIntersection(c,rc,a,b,ra,rb,nullptr);
        top = MinkowskiIntersectionAlgorithm::pointEdgeIntersection(d,rd,a,b,ra,rb,nullptr);

        Interval clipper = Interval(0.0, 1.0);
        //intersection.begin is x-coord of leftmost point, intersection.end is rightmost
        double x1 = (left.is_empty() ? intersection.begin : 0.0);
        double x2 = (right.is_empty() ? intersection.end : 1.0);
        Interval intersection2;
        Interval intersection3;
        MinkowskiIntersectionAlgorithm::pointEdgeIntersection(a+((b-a)*x1),ra+((rb-ra)*x1),c,d,rc,rd,&intersection2);
        //IntersectionAlgorithm::intersection_interval(tba1, tbb1, x1 * tbb2, tbc1, x1 * tbc2, x1 * x1 * tbc3, delta,
        //                                             &intersection2);
        MinkowskiIntersectionAlgorithm::pointEdgeIntersection(a+((b-a)*x2),ra+((rb-ra)*x2),c,d,rc,rd,&intersection3);
        //IntersectionAlgorithm::intersection_interval(tba1, tbb1, x2 * tbb2, tbc1, x2 * tbc2, x2 * x2 * tbc3, delta,
        //                                             &intersection3);
        assert(!intersection2.is_empty());
        assert(!intersection3.is_empty());
        leftPair = {CellPoint(x1, clipper.clip(intersection2.begin)), CellPoint(x1, clipper.clip(intersection2.end))};
        rightPair = {CellPoint(x2, clipper.clip(intersection3.begin)), CellPoint(x2, clipper.clip(intersection3.end))};


        //compute Topmost
        Interval intersection4 = MinkowskiIntersectionAlgorithm::edgeEdgeIntersection(a,b,ra,rb,c,d,rc,rd,nullptr);
        double y1 = (bottom.is_empty() ? intersection4.begin : 0.0);
        double y2 = (top.is_empty() ? intersection4.end : 1.0);
        Interval intersection5, intersection6;
        MinkowskiIntersectionAlgorithm::pointEdgeIntersection(c+((d-c)*y1),rc+((rd-rc)*y1),a,b,ra,rb,&intersection5);
        //IntersectionAlgorithm::intersection_interval(lra1, lrb1, y1 * lrb2, lrc1, y1 * lrc2, y1 * y1 * lrc3, delta,
        //                                             &intersection5);
        MinkowskiIntersectionAlgorithm::pointEdgeIntersection(c+((d-c)*y2),rc+((rd-rc)*y2),a,b,ra,rb,&intersection6);
        //IntersectionAlgorithm::intersection_interval(lra1, lrb1, y2 * lrb2, lrc1, y2 * lrc2, y2 * y2 * lrc3, delta,
        //                                             &intersection6);

        //bottomMost = Point(intersection5.begin,intersection4.begin);
        //topMost = Point(intersection6.end,intersection4.end);
        bottomPair = {CellPoint(clipper.clip(intersection5.begin), y1), CellPoint(clipper.clip(intersection5.end), y1)};
        topPair = {CellPoint(clipper.clip(intersection6.begin), y2), CellPoint(clipper.clip(intersection6.end), y2)};
    }

}

CellPoint MinkowskiCell::leftMostAt(double y, CellPoint *outer) const {
    assert(0<=y && y<=1);
    Interval out;
    //Interval in = IntersectionAlgorithm::intersection_interval(lra1,lrb1,y*lrb2,lrc1,y*lrc2,y*y*lrc3,delta, &out);
    distance_t papa = caca -(2*y*cacd) + (y*y*cdcd);
    distance_t paab = caab - (y*cdab);
    Interval in = MinkowskiIntersectionAlgorithm::pointEdgeIntersectionPrimitive(papa,paab,abab,rc + ((rd-rc)*y),ra,rb,&out);

    assert(!in.is_empty());
    if(outer != nullptr)
        outer->x = out.begin, outer->y = y;
    return {in.begin,y};
}

CellPoint MinkowskiCell::rightMostAt(double y, CellPoint *outer) const {
    assert(0<=y && y<=1);
    Interval out;
    //Interval in = IntersectionAlgorithm::intersection_interval(lra1,lrb1,y*lrb2,lrc1,y*lrc2,y*y*lrc3,delta, &out);
    distance_t papa = caca -(2*y*cacd) + (y*y*cdcd);
    distance_t paab = caab - (y*cdab);
    Interval in = MinkowskiIntersectionAlgorithm::pointEdgeIntersectionPrimitive(papa,paab,abab,rc + ((rd-rc)*y),ra,rb,&out);

    assert(!in.is_empty());
    if(outer != nullptr)
        outer->x = out.end, outer->y = y;
    return {in.end,y};
}

CellPoint MinkowskiCell::topMostAt(double x, CellPoint *outer) const {
    assert(0<=x && x<=1);
    Interval out;
    //Interval in = IntersectionAlgorithm::intersection_interval(tba1,tbb1,x*tbb2,tbc1,x*tbc2,x*x*tbc3,delta,&out);
    distance_t pcpc = acac - (2*x*acab) + (x*x*abab);
    distance_t pccd = accd - (x*abcd);
    Interval in = MinkowskiIntersectionAlgorithm::pointEdgeIntersectionPrimitive(pcpc,pccd,cdcd,ra+((rb-ra)*x),rc,rd,&out);

    assert(!in.is_empty());
    if(outer != nullptr)
        outer->x = x, outer->y = out.end;
    return {x,in.end};
}

CellPoint MinkowskiCell::bottomMostAt(double x, CellPoint *outer) const {
    assert(0<=x && x<=1);
    Interval out;
    //Interval in = IntersectionAlgorithm::intersection_interval(tba1,tbb1,x*tbb2,tbc1,x*tbc2,x*x*tbc3,delta,&out);
    distance_t pcpc = acac - (2*x*acab) + (x*x*abab);
    distance_t pccd = accd - (x*abcd);
    Interval in = MinkowskiIntersectionAlgorithm::pointEdgeIntersectionPrimitive(pcpc,pccd,cdcd,ra+((rb-ra)*x),rc,rd,&out);

    assert(!in.is_empty());
    if(outer != nullptr)
        outer->x = x, outer->y = out.begin;
    return {x,in.begin};
}

/*    Interval intersection = IntersectionAlgorithm::intersection_interval(c,d,delta,a,b);
    if(intersection.is_empty()){
        isEmpty = true;
    }else {
        isEmpty = false;

        Point dc = d - c, ca = c - a, ba = b - a, ac = a - c;

        //topbottom values
        tba1 = (dc).length_sqr();
        tbb1 = (ca).dot(dc);
        tbb2 = (ba).dot(dc);
        tbc1 = (ca).length_sqr();
        tbc2 = (ca).dot(ba);
        tbc3 = (ba).length_sqr();

        //leftright values
        lra1 = (ba).length_sqr();
        lrb1 = (ac).dot(ba);
        lrb2 = (dc).dot(ba);
        lrc1 = (ac).length_sqr();
        lrc2 = (ac).dot(dc);
        lrc3 = (dc).length_sqr();

        //boundary update
        left = IntersectionAlgorithm::intersection_interval(tba1, tbb1, 0.0 * tbb2, tbc1, 0.0 * tbc2, 0.0 * 0.0 * tbc3,
                                                            delta);
        right = IntersectionAlgorithm::intersection_interval(tba1, tbb1, 1.0 * tbb2, tbc1, 1.0 * tbc2, 1.0 * 1.0 * tbc3,
                                                             delta);
        bottom = IntersectionAlgorithm::intersection_interval(lra1, lrb1, 0.0 * lrb2, lrc1, 0.0 * lrc2,
                                                              0.0 * 0.0 * lrc3, delta);
        top = IntersectionAlgorithm::intersection_interval(lra1, lrb1, 1.0 * lrb2, lrc1, 1.0 * lrc2, 1.0 * 1.0 * lrc3,
                                                           delta);


        Interval clipper = Interval(0.0, 1.0);
        //intersection.begin is x-coord of leftmost point, intersection.end is rightmost
        double x1 = (left.is_empty() ? intersection.begin : 0.0);
        double x2 = (right.is_empty() ? intersection.end : 1.0);
        Interval intersection2;
        Interval intersection3;
        IntersectionAlgorithm::intersection_interval(tba1, tbb1, x1 * tbb2, tbc1, x1 * tbc2, x1 * x1 * tbc3, delta,
                                                     &intersection2);
        IntersectionAlgorithm::intersection_interval(tba1, tbb1, x2 * tbb2, tbc1, x2 * tbc2, x2 * x2 * tbc3, delta,
                                                     &intersection3);
        assert(!intersection2.is_empty());
        assert(!intersection3.is_empty());
        leftPair = {CellPoint(x1, clipper.clip(intersection2.begin)), CellPoint(x1, clipper.clip(intersection2.end))};
        rightPair = {CellPoint(x2, clipper.clip(intersection3.begin)), CellPoint(x2, clipper.clip(intersection3.end))};


        //compute Topmost
        Interval intersection4 = IntersectionAlgorithm::intersection_interval(a, b, delta, c, d);
        double y1 = (bottom.is_empty() ? intersection4.begin : 0.0);
        double y2 = (top.is_empty() ? intersection4.end : 1.0);
        Interval intersection5, intersection6;
        IntersectionAlgorithm::intersection_interval(lra1, lrb1, y1 * lrb2, lrc1, y1 * lrc2, y1 * y1 * lrc3, delta,
                                                     &intersection5);
        IntersectionAlgorithm::intersection_interval(lra1, lrb1, y2 * lrb2, lrc1, y2 * lrc2, y2 * y2 * lrc3, delta,
                                                     &intersection6);
        //bottomMost = Point(intersection5.begin,intersection4.begin);
        //topMost = Point(intersection6.end,intersection4.end);
        bottomPair = {CellPoint(clipper.clip(intersection5.begin), y1), CellPoint(clipper.clip(intersection5.end), y1)};
        topPair = {CellPoint(clipper.clip(intersection6.begin), y2), CellPoint(clipper.clip(intersection6.end), y2)};
    }*/