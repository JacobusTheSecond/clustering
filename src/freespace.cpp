//
// Created by Jacobus Conradi on 26.01.22.
//

#include "freespace.h"

// Cell:
// d +-----+
//   |     |
// c +-----+
//   a     b

// FreeSpaceDiagram
//
//    +   +---+-+---+---+
//    |   |   | |   |   |
//    +   +---+-+---+---+
// B  |   |   | |   |   |
//    |   |   | |   |   |
//    +   +---+-+---+---+
//    +   +---+-+---+---+
//
//        +---+-+---+---+
//               T

//TODO: use low dimensional primitives
FreeSpace::FreeSpace(const Curve& B, const Curve& T, double delta, int threadcount) : reseterList(threadcount) {
    this->delta = delta;
    nx = T.size()-1;
    ny = B.size()-1;
    for(auto& list : reseterList)
        list.reserve(10000);

    Interval clipper(0,1);
//#pragma omp parallel for default(none) shared(B,T,delta,nx,ny, clipper)
    for(int i=0;i<nx;++i){
        if((i+1)%100==0){
            std::cout << ".";
        }
        std::vector<Cell> row;
        for(int j=0;j<ny;++j){
            //i,j is current cell
            Cell cell((Point &) T[i],(Point &) T[i + 1],(Point &) B[j],(Point &) B[j + 1], threadcount);

            Point dc=cell.d-cell.c,ca=cell.c-cell.a,ba=cell.b-cell.a,ac=cell.a-cell.c;

            //topbottom values
            cell.tba1 = (dc).length_sqr();
            cell.tbb1 = (ca).dot(dc);
            cell.tbb2 = (ba).dot(dc);
            cell.tbc1 = (ca).length_sqr();
            cell.tbc2 = (ca).dot(ba);
            cell.tbc3 = (ba).length_sqr();

            //leftright values
            cell.lra1 = (ba).length_sqr();
            cell.lrb1 = (ac).dot(ba);
            cell.lrb2 = (dc).dot(ba);
            cell.lrc1 = (ac).length_sqr();
            cell.lrc2 = (ac).dot(dc);
            cell.lrc3 = (dc).length_sqr();

            cell.delta = delta;
            IntersectionAlgorithm::intersection_interval(cell.a,delta,cell.c,cell.d,&cell.left);
            IntersectionAlgorithm::intersection_interval(cell.b,delta,cell.c,cell.d,&cell.right);
            IntersectionAlgorithm::intersection_interval(cell.c,delta,cell.a,cell.b,&cell.bottom);
            IntersectionAlgorithm::intersection_interval(cell.d,delta,cell.a,cell.b,&cell.top);
            cell.left = clipper.clip(cell.left);
            cell.right = clipper.clip(cell.right);
            cell.top = clipper.clip(cell.top);
            cell.bottom = clipper.clip(cell.bottom);
            //cell.leftMarked = false,cell.rightMarked = false;

            Interval intersection = IntersectionAlgorithm::intersection_interval(cell.c,cell.d,delta,cell.a,cell.b);
            if(intersection.is_empty()){
                cell.isEmpty = true;
            }else {
                cell.isEmpty = false;
                //intersection.begin is x-coord of leftmost point, intersection.end is rightmost
                Point p = cell.a + (cell.b - cell.a) * intersection.begin;
                Point q = cell.a + (cell.b - cell.a) * intersection.end;
                Interval intersection2;
                IntersectionAlgorithm::intersection_interval(p, delta, cell.c, cell.d, &intersection2);
                Interval intersection3;
                IntersectionAlgorithm::intersection_interval(q, delta, cell.c, cell.d, &intersection3);
                assert(!intersection2.is_empty());
                assert(!intersection3.is_empty());
                //cell.leftMost = Point(clipper.clip(intersection.begin), clipper.clip(intersection2.begin));
                //cell.rightMost = Point(clipper.clip(intersection.end), clipper.clip(intersection3.end));
                cell.leftPair = {CellPoint(clipper.clip(intersection.begin), clipper.clip(intersection2.begin)),CellPoint(clipper.clip(intersection.begin), clipper.clip(intersection2.end))};
                cell.rightPair = {CellPoint(clipper.clip(intersection.end), clipper.clip(intersection3.begin)),CellPoint(clipper.clip(intersection.end), clipper.clip(intersection3.end))};


                //compute Topmost
                Interval intersection4 = IntersectionAlgorithm::intersection_interval(cell.a,cell.b,delta,cell.c,cell.d);
                Interval intersection5,intersection6;
                Point u = cell.c + (cell.d-cell.c)*intersection4.begin;
                Point v = cell.c + (cell.d-cell.c)*intersection4.end;
                IntersectionAlgorithm::intersection_interval(u,delta,cell.a,cell.b,&intersection5);
                IntersectionAlgorithm::intersection_interval(v,delta,cell.a,cell.b,&intersection6);
                //cell.bottomMost = Point(intersection5.begin,intersection4.begin);
                //cell.topMost = Point(intersection6.end,intersection4.end);
                cell.bottomPair = {CellPoint(clipper.clip(intersection5.begin),clipper.clip(intersection4.begin)),CellPoint(clipper.clip(intersection5.end),clipper.clip(intersection4.begin))};
                cell.topPair = {CellPoint(clipper.clip(intersection6.begin),clipper.clip(intersection4.end)),CellPoint(clipper.clip(intersection6.end),clipper.clip(intersection4.end))};
            }
            row.push_back(cell);
        }
        data.push_back(row);
    }
}

/*void FreeSpace::prettyPrint(int i, int j) {
    std::cout << "-----Cell"<< i << "," << j <<"-----\n";
    std::cout<<"left: " << data[i][j].left << " ";
    std::cout<<"top: " << data[i][j].top << " ";
    std::cout<<"right: " << data[i][j].right << " ";
    std::cout<<"bottom: " << data[i][j].bottom << "\n";
    std::cout << "----------------------\n";
}*/

int FreeSpace::xSize() {
    return nx;
}
int FreeSpace::ySize() {
    return ny;
}

Cell* FreeSpace::getCell(int i, int j) {
    assert(0<=i && i < nx && 0<=j && j<ny);
    return &data[i][j];
}

void FreeSpace::reset(int threadID) {
    for(int x=0;x<nx;x++){
        for(int y=0;y<ny;y++){
            Cell* cell = getCell(x,y);
            //if(alsoMarked)
            //    cell->leftMarked = false,cell->rightMarked=false;
            assert(cell->toAbove(threadID).is_empty());
            assert(cell->toRight(threadID).is_empty());
            assert(cell->toBottom(threadID).is_empty());
            assert(cell->toLeft(threadID).is_empty());
        }
    }
}
/*
void FreeSpace::identifyImportantYs() {


    for(int x=nx-1;x>=0;x--){
        for(int y=ny-1;y>=0;y--) {

            Cell* cell = getCell(x,y);
            if(cell->isEmpty)
                continue;
            FreeSpacePoint start = {{x,cell->rightMost.x()},{y,cell->rightMost.y()}};
            FreeSpacePoint leftMost = findExtremPointInDirectionFromPoint<Direction::Left>(start);

            //check if endcell already contains this Y
            bool contains = false;
            for(std::pair<Point, FreeSpacePoint> startP: getCell(leftMost.x.id, leftMost.y.id)->importantYs)
                if(startP.first.y() == leftMost.y.t) {
                    contains = true;
                    break;
                }
            if(!contains) {
                Point p = {leftMost.x.t, leftMost.y.t};
                std::pair<Point,FreeSpacePoint> pair = {p,start};
                getCell(leftMost.x.id, leftMost.y.id)->importantYs.push_back(pair);
            }
        }
    }
}
*/
/*
void FreeSpace::identifyImportantUpDownYs() {
    int veryshortcount = 0;
    Interval clipper = Interval(0,1);
    for(int x=nx-1;x>=0;x--){
        for(int y=ny-1;y>=0;y--) {

            Cell* cell = getCell(x,y);
            if(cell->isEmpty)
                continue;


            //TODO: rename properly, this is nolonger topright
            CellPoint lefttop = cell->getExtreme<Direction::Left,Direction::Up>();
            CellPoint leftbottom = cell->getExtreme<Direction::Left,Direction::Down>();
            FreeSpacePoint topstart = {{x,lefttop.x},{y,lefttop.y}};
            FreeSpacePoint bottomstart = {{x,leftbottom.x},{y,leftbottom.y}};
            std::vector<FreeSpacePoint> leftDownMost = findExtremePointsFromPoint<Direction::Left,Direction::Down>(topstart);
            std::vector<FreeSpacePoint> leftUpMost = findExtremePointsFromPoint<Direction::Left,Direction::Up>(bottomstart);

            //check if endcell already contains this Y
            for(FreeSpacePoint leftMost : leftDownMost) {
                if(leftMost.y.id == topstart.y.id && std::abs( leftMost.y.t - topstart.y.t) < EPSILON){
                    veryshortcount++;
                    continue;
                }
                if(leftMost.y.id == topstart.y.id && clipper.clip(leftMost.y.t) == clipper.clip(topstart.y.t)) {
                    veryshortcount++;
                    continue;
                }
                bool contains = false;
                for (std::pair<CellPoint, FreeSpacePoint> startP: getCell(leftMost.x.id, leftMost.y.id)->importantDownYs)
                    if (clipper.clip(startP.first.y) == clipper.clip(leftMost.y.t)) {
                        contains = true;
                        break;
                    }
                if (!contains) {
                    CellPoint p = {leftMost.x.t, leftMost.y.t};
                    std::pair<CellPoint, FreeSpacePoint> pair = {p, topstart};
                    getCell(leftMost.x.id, leftMost.y.id)->importantDownYs.push_back(pair);
                }
            }
            for(FreeSpacePoint leftMost : leftUpMost) {
                if(leftMost.y.id == bottomstart.y.id && std::abs(leftMost.y.t- bottomstart.y.t) < EPSILON){
                    veryshortcount++;
                    continue;
                }
                if(leftMost.y.id == bottomstart.y.id && clipper.clip(leftMost.y.t) == clipper.clip(bottomstart.y.t)) {
                    veryshortcount++;
                    continue;
                }
                bool contains = false;
                for (std::pair<CellPoint, FreeSpacePoint> startP: getCell(leftMost.x.id, leftMost.y.id)->importantUpYs)
                    if (clipper.clip(startP.first.y) == clipper.clip(leftMost.y.t)) {
                        contains = true;
                        break;
                    }
                if (!contains) {
                    CellPoint p = {leftMost.x.t, leftMost.y.t};
                    std::pair<CellPoint, FreeSpacePoint> pair = {p, bottomstart};
                    getCell(leftMost.x.id, leftMost.y.id)->importantUpYs.push_back(pair);
                }
            }

            if(!cell->getBoundary<Direction::Left>().is_empty() && cell->getBoundary<Direction::Left>().end != 1.0){
                //try to walk left as far as possible
                double t = cell->getBoundary<Direction::Left>().end;
                int tx = x;
                Cell* curcell = getCell(tx,y);
                while(tx >= 0 && curcell->left.begin <= t && t <= curcell->left.end) {
                    curcell = getCell(tx,y);
                    tx -= 1;
                }
                //offset tx
                tx++;
                if(curcell->leftPair.first.y > t) {
                    CellPoint p = curcell->getExtremePointAt<Direction::Left>(t);
                    FreeSpacePoint upperstart = {{x,0},{y,t}};
                    std::pair<CellPoint, FreeSpacePoint> pair = {p, upperstart};
                    getCell(tx, y)->importantUpYs.push_back(pair);
                }
            }
        }
    }
    std::cout << "Ignoring " << veryshortcount << " same-cell-candidates of length <= " << EPSILON << "\n";
}

void FreeSpace::identifyImportantUpDownYs(int x, int y) {
    int veryshortcount = 0;
    Interval clipper = Interval(0,1);
            Cell* cell = getCell(x,y);
            if(cell->isEmpty)
                return;


            //TODO: rename properly, this is nolonger topright
            CellPoint lefttop = cell->getExtreme<Direction::Left,Direction::Up>();
            CellPoint leftbottom = cell->getExtreme<Direction::Left,Direction::Down>();
            FreeSpacePoint topstart = {{x,lefttop.x},{y,lefttop.y}};
            FreeSpacePoint bottomstart = {{x,leftbottom.x},{y,leftbottom.y}};
            std::vector<FreeSpacePoint> leftDownMost = findExtremePointsFromPoint<Direction::Left,Direction::Down>(topstart);
            std::vector<FreeSpacePoint> leftUpMost = findExtremePointsFromPoint<Direction::Left,Direction::Up>(bottomstart);

            //check if endcell already contains this Y
            for(FreeSpacePoint leftMost : leftDownMost) {
                if(leftMost.y.id == topstart.y.id && leftMost.x.id == topstart.x.id && std::abs( leftMost.y.t - topstart.y.t) < EPSILON){
                    veryshortcount++;
                    continue;
                }
                if(leftMost.y.id == topstart.y.id && clipper.clip(leftMost.y.t) == clipper.clip(topstart.y.t)) {
                    veryshortcount++;
                    continue;
                }
                bool contains = false;
                for (std::pair<CellPoint, FreeSpacePoint> startP: getCell(leftMost.x.id, leftMost.y.id)->importantDownYs)
                    if (clipper.clip(startP.first.y) == clipper.clip(leftMost.y.t)) {
                        contains = true;
                        break;
                    }
                if (!contains) {
                    CellPoint p = {leftMost.x.t, leftMost.y.t};
                    std::pair<CellPoint, FreeSpacePoint> pair = {p, topstart};
                    getCell(leftMost.x.id, leftMost.y.id)->importantDownYs.push_back(pair);
                }
            }
            for(FreeSpacePoint leftMost : leftUpMost) {
                if(leftMost.y.id == bottomstart.y.id && leftMost.x.id == bottomstart.x.id && std::abs(leftMost.y.t- bottomstart.y.t) < EPSILON){
                    veryshortcount++;
                    continue;
                }
                if(leftMost.y.id == bottomstart.y.id && clipper.clip(leftMost.y.t) == clipper.clip(bottomstart.y.t)) {
                    veryshortcount++;
                    continue;
                }
                bool contains = false;
                for (std::pair<CellPoint, FreeSpacePoint> startP: getCell(leftMost.x.id, leftMost.y.id)->importantUpYs)
                    if (clipper.clip(startP.first.y) == clipper.clip(leftMost.y.t)) {
                        contains = true;
                        break;
                    }
                if (!contains) {
                    CellPoint p = {leftMost.x.t, leftMost.y.t};
                    std::pair<CellPoint, FreeSpacePoint> pair = {p, bottomstart};
                    getCell(leftMost.x.id, leftMost.y.id)->importantUpYs.push_back(pair);
                }
            }


            if(!cell->getBoundary<Direction::Left>().is_empty() && cell->getBoundary<Direction::Left>().end != 1.0){
                //try to walk left as far as possible
                double t = cell->getBoundary<Direction::Left>().end;
                int tx = x;
                Cell* curcell = getCell(tx,y);
                while(tx >= 0 && curcell->left.begin <= t && t <= curcell->left.end) {
                    curcell = getCell(tx,y);
                    tx -= 1;
                }
                //offset tx
                tx++;
                if(curcell->leftPair.first.y > t) {
                    CellPoint p = curcell->getExtremePointAt<Direction::Left>(t);
                    FreeSpacePoint upperstart = {{x,0},{y,t}};
                    std::pair<CellPoint, FreeSpacePoint> pair = {p, upperstart};
                    getCell(tx, y)->importantUpYs.push_back(pair);
        }
    }
    std::cout << "Ignoring " << veryshortcount << " same-cell-candidates of length <= " << EPSILON << "\n";
}
*/
/*
bool FreeSpace::updateCellTopRight(int x, int y) {

    Cell* currentCell = getCell(x,y);

    if(currentCell->is_empty())
        return true;

    Interval fromBelow,fromLeft;
    if(x>0)
        fromLeft = getCell(x-1,y)->toRight;
    if(y>0)
        fromBelow = getCell(x,y-1)->toAbove;

    if(fromLeft.is_empty() && fromBelow.is_empty())
        return true;

    if(!fromLeft.is_empty())
        currentCell->toAbove = currentCell->top;
    else
        currentCell->toAbove = currentCell->top.maximum(fromBelow.begin);

    if(!fromBelow.is_empty())
        currentCell->toRight = currentCell->right;
    else
        currentCell->toRight = currentCell->right.maximum(fromLeft.begin);

    return false;
}*/
/*
bool FreeSpace::updateCellBottomLeft(int x, int y)  {

    Cell* currentCell = getCell(x,y);

    if(currentCell->is_empty())
        return true;

    Interval fromAbove,fromRight;
    if(x<nx-1)
        fromRight = getCell(x+1,y)->toLeft;
    if(y<ny-1)
        fromAbove = getCell(x,y+1)->toBottom;

    if(fromRight.is_empty() && fromAbove.is_empty())
        return true;

    if(!fromRight.is_empty())
        currentCell->toBottom = currentCell->bottom;
    else
        currentCell->toBottom = currentCell->bottom.minimum(fromAbove.end);

    if(!fromAbove.is_empty())
        currentCell->toLeft = currentCell->left;
    else
        currentCell->toLeft = currentCell->left.minimum(fromRight.end);

    return false;
}
*/
//TODO: maybe be really fast with trivial cells
template<Direction primary, Direction secondary>
std::vector<FreeSpacePoint> FreeSpace::findExtremePointsFromPoint(FreeSpacePoint start, ParamPoint globalLimit, bool applyToSecondary, int threadID, std::vector<FreeSpacePoint>* localMaxima ) {

    //globalLimit is always in secondary direction

    //startstuff
    bool primaryUpOrRight = (primary == Direction::Right || primary == Direction::Up);
    bool secondaryUpOrRight = (secondary == Direction::Right || secondary == Direction::Up);

    bool primaryRight = (primary == Direction::Right);
    bool primaryLeft = (primary == Direction::Left);
    bool primaryOnX = primaryRight || primaryLeft;
    bool primaryUp = (primary == Direction::Up);
    bool primaryDown = (primary == Direction::Down);
    bool primaryOnY = primaryUp || primaryDown;

    bool secondaryRight = (secondary == Direction::Right);
    bool secondaryLeft = (secondary == Direction::Left);
    bool secondaryOnX = secondaryRight || secondaryLeft;
    bool secondaryUp = (secondary == Direction::Up);
    bool secondaryDown = (secondary == Direction::Down);
    bool secondaryOnY = secondaryUp || secondaryDown;
    assert(secondaryOnX ^ secondaryOnY);

    bool right = (primary == Direction::Right || secondary == Direction::Right);
    bool up = (primary == Direction::Up || secondary == Direction::Up);

    bool applyToPrimary = !applyToSecondary;

    if(applyToSecondary) {
        if (secondaryUpOrRight && globalLimit.id == -1)
            globalLimit = {std::numeric_limits<int>::max(), 1.0};
    }else{
        if (primaryUpOrRight && globalLimit.id == -1)
            globalLimit = {std::numeric_limits<int>::max(), 1.0};
    }

    Interval localPrimaryLimit = Interval(),localSecondaryLimit = Interval();

    //initialize limits
    if(applyToSecondary){
        if(secondaryUpOrRight)
            localSecondaryLimit = Interval(0,globalLimit.t);
        else
            localSecondaryLimit = Interval(globalLimit.t,1);
    }else{
        if(primaryUpOrRight)
            localPrimaryLimit = Interval(0,globalLimit.t);
        else
            localPrimaryLimit = Interval(globalLimit.t,1);
    }

    Interval primaryStart, secondaryStart;
    switch(primary){
        case Direction::Right:
        case Direction::Left:
            primaryStart = secondaryUp?Interval(start.y.t,1):Interval(0,start.y.t);
            break;
        case Direction::Up:
        case Direction::Down:
            primaryStart = secondaryRight?Interval(start.x.t,1):Interval(0,start.x.t);
            break;
    }


    switch(secondary){
        case Direction::Right:
        case Direction::Left:
            secondaryStart = primaryUp?Interval(start.y.t,1):Interval(0,start.y.t);
            break;
        case Direction::Up:
        case Direction::Down:
            secondaryStart = primaryRight?Interval(start.x.t,1):Interval(0,start.x.t);
            break;
    }

    bool limitX = (applyToSecondary && secondaryOnX) || (applyToPrimary && primaryOnX);
    bool limitXRight = (applyToSecondary && secondaryRight) || (applyToPrimary && primaryRight);

    bool limitY = (applyToSecondary && secondaryOnY) || (applyToPrimary && primaryOnY);
    bool limitYUp = (applyToSecondary && secondaryUp) || (applyToPrimary && primaryUp);


    //initialize startcell
    Cell* startCell = getCell(start.x.id,start.y.id);

    if((limitY && start.y.id == globalLimit.id) ||(limitX && start.x.id == globalLimit.id)){
        startCell->getTo<primary>(threadID) = startCell->getBoundary<primary>().intersect(primaryStart).intersect(localSecondaryLimit);
        startCell->getTo<secondary>(threadID) = startCell->getBoundary<secondary>().intersect(secondaryStart).intersect(localPrimaryLimit);
    }else{
        startCell->getTo<primary>(threadID) = startCell->getBoundary<primary>().intersect(primaryStart);
        startCell->getTo<secondary>(threadID) = startCell->getBoundary<secondary>().intersect(secondaryStart);
    }
    reseterList[threadID].emplace_back(start.x.id,start.y.id);

    //this (above) needs to restrict the global limit...

    int xstep = right?1:-1;
    int ystep = up?1:-1;

    std::vector<FreeSpacePoint> extremePoints;

    //find initial extremePoint
    CellPoint p(-1,-1);
    if(primaryOnX) {
        p = startCell->getExtremeWithLimit<primary, secondary>(start.y.t);
    }else {
        p = startCell->getExtremeWithLimit<primary,secondary>(start.x.t);
    }
    FreeSpacePoint fp = {{start.x.id,p.x},{start.y.id,p.y}};
    if(applyToPrimary){
        if(fp.isFurtherThan<primary>(globalLimit)){
            p = startCell->getExtremePointAt<secondary>(globalLimit.t);
            fp = {{start.x.id,p.x},{start.y.id,p.y}};
        }
    }else{
        if(fp.isFurtherThan<secondary>(globalLimit)){
            p = startCell->getExtremePointAt<primary>(globalLimit.t);
            fp = {{start.x.id,p.x},{start.y.id,p.y}};
        }
    }
    extremePoints.emplace_back(fp);

    if(localMaxima != nullptr){
        //check if fp is a local maximum, by checking if it is on the boundary or whatever

        //fuck applyToPrimary
        assert(applyToSecondary);

        FreeSpacePoint plm = fp; //potentiallocalmaximum

        switch (primary) {
            case Direction::Up:
                if(plm.y.t < 1 && (plm.x == globalLimit || ((secondaryRight && plm.x.t<1) || (secondaryLeft && plm.x.t > 0))))
                    localMaxima->push_back(plm);
                break;
            case Direction::Down:
                if(plm.y.t > 0 && (plm.x == globalLimit || ((secondaryRight && plm.x.t<1) || (secondaryLeft && plm.x.t > 0))))
                    localMaxima->push_back(plm);
                break;
            case Direction::Left:
                if(plm.x.t > 0 && (plm.y == globalLimit || ((secondaryUp && plm.y.t<1) || (secondaryDown && plm.y.t > 0))))
                    localMaxima->push_back(plm);
                break;
            case Direction::Right:
                if(plm.x.t < 1 && (plm.y == globalLimit || ((secondaryUp && plm.y.t<1) || (secondaryDown && plm.y.t > 0))))
                    localMaxima->push_back(plm);
                break;
        }
    }

    //first handle the row
    //stopping condition has to respect the global limit

    int y = start.y.id;

    int tempstart = -1;
    int tempend = -1;
    int startIdx = -1;
    int endIdx = -1;

    if((up && !startCell->getTo<Direction::Up>(threadID).is_empty()) || (!up && !startCell->getTo<Direction::Down>(threadID).is_empty())) {
        tempstart = start.x.id;
        tempend = start.x.id;
    }
    //now step through the rest
    while(((y == start.y.id) || (startIdx != -1)) &&( (limitY)?(limitYUp?(y<=std::min(globalLimit.id,ny-1)):(y>=std::max(globalLimit.id,0))):isOk(0,y))){

        for(int x = (startIdx==-1)?start.x.id+xstep:startIdx;(limitX)?((limitXRight)?x<=std::min(globalLimit.id,nx-1):x>=std::max(globalLimit.id,0)):isOk(x,0);x+=xstep){
            Cell* cell = getCell(x,y);
            if(cell->is_empty()) {
                if(startIdx == -1)
                    break;
                if( (right)?x>endIdx:x<endIdx )
                    break;
                continue;
            }

            //TODO:update cell in here for optimizations

            //limitPair.snd encodes the lower limit (in secondary direction) on which we can walk in primary direction
            std::pair<bool,double> limitPair;
            if((limitY && y == globalLimit.id) ||(limitX && x == globalLimit.id))
                limitPair = updateCellAndGetLimit<primary,secondary>(x,y,threadID,localPrimaryLimit,localSecondaryLimit);
            else
                limitPair = updateCellAndGetLimit<primary,secondary>(x,y,threadID);
            reseterList[threadID].emplace_back(x,y);
            if(!limitPair.first){
                //then there is nothing to do in this cell as both "from"intervals are empty
                if( (right)?x>endIdx:x<endIdx )
                    break;
                continue;
            }



            if((up && !cell->getTo<Direction::Up>(threadID).is_empty()) || (!up && !cell->getTo<Direction::Down>(threadID).is_empty())) {
                if (tempstart == -1) {
                    tempstart = x;
                    tempend = x;
                } else {
                    tempend = x;
                }
            }

            //check if we stepped to far
            CellPoint primaryMostPoint(-1,-1);
            if(limitPair.second == -1)
                primaryMostPoint = cell->getExtreme<primary,secondary>();
            else
                primaryMostPoint = cell->getExtremeWithLimit<primary,secondary>(limitPair.second);
            FreeSpacePoint primaryMostFSP = {{x,primaryMostPoint.x},{y,primaryMostPoint.y}};
            if(applyToSecondary) {
                if (primaryMostFSP.isFurtherThan<secondary>(globalLimit)) {
                    CellPoint limit = cell->getExtremePointAt<primary>(globalLimit.t);
                    primaryMostFSP = {{x, limit.x},
                                      {y, limit.y}};
                }
            }else{
                if (primaryMostFSP.isFurtherThan<primary>(globalLimit)) {
                    CellPoint limit = cell->getExtremePointAt<secondary>(globalLimit.t);
                    primaryMostFSP = {{x, limit.x},
                                      {y, limit.y}};
                }
            }
            //update extremepoints
            bool canstepFurther = false;
            switch(secondary){
                case Direction::Up:
                    if(!(applyToSecondary && primaryMostFSP.y == globalLimit) && primaryMostFSP.y.t>=1.0 && y != ny - 1)
                        canstepFurther = true;
                    break;
                case Direction::Down:
                    if(!(applyToSecondary && primaryMostFSP.y == globalLimit) && primaryMostFSP.y.t<=0.0 && y != 0)
                        canstepFurther = true;
                    break;
                case Direction::Left:
                    if(!(applyToSecondary && primaryMostFSP.x == globalLimit) && primaryMostFSP.x.t<=0.0 && x != 0)
                        canstepFurther = true;
                    break;
                case Direction::Right:
                    if(!(applyToSecondary && primaryMostFSP.x == globalLimit) && primaryMostFSP.x.t>=1.0 && x!=nx-1)
                        canstepFurther = true;
                    break;
            }
            if(!canstepFurther) {
                if (!(extremePoints.empty() ) && primaryMostFSP.isFurtherThan<primary>(extremePoints[0]))
                    extremePoints.clear();
                if(extremePoints.empty() || (primaryOnX?(extremePoints[0].x == primaryMostFSP.x):(extremePoints[0].y == primaryMostFSP.y)))
                    extremePoints.emplace_back(primaryMostFSP);

                if(localMaxima != nullptr){
                    //check if fp is a local maximum, by checking if it is on the boundary or whatever

                    //fuck applyToPrimary
                    assert(applyToSecondary);

                    FreeSpacePoint plm = primaryMostFSP; //potentiallocalmaximum

                    switch (primary) {
                        case Direction::Up:
                            if(plm.y.t < 1 && (plm.x == globalLimit || ((secondaryRight && plm.x.t<1) || (secondaryLeft && plm.x.t > 0))))
                                localMaxima->push_back(plm);
                            break;
                        case Direction::Down:
                            if(plm.y.t > 0 && (plm.x == globalLimit || ((secondaryRight && plm.x.t<1) || (secondaryLeft && plm.x.t > 0))))
                                localMaxima->push_back(plm);
                            break;
                        case Direction::Left:
                            if(plm.x.t > 0 && (plm.y == globalLimit || ((secondaryUp && plm.y.t<1) || (secondaryDown && plm.y.t > 0))))
                                localMaxima->push_back(plm);
                            break;
                        case Direction::Right:
                            if(plm.x.t < 1 && (plm.y == globalLimit || ((secondaryUp && plm.y.t<1) || (secondaryDown && plm.y.t > 0))))
                                localMaxima->push_back(plm);
                            break;
                    }
                }
            }
        }


        //end of loop

        startIdx = tempstart;
        endIdx = tempend;
        tempstart = -1;
        tempend = -1;
        y += ystep;
    }
    for(auto pr : reseterList[threadID]){
        Cell* c = getCell(pr.first,pr.second);
        c->getTo<primary>(threadID) = Interval();
        c->getTo<secondary>(threadID) = Interval();
    }
    reseterList[threadID].clear();
    //reset(threadID);
    return extremePoints;
}

template<Direction primary, Direction secondary>
std::pair<bool,double> FreeSpace::updateCellAndGetLimit(int x, int y,int threadID, Interval primaryLimit, Interval secondaryLimit) {
    assert(!getCell(x,y)->isEmpty);
    int primaryFromX, primaryFromY;
    int secondaryFromX, secondaryFromY;

    //maybe put this in its own function...
    switch (primary) {
        case Direction::Up:
            primaryFromX = x;
            primaryFromY = y-1;
            break;
        case Direction::Left:
            primaryFromX = x+1;
            primaryFromY = y;
            break;
        case Direction::Right:
            primaryFromX = x-1;
            primaryFromY = y;
            break;
        case Direction::Down:
            primaryFromX = x;
            primaryFromY = y+1;
            break;
    }

    switch (secondary) {
        case Direction::Up:
            secondaryFromX = x;
            secondaryFromY = y-1;
            break;
        case Direction::Left:
            secondaryFromX = x+1;
            secondaryFromY = y;
            break;
        case Direction::Right:
            secondaryFromX = x-1;
            secondaryFromY = y;
            break;
        case Direction::Down:
            secondaryFromX = x;
            secondaryFromY = y+1;
            break;
    }
    Interval fromPrimary, fromSecondary;
    if(isOk(primaryFromX,primaryFromY))
        fromPrimary = getCell(primaryFromX,primaryFromY)->getTo<primary>(threadID);
    if(isOk(secondaryFromX,secondaryFromY))
        fromSecondary = getCell(secondaryFromX,secondaryFromY)->getTo<secondary>(threadID);

    if(fromPrimary.is_empty() && fromSecondary.is_empty())
        return {false,-1.0};

    double limit;

    if(!fromSecondary.is_empty()) {
        limit = -1;
    }else{
        //this is awful. this should also be templated
        switch (primary) {
            case Direction::Up:
            case Direction::Down:
                if(secondary == Direction::Left){
                    limit = fromPrimary.end;
                }else if(secondary == Direction::Right){
                    limit = fromPrimary.begin;
                }else{
                    assert(false);
                }
                break;
            case Direction::Left:
            case Direction::Right:
                if(secondary == Direction::Up){
                    limit = fromPrimary.begin;
                }else if(secondary == Direction::Down){
                    limit = fromPrimary.end;
                }else{
                    assert(false);
                }
                break;
        }
    }

    //update getTo's
    Cell* cell = getCell(x,y);
    if(!fromSecondary.is_empty()) {
        cell->getTo<primary>(threadID) = cell->getBoundary<primary>().intersect(secondaryLimit);
    }else {
        switch (primary) {
            case Direction::Up:
            case Direction::Down:
                    if(secondary == Direction::Left){
                        cell->getTo<primary>(threadID) = cell->getBoundary<primary>().minimum(fromPrimary.end).intersect(secondaryLimit);
                    }else if(secondary == Direction::Right){
                        cell->getTo<primary>(threadID) = cell->getBoundary<primary>().maximum(fromPrimary.begin).intersect(secondaryLimit);
                    }else{
                        assert(false);
                    }
                break;
            case Direction::Left:
            case Direction::Right:
                if(secondary == Direction::Down){
                    cell->getTo<primary>(threadID) = cell->getBoundary<primary>().minimum(fromPrimary.end).intersect(secondaryLimit);
                }else if(secondary == Direction::Up){
                    cell->getTo<primary>(threadID) = cell->getBoundary<primary>().maximum(fromPrimary.begin).intersect(secondaryLimit);
                }else{
                    assert(false);
                }
                break;
        }
    }


    if(!fromPrimary.is_empty()) {
        cell->getTo<secondary>(threadID) = cell->getBoundary<secondary>().intersect(primaryLimit);
    }else {
        switch (secondary) {
            case Direction::Up:
            case Direction::Down:
                if(primary == Direction::Left){
                    cell->getTo<secondary>(threadID) = cell->getBoundary<secondary>().minimum(fromSecondary.end).intersect(primaryLimit);
                }else if(primary == Direction::Right){
                    cell->getTo<secondary>(threadID) = cell->getBoundary<secondary>().maximum(fromSecondary.begin).intersect(primaryLimit);
                }else{
                    assert(false);
                }
                break;
            case Direction::Left:
            case Direction::Right:
                if(primary == Direction::Down){
                    cell->getTo<secondary>(threadID) = cell->getBoundary<secondary>().minimum(fromSecondary.end).intersect(primaryLimit);
                }else if(primary == Direction::Up){
                    cell->getTo<secondary>(threadID) = cell->getBoundary<secondary>().maximum(fromSecondary.begin).intersect(primaryLimit);
                }else{
                    assert(false);
                }
                break;
        }
    }

    return {true,limit};

}

bool FreeSpace::isOk(int x, int y) {
    return (0<=x && x < nx && 0<=y && y<ny);
}

std::vector<FreeSpacePoint> FreeSpace::findExtremePointsFromPointWrapper(Direction primary,Direction secondary,FreeSpacePoint start, ParamPoint globalLimit, bool applyToSecondary, int threadID, std::vector<FreeSpacePoint>* localMaxima) {
    switch (primary) {

        case Direction::Up:
            switch(secondary){
                case Direction::Up:
                case Direction::Down:
                    assert(false);
                case Direction::Left:
                    return findExtremePointsFromPoint<Direction::Up,Direction::Left>(start,globalLimit, applyToSecondary,threadID,localMaxima);
                case Direction::Right:
                    return findExtremePointsFromPoint<Direction::Up,Direction::Right>(start,globalLimit, applyToSecondary,threadID,localMaxima);
            }
            break;
        case Direction::Down:
            switch(secondary){
                case Direction::Up:
                case Direction::Down:
                    assert(false);
                case Direction::Left:
                    return findExtremePointsFromPoint<Direction::Down,Direction::Left>(start,globalLimit, applyToSecondary,threadID,localMaxima);
                case Direction::Right:
                    return findExtremePointsFromPoint<Direction::Down,Direction::Right>(start,globalLimit, applyToSecondary,threadID,localMaxima);
            }
            break;
        case Direction::Left:
            switch (secondary) {

                case Direction::Up:
                    return findExtremePointsFromPoint<Direction::Left,Direction::Up>(start,globalLimit, applyToSecondary,threadID,localMaxima);
                case Direction::Down:
                    return findExtremePointsFromPoint<Direction::Left,Direction::Down>(start,globalLimit, applyToSecondary,threadID,localMaxima);
                case Direction::Left:
                case Direction::Right:
                    assert(false);
            }
            break;
        case Direction::Right:
            switch (secondary) {

                case Direction::Up:
                    return findExtremePointsFromPoint<Direction::Right,Direction::Up>(start,globalLimit, applyToSecondary,threadID,localMaxima);
                case Direction::Down:
                    return findExtremePointsFromPoint<Direction::Right,Direction::Down>(start,globalLimit, applyToSecondary,threadID,localMaxima);
                case Direction::Left:
                case Direction::Right:
                    assert(false);
            }
            break;
    }
}
/*
void FreeSpace::identifyImportantUpDownYs(int l) {
    //ultrafastIdentify(l);
    //Length-limit: l
        int veryshortcount = 0;
        Interval clipper = Interval(0,1);
        for(int x=nx-1;x>=0;x--){
            if((nx-x)%10==0)
                std::cout << "Down to column " << x << std::endl;
            for(int y=ny-1;y>=0;y--) {

                Cell* cell = getCell(x,y);
                if(cell->isEmpty)
                    continue;


                //TODO: rename properly, this is nolonger topright
                CellPoint lefttop = cell->getExtreme<Direction::Left,Direction::Up>();
                CellPoint leftbottom = cell->getExtreme<Direction::Left,Direction::Down>();
                FreeSpacePoint topstart = {{x,lefttop.x},{y,lefttop.y}};
                FreeSpacePoint bottomstart = {{x,leftbottom.x},{y,leftbottom.y}};
                ParamPoint downlimit = {y-l+1,0.0};
                ParamPoint uplimit = {y+l-1,1.0};
                std::vector<FreeSpacePoint> leftDownMost = findExtremePointsFromPoint<Direction::Left,Direction::Down>(topstart,downlimit,true);
                std::vector<FreeSpacePoint> leftUpMost = findExtremePointsFromPoint<Direction::Left,Direction::Up>(bottomstart,uplimit,true);

                //check if endcell already contains this Y
                for(FreeSpacePoint leftMost : leftDownMost) {
                    if(leftMost.y.id == topstart.y.id && std::abs( leftMost.y.t - topstart.y.t) < EPSILON){
                        veryshortcount++;
                        continue;
                    }
                    if(leftMost.y.id == topstart.y.id && clipper.clip(leftMost.y.t) == clipper.clip(topstart.y.t)) {
                        veryshortcount++;
                        continue;
                    }
                    bool contains = false;
                    for (std::pair<CellPoint, FreeSpacePoint> startP: getCell(leftMost.x.id, leftMost.y.id)->importantDownYs)
                        if (clipper.clip(startP.first.y) == clipper.clip(leftMost.y.t)) {
                            contains = true;
                            break;
                        }
                    if (!contains) {
                        CellPoint p = {leftMost.x.t, leftMost.y.t};
                        std::pair<CellPoint, FreeSpacePoint> pair = {p, topstart};
                        getCell(leftMost.x.id, leftMost.y.id)->importantDownYs.push_back(pair);
                    }
                }
                for(FreeSpacePoint leftMost : leftUpMost) {
                    if(leftMost.y.id == bottomstart.y.id && std::abs(leftMost.y.t- bottomstart.y.t) < EPSILON){
                        veryshortcount++;
                        continue;
                    }
                    if(leftMost.y.id == bottomstart.y.id && clipper.clip(leftMost.y.t) == clipper.clip(bottomstart.y.t)) {
                        veryshortcount++;
                        continue;
                    }
                    bool contains = false;
                    for (std::pair<CellPoint, FreeSpacePoint> startP: getCell(leftMost.x.id, leftMost.y.id)->importantUpYs)
                        if (clipper.clip(startP.first.y) == clipper.clip(leftMost.y.t)) {
                            contains = true;
                            break;
                        }
                    if (!contains) {
                        CellPoint p = {leftMost.x.t, leftMost.y.t};
                        std::pair<CellPoint, FreeSpacePoint> pair = {p, bottomstart};
                        getCell(leftMost.x.id, leftMost.y.id)->importantUpYs.push_back(pair);
                    }
                }

                //TODO: wtf is  i probably have to walk left as long as the left-interval is contained in what i have checked so far... I dont know, what this trash if condition is
                if(!cell->getBoundary<Direction::Left>().is_empty() && cell->getBoundary<Direction::Left>().end != 1.0){
                    //try to walk left as far as possible
                    double t = cell->getBoundary<Direction::Left>().end;
                    int tx = x;
                    Cell* curcell = getCell(tx,y);
                    while(tx >= 0 && curcell->left.begin <= t && t <= curcell->left.end) {
                        curcell = getCell(tx,y);
                        tx -= 1;
                    }
                    //offset tx
                    tx++;
                    if(curcell->leftPair.first.y > t) {
                        CellPoint p = curcell->getExtremePointAt<Direction::Left>(t);
                        FreeSpacePoint upperstart = {{x,0},{y,t}};
                        std::pair<CellPoint, FreeSpacePoint> pair = {p, upperstart};
                        getCell(tx, y)->importantUpYs.push_back(pair);
                    }
                }
            }
        }
        std::cout << "Ignoring " << veryshortcount << " same-cell-candidates of length <= " << EPSILON << "\n";
}
*/
void FreeSpace::ultrafastIdentify(int l) {
    assert(l>0);
    upStarts.clear();
    downStarts.clear();
    for(int y=0;y<ny;y++){
        //will handle trivial starts later
        //upStarts.push_back({y, 0});
        //downStarts.push_back({y, 1});
        for(int x=nx-1;x>=0;x--){
            if(x==30 && y==15){
                int t = 0;
            }
            Cell* cell = getCell(x,y);

            if(cell->is_empty())
                continue;

            //we are now at an interesting cell

            if(cell->left.is_empty() || x == 0){
                CellPoint leftdown = cell->getExtreme<Direction::Left,Direction::Down>();
                CellPoint leftup = cell->getExtreme<Direction::Left,Direction::Up>();

                //check if we are at 0 or 1
                if(x == 0){
                    if(leftdown.x <= 0 && leftdown.y > 0){
                        upStarts.push_back({{x,leftdown.x},{y, leftdown.y}});
                    }
                    if(leftup.x <= 0 && leftup.y < 1){
                        downStarts.push_back({{x,leftup.x},{y,leftup.y}});
                    }
                }else{
                    if(leftdown.y > 0 && leftup.y<1){
                        upStarts.push_back({{x,leftdown.x},{y,leftdown.y}});
                        downStarts.push_back({{x,leftup.x},{y,leftup.y}});
                    }
                }

                //trivialstartcounter++;
            }else{
                bool down=false,up=false;
                if(!cell->bottom.is_empty() || //we can escape to bottom
                    cell->getExtreme<Direction::Right,Direction::Down>().y < cell->left.begin ) //down-extreme is below monotonicity
                    down = true;
                if(!cell->top.is_empty() ||
                cell->getExtreme<Direction::Right,Direction::Up>().y > cell->left.end)
                    up = true;
                if(!up && !down)
                    continue;
                //step left as long as possible, and then pushback monotonicity events
                Interval leftmonos = cell->left;
                int tx = x-1;
                while(tx>=1 && leftmonos.intersects(getCell(tx, y)->left)){
                    leftmonos = leftmonos.intersect(getCell(tx, y)->left);
                    tx -=1;
                }
                //tx,y is now the cell we want to know the monotonicities of
                Cell* endcell = getCell(tx,y);
                if(up && endcell->getExtreme<Direction::Left,Direction::Down>().y > leftmonos.end && leftmonos.end < 1.0) {
                    auto cp = endcell->getExtremePointAt<Direction::Left>(leftmonos.end);
                    upStarts.push_back({{tx,cp.x},{y,cp.y}});
                }
                if(down && endcell->getExtreme<Direction::Left,Direction::Up>().y < leftmonos.begin && leftmonos.begin > 0) {
                    auto cp = endcell->getExtremePointAt<Direction::Left>(leftmonos.begin);
                    downStarts.push_back({{tx,cp.x},{y,cp.y}});
                }
            }
        }
    }

    //now for the trivial up ones
    for(int id = ySize()-1;id > l-1;id--){
        for(int x=0;x<xSize();++x){
            Cell* cell = getCell(x,id);
            if(!cell->bottom.is_empty() && (x==0 || cell->bottom.begin > 0)){
                upStarts.push_back({{x,cell->bottom.begin},{id,0}});
            }
        }
    }
    for(int id = ySize()-1;id > l-1;id--){
        for(int x=0;x<xSize();++x){
            Cell* cell = getCell(x,id);
            if(!cell->top.is_empty() && (x==0 || cell->top.begin > 0)){
                downStarts.push_back({{x,cell->top.begin},{id,1}});
            }
        }
    }

    sort(upStarts.begin(),upStarts.end());
    upStarts.erase(std::unique(upStarts.begin(),upStarts.end()),upStarts.end());

    sort(downStarts.begin(),downStarts.end());
    downStarts.erase(std::unique(downStarts.begin(),downStarts.end()),downStarts.end());

    //std::cout << upStarts.size() << " " << downStarts.size() << std::endl;

    //ends need to be done in a more output sensitive way.
    //for every possible start at a startY i need to find all possible local maxima respecting the limit to the right

    /*Case 1: non trivial starts
     * the limit is determined by the y-coordinate of the rightmost maximum
     *
     * Case 2: trivial starts
     * the limit is determined by the y-coordinate of the rightmost maximum iff the length is exactly l, otherwise it is the lth vertex
     * */

    //FreeSpaceVisualizer fsv(*this);
    //fsv.show(true);
}


template<Direction secondary> std::vector<FreeSpacePoint> FreeSpace::findExtremePointsFromPointUntilEnds(FreeSpacePoint start,
                                                                           const std::vector<ParamPoint>& limits, int threadID) {

    //TODO: invalidate wrong ones
    assert(secondary == Direction::Up || secondary == Direction::Down);
    //globalLimit is always in secondary direction

    assert(limits.size() > 0);
    //assume it is sorted if Up, and reverse sorted if Down
    //assert(std::is_sorted(limits.begin(),limits.end()));

    ParamPoint globalLimit = limits.back();

    std::vector<FreeSpacePoint> result(limits.size());
    for(int i=0;i<result.size();i++){
        result[i].y = limits[i];
        result[i].x = {-1,-1};
    }

    Interval localPrimaryLimit = Interval(),localSecondaryLimit = Interval();

    //initialize limits
    if(secondary == Direction::Up)
        localSecondaryLimit = Interval(0,globalLimit.t);
    else
        localSecondaryLimit = Interval(globalLimit.t,1);

    Interval primaryStart, secondaryStart;
    primaryStart = (secondary==Direction::Up)?Interval(start.y.t,1):Interval(0,start.y.t);
    secondaryStart = Interval(start.x.t,1);

    //initialize startcell
    Cell* startCell = getCell(start.x.id,start.y.id);

    if(start.y.id == globalLimit.id){
        startCell->getTo<Direction::Right>(threadID) = startCell->getBoundary<Direction::Right>().intersect(primaryStart).intersect(localSecondaryLimit);
        startCell->getTo<secondary>(threadID) = startCell->getBoundary<secondary>().intersect(secondaryStart).intersect(localPrimaryLimit);
    }else{
        startCell->getTo<Direction::Right>(threadID) = startCell->getBoundary<Direction::Right>().intersect(primaryStart);
        startCell->getTo<secondary>(threadID) = startCell->getBoundary<secondary>().intersect(secondaryStart);
    }

    reseterList[threadID].emplace_back(start.x.id,start.y.id);

    //this (above) needs to restrict the global limit...

    int xstep = 1;
    int ystep = secondary==Direction::Up?1:-1;

    CellPoint startInStartCell = {start.x.t,start.y.t};

    CellPoint mostSecondaryYPossible = (startCell->getExtreme<secondary,Direction::Right>().x < start.x.t)?startCell->getExtremePointAt<secondary>(startInStartCell.x):startCell->getExtreme<secondary,Direction::Right>();
    for(int i=0;i<limits.size();++i){
        ParamPoint pp = limits[i];
        if(start.y.id == pp.id){
            if((secondary==Direction::Up && pp.t<=mostSecondaryYPossible.y) ||
                secondary==Direction::Down && pp.t >= mostSecondaryYPossible.y){
                CellPoint potentialNew = startCell->getExtremePointAt<Direction::Right>(pp.t);
                FreeSpacePoint potentialNewFSP = {{start.x.id,potentialNew.x},{start.y.id,potentialNew.y}};
                if(result[i].x < potentialNewFSP.x){
                    result[i] = potentialNewFSP;
                }
            }
        }
    }

    //first handle the row
    //stopping condition has to respect the global limit

    int y = start.y.id;

    int tempstart = -1;
    int tempend = -1;
    int startIdx = -1;
    int endIdx = -1;

    if((secondary == Direction::Up && !startCell->getTo<Direction::Up>(threadID).is_empty()) || (secondary == Direction::Down && !startCell->getTo<Direction::Down>(threadID).is_empty())) {
        tempstart = start.x.id;
        tempend = start.x.id;
    }
    //now step through the rest
    int baselimitsindex = 0;
    while(((y == start.y.id) || (startIdx != -1)) && ((secondary==Direction::Up)?(y<=std::min(globalLimit.id,ny-1)):(y>=std::max(globalLimit.id,0)))){

        //increase baselimitindex
        if(secondary == Direction::Up){
            while(limits[baselimitsindex].id < y){
                baselimitsindex++;
            }
        }else{
            while(limits[baselimitsindex].id > y){
                baselimitsindex++;
            }
        }

        for(int x = (startIdx==-1)?start.x.id+xstep:startIdx;isOk(x,0);x+=xstep){
            Cell* cell = getCell(x,y);
            if(cell->is_empty()) {
                if(startIdx == -1)
                    break;
                if(x>endIdx)
                    break;
                continue;
            }

            //TODO:update cell in here for optimizations

            //limitPair.snd encodes the lower limit (in secondary direction) on which we can walk in primary direction
            std::pair<bool,double> limitPair;
            if(y == globalLimit.id)
                limitPair = updateCellAndGetLimit<Direction::Right,secondary>(x,y,threadID,localPrimaryLimit,localSecondaryLimit);
            else
                limitPair = updateCellAndGetLimit<Direction::Right,secondary>(x,y,threadID);
            reseterList[threadID].emplace_back(x,y);
            if(!limitPair.first){
                //then there is nothing to do in this cell as both "from"intervals are empty
                if(x>endIdx)
                    break;
                continue;
            }



            //TODO: this is not exhaustive
            if((secondary==Direction::Up && !cell->getTo<Direction::Up>(threadID).is_empty()) || (secondary==Direction::Down && !cell->getTo<Direction::Down>(threadID).is_empty())) {
                if (tempstart == -1) {
                    tempstart = x;
                    tempend = x;
                } else {
                    tempend = x;
                }
            }

            //TODO: i am here

            if(limits[baselimitsindex].id == y) {

                //first from below

                Interval fromPrimary, fromSecondary;
                if (isOk(x - xstep, y)) {
                    fromPrimary = getCell(x - xstep, y)->getTo<Direction::Right>(threadID);
                }
                if (isOk(x, y - ystep)) {
                    fromSecondary = getCell(x, y - ystep)->getTo<secondary>(threadID);
                }
                if(fromPrimary.is_empty() && fromSecondary.is_empty())
                    continue;

                Interval primaryOKs, secondaryOKs;
                if (!fromPrimary.is_empty()) {
                    if (secondary == Direction::Up)
                        primaryOKs = Interval(fromPrimary.begin, cell->getExtreme<Direction::Up, Direction::Right>().y);
                    else
                        primaryOKs = Interval(cell->getExtreme<Direction::Down, Direction::Right>().y, fromPrimary.end);
                }
                if (!fromSecondary.is_empty()) {
                    if (secondary == Direction::Up) {
                        if (cell->getExtreme<secondary, Direction::Right>().x < fromSecondary.begin)
                            secondaryOKs = Interval(0, cell->getExtremePointAt<secondary>(fromSecondary.begin).y);
                        else
                            secondaryOKs = Interval(0, cell->getExtreme<secondary, Direction::Right>().y);
                    } else {
                        if (cell->getExtreme<secondary, Direction::Right>().x < fromSecondary.begin)
                            secondaryOKs = Interval(cell->getExtremePointAt<secondary>(fromSecondary.begin).y, 1);
                        else
                            secondaryOKs = Interval(cell->getExtreme<secondary, Direction::Right>().y, 1);
                    }
                }

                //now iterate over all limits, with same id as y
                for (int i = 0; baselimitsindex + i < limits.size() && limits[baselimitsindex + i].id == y; ++i) {
                    ParamPoint pp = limits[baselimitsindex + i];
                    if ((!primaryOKs.is_empty() && primaryOKs.begin <= pp.t && pp.t <= primaryOKs.end) ||
                        (!secondaryOKs.is_empty() && secondaryOKs.begin <= pp.t && pp.t <= secondaryOKs.end)) {
                        CellPoint potentialNew = cell->getExtremePointAt<Direction::Right>(pp.t);
                        FreeSpacePoint potentialNewFSP = {{x, potentialNew.x},
                                                          {y, potentialNew.y}};
                        if (result[baselimitsindex + i].x < potentialNewFSP.x) {
                            result[baselimitsindex + i] = potentialNewFSP;
                        }
                    }
                }
            }
        }

        //end of loop

        startIdx = tempstart;
        endIdx = tempend;
        tempstart = -1;
        tempend = -1;
        y += ystep;
    }
    //more aimed reset
    for(auto p : reseterList[threadID]){
        Cell* c = getCell(p.first,p.second);
        c->getTo<Direction::Right>(threadID) = Interval();
        c->getTo<secondary>(threadID) = Interval();
    }
    reseterList[threadID].clear();
    //reset(threadID);
    return result;
}

std::vector<FreeSpacePoint> FreeSpace::findExtremePointsFromPointUntilEndsWrapper(Direction dir, FreeSpacePoint start,
                                                                                  const std::vector<ParamPoint> &limits, int threadID) {
    switch (dir) {

        case Direction::Up:
            return findExtremePointsFromPointUntilEnds<Direction::Up>(start,limits,threadID);
            break;
        case Direction::Down:
            return findExtremePointsFromPointUntilEnds<Direction::Down>(start,limits,threadID);
            break;
        case Direction::Left:
        case Direction::Right:
            assert(false);
            break;
    }
    assert(false);
}

void FreeSpace::resetXSlice(int startX, int endX, int threadID) {
    for(int x = startX; x <= endX; ++x){
        for (int y=0;y<ny;++y){
            Cell* cell = getCell(x,y);
            cell->toAbove(threadID) = Interval(),cell->toRight(threadID)=Interval();
            cell->toBottom(threadID) = Interval(),cell->toLeft(threadID)=Interval();
        }
    }
}

template<Direction>
std::vector<std::vector<FreeSpacePoint>>
FreeSpace::findExtremePointsFromYUntilEnds(ParamPoint start, const std::vector<ParamPoint> &limits, int threadID) {
    //TODO:
    //first handle startRow initializing all possible start cells
    //Then do everything as before
    return std::vector<std::vector<FreeSpacePoint>>();
}

bool Cell::is_empty() const {
    return isEmpty;
}
/*
Point Cell::rightMostAtOrAbove(double y) const {
    assert(0<=y && y<=1);
    if(y <= rightMost.y())
        return rightMost;
    return rightMostAt(y);
}*/

CellPoint Cell::rightMostAt(double y) const {
    assert(0<=y && y<=1);
    //Point p = c + (d-c)*y;
    Interval intersection;
    //IntersectionAlgorithm::intersection_interval(p,delta,a,b,&intersection);
    //out with the old, in with the new:
    IntersectionAlgorithm::intersection_interval(lra1,lrb1,y*lrb2,lrc1,y*lrc2,y*y*lrc3,delta, &intersection);

    assert(!intersection.is_empty());
    return CellPoint(intersection.end,y);
}
/*
Point Cell::leftMostAtOrBelow(double y) const {
    assert(0<=y && y<=1);
    if(y >= leftMost.y())
        return leftMost;
    return leftMostAt(y);
}*/

CellPoint Cell::leftMostAt(double y) const {
    assert(0<=y && y<=1);
    //Point p = c + (d-c)*y;
    Interval intersection;
    //IntersectionAlgorithm::intersection_interval(p,delta,a,b,&intersection);
    //out with the old, in with the new:
    IntersectionAlgorithm::intersection_interval(lra1,lrb1,y*lrb2,lrc1,y*lrc2,y*y*lrc3,delta, &intersection);

    assert(!intersection.is_empty());
    return CellPoint(intersection.begin,y);
}
/*
Point Cell::topMostAtOrRight(double x) const {
    assert(0<=x && x<=1);
    if(x<=topMost.x())
        return topMost;
    return topMostAt(x);
}*/

CellPoint Cell::topMostAt(double x) const {
    assert(0<=x && x<=1);
    //Point p = a + (b-a)*x;
    Interval intersection;
    //IntersectionAlgorithm::intersection_interval(p,delta,c,d,&intersection);
    IntersectionAlgorithm::intersection_interval(tba1,tbb1,x*tbb2,tbc1,x*tbc2,x*x*tbc3,delta,&intersection);
    assert(!intersection.is_empty());
    return {x,intersection.end};
}
/*
Point Cell::bottomMostAtOrLeft(double x) const {
    assert(0<=x && x<=1);
    if(x>=bottomMost.x())
        return bottomMost;
    return bottomMostAt(x);
}*/
CellPoint Cell::bottomMostAt(double x) const {
    assert(0<=x && x<=1);
    //Point p = a + (b-a)*x;
    Interval intersection;
    //IntersectionAlgorithm::intersection_interval(p,delta,c,d,&intersection);
    IntersectionAlgorithm::intersection_interval(tba1,tbb1,x*tbb2,tbc1,x*tbc2,x*x*tbc3,delta,&intersection);
    assert(!intersection.is_empty());
    return {x,intersection.begin};
}


template<Direction dir>
Interval &Cell::getBoundary() {
    switch(dir){
        case Direction::Right:
            return right;
        case Direction::Left:
            return left;
        case Direction::Up:
            return top;
        case Direction::Down:
            return bottom;
    }
}

template<Direction dir>
Interval &Cell::getTo(int threadID) {
    switch(dir){
        case Direction::Right:
            return toRight(threadID);
        case Direction::Left:
            return toLeft(threadID);
        case Direction::Up:
            return toAbove(threadID);
        case Direction::Down:
            return toBottom(threadID);
    }
}

Cell::Cell(Point &a,Point &b,Point &c,Point &d,int tc) : a(a),b(b),c(c),d(d),toAboves(tc), toRights(tc), toBottoms(tc), toLefts(tc) {
}


bool FreeSpacePoint::isLeftOf(const FreeSpacePoint &other) {
    return x<=other.x;
}

bool FreeSpacePoint::isAbove(const FreeSpacePoint &other) {
    return y>=other.y;
}

bool FreeSpacePoint::isRightOf(const FreeSpacePoint &other) const {
    return x>=other.x;
}
bool FreeSpacePoint::isBelow(const FreeSpacePoint &other) const {
    return y<=other.y;
}

CellPoint::CellPoint(double _x, double _y) {
x=_x;y=_y;
}

CellPoint::CellPoint() {
x=-1;y=-1;
}
