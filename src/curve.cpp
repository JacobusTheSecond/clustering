#include <fstream>
#include "curve.h"
#include "geometry_basics.h"

Curve::Curve(const Points& points)
	: points(points), prefix_length(points.size())
{
	if (points.empty()) { return; }

	auto const& front = points.front();
	extreme_points = {front,front};
	prefix_length[0] = 0;

	for (PointID i = 1; i < points.size(); ++i)
	{
		auto segment_distance = points[i - 1].dist(points[i]);
		prefix_length[i] = prefix_length[i - 1] + segment_distance;

        extreme_points.min = extreme_points.min.min(points[i]);
        extreme_points.max = extreme_points.max.max(points[i]);
	}
}

void Curve::push_back(Point const& point)
{
    if(size() == 0)
        extreme_points = {point,point};
	if (prefix_length.size()) {
		auto segment_distance = points.back().dist(point);
		prefix_length.push_back(prefix_length.back() + segment_distance);
	}
	else {
		prefix_length.push_back(0);
	}

    extreme_points.min = extreme_points.min.min(point);
    extreme_points.max = extreme_points.max.max(point);

	points.push_back(point*1.0);
}

auto Curve::getExtremePoints() const -> ExtremePoints const&
{
	return extreme_points;
}

distance_t Curve::getUpperBoundDistance(Curve const& other) const
{
	auto const& extreme1 = this->getExtremePoints();
	auto const& extreme2 = other.getExtremePoints();

	Point min_point = extreme1.min.min(extreme2.min);
	Point max_point = extreme1.max.max(extreme2.max);

	return min_point.dist(max_point);
}

Curve::Curve(std::string filename,dimensions_t d) {
    std::ifstream file(filename);
    std::string line;
    while(file){
        double number;
        Point p(0);
        for (int i=0;i<d;i++){
            file >> number;
            p.push_back(number);
        }
        //somehow the last point is added 96 times without this... very weird
        if(file)
            push_back(p);
    }

    prefix_length = std::vector<distance_t>(points.size());
    prefix_length[0] = 0;
    for(int i=1;i<points.size();++i){
        prefix_length[i] = prefix_length[i-1] + (points[i-1].dist(points[i]));
    }
}
/*
std::vector<BreakpointCandidate> Curve::computeBreakpointCandidates(double delta) {
    assert(this->size()>=2);
    assert(delta>0);
    std::vector<BreakpointCandidate> candidates;
    candidates.push_back({0,this->operator[](0),0,this->operator[](0)});
    for (int i = 0; i < this->size()-2; ++i) {
        Point x = this->operator[](i);
        Point y = this->operator[](i+1);
        Point v = y-x;
        double length = sqrt(x.dist_sqr(y));
        double stepsize = delta/length;
        int ints = 0;
        for (double s = 0;s<1;ints++,s=stepsize*ints) {
            int intt = ints+1;
            for (double t = stepsize*intt; t < 1; intt++,t = stepsize*intt) {
                candidates.push_back({i, x + (v * s),i, x + (v * t)});
            }
            candidates.push_back({i, x + (v * s),i, y});
            for(int j=i+1;j<size()-2;++j){
                intt = 1;
                double slength = sqrt(operator[](j).dist_sqr(operator[](j+1)));
                double sstepsize = delta/slength;
                for(double t=sstepsize*intt;t<1;intt++,t=sstepsize*intt){
                    candidates.push_back({i,x+(v*s),j, operator[](j)*(1-t)+ operator[](j+1)*t});
                }
                candidates.push_back({i,x+(v*s),j, operator[](j+1)});

            }
        }
        candidates.push_back({i,y,i,y});
    }
    return candidates;
}

Curves Curve::computeCandidates(double delta) {
    Curves curves;
    for (int startEdge = 0;startEdge<size()-1;++startEdge){
        std::cout << startEdge << "/" << size()-2 << "\n";
        for (int endEdge = startEdge;endEdge<size()-1;++endEdge){
            Points subpoints(&(operator[](startEdge+1)),&(operator[](endEdge+1)));
            double startLength = operator[](startEdge).dist(operator[](startEdge+1));
            double endLength = operator[](endEdge).dist(operator[](endEdge+1));
            double startStepSize = delta/startLength;
            double endStepSize = delta/endLength;
            Point ps = operator[](startEdge);
            Point pt = operator[](endEdge);
            Point vs = operator[](startEdge+1)- operator[](startEdge);
            Point vt = operator[](endEdge+1) - operator[](endEdge);
            int scount=0;
            for (double s = scount*startStepSize;s<=1;scount++,s=scount*startStepSize){
                int tcount=0;
                for (double t = tcount*endStepSize;t<=1;tcount++,t=tcount*endStepSize){
                    Curve candidate;
                    candidate.push_back(ps + vs * s);
                    for (auto subpoint : subpoints)
                        candidate.push_back(subpoint);
                    candidate.push_back(pt+vt*t);
                    curves.push_back(candidate);
                }
            }
        }
    }
    return curves;
}
*/
//TODO: add length restriction capability
//what does it mean to be marked...
//we dont need to step into marked cells.
//TODO: mark cells if we reached them from below
//TODO: mark cells, if we can reach the rightmost point from left
/*
std::vector<Curve> Curve::canonicalPointCandidates(const Curve &c, double delta) {
    std::vector<Curve> candidates;
    //for every cell
    //step as far as possible
    //after updating every cell, go over every cell, and compute the candidate(s?)
    FreeSpace fs = FreeSpace(c,*this,delta);
    std::cout << "FreeSpaceSize:" << fs.xSize()<<","<<fs.ySize()<<"\n";
    for(int j=0;j<fs.ySize();j++){
        std::cout << "Handling row "<<j<<"/"<<fs.ySize()<<"\n";
        int i=0;
        while (i < fs.xSize() && (fs.getCell(i,j)->leftMarked || fs.getCell(i,j)->is_empty()) ) {
            i++;
        }
        while(i<fs.xSize()){
            Cell* cell = fs.getCell(i,j);
            assert(!cell->is_empty());
            Point leftMost = cell->leftMost;
            Point p = cell->c + (cell->d-cell->c)*cell->leftMost.y();
            Point q = cell->c + (cell->d - cell->c) * cell->rightMost.y();

            //handle current cell
            Curve bestCandidate;
            bestCandidate.push_back(p);
            bestCandidate.push_back(q);
            int startIndex = i;
            int startParam = cell->leftMost.x();
            int endIndex = i;
            double endParam = cell->rightMost.x();
            cell->leftMarked = true;
            if(cell->leftMost.y()<= cell->right.begin)
                cell->rightMarked = true;

            //start search
            cell->toAbove = cell->top.maximum(leftMost.x());
            cell->toRight = cell->right.maximum(leftMost.y());

            //now run wild
            for (int y = j;y < fs.ySize(); ++y) {
                Points subpoints;
                if(y>j)
                    subpoints = Points(&(operator[](j + 1)), &(operator[](y)));

                for (int x = ((y == j) ? i + 1 : i); x < fs.xSize(); ++x) {
                    Cell *endcell = fs.getCell(x, y);
                    //if the current cell is rightMarked, then there is no need, to process any further
                    if(endcell->rightMarked)
                        break;
                    Interval fromLeft, fromBelow;
                    if (y > 0)
                        fromBelow = fs.getCell(x, y - 1)->toAbove;
                    if (x > 0)
                        fromLeft = fs.getCell(x - 1, y)->toRight;

                    if (fromBelow.is_empty() && fromLeft.is_empty())
                        continue;

                    bool entireTopReachable = !fromLeft.is_empty() || fromBelow.begin <= endcell->top.begin;
                    bool entireRightReachable = !fromBelow.is_empty() || fromLeft.begin <= endcell->right.begin;

                    //a cell is rightmarked, if we have reached everything already
                    if(entireTopReachable && entireRightReachable)
                        endcell->rightMarked = true;

                    bool startLeftOfLeftMost = startIndex<x || (startIndex == x && startParam <= endcell->leftMost.x());
                    bool entireLeftReached = !fromLeft.is_empty() && fromLeft.begin <= endcell->left.begin;
                    bool leftMostPointAtBottom = endcell->leftMost.y() <= 0;

                    //a cell is leftmarked, if there is no need, to start from its leftmost point
                    if( (endcell->rightMarked && startLeftOfLeftMost) ||
                            (entireLeftReached) ||
                            (leftMostPointAtBottom))
                        endcell->leftMarked = true;

                    //first handle outgoing
                    if (!fromLeft.is_empty()) {
                        endcell->toAbove = endcell->top;
                    } else {
                        endcell->toAbove = endcell->top.maximum(fromBelow.begin);
                    }

                    if (!fromBelow.is_empty()) {
                        endcell->toRight = endcell->right;
                    } else {
                        endcell->toRight = endcell->right.maximum(fromLeft.begin);
                    }

                    //only update candidate if conditions are met
                    if ((endcell->right.is_empty() || x == fs.xSize()-1)) {
                        double t = fromBelow.is_empty() ? std::max(fromLeft.begin,endcell->rightMost.y())
                                                            : endcell->rightMost.y();
                        q = endcell->c + (endcell->d - endcell->c) * t;
                        Interval intersection = IntersectionAlgorithm::intersection_interval(q,delta,endcell->a,endcell->b);
                        if(x>endIndex || (x==endIndex && intersection.end > endParam))
                            endIndex = x;
                        endParam = intersection.end;
                        bestCandidate = Curve();
                        bestCandidate.push_back(p);
                        for (const auto &subpoint: subpoints)
                            bestCandidate.push_back(subpoint);
                        bestCandidate.push_back(q);
                    }
                }
            }
            candidates.push_back(bestCandidate);

            fs.reset();

            //go to next interesting cell
            i++;
            while (i < fs.xSize() && (fs.getCell(i,j)->leftMarked || fs.getCell(i,j)->is_empty()) ) {
                i++;
            }
        }
    }
    return candidates;
}*/

void Curve::clear() {
    points.clear();
}

Curve::Curve(const Curve& c, std::vector<int> times) {
    for(int startindex = 0; startindex<c.size()-1;++startindex){
        Point s = c[startindex];
        Point t = c[startindex+1];
        int stime = times[startindex];
        int ttime = times[startindex+1];
        for(int i=stime;i<ttime;++i){
            push_back(s + (t-s)*((double)(i-stime)/(double)(ttime-stime)));
        }
    }
    push_back(c.back());
}
/*
std::vector<Curve> Curve::canonicalPointCandidates(const Curve &T1, const Curve &T2, double delta) {


    //if(T1 == T2){
    //    return canonicalPointCandidates(T1,delta);
    //}

    std::cout << "Basecurve complexity: "<<size()<<"\n";
    std::cout << "T1 complexity: "<<T1.size()<<"\n";
    std::cout << "T2 complexity: "<<T2.size()<<"\n";


    FreeSpace ff1(*this,T1,delta);
    FreeSpace ff2(*this,T2,delta);

    ff1.identifyImportantUpDownYs();
    ff2.identifyImportantUpDownYs();

    FreeSpaceVisualizer fsv1 = FreeSpaceVisualizer(ff1);
    fsv1.show();
    FreeSpaceVisualizer fsv2 = FreeSpaceVisualizer(ff2);
    fsv2.show();

    //provisory candidatecount
    int count = 0,preprunecount = 0;

    //iterate over edges on B
    for(int edgeIdxB = 0; edgeIdxB < size()-1;edgeIdxB++){
        //identify important Y (and corresponding cell)
        for(int edgeIdxT1 = 0; edgeIdxT1 < T1.size()-1;edgeIdxT1++){

            ///                             ///
            /// FIRST OF GO TOWARDS THE TOP ///
            ///                             ///

            if(!ff1.getCell(edgeIdxT1, edgeIdxB)->importantDownYs.empty()){

                std::vector<std::pair<ParamPoint, ParamPoint>> endpoints;

                //now FOR EVERY START

                Cell *startCellInFF1 = ff1.getCell(edgeIdxT1, edgeIdxB);

                //we assume, that importantYs is sorted

                for (std::pair<CellPoint, FreeSpacePoint> paramPointPair: startCellInFF1->importantDownYs) {

                    CellPoint paramPoint = paramPointPair.first;
                    ParamPoint ystart = {edgeIdxB, paramPoint.y};

                    //Identify limits (in y-direction) for B
                    //first upwards:
                    FreeSpacePoint importantStart = {{edgeIdxT1, paramPoint.x}, ystart};
                    std::vector<FreeSpacePoint> up_extremes = ff1.findExtremePointsFromPointWrapper(Direction::Up,
                                                                                                    Direction::Right,
                                                                                                    importantStart);

                    ParamPoint limit = up_extremes[0].y;

                    //now we step through T2, finding the first cell, that has a point between ystart and limit

                    int startOnT2 = 0;

                    while (startOnT2 < ff2.xSize()) {

                        //first zoom to first nonempty cell
                        Cell *cell = ff2.getCell(startOnT2, edgeIdxB);

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



                        CellPoint leftMostInT2Cell = cell->getExtremePointAt<Direction::Left>(importantStart.y.t);
                        FreeSpacePoint leftMostInCellFSP = {{startOnT2, leftMostInT2Cell.x},
                                                            {edgeIdxB,  leftMostInT2Cell.y}};

                        assert(leftMostInCellFSP.y >= importantStart.y && leftMostInCellFSP.y <= limit);


                        std::vector<FreeSpacePoint> localExtremes = ff2.findExtremePointsFromPointWrapper(
                                Direction::Right, Direction::Up, leftMostInCellFSP, limit);

                        ParamPoint maximumOnT2 = localExtremes[0].x;

                        //add maximum to endpoints
                        preprunecount++;
                        if (endpoints.empty() || endpoints.back().second < localExtremes[0].x) {
                            endpoints.emplace_back(localExtremes[0].y, localExtremes[0].x);
                            count++;
                        }
                        startOnT2++;
                    }
                }
            }///end of topscope

            ///                             ///
            /// THEN GO TOWARDS THE BOTTOM  ///
            ///                             ///

            if(!ff1.getCell(edgeIdxT1, edgeIdxB)->importantUpYs.empty()){

                std::vector<std::pair<ParamPoint, ParamPoint>> endpoints;

                //now FOR EVERY START

                Cell *startCellInFF1 = ff1.getCell(edgeIdxT1, edgeIdxB);

                //we assume, that importantYs is sorted

                for (std::pair<CellPoint, FreeSpacePoint> paramPointPair: startCellInFF1->importantUpYs) {

                    CellPoint paramPoint = paramPointPair.first;
                    ParamPoint ystart = {edgeIdxB, paramPoint.y};

                    //Identify limits (in y-direction) for B
                    //first upwards:
                    FreeSpacePoint importantStart = {{edgeIdxT1, paramPoint.x}, ystart};
                    std::vector<FreeSpacePoint> down_extremes = ff1.findExtremePointsFromPointWrapper(Direction::Down,
                                                                                                    Direction::Right,
                                                                                                    importantStart);

                    ParamPoint limit = down_extremes[0].y;

                    //now we step through T2, finding the first cell, that has a point between ystart and limit

                    int startOnT2 = 0;

                    while (startOnT2 < ff2.xSize()) {

                        //first zoom to first nonempty cell
                        Cell *cell = ff2.getCell(startOnT2, edgeIdxB);

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



                        CellPoint leftMostInT2Cell = cell->getExtremePointAt<Direction::Left>(importantStart.y.t);
                        FreeSpacePoint leftMostInCellFSP = {{startOnT2, leftMostInT2Cell.x},
                                                            {edgeIdxB,  leftMostInT2Cell.y}};

                        assert(leftMostInCellFSP.y <= importantStart.y && leftMostInCellFSP.y >= limit);


                        std::vector<FreeSpacePoint> localExtremes = ff2.findExtremePointsFromPointWrapper(
                                Direction::Right, Direction::Down, leftMostInCellFSP, limit);

                        ParamPoint maximumOnT2 = localExtremes[0].x;

                        //add maximum to endpoints
                        preprunecount++;
                        if (endpoints.empty() || endpoints.back().second < localExtremes[0].x) {
                            endpoints.emplace_back(localExtremes[0].y, localExtremes[0].x);
                            count++;
                        }
                        startOnT2++;
                    }
                }
            }///end of topscope

        }
    }
    std::cout << "Candidates generated: " << count << "/" << preprunecount << "\n";
}
*/
Point Curve::eval(ParamPoint t) {
    Point ps = operator[](t.id);
    Point pt = operator[](t.id+1);
    return ps + (pt - ps)*t.t;
}

Curve Curve::constructSubcurve(ParamPoint s, ParamPoint t) {
    Curve result;

    result.push_back(eval(s));

    int step = (s>t)?-1:1;

    for(int i = s.id+std::max(step,0);(s>t)?i>t.id:i<=t.id;i+=step){
        result.push_back(operator[](i));
    }
    result.push_back(eval(t));
    return result;
}
