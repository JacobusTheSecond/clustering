//
// Created by Jacobus Conradi on 26.01.22.
//

#ifndef CLUSTERING_FREESPACE_H
#define CLUSTERING_FREESPACE_H

#include "geometry_basics.h"
#include "curve.h"

struct CellPoint{
    CellPoint();
    CellPoint(double _x, double _y);

    double x;double y;};

struct FreeSpacePoint{
    ParamPoint x;
    ParamPoint y;

    bool isLeftOf(const FreeSpacePoint& other);
    bool isAbove(const FreeSpacePoint& other);
    bool isRightOf(const FreeSpacePoint& other) const;
    bool isBelow(const FreeSpacePoint& other) const;



    template<Direction dir>
    bool isFurtherThan(const FreeSpacePoint &other) {
        switch (dir) {
            case Direction::Up:
                return isAbove(other);
            case Direction::Right:
                return isRightOf(other);
            case Direction::Down:
                return isBelow(other);
            case Direction::Left:
                return isLeftOf(other);
        }
    }

    template<Direction dir>
    bool isFurtherThan(const ParamPoint &other) {
        switch (dir) {
            case Direction::Up:
                return y > other;
            case Direction::Right:
                return x>other;
            case Direction::Down:
                return y<other;
            case Direction::Left:
                return x<other;
        }
    }

    bool operator==(FreeSpacePoint const& other) const {
        return x == other.x && y == other.y;
    }

    bool operator<(FreeSpacePoint const& other) const {
        return x < other.x;
    }
};
//TODO: multithreading this bitch
struct Cell{
    Cell(Point &a, Point &b, Point &c, Point &d,int threadcount = 1);

    //topbottom values
    double tba1;// = (d-c).length_sqr();
    double tbb1;// = (c-a).dot(d-c);
    double tbb2;// = (b-a).dot(d-c);
    double tbc1;// = (c-a).length_sqr();
    double tbc2;// = (c-a).dot(b-a);
    double tbc3;// = (b-a).length_sqr();

    //leftright values
    double lra1;// = (b-a).length_sqr();
    double lrb1;// = (a-c).dot(b-a);
    double lrb2;// = (d-c).dot(b-a);
    double lrc1;// = (a-c).length_sqr();
    double lrc2;// = (a-c).dot(d-c);
    double lrc3;// = (d-c).length_sqr();

    bool isEmpty;
    Point& a,b,c,d;
    Interval left,top,right,bottom;

    std::pair<CellPoint,CellPoint> leftPair, rightPair, topPair, bottomPair;

    //Interval toAbove,toRight,toLeft,toBottom;

    Interval& toAbove(int threadID = 0){return toAboves[threadID];}
    Interval& toRight(int threadID = 0){return toRights[threadID];}
    Interval& toBottom(int threadID = 0){return toBottoms[threadID];}
    Interval& toLeft(int threadID = 0){return toLefts[threadID];}

    std::vector<Interval> toAboves;
    std::vector<Interval> toRights;
    std::vector<Interval> toLefts;
    std::vector<Interval> toBottoms;

    //std::vector<std::pair<CellPoint,FreeSpacePoint>> importantUpYs;
    //std::vector<std::pair<CellPoint,FreeSpacePoint>> importantDownYs;
    double delta;

    [[nodiscard]] bool is_empty() const;

    template<Direction> Interval& getBoundary();
    template<Direction> Interval& getTo(int threadID=0);

    CellPoint leftMostAt(double y) const;
    CellPoint rightMostAt(double y) const;
    CellPoint topMostAt(double x) const;
    CellPoint bottomMostAt(double x) const;

    template<Direction dir> CellPoint getExtremePointAt(double limit){
        switch(dir){
            case Direction::Up:
                return topMostAt(limit);
            case Direction::Right:
                return rightMostAt(limit);
            case Direction::Down:
                return bottomMostAt(limit);
            case Direction::Left:
                return leftMostAt(limit);
        }
    }
    template<Direction primary, Direction secondary>
    CellPoint getExtreme() {
        switch (primary) {
            case Direction::Up:
                if(secondary == Direction::Left)
                    return topPair.first;
                if(secondary == Direction::Right)
                    return topPair.second;
                assert(false);
            case Direction::Down:
                if(secondary == Direction::Left)
                    return bottomPair.first;
                if(secondary == Direction::Right)
                    return bottomPair.second;
                assert(false);
            case Direction::Left:
                if(secondary == Direction::Down)
                    return leftPair.first;
                if(secondary == Direction::Up)
                    return leftPair.second;
                assert(false);
            case Direction::Right:
                if(secondary == Direction::Down)
                    return rightPair.first;
                if(secondary == Direction::Up)
                    return rightPair.second;
                assert(false);
        }
    }
    template<Direction primary, Direction secondary>
    CellPoint getExtremeWithLimit(double limit) {

        //limit is always in secondary direction

        switch (primary) {
            case Direction::Up:
                //maximize upwards
                if(secondary == Direction::Left){
                    if(limit > topPair.first.x){
                        return topPair.first;
                    }else{
                        return topMostAt(limit);
                    }
                }
                if(secondary == Direction::Right){
                    if(limit < topPair.second.x){
                        return topPair.second;
                    }else{
                        return topMostAt(limit);
                    }
                }
                assert(false);
            case Direction::Down:
                if(secondary == Direction::Left){
                    if(limit > bottomPair.first.x){
                        return bottomPair.first;
                    }else{
                        return bottomMostAt(limit);
                    }
                }
                if(secondary == Direction::Right){
                    if(limit < bottomPair.second.x){
                        return bottomPair.second;
                    }else{
                        return bottomMostAt(limit);
                    }
                }
                assert(false);
            case Direction::Left:
                if(secondary == Direction::Down) {
                    if (limit > leftPair.first.y) {
                        return leftPair.first;
                    } else {
                        return leftMostAt(limit);
                    }
                }
                if(secondary == Direction::Up){
                    if(limit < leftPair.second.y){
                        return leftPair.second;
                    }else{
                        return leftMostAt(limit);
                    }
                }
                assert(false);
            case Direction::Right:
                if(secondary == Direction::Down) {
                    if (limit > rightPair.first.y) {
                        return rightPair.first;
                    } else {
                        return rightMostAt(limit);
                    }
                }
                if(secondary == Direction::Up){
                    if(limit < rightPair.second.y){
                        return rightPair.second;
                    }else{
                        return rightMostAt(limit);
                    }
                }
                assert(false);
        }
    }
};


class FreeSpace{

private:
    int nx,ny;
    std::vector<std::vector<Cell>> data;
    double delta;

    //bool updateCellTopRight(int x, int y);
    //bool updateCellBottomLeft(int x, int y);
    template<Direction,Direction> std::pair<bool,double> updateCellAndGetLimit(int x, int y,int threadID, Interval primaryLimit = Interval(0,1), Interval secondaryLimit = Interval(0,1));
    bool isOk(int x, int y);
    //template<Direction,Direction> bool updateCell(int x, int y, int threadID);

    std::vector<std::vector<std::pair<int,int>>> reseterList;


public:
    std::vector<FreeSpacePoint> upStarts, downStarts;
    std::vector<FreeSpacePoint> pointvector;
    FreeSpace(FreeSpace const &) = default;
    FreeSpace(const Curve& B, const Curve& T,double delta, int threadcout = 1);
    void prettyPrint(int i,int j);
    int xSize();
    int ySize();
    Cell* getCell(int,int);
    void ultrafastIdentify(int l);
    //void identifyImportantUpDownYs();
    //void identifyImportantUpDownYs(int l);
    //void identifyImportantYs();
    void reset(int threadID = 0);
    void resetXSlice(int startX, int endX, int threadID = 0);
    void resetYSlice(int startY, int endY, int threadID = 0);
    //template<Direction> FreeSpacePoint findExtremPointInDirectionFromPoint(FreeSpacePoint start, ParamPoint globalLimit = {-1,0.});
    template<Direction,Direction> std::vector<FreeSpacePoint> findExtremePointsFromPoint(FreeSpacePoint start, ParamPoint globalLimit = {-1,0.}, bool applyToSecondary = true, int threadID=0, std::vector<FreeSpacePoint>* localMaxima = nullptr);
    std::vector<FreeSpacePoint> findExtremePointsFromPointWrapper(Direction,Direction,FreeSpacePoint start, ParamPoint globalLimit = {-1,0.}, bool applyToSecondary = true, int threadID=0, std::vector<FreeSpacePoint>* localMaxima = nullptr);

    std::vector<FreeSpacePoint> findExtremePointsFromPointUntilEndsWrapper(Direction,FreeSpacePoint start, const std::vector<ParamPoint>& limits,int threadID=0);
    template<Direction> std::vector<FreeSpacePoint> findExtremePointsFromPointUntilEnds(FreeSpacePoint start, const std::vector<ParamPoint>& limits,int threadID=0);

    template<Direction> std::vector<std::vector<FreeSpacePoint>> findExtremePointsFromYUntilEnds(ParamPoint start, const std::vector<ParamPoint>& limits, int threadID=0);

    //void identifyImportantUpDownYs(int x, int y);
};

#endif //CLUSTERING_FREESPACE_H
