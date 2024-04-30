//
// Created by Jacobus Conradi on 17.04.23.
//

#ifndef CLUSTERING_SPARSEFREESPACE_H
#define CLUSTERING_SPARSEFREESPACE_H

#include "geometry_types.h"
#include "Curve.h"


struct CellPoint{
    CellPoint(){};
    CellPoint(distance_t _x, distance_t _y):x(_x),y(_y){};
    distance_t x{};
    distance_t y{};
};


inline std::ostream &operator<<(std::ostream &os, CellPoint const &m) {
    return os << "( " << m.x << " , " << m.y << " )";
}

class Cell{
public:
    std::pair<CellPoint,CellPoint> leftPair, rightPair, topPair, bottomPair;
    Interval left,top,right,bottom;

    //this isnt exactly clean
    Interval& toAbove(int threadID = 0){return toAboves[threadID];}
    Interval& toRight(int threadID = 0){return toRights[threadID];}
    Interval& toBottom(int threadID = 0){return toBottoms[threadID];}
    Interval& toLeft(int threadID = 0){return toLefts[threadID];}
    std::vector<Interval> toAboves;
    std::vector<Interval> toRights;
    std::vector<Interval> toLefts;
    std::vector<Interval> toBottoms;

    virtual CellPoint leftMostAt(double y, CellPoint* outer = nullptr) const = 0;
    virtual CellPoint rightMostAt(double y, CellPoint* outer = nullptr) const = 0;
    virtual CellPoint topMostAt(double x, CellPoint* outer = nullptr) const = 0;
    virtual CellPoint bottomMostAt(double x, CellPoint* outer = nullptr) const = 0;
    virtual bool is_empty()=0;
    Cell(int tc):toAboves(tc),toRights(tc),toLefts(tc),toBottoms(tc){};
    virtual ~Cell(){};
    void dump(int threadID = 0){
        std::cout << "left: " << left  << std::endl;
        std::cout << "top: " << top << std::endl;
        std::cout << "right: " << right << std::endl;
        std::cout << "bottom: " << bottom << std::endl;

        std::cout << "leftPair: " << leftPair.first << " , " << leftPair.second << std::endl;
        std::cout << "topPair: " << topPair.first << " , " << topPair.second << std::endl;
        std::cout << "rightPair: " << rightPair.first << " , " << rightPair.second << std::endl;
        std::cout << "bottomPair: " << bottomPair.first << " , " << bottomPair.second << std::endl;

        std::cout << "toAbove:" << toAbove(threadID) << std::endl;
        std::cout << "toRight:" << toRight(threadID) << std::endl;
    };
};

class SparseCell : public Cell{
    //void updateDelta(distance_t delta);
private:
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

    //Point& a,b,c,d;

    bool isEmpty;
public:
    double delta=0;

    CellPoint leftMostAt(double y, CellPoint* outer = nullptr) const;
    CellPoint rightMostAt(double y, CellPoint* outer = nullptr) const;
    CellPoint topMostAt(double x, CellPoint* outer = nullptr) const;
    CellPoint bottomMostAt(double x, CellPoint* outer = nullptr) const;
    //SparseCell(Point &a, Point &b, Point &c, Point &d,int threadcount = 1);
    SparseCell(Point &a, Point &b, Point &c, Point &d,distance_t delta, int threadcount = 1);
    bool is_empty(){return isEmpty;};
};

class MinkowskiCell : public Cell{

private:
    Point a,b,c,d;
    distance_t ra,rb,rc,rd;
    bool isEmpty;

    distance_t acac,acab,accd,abab,abcd;
    distance_t caca,cacd,caab,cdcd,cdab;
public:
    CellPoint leftMostAt(double y, CellPoint* outer = nullptr) const;
    CellPoint rightMostAt(double y, CellPoint* outer = nullptr) const;
    CellPoint topMostAt(double x, CellPoint* outer = nullptr) const;
    CellPoint bottomMostAt(double x, CellPoint* outer = nullptr) const;

    MinkowskiCell(Point &a, Point &b, Point &c, Point &d, distance_t ra, distance_t rb, distance_t rc, distance_t rd, int tc = 1);
    bool is_empty(){return isEmpty;};
};

template <typename T> class SparseGridCell{
public:
    T data;
    int x;
    int y;
    int upId = -1;
    int downId = -1;
    int leftId = -1;
    int rightId = -1;
    //SparseGridCell<T>* up=nullptr;
    //SparseGridCell<T>* down=nullptr;
    //SparseGridCell<T>* left=nullptr;
    //SparseGridCell<T>* right=nullptr;
    SparseGridCell(T _data, int _x, int _y):data(std::move(_data)),x(_x),y(_y){}

    void dump(int threadID = 0){
        std::cout << "x:" << x << std::endl;
        std::cout << "y:" << y << std::endl;
        std::cout << "upId:" << upId << std::endl;
        std::cout << "downId:" << downId << std::endl;
        std::cout << "leftId:" << leftId << std::endl;
        std::cout << "rightId:" << rightId << std::endl;
        std::cout << "----data----" << std::endl;
        data->dump(threadID);
        std::cout << "------------" << std::endl;
    };
};

class SparseFreespace :std::vector<std::vector<SparseGridCell<std::unique_ptr<Cell>>>> {
private:
    using Parent = std::vector<std::vector<SparseGridCell<std::unique_ptr<Cell>>>>;
    const unsigned int nx,ny;
public:
    std::vector<CPoint> upStarts, downStarts, upEnds, downEnds;
    const CurveID BID,TID;
    const distance_t delta;
    SparseFreespace(Curve& B, Curve& T, distance_t _delta, int threadcount = 1, CurveID BID = -1, CurveID TID = -1);
    std::vector<SparseGridCell<std::unique_ptr<Cell>>>& row(PointID y){return operator[](y);}
    SparseGridCell<std::unique_ptr<Cell>>* cell(PointID y, PointID xidx){return &(operator[](y).operator[](xidx));}
    int xSize();
    int ySize();
    void identifyStarts();
    void identifyEnds();
    void identifyStartsAndEnds(){identifyStarts();/*identifyEnds();*/}
    //SparseFreespace(SparseFreespace& sf)=default;
};


#endif //CLUSTERING_SPARSEFREESPACE_H
