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

class SparseCell{
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

    //trying numerical stability stuff. they are always orderd bottom to top, and left to right
    std::pair<CellPoint,CellPoint> leftPair, rightPair, topPair, bottomPair;
    Interval left,top,right,bottom;

    Interval& toAbove(int threadID = 0){return toAboves[threadID];}
    Interval& toRight(int threadID = 0){return toRights[threadID];}
    Interval& toBottom(int threadID = 0){return toBottoms[threadID];}
    Interval& toLeft(int threadID = 0){return toLefts[threadID];}

    std::vector<Interval> toAboves;
    std::vector<Interval> toRights;
    std::vector<Interval> toLefts;
    std::vector<Interval> toBottoms;

    double delta=0;

    CellPoint leftMostAt(double y, CellPoint* outer = nullptr) const;
    CellPoint rightMostAt(double y, CellPoint* outer = nullptr) const;
    CellPoint topMostAt(double x, CellPoint* outer = nullptr) const;
    CellPoint bottomMostAt(double x, CellPoint* outer = nullptr) const;

public:
//SparseCell(Point &a, Point &b, Point &c, Point &d,int threadcount = 1);
SparseCell(Point &a, Point &b, Point &c, Point &d,distance_t delta, int threadcount = 1);

    bool is_empty();
};

template <typename T> class SparseGridCell{
public:
    T data;
    int x;
    int y;
    SparseGridCell<T>* up=nullptr;
    SparseGridCell<T>* down=nullptr;
    SparseGridCell<T>* left=nullptr;
    SparseGridCell<T>* right=nullptr;
    SparseGridCell(T& _data, int _x, int _y):data(_data),x(_x),y(_y){}
};

class SparseFreespace :std::vector<std::vector<SparseGridCell<SparseCell>>> {
private:
    using Parent = std::vector<std::vector<SparseGridCell<SparseCell>>>;
    const unsigned int nx,ny;
public:
    std::vector<CPoint> upStarts, downStarts, upEnds, downEnds;
    const CurveID BID,TID;
    const distance_t delta;
    SparseFreespace(Curve& B, Curve& T, distance_t _delta, int threadcount = 1, CurveID BID = -1, CurveID TID = -1);
    std::vector<SparseGridCell<SparseCell>>& row(PointID y){return operator[](y);}
    SparseGridCell<SparseCell>* cell(PointID y, PointID xidx){return &(operator[](y).operator[](xidx));}
    int xSize();
    int ySize();
    void identifyStarts();
    void identifyEnds();
    void identifyStartsAndEnds(){identifyStarts();/*identifyEnds();*/}
    //SparseFreespace(SparseFreespace& sf)=default;
};


#endif //CLUSTERING_SPARSEFREESPACE_H
