//
// Created by Jacobus Conradi on 17.04.23.
//

#ifndef CLUSTERING_FREESPACEVISUALIZER_H
#define CLUSTERING_FREESPACEVISUALIZER_H


#include "SparseFreespace.h"
#include "Candidate.h"

class FreeSpaceVisualizer{
private:
    SparseFreespace& freespace;
public:
    explicit FreeSpaceVisualizer(SparseFreespace &freespace) : freespace(freespace){}
    void show(bool withPoints=false);
    //void show(FreeSpacePoint& a, FreeSpacePoint& b);
    //void show(FreeSpacePoint& a, FreeSpacePoint& b, ParamPoint x1, ParamPoint x2);
    //void showCandidate(Candidate & c, int curveId);
    //void showCandidates(std::vector<Candidate> candidates, int curveId);
};

class FreeSpacesVisualizer{
private:
    std::vector<std::vector<SparseFreespace>> & freespaces;
public:
    explicit FreeSpacesVisualizer(std::vector<std::vector<SparseFreespace>> & fs) : freespaces(fs){};
    void show();
    //void showCandidates(CandidateSet& candidates, std::vector<std::pair<int,int>> & indices);
    void showCandidates(std::vector<Candidate> candidates);
};

class SparseFreeSpacesVisualizer{
private:
    SparseFreeSpaces & freespaces;
public:
    explicit SparseFreeSpacesVisualizer(SparseFreeSpaces& fs) : freespaces(fs){};
    void show();
    //void showCandidates(CandidateSet& candidates, std::vector<std::pair<int,int>> & indices);
    void showCandidates(std::vector<Candidate> candidates);
};

class ClusteringVisulaizer{
public:
    bool withAutocoloring;
    void showClusteringStretched(Curves c, std::vector<std::vector<std::pair<Label,CPoint>>> groundthruth, std::vector<Candidate> candidates);
};

#endif //CLUSTERING_FREESPACEVISUALIZER_H
