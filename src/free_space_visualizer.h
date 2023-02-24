//
// Created by Jacobus Conradi on 04.02.22.
//

#ifndef CLUSTERING_FREE_SPACE_VISUALIZER_H
#define CLUSTERING_FREE_SPACE_VISUALIZER_H

#include "freespace.h"
#include "candidate.h"
#include "io.h"


class FreeSpaceVisualizer{
private:
    FreeSpace& freespace;
public:
    explicit FreeSpaceVisualizer(FreeSpace &freespace);
    void show(bool withPoints=false);
    void show(FreeSpacePoint& a, FreeSpacePoint& b);
    void show(FreeSpacePoint& a, FreeSpacePoint& b, ParamPoint x1, ParamPoint x2);
    void showCandidate(Candidate & c, int curveId);
    void showCandidates(std::vector<Candidate> candidates, int curveId);
};

class FreeSpacesVisualizer{
private:
    std::vector<std::vector<FreeSpace>> & freespaces;
public:
    explicit FreeSpacesVisualizer(std::vector<std::vector<FreeSpace>> & fs) : freespaces(fs){};
    void show();
    void showCandidates(CandidateSet& candidates, std::vector<std::pair<int,int>> & indices);
    void showCandidates(std::vector<Candidate> candidates);
};



class ClusteringVisulaizer{
public:
    void showClustering(Curves c, std::vector<std::vector<std::pair<Label,ParamPoint>>> groundthruth, std::vector<Candidate> candidates);
    void showClusteringStretched(Curves c, std::vector<std::vector<std::pair<Label,ParamPoint>>> groundthruth, std::vector<Candidate> candidates);
};

#endif //CLUSTERING_FREE_SPACE_VISUALIZER_H
