#include "center_clustering_algs.h"
#include "io.h"
#include "curve_simplification.h"
#include "frechet_light.h"
#include <random>
#include "SWatch.h"
#include <map>
#include "free_space_visualizer.h"

#include <iostream>

//void printUsage()
//{
//	std::cout << "USAGE: ./main <dataset_file> <k> <l> [<header_size>]\n"
//	          << "The dataset_file should contain a newline separated list of filenames of curve files.\n"
//	          << "The header size (default: 1) gives the number of lines which are ignored at the beginning of each curve file.\n";
//}

void experiments(){
    Curve c1 = Curve("../data/86_1.txt",93);
    Curve c2 = Curve("../data/86_2.txt",93);
    Curve c4 = Curve("../data/86_4.txt",93);
    Curve c6 = Curve("../data/86_6.txt",93);

    //groundtruths are 1-indexed, as they are from a mathlab file
    FrameLabeling gt86_01 = {{walk,500},{transition,600},
                                                 {jump,1100},{transition,1200},
                                                 {walk,1930},{transition,2030},
                                                 {punch,2450},{transition,2550},
                                                 {walk,3150},{transition,3250},
                                                 {leg_kick,4015},{transition,4115},
                                                 {punch,4579}};
    FrameLabeling gt86_02 = {{walk,980},{transition,1100},
                                                 {squat,1850},{transition,1950},
                                                 {run,2580},{transition,2780},
                                                 {stand,3080},{transition,3180},
                                                 {arm_up,4670},{transition,4750},
                                                 {walk,5800},{transition,6000},
                                                 {jump,7250},{transition,7420},
                                                 {drink,8725},{transition,8880},
                                                 {punch,9570},{transition,9660},
                                                 {walk,10617}};
    FrameLabeling gt86_04 = {{walk,970},{transition,1100},
                                                 {stretch,2180},{transition,2295},
                                                 {punch,3368},{transition,3500},
                                                 {stand,4000},{transition,4130},
                                                 {walk,5050},{transition,5090},
                                                 {slap,5830},{transition,5930},
                                                 {turn,6600},{transition,6750},
                                                 {drink,8000},{transition,8150},
                                                 {punch,9080},{transition,9200},
                                                 {walk,10078}};
    // 86_6 Labels DONT LINE UP WITH ABOVE LABELS, ONLY CLUSTER 86_6 ALONE!!!!!
    FrameLabeling gt86_06 = {{walk,1050},{transition,1150},
                            {stand,1540},{transition,1628},
                            {run,2550},{transition,2650},
                            {stand,3100},{transition,3220},
                            {leg_kick,3900},{transition,3980},
                            {punch,4570},{transition,4670},
                            {slap,5450},{transition,5530},
                            {punch,6170},{transition,6250},
                            {drink,6960},{transition,7050}, //knee_shot
                            {jump,7930},{transition,8060}, //cheer
                            {arm_up,8870},{transition,8950}, //raise_arms
                            {walk,9939},{transition,9939}};

    double delta = 1.25;

    CurveClusterer cc(25,true);
    Curves allCurves = {c1,c2,c4,c6};
    std::vector<FrameLabeling> allLabelings = {gt86_01,gt86_02,gt86_04,gt86_06};

    for(int i=0;i<4;i++) {
        Curves curves = {allCurves[i]};
        std::vector<FrameLabeling> labelings = {allLabelings[i]};
        cc.initCurves(curves, delta, labelings);

        //custom filter which accesses ccs curves (not very cool tbh)
        int length = 0;
        for (const Curve &c: cc.simplifiedCurves) {
            length += c.size();
        }
        int complexity = (int) ((length) / (25 * cc.simplifiedCurves.size()));
        double guarantee = 1.0 + 1.0 + 2 * (7.0 / 3.0);
        auto filter = [=](const Candidate &a) {
            return cc.simplifiedCurves[a.getIndex()].subcurve_length(a.getStart(), a.getEnd()) >
                   2 * guarantee * delta && a.getEnd().id - a.getStart().id > complexity / 3;
        };

        //cover algorithm
        auto result = cc.greedyCover(complexity, 10, filter);

        //visulaization
        std::cout << "Cutoff length: " << 2 * guarantee * delta << " & " << complexity / 3 << std::endl;
        for (auto r: result) {
            std::cout << "(" << r.getEnd().id - r.getStart().id << ","
                      << cc.simplifiedCurves[r.getIndex()].subcurve_length(r.getStart(), r.getEnd()) << ")"
                      << std::endl;
        }

        ClusteringVisulaizer cv;
        cv.showClusteringStretched(cc.simplifiedCurves, cc.simplifiedGTs, result);
    }
}

int main(int argc, char* argv[])
{
#ifdef HASVISUAL
    std::cout << "COMPILED WITH OPENCV\n";
#else
    std::cout << "NOT COMPILED WITH OPENCV\n";
#endif

    experiments();

    return 0;
}

