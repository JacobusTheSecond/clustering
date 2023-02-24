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

void localData(){
    //get curves
    Curve hm36c1 = Curve("/Users/styx/data/curveclustering/anchoredlegseq3d/S11_eating_1.txt",33);
    Curve hm36c2 = Curve("/Users/styx/data/curveclustering/anchoredlegseq3d/S11_eating_2.txt",33);
    Curve hm36c3 = Curve("/Users/styx/data/curveclustering/anchoredlegseq3d/S1_eating_1.txt",33);
    Curve hm36c4 = Curve("/Users/styx/data/curveclustering/anchoredlegseq3d/S1_eating_2.txt",33);
    Curve hm36c5 = Curve("/Users/styx/data/curveclustering/anchoredlegseq3d/S5_eating_1.txt",33);
    Curve hm36c6 = Curve("/Users/styx/data/curveclustering/anchoredlegseq3d/S5_eating_2.txt",33);

    FrameLabeling hm36gt1 = io::readHM36GroundTruth("/Users/styx/data/curveclustering/label11/label11_S11_eating_1.txt",11);
    FrameLabeling hm36gt2 = io::readHM36GroundTruth("/Users/styx/data/curveclustering/label11/label11_S11_eating_2.txt",11);
    FrameLabeling hm36gt3 = io::readHM36GroundTruth("/Users/styx/data/curveclustering/label11/label11_S1_eating_1.txt",11);
    FrameLabeling hm36gt4 = io::readHM36GroundTruth("/Users/styx/data/curveclustering/label11/label11_S1_eating_2.txt",11);
    FrameLabeling hm36gt5 = io::readHM36GroundTruth("/Users/styx/data/curveclustering/label11/label11_S5_eating_1.txt",11);
    FrameLabeling hm36gt6 = io::readHM36GroundTruth("/Users/styx/data/curveclustering/label11/label11_S5_eating_2.txt",11);

    double delta = 0.04;

    CurveClusterer cc(25,true);
    Curves hm36cs = {hm36c1,hm36c2,hm36c3,hm36c4,hm36c5,hm36c6};
    std::vector<FrameLabeling> hm36gts = {hm36gt1,hm36gt2,hm36gt3,hm36gt4,hm36gt5,hm36gt6};
    cc.initCurves(hm36cs,delta,hm36gts);

    //custom filter which accesses ccs curves (not very cool tbh)
    int length = 0;
    for(const Curve& c : cc.simplifiedCurves){
        length += c.size();
    }
    int complexity = (int) ((length) / (25 * cc.simplifiedCurves.size()));
    double guarantee = 1.0+1.0+2*(7.0/3.0);
    auto filter = [=](const Candidate& a){return cc.simplifiedCurves[a.getIndex()].subcurve_length(a.getStart(),a.getEnd())>2*guarantee*delta && a.getEnd().id - a.getStart().id > complexity/3;};

    //cover algorithm
    auto result = cc.greedyCover(complexity,10,filter);

    //visulaization
    std::cout << "Cutoff length: " << 2*guarantee*delta  << " & " << complexity/3 << std::endl;
    for (auto r : result){
        std::cout << "(" << r.getEnd().id - r.getStart().id << "," << cc.simplifiedCurves[r.getIndex()].subcurve_length(r.getStart(),r.getEnd()) << ")"<<std::endl;
    }

    ClusteringVisulaizer cv;
    cv.showClustering(cc.simplifiedCurves,cc.simplifiedGTs,result);
}

int main(int argc, char* argv[])
{
#ifdef HASVISUAL
    std::cout << "COMPILED WITH OPENCV\n";
#else
    std::cout << "NOT COMPILED WITH OPENCV\n";
#endif

    localData();

    return 0;
}

//OLD STUFF THAT MIGHT STILL BE NEEDED



/**
    std::vector<std::pair<Label,ParamPoint>> shm36gt1,shm36gt2,shm36gt3,shm36gt4,shm36gt5,shm36gt6;

    Curve simpHM36C1 = good_simplification(hm36c1,delta,hm36gt1,&shm36gt1);
    Curve simpHM36C2 = good_simplification(hm36c2,delta,hm36gt2,&shm36gt2,25);
    Curve simpHM36C3 = good_simplification(hm36c3,delta,hm36gt3,&shm36gt3,25);
    Curve simpHM36C4 = good_simplification(hm36c4,delta,hm36gt4,&shm36gt4,25);
    Curve simpHM36C5 = good_simplification(hm36c5,delta,hm36gt5,&shm36gt5,25);
    Curve simpHM36C6 = good_simplification(hm36c6,delta,hm36gt6,&shm36gt6,25);

    Curves curvesHM36 = {simpHM36C1};//,simpHM36C2,simpHM36C3,simpHM36C4,simpHM36C5,simpHM36C6};
    auto simplifiedHM36gts = {shm36gt1};//,shm36gt2,shm36gt3,shm36gt4,shm36gt5,shm36gt6};



    //Curve unsimpC = Curve(simpC1,times);
    //io::exportCurve("/Users/styx/data/curveclustering/results/testcurve.txt",unsimpC);
    //io::exportCurve("/Users/styx/data/curveclustering/results/testcurveorigin.txt",c1);
    Curves curvesCMU = {simpC1,simpC2,simpC4};//,simpC6};
    std::vector<std::vector<std::pair<Label,ParamPoint>>> simplifiedGTs = {simplifiedGT86_01,simplifiedGT86_02,simplifiedGT86_04};

    double guaranteeCMU = 1.0+1.0+2*(7.0/3.0);
    int length = 0;
    for(const Curve& c : curvesHM36){
        length += c.size();
    }
    int complexity = (int) ((length) / (25 * curvesHM36.size()));
    auto filter = [=](const Candidate& a){return curvesHM36[a.getIndex()].subcurve_length(a.getStart(),a.getEnd())>2*guaranteeCMU*deltaprime && a.getEnd().id - a.getStart().id > complexity/3;};
    //auto filter2 = [](const std::pair<int,Candidate>&){return true;};
    auto result = greedyCoverUnsanitizedOutput(curvesHM36,deltaprime, complexity, 1, true, filter);

    std::cout << "Cutoff length: " << 2*guaranteeCMU*deltaprime  << " & " << complexity/3 << std::endl;
    for (auto r : result){
        std::cout << "(" << r.getEnd().id - r.getStart().id << "," << curvesHM36[r.getIndex()].subcurve_length(r.getStart(),r.getEnd()) << ")"<<std::endl;
    }

    ClusteringVisulaizer cv;
    cv.showClustering(curvesHM36,simplifiedHM36gts,result);
    return 0;
    */


//std::cout << "Hello world!\n";
//FrechetLight fl;
/*
    Curve c1 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S1_directions_1.txt",96);
    Curve c2 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S1_directions_2.txt",96);
    Curve c3 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S5_directions_1.txt",96);
    Curve c4 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S5_directions_2.txt",96);
    Curve c5 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S7_directions_1.txt",96);
    Curve c6 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S7_directions_2.txt",96);
    Curve c7 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S8_directions_1.txt",96);
    Curve c8 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S8_directions_2.txt",96);
    */

//Curve c1 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_directions_1.txt",96);
//Curve c2 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_directions_2.txt",96);
//Curve c3 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_discussion_1.txt",96);
//Curve c4 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_discussion_2.txt",96);
//Curve c5 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_eating_1.txt",96);
//Curve c6 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_eating_2.txt",96);
//Curve c7 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_greeting_1.txt",96);
//Curve c8 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_greeting_2.txt",96);

//andreas weber heuristik
//visualize candidates DONE
//normalize direction DONE
/**
    Curve c1 = Curve("../data/86_1.txt",93);
    Curve c2 = Curve("../data/86_2.txt",93);
    Curve c4 = Curve("../data/86_4.txt",93);
    Curve c6 = Curve("../data/86_6.txt",93);

    //groundtruths are 1-indexed, as they are from a mathlab file
    std::vector<std::pair<Label,int>> gt86_01 = {{walk,500},{transition,600},
                                                 {jump,1100},{transition,1200},
                                                 {walk,1930},{transition,2030},
                                                 {punch,2450},{transition,2550},
                                                 {walk,3150},{transition,3250},
                                                 {leg_kick,4015},{transition,4115},
                                                 {punch,4579}};
    std::vector<std::pair<Label,int>> gt86_02 = {{walk,980},{transition,1100},
                                                 {squat,1850},{transition,1950},
                                                 {run,2580},{transition,2780},
                                                 {stand,3080},{transition,3180},
                                                 {arm_up,4670},{transition,4750},
                                                 {walk,5800},{transition,6000},
                                                 {jump,7250},{transition,7420},
                                                 {drink,8725},{transition,8880},
                                                 {punch,9570},{transition,9660},
                                                 {walk,10617}};
    std::vector<std::pair<Label,int>> gt86_04 = {{walk,970},{transition,1100},
                                                 {stretch,2180},{transition,2295},
                                                 {punch,3368},{transition,3500},
                                                 {stand,4000},{transition,4130},
                                                 {walk,5050},{transition,5090},
                                                 {slap,5830},{transition,5930},
                                                 {turn,6600},{transition,6750},
                                                 {drink,8000},{transition,8150},
                                                 {punch,9080},{transition,9200},
                                                 {walk,10078}};

    //std::cout << "Input complexities: " << c1.size() << " " << c2.size() << std::endl;

    double delta = 0.04;//1.0;//1.0;
    double deltaprime = 0.04;//1.25;//1.25;

    std::vector<std::pair<Label,ParamPoint>> simplifiedGT86_01,simplifiedGT86_02,simplifiedGT86_04;

    Curve simpC1 = good_simplification(c1,delta,gt86_01,&simplifiedGT86_01);
    Curve simpC2 = good_simplification(c2,delta,gt86_02,&simplifiedGT86_02);
    Curve simpC4 = good_simplification(c4,delta,gt86_04,&simplifiedGT86_04);
    Curve simpC6 = good_simplification(c6,delta);
*/



/*
for(int i=0;i<result.size();++i){
    io::exportCurve("/Users/styx/data/curveclustering/results/resultcurve" + std::to_string(i) + ".txt",result[i]);
}*/

/*
    std::vector<std::vector<int>> sizes;


    for(double sd=0.25;sd<5;sd+=0.25) {
        std::cout << "START";
        sizes.emplace_back();
        Curve sc1 = good_simplification(c1,sd);
        Curve sc2 = good_simplification(c2,sd);
        Curves scs = {sc1,sc2};
        int l = 0;
        for(const Curve& c : scs){
            l += c.size();
        }
        stdc::SWatch watch;
        watch.start();
        Curves result = greedyCover(scs, 2.25, (int) ((l) / (25 * scs.size())), 2);
        watch.stop();
        sizes.back().push_back(std::chrono::duration_cast<std::chrono::milliseconds>(watch.elapsed()).count());
        watch.reset();
        std::cout << "--------\n";
        std::cout << sizes.back().back() << std::endl;
        //for (double d = 2.25; d <= 5; d += 0.25) {
        //    std::cout << "\nSIMPLIFICATION DELTA = " << sd << " FREESPACE DELTA = " << d << "\n";
        //    Curves result = greedyCoverAlreadySimplified(scs, d, (int) ((l) / (25 * scs.size())), 25);
        //    sizes.back().push_back(result.size());
        //}
    }
    std::cout << "\n\nResuls:\n";
    double sd = 0.25;
    for(const auto& subsizes : sizes) {
        std::cout << sd << "  : ";
        for (int size: subsizes) {
            std::cout << size << " ";
        }
        std::cout << std::endl;
        sd += 0.25;
    }
*/
/**
 *
Resuls:
0.25  : 33 33 35 28 20 17 13 10 9 6 6 5 5 6 4 4 4 4 3 3
0.5   : 35 34 33 29 21 16 12 9 8 6 6 6 5 5 4 5 4 3 3 3
0.75  : 34 31 30 26 21 15 11 9 8 6 6 5 5 5 5 4 4 3 3 3
1.0   : 35 34 33 29 22 15 12 9 7 6 6 5 5 4 4 5 4 3 3 3
1.25  : 34 32 33 27 22 15 13 9 7 7 5 5 5 4 4 4 4 4 3 3
1.5   : 33 33 29 28 22 14 11 8 8 7 6 5 5 4 4 4 4 4 3 3
1.75  : 32 32 31 27 21 14 11 9 7 7 6 5 5 4 4 4 4 3 3 3
2.0   : 39 38 31 29 20 15 11 11 8 6 7 5 5 4 4 5 4 4 3 3
2.25  : 34 33 31 26 21 15 11 11 7 7 6 5 5 4 4 5 4 3 3 3
2.5   : 33 30 30 23 19 13 11 9 8 7 6 5 5 5 4 5 4 4 3 3
2.75  : 37 39 35 26 20 15 13 10 8 7 6 5 5 5 4 4 3 3 3 3
3.0   : 37 38 32 23 19 15 10 10 8 9 6 6 5 4 4 4 4 3 3 3
3.25  : 35 37 31 21 19 14 11 9 9 8 6 5 5 4 4 4 4 4 3 3
3.5   : 35 34 29 22 16 14 12 8 8 7 5 5 5 4 4 4 3 3 2 3
3.75  : 33 33 29 21 17 13 10 8 7 6 5 5 5 4 4 4 2 2 2 3
4.0   : 39 37 35 24 19 13 12 9 8 7 6 5 5 4 4 4 3 3 3 2
4.25  : 37 38 34 26 18 15 13 11 8 7 6 5 4 4 5 4 3 3 3 2
4.5   : 37 37 33 26 19 16 13 9 7 8 6 6 5 5 4 3 2 3 3 3
4.75  : 35 33 30 23 17 18 12 9 8 6 6 6 4 4 3 2 2 3 3 2
 * */


//0.25  : 64775
//0.5  : 37201
//0.75  : 15123
//1  : 22767
//1.25  : 17945
//1.5  : 10853
//1.75  : 7876
//2  : 5176
//2.25  : 6202
//2.5  : 4330
/*
    for(int i=0;i<result.size();++i){
        io::exportCurve("/Users/styx/data/curveclustering/results/resultcurve" + std::to_string(i) + ".txt",result[i]);
    }
*/
/*

Curve c1 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S1_eating_1.txt",96);
Curve c2 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S1_eating_2.txt",96);
Curve c3 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S5_eating_1.txt",96);
Curve c4 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S5_eating_2.txt",96);
Curve c5 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S6_eating_1.txt",96);
Curve c6 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S6_eating_2.txt",96);
Curve c7 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S7_eating_1.txt",96);
Curve c8 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S7_eating_2.txt",96);
Curve c9 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S8_eating_1.txt",96);
Curve c10 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S8_eating_2.txt",96);
Curve c11 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S9_eating_1.txt",96);
Curve c12 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S9_eating_2.txt",96);
Curve c13 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_eating_1.txt",96);
Curve c14 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_eating_2.txt",96);
//std::vector<int> times;
double delta = 0.1;
double guarantee = 1.0+1.0+2*(7.0/3.0);
Curve simp1 = good_simplification(c1,delta);
Curve simp2 = good_simplification(c2,delta,&times);
Curve simp3 = good_simplification(c3,delta);
Curve simp4 = good_simplification(c4,delta);
Curve simp5 = good_simplification(c5,delta);
Curve simp6 = good_simplification(c6,delta);
Curve simp7 = good_simplification(c7,delta);
Curve simp8 = good_simplification(c8,delta);
Curve simp9 = good_simplification(c9,delta);
Curve simp10 = good_simplification(c10,delta);
Curve simp11 = good_simplification(c11,delta);
Curve simp12 = good_simplification(c12,delta);
Curve simp13 = good_simplification(c13,delta);
Curve simp14 = good_simplification(c14,delta);
Curves curves = {simp1,simp2,simp3,simp4,simp5,simp6,simp7,simp8,simp9,simp10,simp11,simp12,simp13,simp14};
//CandidateSet cs = CandidateSet(curves,guarantee*delta);
//cs.computeCandidates(20);

//Gudmunsen Vahrenhold GPU frechet distance

//visualize
Curve unsimp2 = Curve(simp2,times);
io::exportCurve("/Users/styx/data/curveclustering/results/testcurve.txt",unsimp2);

greedyCoverAlreadySimplified(curves,delta,10);
 */
/*
    std::vector<std::pair<int,int>> indices;

    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    mpz_t r;
    mpz_init(r);

    int generations = 10000;

    for(int i=0;i<generations;i++) {
        std::cout << "Round " << i << "\n";
        std::cout << "Total accumulated weight so far: " << cs.tW << "\n";
        indices.reserve(100);
        for (int n = 0; n < 100; n++) {
            //TODO: dont draw multiple elements
            mpz_urandomm(r,rs,cs.tW);
            indices.push_back(cs.draw(r));
        }
        std::pair<int,ParamPoint> t = cs.findNonCovered(indices);
        if(t.first == -1){
            for(auto ind : indices){
                std::cout << "(" << ind.first << "," << ind.second << ") ";
            }
            std::cout << "\nFound one\n";
            cs.showCovering(indices);
            continue;
        }
        std::cout <<"Non-covered point: " << t.first << " " << t.second.id << " " << t.second.t << "\n";
        cs.reweight(t);
        indices.clear();
    }
    */

//Some old code snippets that i maybe want to reuse

//cs.showCandidates();
/*simp1.canonicalPointCandidates(simp1,simp1,0.5);
Curves candidates = simp2.canonicalPointCandidates(simp1,1.5);
std::cout << candidates.size() << "\n";
std::cout << candidates[0].size() << "\n";
Clustering clusters = computeCenterClustering(candidates,100,simp2.size()/2, ClusterAlg::Gonzalez, CenterAlg::FSA,1);
for (auto cluster : clusters)
    std::cout << cluster.cost << " ";
io::exportClustering("/Users/styx/data/curveclustering/results/testclusters","/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_directions_2.txt",clusters,candidates);
*/

	/*int header_size = 1;

	if (argc == 5) {
		header_size = std::stoi(argv[4]);
	}
	else if (argc != 4) {
		printUsage();
		ERROR("Wrong number of arguments.");
	}

	std::string base_path = argv[1];
	int k = std::stoi(argv[2]);
	int l = std::stoi(argv[3]);

	auto curves = io::readCurves(base_path, header_size);
	auto clustering = computeCenterClustering(curves, k, l, ClusterAlg::Gonzalez, CenterAlg::FSA);

	io::exportClustering("out.clustering", base_path, clustering, curves);
	io::exportCentersGPX("centers.gpx", clustering);*/
