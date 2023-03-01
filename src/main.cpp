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


void experiment1(const std::string& root_dir) {
    Curve c1 = Curve("../data/86_1.txt", 93);

    //groundtruths are 1-indexed, as they are from a mathlab file
    FrameLabeling gt86_01 = {{walk,       500},
                             {transition, 600},
                             {jump,       1100},
                             {transition, 1200},
                             {walk,       1930},
                             {transition, 2030},
                             {punch,      2450},
                             {transition, 2550},
                             {walk,       3150},
                             {transition, 3250},
                             {leg_kick,   4015},
                             {transition, 4115},
                             {punch,      4579}};

    double delta = 1.25;

    CurveClusterer cc(-1, false);
    Curves curves = {c1};
    std::vector<FrameLabeling> labelings = {gt86_01};

    cc.initCurves(curves, delta, labelings);
    cc.setWithSort(false);

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
    cv.withAutocoloring = false;
    cv.showClusteringStretched(cc.simplifiedCurves, cc.simplifiedGTs, result);

    for (int i = 0; i < std::min(20, (int) result.size()); ++i) {
        Candidate c = result[i];
        for (int j = 0; j < c.visualMatchings.size(); ++j) {
            auto m = c.visualMatchings[j];
            io::exportSubcurve(
                    root_dir + "/matching" + std::to_string(i) + "/interval" +
                    std::to_string(j) + ".txt", cc.simplifiedCurves[m.getIndex()], m.getStart(), m.getEnd(), 100);
        }
    }
}

void experiments() {
    Curve c1 = Curve("../data/86_1.txt", 93);
    Curve c2 = Curve("../data/86_2.txt", 93);
    Curve c3 = Curve("../data/86_3.txt", 93);
    Curve c4 = Curve("../data/86_4.txt", 93);
    Curve c5 = Curve("../data/86_5.txt", 93);
    Curve c6 = Curve("../data/86_6.txt", 93);
    Curve c7 = Curve("../data/86_7.txt", 93);
    Curve c8 = Curve("../data/86_8.txt", 93);
    Curve c9 = Curve("../data/86_9.txt", 93);
    Curve c10 = Curve("../data/86_10.txt", 93);
    Curve c11 = Curve("../data/86_11.txt", 93);
    Curve c12 = Curve("../data/86_12.txt", 93);
    Curve c13 = Curve("../data/86_13.txt", 93);
    Curve c14 = Curve("../data/86_14.txt", 93);

    //groundtruths are 1-indexed, as they are from a mathlab file
    FrameLabeling gt86_01 = {{walk,       500},
                             {transition, 600},
                             {jump,       1100},
                             {transition, 1200},
                             {walk,       1930},
                             {transition, 2030},
                             {punch,      2450},
                             {transition, 2550},
                             {walk,       3150},
                             {transition, 3250},
                             {leg_kick,   4015},
                             {transition, 4115},
                             {punch,      4579}};
    FrameLabeling gt86_02 = {{walk,       980},
                             {transition, 1100},
                             {squat,      1850},
                             {transition, 1950},
                             {run,        2580},
                             {transition, 2780},
                             {stand,      3080},
                             {transition, 3180},
                             {arm_up,     4670},
                             {transition, 4750},
                             {walk,       5800},
                             {transition, 6000},
                             {jump,       7250},
                             {transition, 7420},
                             {drink,      8725},
                             {transition, 8880},
                             {punch,      9570},
                             {transition, 9660},
                             {walk,       10617}};
    FrameLabeling gt86_03 = {{walk,       900},
                             {transition, 1050},
                             {run,        1800},
                             {transition, 1970},
                             {jump,       2351},
                             {transition, 2500},
                             {walk,       3414},
                             {transition, 3570},
                             {leg_kick,   4612},
                             {transition, 4730},
                             {drink,      5350},//jump_un_left_leg
                             {transition, 5500},
                             {punch,      6150},//jump_on_right_leg
                             {transition, 6280},
                             {arm_up,     6900},//arm_circle
                             {transition, 7070},
                             {walk,       8401}};
    FrameLabeling gt86_04 = {{walk,       970},
                             {transition, 1100},
                             {stretch,    2180},
                             {transition, 2295},
                             {punch,      3368},
                             {transition, 3500},
                             {stand,      4000},
                             {transition, 4130},
                             {walk,       5050},
                             {transition, 5090},
                             {slap,       5830},
                             {transition, 5930},
                             {turn,       6600},
                             {transition, 6750},
                             {drink,      8000},
                             {transition, 8150},
                             {punch,      9080},
                             {transition, 9200},
                             {walk,       10078}};
    FrameLabeling gt86_05 = {{walk,       760},
                             {transition, 850},
                             {jump,       1480},
                             {transition, 1675},
                             {drink,      2270}, // jump_jack
                             {transition, 2400},
                             {jump,       3881},
                             {transition, 4000},
                             {walk,       4480},
                             {transition, 4590},
                             {punch,      5150},
                             {transition, 5250},
                             {leg_kick,   5800},//cheer
                             {transition, 5940},
                             {arm_up,     6510},
                             {transition, 6630},
                             {slap,       7300},
                             {transition, 7400},
                             {walk,       8340},};
    FrameLabeling gt86_06 = {{walk,       1050},
                             {transition, 1150},
                             {stand,      1540},
                             {transition, 1628},
                             {run,        2550},
                             {transition, 2650},
                             {stand,      3100},
                             {transition, 3220},
                             {leg_kick,   3900},
                             {transition, 3980},
                             {punch,      4570},
                             {transition, 4670},
                             {slap,       5450},
                             {transition, 5530},
                             {punch,      6170},
                             {transition, 6250},
                             {drink,      6960},
                             {transition, 7050}, //knee_shot
                             {jump,       7930},
                             {transition, 8060}, //cheer
                             {arm_up,     8870},
                             {transition, 8950}, //raise_arms
                             {walk,       9939},
                             {transition, 9939}};
    FrameLabeling gt86_07 = {{walk,       1040},
                             {transition, 1140},
                             {drink,      1870},//rotate_body
                             {transition, 1950},
                             {arm_up,     2525},//rotate_arms
                             {transition, 2590},
                             {drink,      3640},//rotate_body
                             {transition, 3720},
                             {arm_up,     4420},//rotate_arms
                             {transition, 4500},
                             {jump,       5040},
                             {transition, 5200},
                             {punch,      5700}, //jump_2
                             {transition, 5840},
                             {walk,       6910},
                             {transition, 7040},
                             {run,        7750},
                             {transition, 7850},
                             {walk,       8702}};
    FrameLabeling gt86_08 = {{walk,       900},
                             {transition, 1130},
                             {drink,      1845},//sit_up
                             {transition, 1920},
                             {slap,       2645},//body_rotation
                             {transition, 2740},
                             {arm_up,     3300},//arm_rotation
                             {transition, 3360},
                             {stand,      3915},
                             {transition, 4010},
                             {leg_kick,   4700},
                             {transition, 4850},
                             {run,        5580},
                             {transition, 5775},
                             {stand,      6330},
                             {transition, 6420},
                             {stretch,    7150},
                             {transition, 7180},
                             {punch,      8120},
                             {transition, 8230},
                             {walk,       9206}};
    FrameLabeling gt86_09 = {{walk,       880},
                             {transition, 1045},
                             {arm_up,     2060},//look_around
                             {transition, 2170},
                             {slap,       2825},//clap
                             {transition, 2880},
                             {drink,      3600},//stand_and_clap
                             {transition, 3700},
                             {walk,       4794}};
    FrameLabeling gt86_10 = {{walk,       1900},
                             {transition, 2010},
                             {arm_up,     3780},//sit
                             {transition, 3850},
                             {walk,       4930},
                             {transition, 5150},
                             {stand,      5440},//stand_in_inclined_position
                             {transition, 5680},
                             {run,        6600},
                             {transition, 6720},
                             {walk,       7583}};
    FrameLabeling gt86_11 = {{walk,       1020},
                             {transition, 1190},
                             {arm_up,     1685},//both_arms_rotation
                             {transition, 1730},
                             {drink,      2330},//right_arm_rotation
                             {transition, 2365},
                             {arm_up,     2730},//both_arms_rotation
                             {transition, 2766},
                             {slap,       3300},//left_arm_rotation
                             {transition, 3370},
                             {drink,      4020},//right_arm_rotation
                             {transition, 4050},
                             {arm_up,     4600},//both_arms_rotation
                             {transition, 4720},
                             {walk,       5674}};
    FrameLabeling gt86_12 = {{walk,       910},
                             {transition, 1120},
                             {punch,      1560},//drag
                             {transition, 1790},
                             {leg_kick,   3210},//sweep_floor
                             {transition, 3340},
                             {squat,      4115},//collect_dirt
                             {transition, 4200},
                             {run,        4935},//throw_it_away
                             {transition, 5020},
                             {drink,      5300},//get_up_and_short_walk
                             {transition, 5320},
                             {arm_up,     7500},//wash_window
                             {transition, 7680},
                             {walk,       8856}};
    FrameLabeling gt86_13 = {{walk,       960},
                             {transition, 1105},
                             {squat,      1460},//climp_up
                             {transition, 1650},
                             {drink,      2280},//stand_move_hands_a_bit
                             {transition, 2380},
                             {arm_up,     2580},//climp_down
                             {transition, 2710},
                             {stand,      3010},
                             {transition, 3150},
                             {squat,      3560},//climb_up
                             {transition, 3770},
                             {punch,      4750},//look_around
                             {transition, 4790},
                             {arm_up,     5170},//climp_down
                             {transition, 5360},
                             {walk,       6221}};
    FrameLabeling gt86_14 = {{squat,      630},//walk_lead_ball
                             {transition, 700},
                             {arm_up,     1890},//throw_ball
                             {transition, 1980},
                             {squat,      2900},//walk_lead_ball
                             {transition, 3000},
                             {drink,      4090},//lead_ball_both_hands
                             {transition, 4240},
                             {squat,      5030},//walk_lead_ball
                             {transition, 5105},
                             {arm_up,     5240},//throw_ball
                             {transition, 5320},
                             {walk,       6055}};
    double delta = 1.25;

    CurveClusterer cc(25, false);
    Curves allCurves = {c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14};
    std::vector<FrameLabeling> allLabelings = {gt86_01,gt86_02,gt86_03,gt86_04,gt86_05,gt86_06,gt86_07,gt86_08,gt86_09,gt86_10,gt86_11,gt86_12,gt86_13,gt86_14};

    for (int i = 0; i < allCurves.size(); i++) {
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

        ClusteringVisulaizer cv{true};
        cv.showClusteringStretched(cc.simplifiedCurves, cc.simplifiedGTs, result);
    }
}

int main(int argc, char *argv[]) {
#ifdef HASVISUAL
    std::cout << "COMPILED WITH OPENCV\n";
#else
    std::cout << "NOT COMPILED WITH OPENCV\n";
#endif

    experiment1("/Users/styx/data/curveclustering/results/cluster");
    experiments();

    return 0;
}

