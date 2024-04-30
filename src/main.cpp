#include <iostream>
#include "Curve.h"
#include "CurveSimplification.h"
#include "SparseFreespace.h"
#include "FreespaceVisualizer.h"
#include "Candidate.h"
#include "center_clustering_algs.h"
#include "io.h"

//#include <opencv2/opencv.hpp>
//#include <opencv2/highgui.hpp>
#include <pybind11/pybind11.h>


//void printUsage()
//{
//	std::cout << "USAGE: ./main <dataset_file> <k> <l> [<header_size>]\n"
//	          << "The dataset_file should contain a newline separated list of filenames of curve files.\n"
//	          << "The header size (default: 1) gives the number of lines which are ignored at the beginning of each curve file.\n";
//}
/*

void experiment1(const std::string& root_dir) {
    Curve c1 = Curve("../data/86_1.txt", 93);

    //groundtruths are 1-indexed, as they are from a mathlab file
    FrameLabeling gt86_01 = {{2,       500},
                             {1, 600},
                             {3,       1100},
                             {1, 1200},
                             {2,       1930},
                             {1, 2030},
                             {4,      2450},
                             {1, 2550},
                             {2,       3150},
                             {1, 3250},
                             {5,   4015},
                             {1, 4115},
                             {4,      4579}};

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

    //print parameters
    for(auto r : result){
        for(auto m : r.matchings){
            ParamPoint s = cc.mapSimplificationToBase(m.getIndex(),m.getStart());
            ParamPoint t = cc.mapSimplificationToBase(m.getIndex(),m.getEnd());
            std::cout << "(" << (double)s.id + s.t << "," << (double)t.id + t.t << ") ";
        }
        std::cout << std::endl;
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
*/
//enum Label {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};

void experiments() {

    //groundtruths are 1-indexed, as they are from a mathlab file
    FrameLabeling gt86_01 = {{2,       500},
                             {1, 600},
                             {3,       1100},
                             {1, 1200},
                             {2,       1930},
                             {1, 2030},
                             {4,      2450},
                             {1, 2550},
                             {2,       3150},
                             {1, 3250},
                             {5,   4015},
                             {1, 4115},
                             {4,      4579}};
    FrameLabeling gt86_02 = {{2,       980},
                             {1, 1100},
                             {6,      1850},
                             {1, 1950},
                             {7,        2580},
                             {1, 2780},
                             {8,      3080},
                             {1, 3180},
                             {9,     4670},
                             {1, 4750},
                             {2,       5800},
                             {1, 6000},
                             {3,       7250},
                             {1, 7420},
                             {10,      8725},
                             {1, 8880},
                             {4,      9570},
                             {1, 9660},
                             {2,       10617}};
    FrameLabeling gt86_03 = {{2,       900},
                             {1, 1050},
                             {7,        1800},
                             {1, 1970},
                             {3,       2351},
                             {1, 2500},
                             {2,       3414},
                             {1, 3570},
                             {5,   4612},
                             {1, 4730},
                             {10,      5350},//3_un_left_leg
                             {1, 5500},
                             {4,      6150},//3_on_right_leg
                             {1, 6280},
                             {9,     6900},//arm_circle
                             {1, 7070},
                             {2,       8401}};
    FrameLabeling gt86_04 = {{2,       970},
                             {1, 1100},
                             {11,    2180},
                             {1, 2295},
                             {4,      3368},
                             {1, 3500},
                             {8,      4000},
                             {1, 4130},
                             {2,       5050},
                             {1, 5090},
                             {12,       5830},
                             {1, 5930},
                             {13,       6600},
                             {1, 6750},
                             {10,      8000},
                             {1, 8150},
                             {4,      9080},
                             {1, 9200},
                             {2,       10078}};
    FrameLabeling gt86_05 = {{2,       760},
                             {1, 850},
                             {3,       1480},
                             {1, 1675},
                             {10,      2270}, // 3_jack
                             {1, 2400},
                             {3,       3881},
                             {1, 4000},
                             {2,       4480},
                             {1, 4590},
                             {4,      5150},
                             {1, 5250},
                             {5,   5800},//cheer
                             {1, 5940},
                             {9,     6510},
                             {1, 6630},
                             {12,       7300},
                             {1, 7400},
                             {2,       8340},};
    FrameLabeling gt86_06 = {{2,       1050},
                             {1, 1150},
                             {8,      1540},
                             {1, 1628},
                             {7,        2550},
                             {1, 2650},
                             {8,      3100},
                             {1, 3220},
                             {5,   3900},
                             {1, 3980},
                             {4,      4570},
                             {1, 4670},
                             {12,       5450},
                             {1, 5530},
                             {4,      6170},
                             {1, 6250},
                             {10,      6960},
                             {1, 7050}, //knee_shot
                             {3,       7930},
                             {1, 8060}, //cheer
                             {9,     8870},
                             {1, 8950}, //raise_arms
                             {2,       9939},
                             {1, 9939}};
    FrameLabeling gt86_07 = {{2,       1040},
                             {1, 1140},
                             {10,      1870},//rotate_body
                             {1, 1950},
                             {9,     2525},//rotate_arms
                             {1, 2590},
                             {10,      3640},//rotate_body
                             {1, 3720},
                             {9,     4420},//rotate_arms
                             {1, 4500},
                             {3,       5040},
                             {1, 5200},
                             {4,      5700}, //3_2
                             {1, 5840},
                             {2,       6910},
                             {1, 7040},
                             {7,        7750},
                             {1, 7850},
                             {2,       8702}};
    FrameLabeling gt86_08 = {{2,       900},
                             {1, 1130},
                             {10,      1845},//sit_up
                             {1, 1920},
                             {12,       2645},//body_rotation
                             {1, 2740},
                             {9,     3300},//arm_rotation
                             {1, 3360},
                             {8,      3915},
                             {1, 4010},
                             {5,   4700},
                             {1, 4850},
                             {7,        5580},
                             {1, 5775},
                             {8,      6330},
                             {1, 6420},
                             {11,    7150},
                             {1, 7180},
                             {4,      8120},
                             {1, 8230},
                             {2,       9206}};
    FrameLabeling gt86_09 = {{2,       880},
                             {1, 1045},
                             {9,     2060},//look_around
                             {1, 2170},
                             {12,       2825},//clap
                             {1, 2880},
                             {10,      3600},//8_and_clap
                             {1, 3700},
                             {2,       4794}};
    FrameLabeling gt86_10 = {{2,       1900},
                             {1, 2010},
                             {9,     3780},//sit
                             {1, 3850},
                             {2,       4930},
                             {1, 5150},
                             {8,      5440},//8_in_inclined_position
                             {1, 5680},
                             {7,        6600},
                             {1, 6720},
                             {2,       7583}};
    FrameLabeling gt86_11 = {{2,       1020},
                             {1, 1190},
                             {9,     1685},//both_arms_rotation
                             {1, 1730},
                             {10,      2330},//right_arm_rotation
                             {1, 2365},
                             {9,     2730},//both_arms_rotation
                             {1, 2766},
                             {12,       3300},//left_arm_rotation
                             {1, 3370},
                             {10,      4020},//right_arm_rotation
                             {1, 4050},
                             {9,     4600},//both_arms_rotation
                             {1, 4720},
                             {2,       5674}};
    FrameLabeling gt86_12 = {{2,       910},
                             {1, 1120},
                             {4,      1560},//drag
                             {1, 1790},
                             {5,   3210},//sweep_floor
                             {1, 3340},
                             {6,      4115},//collect_dirt
                             {1, 4200},
                             {7,        4935},//throw_it_away
                             {1, 5020},
                             {10,      5300},//get_up_and_short_2
                             {1, 5320},
                             {9,     7500},//wash_window
                             {1, 7680},
                             {2,       8856}};
    FrameLabeling gt86_13 = {{2,       960},
                             {1, 1105},
                             {6,      1460},//climp_up
                             {1, 1650},
                             {10,      2280},//8_move_hands_a_bit
                             {1, 2380},
                             {9,     2580},//climp_down
                             {1, 2710},
                             {8,      3010},
                             {1, 3150},
                             {6,      3560},//climb_up
                             {1, 3770},
                             {4,      4750},//look_around
                             {1, 4790},
                             {9,     5170},//climp_down
                             {1, 5360},
                             {2,       6221}};
    FrameLabeling gt86_14 = {{6,      630},//2_lead_ball
                             {1, 700},
                             {9,     1890},//throw_ball
                             {1, 1980},
                             {6,      2900},//2_lead_ball
                             {1, 3000},
                             {10,      4090},//lead_ball_both_hands
                             {1, 4240},
                             {6,      5030},//2_lead_ball
                             {1, 5105},
                             {9,     5240},//throw_ball
                             {1, 5320},
                             {2,       6055}};

    FrameLabeling haca_01 = {{2,239},
                             {2,574},
                             {3,842},
                             {3,1130},
                             {2,1513},
                             {2,1896},
                             {4,2279},
                             {4,2615},
                             {2,2902},
                             {2,3285},
                             {5,3668},
                             {5,4004},
                             {4,4291},
                             {4,4579}};
    FrameLabeling  haca_02 = {{2,241},
                              {2,435},
                              {2,676},
                              {2,966},
                              {6,1257},
                              {6,1547},
                              {6,1837},
                              {7,2127},
                              {7,2417},
                              {7,2707},
                              {2,2997},
                              {9,3287},
                              {9,3577},
                              {9,3867},
                              {9,4157},
                              {9,4447},
                              {9,4738},
                              {2,5028},
                              {2,5318},
                              {2,5511},
                              {2,5695},
                              {2,5985},
                              {8,6265},
                              {3,6555},
                              {3,6845},
                              {3,7136},
                              {3,7426},
                              {10,7706},
                              {10,7991},
                              {10,8281},
                              {10,8571},
                              {10,8862},
                              {4,9152},
                              {4,9359},
                              {4,9553},
                              {2,9843},
                              {2,10036},
                              {2,10326},
                              {2,10617}};

    FrameLabeling haca_03 = {{2,286},
                             {2,573},
                             {2,859},
                             {7,1098},
                             {7,1385},
                             {7,1671},
                             {7,1958},
                             {3,2244},
                             {3,2531},
                             {2,2817},
                             {2,3056},
                             {2,3295},
                             {2,3582},
                             {5,3868},
                             {5,4155},
                             {5,4441},
                             {5,4728},
                             {10,5014},
                             {10,5205},
                             {10,5396},
                             {4,5683},
                             {4,5965},
                             {4,6251},
                             {9,6538},
                             {9,6777},
                             {9,7063},
                             {2,7350},
                             {2,7589},
                             {2,7875},
                             {2,8114},
                             {2,8401}};

    FrameLabeling haca_04 = {{2,287},
                             {2,575},
                             {2,814},
                             {2,1102},
                             {11,1389},
                             {11,1677},
                             {11,1964},
                             {11,2252},
                             {4,2525},
                             {4,2808},
                             {4,3081},
                             {4,3368},
                             {5,3656},
                             {5,3943},
                             {5,4231},
                             {2,4519},
                             {2,4806},
                             {2,5046},
                             {12,5333},
                             {12,5621},
                             {12,5908},
                             {13,6196},
                             {13,6435},
                             {13,6723},
                             {10,7010},
                             {10,7236},
                             {5,7523},
                             {10,7792},
                             {10,8079},
                             {4,8367},
                             {4,8640},
                             {4,8927},
                             {4,9215},
                             {2,9502},
                             {2,9790},
                             {2,10078}};

    FrameLabeling haca_05 = {{2,287},
                             {2,575},
                             {2,862},
                             {3,1150},
                             {3,1437},
                             {3,1725},
                             {10,2008},
                             {10,2291},
                             {3,2578},
                             {3,2866},
                             {3,3153},
                             {3,3441},
                             {3,3729},
                             {3,4016},
                             {2,4304},
                             {2,4591},
                             {4,4783},
                             {4,4975},
                             {4,5166},
                             {5,5454},
                             {5,5742},
                             {5,5957},
                             {9,6245},
                             {9,6523},
                             {12,6811},
                             {12,7093},
                             {12,7381},
                             {2,7668},
                             {2,7908},
                             {2,8148},
                             {2,8340}};

    FrameLabeling haca_06 = {{2,288},
                             {2,572},
                             {2,861},
                             {2,1149},
                             {8,1438},
                             {8,1630},
                             {7,1919},
                             {7,2208},
                             {7,2496},
                             {8,2785},
                             {8,3074},
                             {5,3362},
                             {5,3651},
                             {5,3940},
                             {4,4180},
                             {4,4372},
                             {4,4565},
                             {12,4839},
                             {12,5080},
                             {12,5272},
                             {12,5465},
                             {4,5753},
                             {4,5946},
                             {4,6234},
                             {10,6523},
                             {10,6763},
                             {10,7052},
                             {3,7341},
                             {3,7629},
                             {3,7918},
                             {9,8207},
                             {9,8495},
                             {9,8784},
                             {9,9073},
                             {2,9361},
                             {2,9650},
                             {2,9939}};

    FrameLabeling haca_07 ={{2,286},
                            {2,572},
                            {2,849},
                            {2,1135},
                            {10,1421},
                            {10,1707},
                            {9,1989},
                            {9,2275},
                            {9,2561},
                            {10,2848},
                            {10,3086},
                            {10,3372},
                            {10,3659},
                            {9,3940},
                            {9,4226},
                            {9,4513},
                            {3,4704},
                            {3,4990},
                            {3,5276},
                            {4,5562},
                            {4,5849},
                            {2,6135},
                            {2,6421},
                            {2,6660},
                            {2,6946},
                            {7,7232},
                            {7,7514},
                            {7,7800},
                            {2,8086},
                            {2,8325},
                            {2,8515},
                            {2,8702}};

    FrameLabeling haca_08 = {{2,290},
                             {2,575},
                             {2,817},
                             {2,1107},
                             {10,1300},
                             {10,1494},
                             {10,1784},
                             {12,2059},
                             {12,2349},
                             {12,2639},
                             {9,2930},
                             {9,3220},
                             {2,3510},
                             {2,3703},
                             {2,3993},
                             {5,4283},
                             {5,4573},
                             {8,4864},
                             {7,5154},
                             {7,5444},
                             {7,5734},
                             {2,6024},
                             {2,6314},
                             {11,6595},
                             {11,6885},
                             {11,7175},
                             {4,7465},
                             {4,7755},
                             {4,8045},
                             {2,8335},
                             {2,8625},
                             {2,8915},
                             {2,9206}};

    FrameLabeling haca_09 = {{2,222},
                             {2,458},
                             {2,728},
                             {2,1011},
                             {9,1295},
                             {9,1579},
                             {9,1862},
                             {9,2146},
                             {12,2382},
                             {12,2666},
                             {12,2950},
                             {10,3143},
                             {10,3422},
                             {10,3659},
                             {2,3942},
                             {2,4226},
                             {2,4510},
                             {2,4794}};

    FrameLabeling haca_10 = {{2,282},
                             {2,564},
                             {2,828},
                             {2,1007},
                             {2,1289},
                             {2,1572},
                             {2,1807},
                             {9,2089},
                             {9,2372},
                             {9,2654},
                             {9,2937},
                             {9,3219},
                             {9,3497},
                             {9,3779},
                             {2,4057},
                             {2,4335},
                             {2,4617},
                             {2,4900},
                             {8,5182},
                             {8,5323},
                             {8,5606},
                             {7,5888},
                             {7,6170},
                             {7,6453},
                             {7,6735},
                             {2,7018},
                             {2,7300},
                             {2,7583}};
    FrameLabeling haca_11 = {{2,167},
                             {2,334},
                             {2,497},
                             {2,665},
                             {2,827},
                             {2,995},
                             {9,1162},
                             {9,1320},
                             {9,1488},
                             {9,1655},
                             {10,1823},
                             {10,1990},
                             {10,2157},
                             {10,2325},
                             {9,2492},
                             {9,2618},
                             {9,2785},
                             {12,2911},
                             {12,3036},
                             {12,3162},
                             {12,3288},
                             {12,3413},
                             {10,3539},
                             {10,3706},
                             {10,3874},
                             {10,4041},
                             {9,4208},
                             {9,4376},
                             {9,4543},
                             {9,4711},
                             {2,4878},
                             {2,5041},
                             {2,5208},
                             {2,5357},
                             {2,5525},
                             {2,5674}};

    FrameLabeling haca_12 = {{2,289},
                             {2,481},
                             {2,722},
                             {2,1011},
                             {4,1296},
                             {4,1585},
                             {4,1874},
                             {5,2158},
                             {5,2447},
                             {5,2736},
                             {5,3025},
                             {5,3314},
                             {6,3604},
                             {6,3893},
                             {6,4182},
                             {7,4471},
                             {7,4760},
                             {10,5049},
                             {10,5338},
                             {9,5627},
                             {9,5868},
                             {9,6157},
                             {9,6446},
                             {9,6735},
                             {9,7025},
                             {9,7314},
                             {9,7603},
                             {2,7892},
                             {2,8085},
                             {2,8277},
                             {2,8566},
                             {2,8856}};

    FrameLabeling haca_13 = {{2,285},
                             {2,479},
                             {2,770},
                             {2,1061},
                             {6,1351},
                             {6,1594},
                             {10,1865},
                             {10,2156},
                             {10,2446},
                             {8,2737},
                             {8,3028},
                             {8,3318},
                             {6,3609},
                             {4,3900},
                             {4,4190},
                             {4,4481},
                             {4,4772},
                             {9,5058},
                             {9,5348},
                             {2,5639},
                             {2,5930},
                             {2,6221}};

    FrameLabeling haca_14 = {{6,279},
                             {6,558},
                             {9,838},
                             {9,1117},
                             {9,1397},
                             {9,1676},
                             {9,1956},
                             {6,2235},
                             {6,2515},
                             {6,2794},
                             {10,3074},
                             {10,3260},
                             {10,3539},
                             {10,3819},
                             {10,4098},
                             {6,4378},
                             {6,4657},
                             {6,4937},
                             {6,5216},
                             {2,5496},
                             {2,5775},
                             {2,6055}};

    stdc::SWatch swatch;

    swatch.start();
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
    CurveClusterer cc(-1, false);
    Curves allCurves={c1};//,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14}};
    std::vector<FrameLabeling> hacaAllLabelings = {haca_01,haca_02,haca_03,haca_04,haca_05,haca_06,haca_07,haca_08,haca_09,haca_10,haca_11,haca_12,haca_13,haca_14};
    std::vector<FrameLabeling> allLabelings = {gt86_01};//,gt86_02,gt86_03,gt86_04,gt86_05,gt86_06,gt86_07,gt86_08,gt86_09,gt86_10,gt86_11,gt86_12,gt86_13,gt86_14};
    swatch.stop();
    std::cout << "----------\nElapsed time loading curves: " << std::chrono::duration<double>(swatch.elapsed()).count() << "\n----------\n";

    std::vector<std::vector<double>> times;
    std::vector<std::vector<int>> lengths;
    std::vector<std::vector<int>> results;
    //for(double delta = 0.25;delta <= 3.0; delta += 0.25) {
    double delta = 1.25;
        times.emplace_back();
        lengths.emplace_back();
        results.emplace_back();
        for (int i = 0; i < allCurves.size(); i++) {
            swatch.reset();
            swatch.start();
            Curves curves = allCurves;
            std::vector<FrameLabeling> labelings = allLabelings;
            cc.initCurvesWithGT(curves, delta, labelings);

            lengths.back().push_back(cc.simplifiedCurves[0].size());

            //custom filter which accesses ccs curves (not very cool tbh)
            int length = 0;
            for (const Curve &c: cc.simplifiedCurves) {
                length += c.size();
            }
            int complexity = 10;//(int) ((length) / (25 * cc.simplifiedCurves.size()));
            std::cout << "Complexity: " << complexity << std::endl;
            double guarantee = 1.0 + 1.0 + 2 * (7.0 / 3.0);
            auto filter = [=](const Candidate &a) {
                bool istrivial = complexity == 1;
                bool isdown = a.getEnd() < a.getBegin();
                bool nontrivial_length =
                        cc.simplifiedCurves[a.getCurveIndex()].subcurve_length(a.getBegin(), a.getEnd()) >
                        2 * guarantee * delta;
                bool nontrivial_complexity = a.getEnd().getPoint() > a.getBegin().getPoint() + complexity / 3;
                return istrivial || isdown || (nontrivial_length && nontrivial_complexity);
            };

            auto trivialFilter = [](const Candidate &a) { return true; };

            //cover algorithm
            auto result = cc.greedyCover(complexity, 25, trivialFilter);
            swatch.stop();

            for (int i = 0; i < std::min(20, (int) result.size()); ++i) {
                Cluster c = result[i];
                for (int j = 0; j < c.getMatching().size(); ++j) {
                    auto m = c.getMatching()[j];
                    io::exportSubcurve(
                            "/Users/styx/data/curveclustering/results/cluster/matching" + std::to_string(i) + "/interval" +
                            std::to_string(j) + ".txt", cc.simplifiedCurves[m.getCurveIndex()], m.getBegin(), m.getEnd(), 100);
                }
            }

            //std::cout << "Solution of size " << result.size() << " computed!\n";
            //std::cout << "----------\nElapsed time computing covering: " << std::chrono::duration<double>(swatch.elapsed()).count() << "\n----------\n";
            times.back().push_back(std::chrono::duration<double>(swatch.elapsed()).count());
            results.back().push_back(result.size());
            //ClusteringVisulaizer cv{true};
            //cv.showClusteringStretched(cc.simplifiedCurves, cc.simplifiedGTs, result);
        }
    //}
    /*
    double delta = 0.25;
    for(auto & i : times) {
        std::cout << delta << ":  ";
        for (auto time: i) {
            std::cout << time << " ";
        }
        std::cout << std::endl;
        delta += 0.25;
    }
    std::cout << std::endl;
    delta = 0.25;
    for(auto & i : lengths) {
        std::cout << delta << ":  ";
        for (auto length: i) {
            std::cout << length << " ";
        }
        std::cout << std::endl;
        delta += 0.25;
    }

    std::cout << std::endl;
    delta = 0.25;
    for(auto & i : results) {
        std::cout << delta << ":  ";
        for (auto length: i) {
            std::cout << length << " ";
        }
        std::cout << std::endl;
        delta += 0.25;
    }
     */
}

void experiments2(){
    Curves curves;
    for(int i=1;i<210/*9*/;i++){
        std::string name = "/Users/styx/data/gdac2/world3d_txt/"+std::to_string(i)+"_drifter.txt";
        curves.push_back(Curve(name,3));
        if(curves.back().size() <= 1){
            curves.pop_back();
            continue;
        }
        std::string delimiter = "world3d_txt";
        std::string token1 = name.substr(0, name.find(delimiter));
        std::string token2 = name.substr(name.find(delimiter) + delimiter.length());
        curves.back().set_name( token1 + "simp" + token2);
    }
    auto f = [](Point& p){
        double x = p[0];
        double y = p[1];
        double z = p[2];
        double tan = z/(sqrt((x*x)+(y*y)));
        double c1 = sqrt(9.81*3600)/100; //Magic 100
        double sin = tan/sqrt(1+(tan*tan));
        double cos = 1/sqrt(1+(tan*tan));
        double twoOmega = 2*7.2921*0.00001;
        double L = c1/abs(twoOmega*sin);
        double Ls = sqrt(c1/(2*twoOmega*(1.0/6371000.0)*cos));

        double phi = atan(tan)*180/M_PI;
        double r = fmin(L,Ls)/200000; //Normalize with 200km
        //double f = 2*7.2921*0.00001*(tan/(sqrt(1+(tan*tan))));
        return r;//sqrt(9.81*3600)/abs(f);
    };

    //for(auto& c : curves){
    //    c.assignWeights(f);
    //}

    //int complexity = 10;
    double guarantee = 1.0 + 1.0 + 2 * (7.0 / 3.0);
    //double delta = 10000;

    std::vector<double> deltas {5000,10000,25000,50000};
    std::vector<int> complexities{1,5,10};
    std::vector<int> lengths{(int)(curves.size())/50,(int)(curves.size())/20,(int)(curves.size())/10,(int)(curves.size())/5,(int)(curves.size())/3,(int)(curves.size())/2,(int)(curves.size())};
    //std::vector<double> deltas {10000};
    //std::vector<int> complexities{10};
    //std::vector<int> lengths{((int)(curves.size())/10)};

    std::vector<std::vector<std::vector<long long>>> results;
    std::vector<std::vector<std::vector<int>>> sizes;
    std::vector<int> realLengths;
    std::vector<std::vector<std::vector<double>>> times;

    //std::cout << delta << std::endl;

    for(auto delta:deltas){
        results.emplace_back();
        sizes.emplace_back();
        times.emplace_back();
        for(auto complexity : complexities){
            results.back().emplace_back();
            sizes.back().emplace_back();
            times.back().emplace_back();
            for(auto length : lengths){
                //construct subset of curves
                Curves subset;
                int l = 0;
                for(int i=1;i<length;i++) {
                    subset.emplace_back(curves[i]);
                    l += curves[i].size();
                }
                realLengths.emplace_back(l);

                stdc::SWatch swatch;
                swatch.start();
                CurveClusterer cc(-1, false);
                cc.initCurves(subset,delta);
                auto filter = [=](Candidate &a) {
                    bool withIsTrivial = false;
                    bool withIsDown = false;
                    bool istrivial = complexity == 1;
                    bool isdown = a.getEnd() < a.getBegin();
                    bool nontrivial_length =
                            cc.simplifiedCurves[a.getCurveIndex()].subcurve_length(a.getBegin(), a.getEnd()) >
                            2 * guarantee * delta;
                    bool nontrivial_complexity = a.getEnd().getPoint() > a.getBegin().getPoint() + complexity / 4;
                    return (withIsTrivial && istrivial) || (withIsDown && isdown) || (nontrivial_length && nontrivial_complexity);
                };
                auto trivialFilter = [&](Candidate &c){return true;};
                long long candidateSetSize;
                auto result = cc.greedyCover(complexity,1,trivialFilter,&candidateSetSize);
                swatch.stop();
                double time = std::chrono::duration<double>(swatch.elapsed()).count();
                std::cout << candidateSetSize << " " <<time<< std::endl;
                sizes.back().back().emplace_back(result.size());
                results.back().back().emplace_back(candidateSetSize);
                times.back().back().emplace_back(time);

            }
        }
    }

    for(int i=0;i<deltas.size();i++) {
        for (int j = 0; j < complexities.size(); j++) {
            for (int k = 0; k < lengths.size(); k++) {
                std::cout << deltas[i] << " " << complexities[j] << " " << lengths[k] << " " << realLengths[k] << " " << results[i][j][k] << " " << times[i][j][k] << std::endl;
            }
        }
    }

/*
    for (int i = 0; i < std::min(50, (int) result.size()); ++i) {
        Candidate c = result[i];
        io::exportSubcurve("/Users/styx/data/gdac2/resultcenters/candidate"+ std::to_string(i)+".txt",
                           cc.simplifiedCurves[c.getCurveIndex()],c.getBegin(),c.getEnd());
        for (int j = 0; j < c.visualMatching.size(); ++j) {
            auto m = c.visualMatching[j];
            CurveID originalID = cc.simpIDtoOriginID[m.getCurveIndex()];
            Curve& originalCurve = curves[originalID];
            auto s = cc.mapSimplificationToBase(m.getCurveIndex(),m.getBegin());
            auto t = cc.mapSimplificationToBase(m.getCurveIndex(),m.getEnd());
            io::exportSubcurve(
                    "/Users/styx/data/gdac2/results/matching" + std::to_string(i) + "/interval" +
                    std::to_string(j) + ".txt", originalCurve, s, t);
            io::exportSubcurve(
                    "/Users/styx/data/gdac2/resultsimp/matching" + std::to_string(i) + "/interval" +
                    std::to_string(j) + ".txt", cc.simplifiedCurves[m.getCurveIndex()], m.getBegin(), m.getEnd());
        }
    }
*/
    std::cout << "Done clustering!\n";
}

void experiments3(){
    Curves curves;
    for(int i=1;i<2109;i++){
        std::string name = "/Users/styx/data/gdac2/world3d_txt/"+std::to_string(i)+"_drifter.txt";
        curves.emplace_back(name,3);
        if(curves.back().size() <= 1){
            curves.pop_back();
            continue;
        }
        std::string delimiter = "world3d_txt";
        std::string token1 = name.substr(0, name.find(delimiter));
        std::string token2 = name.substr(name.find(delimiter) + delimiter.length());
        curves.back().set_name( token1 + "simp" + token2);
    }

    auto f = [](Point& p){
        double x = p[0];
        double y = p[1];
        double z = p[2];
        double tan = z/(sqrt((x*x)+(y*y)));
        double c1 = sqrt(9.81*3600)/100; //Magic 100
        double sin = tan/sqrt(1+(tan*tan));
        double cos = 1/sqrt(1+(tan*tan));
        double twoOmega = 2*7.2921*0.00001;
        double L = c1/abs(twoOmega*sin);
        double Ls = sqrt(c1/(2*twoOmega*(1.0/6371000.0)*cos));

        double phi = atan(tan)*180/M_PI;
        double r = fmin(L,Ls)/200000; //Normalize with 200km
        //double f = 2*7.2921*0.00001*(tan/(sqrt(1+(tan*tan))));
        return r;//sqrt(9.81*3600)/abs(f);
    };

    for(auto& c : curves){
        c.assignWeights(f);
    }

    int complexity = 10;
    double guarantee = 1.0 + 1.0 + 2 * (7.0 / 3.0);
    double delta = 200000/guarantee;

    CurveClusterer cc(-1, false);
    cc.initCurves(curves,delta);
    auto filter = [=](const Candidate &a) {
        bool withIsTrivial = false;
        bool withIsDown = false;
        bool istrivial = complexity == 1;
        bool isdown = a.getEnd() < a.getBegin();
        bool nontrivial_length =
                cc.simplifiedCurves[a.getCurveIndex()].subcurve_length(a.getBegin(), a.getEnd()) >
                2 * guarantee * delta;
        bool nontrivial_complexity = a.getEnd().getPoint() > a.getBegin().getPoint() + complexity / 4;
        return (withIsTrivial && istrivial) || (withIsDown && isdown) || (nontrivial_length && nontrivial_complexity);
    };
    auto trivialFilter = [=](const Candidate &c){return true;};
    auto result = cc.greedyCover(complexity,1,trivialFilter);
    std::cout << "Done clustering!\n";


    for (int i = 0; i < std::min(500, (int) result.size()); ++i) {
        auto c = result[i];
        io::exportSubcurve("/Users/styx/data/gdac2/result/resultcenters/candidate"+ std::to_string(i)+".txt",
                           cc.simplifiedCurves[c.getCenter().getCurveIndex()],c.getCenter().getBegin(),c.getCenter().getEnd());
        for (int j = 0; j < c.getMatching().size(); ++j) {
            auto m = c.getMatching()[j];
            CurveID originalID = cc.simpIDtoOriginID[m.getCurveIndex()];
            Curve& originalCurve = curves[originalID];
            auto s = cc.mapSimplificationToBase(m.getCurveIndex(),m.getBegin());
            auto t = cc.mapSimplificationToBase(m.getCurveIndex(),m.getEnd());
            io::exportSubcurve(
                    "/Users/styx/data/gdac2/result/results/matching" + std::to_string(i) + "/interval" +
                    std::to_string(j) + ".txt", originalCurve, s, t);
            io::exportSubcurve(
                    "/Users/styx/data/gdac2/result/resultsimp/matching" + std::to_string(i) + "/interval" +
                    std::to_string(j) + ".txt", cc.simplifiedCurves[m.getCurveIndex()], m.getBegin(), m.getEnd());
        }
    }
    std::cout << "Done writing!\n";
}

void experiments4(){
    Curves curves;
    for(int i=1;i<5000;i++){
        std::string name = "/Users/styx/data/gdac3/world3d_txt/"+std::to_string(i)+"_drifter.txt";
        curves.emplace_back(name,3);
        if(curves.back().size() <= 1){
            curves.pop_back();
            continue;
        }
        std::string delimiter = "world3d_txt";
        std::string token1 = name.substr(0, name.find(delimiter));
        std::string token2 = name.substr(name.find(delimiter) + delimiter.length());
        curves.back().set_name( token1 + "simp" + token2);
    }

    auto f = [](Point& p){
        double x = p[0];
        double y = p[1];
        double z = p[2];
        double tan = z/(sqrt((x*x)+(y*y)));
        double c1 = sqrt(9.81*3600)/100; //Magic 100
        double sin = tan/sqrt(1+(tan*tan));
        double cos = 1/sqrt(1+(tan*tan));
        double twoOmega = 2*7.2921*0.00001;
        double L = c1/abs(twoOmega*sin);
        double Ls = sqrt(c1/(2*twoOmega*(1.0/6371000.0)*cos));

        double phi = atan(tan)*180/M_PI;
        double r = fmin(L,Ls)/200000; //Normalize with 200km
        //double f = 2*7.2921*0.00001*(tan/(sqrt(1+(tan*tan))));
        return r;//sqrt(9.81*3600)/abs(f);
    };

    for(auto& c : curves){
        c.assignWeights(f);
    }

    int complexity = 10;
    double guarantee = 1.0 + 1.0 + 2 * (7.0 / 3.0);
    double delta = 100000/guarantee;

    CurveClusterer cc(-1, false);
    cc.initCurves(curves,delta);
    auto filter = [=](const Candidate &a) {
        bool withIsTrivial = false;
        bool withIsDown = false;
        bool istrivial = complexity == 1;
        bool isdown = a.getEnd() < a.getBegin();
        bool nontrivial_length =
                cc.simplifiedCurves[a.getCurveIndex()].subcurve_length(a.getBegin(), a.getEnd()) >
                2 * guarantee * delta;
        bool nontrivial_complexity = a.getEnd().getPoint() > a.getBegin().getPoint() + complexity / 4;
        return (withIsTrivial && istrivial) || (withIsDown && isdown) || (nontrivial_length && nontrivial_complexity);
    };
    auto trivialFilter = [=](const Candidate &c){return true;};
    auto result = cc.greedyCover(complexity,1,trivialFilter);
    std::cout << "Done clustering!\n";

    //rossby radius
    //hierarchische zeit

    for (int i = 0; i < std::min(500, (int) result.size()); ++i) {
        auto c = result[i];
        //std::cout <<"[" << (c.getEnd().getPoint() - c.getBegin().getPoint()) << "," << c.visualMatching.size() << "],";
        io::exportSubcurve("/Users/styx/data/gdac3/result/resultcentersl10/candidate"+ std::to_string(i)+".txt",
                           cc.simplifiedCurves[c.getCenter().getCurveIndex()],c.getCenter().getBegin(),c.getCenter().getEnd());
    }
    std::cout << "Done writing!\n";
}

int main(int argc, char *argv[]) {
#ifdef HASVISUAL
    std::cout << "COMPILED WITH OPENCV\n";
#else
    std::cout << "NOT COMPILED WITH OPENCV\n";
#endif

    experiments();
    //experiments2();
    //experiments3();
    //experiments2();
    //intersectionprimitivetest();
/*
    Curve c1 = Curve("../data/86_1.txt", 93);
    Curve c2 = Curve("../data/86_2.txt", 93);

    double delta = 1.25;
    CurveClusterer cc(-1,true);
    Curves curves ={c1,c2};
    cc.initCurves(curves,delta);
    int complexity = 10;
    cc.greedyCover(complexity,10,[&complexity, &delta, &cc](const Candidate& c){return (c.getBegin()>c.getEnd() || (c.getEnd().getPoint() > c.getBegin().getPoint() + complexity/4)) && cc.simplifiedCurves[c.getCurveIndex()].subcurve_length(c.getBegin(),c.getEnd())>delta*5;});
    */
    //experiment1("/Users/styx/data/curveclustering/results/cluster");
    //experiments();

    return 0;
}

