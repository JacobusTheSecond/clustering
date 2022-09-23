#include "center_clustering_algs.h"
#include "io.h"
#include "curve_simplification.h"
#include "frechet_light.h"
#include <random>

#include <iostream>

void printUsage()
{
	std::cout << "USAGE: ./main <dataset_file> <k> <l> [<header_size>]\n"
	          << "The dataset_file should contain a newline separated list of filenames of curve files.\n"
	          << "The header size (default: 1) gives the number of lines which are ignored at the beginning of each curve file.\n";
}

int main(int argc, char* argv[])
{
#ifdef HASVISUAL
    std::cout << "COMPILED WITH OPENCV\n";
#elif
    std::cout << "NOT COMPILED WITH OPENCV\n";
#endif
    //std::cout << "Hello world!\n";
    FrechetLight fl;
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
/*
    Curve c1 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_directions_1.txt",96);
    Curve c2 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_directions_2.txt",96);
    Curve c3 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_discussion_1.txt",96);
    Curve c4 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_discussion_2.txt",96);
    Curve c5 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_eating_1.txt",96);
    Curve c6 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_eating_2.txt",96);
    Curve c7 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_greeting_1.txt",96);
    Curve c8 = Curve("/Users/styx/data/curveclustering/fixedseq3d/fixedseq3d_S11_greeting_2.txt",96);
*/


    //andreas weber heuristik
    //visualize candidates DONE
    //normalize direction DONE
    std::vector<int> times;

    Curve c1 = Curve("/Users/styx/data/curveclustering/cmu_fixedseq3d/86_1.txt",93);
    Curve c2 = Curve("/Users/styx/data/curveclustering/cmu_fixedseq3d/86_6.txt",93);

    std::cout << "Input complexities: " << c1.size() << " " << c2.size() << std::endl;

    double delta = 1.0;
    double deltaprime = 1.5;

    Curve simpC1 = good_simplification(c1,delta,&times);
    Curve simpC2 = good_simplification(c2,delta);


    Curve unsimpC = Curve(simpC1,times);
    io::exportCurve("/Users/styx/data/curveclustering/results/testcurve.txt",unsimpC);
    io::exportCurve("/Users/styx/data/curveclustering/results/testcurveorigin.txt",c1);
    Curves curvesCMU = {simpC1,simpC2};
    double guaranteeCMU = 1.0+1.0+2*(7.0/3.0);
    int length = 0;
    for(const Curve& c : curvesCMU){
        length += c.size();
    }
/*
    Curves result = greedyCover(curvesCMU,deltaprime, (int) ((length) / (25 * curvesCMU.size())), 25,true);

    for(int i=0;i<result.size();++i){
        io::exportCurve("/Users/styx/data/curveclustering/results/resultcurve" + std::to_string(i) + ".txt",result[i]);
    }
*/

    std::vector<std::vector<int>> sizes;

    for(double sd=0.25;sd<5;sd+=0.25) {
        sizes.emplace_back();
        Curve sc1 = good_simplification(c1,sd);
        Curve sc2 = good_simplification(c2,sd);
        Curves scs = {sc1,sc2};
        int l = 0;
        for(const Curve& c : scs){
            l += c.size();
        }
        for (double d = 0.25; d <= 5; d += 0.25) {
            std::cout << "\nSIMPLIFICATION DELTA = " << sd << " FREESPACE DELTA = " << d << "\n";
            Curves result = greedyCover(scs, d, (int) ((l) / (25 * scs.size())), 25);
            sizes.back().push_back(result.size());
        }
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

    greedyCover(curves,delta,10);
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
    return 0;
}

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
