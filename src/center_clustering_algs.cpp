#include "center_clustering_algs.h"
//#include <gmpxx.h>
#include "io.h"
#include <algorithm>
#include "SWatch.h"

/*
namespace
{

Clustering computeCenterClusteringRound(
	Curves const& curves, int k, int l, ClusterAlg cluster_alg, CenterAlg center_alg)
{
	auto clustering = computeClustering(curves, k, l, cluster_alg);
	updateClustering(curves, clustering);

	// iterate as long as there are new centers
	int count = 1;
	int const max_count = 10;
	while (count <= max_count && computerCenters(curves, clustering, l, center_alg)) {
		updateClustering(curves, clustering);
		++count;
	}
	// std::cout << "Number of iterations: " << count << std::endl;

	return clustering;
}

} // end anonymous namespace

Clustering computeCenterClustering(
	Curves const& curves, int k, int l, ClusterAlg cluster_alg, CenterAlg center_alg, int max_rounds)
{
	Clustering min_clustering;
	distance_t min_cost = std::numeric_limits<distance_t>::max();

	for (int round = 0; round < max_rounds; ++round) {
        std::cout << "Round " << round << "/" << max_rounds << "\n";
		auto clustering = computeCenterClusteringRound(curves, k, l, cluster_alg, center_alg);

		distance_t cost_sum = 0.;
		for (auto const& cluster: clustering) {
			cost_sum += cluster.cost;
		}
		if (cost_sum < min_cost) {
			min_clustering = std::move(clustering);
			min_cost = cost_sum;
		}
	}

	// remove empty clusters
	for (ClusterID cluster_id = 0; cluster_id < min_clustering.size(); ++cluster_id) {
		if (min_clustering[cluster_id].curve_ids.empty()) {
			std::swap(min_clustering[cluster_id], min_clustering.back());
			min_clustering.pop_back();
		}
	}

	return min_clustering;
}
 */
/*
std::vector<std::pair<int,int>> kApproxCover(Curves& curves, CandidateSet& cs, double delta,int r, int kL,int kU, int imax){
    //1. iterate up to k
    //sample k'
    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    mpz_t ind;
    mpz_init(ind);

    mpz_t indicesW;
    mpz_init(indicesW);

    cs.resetWeights();

    std::vector<std::pair<int,int>> indices;

    int i=1;
    int ks = kL;
    while(i<imax){
        if(i%(imax/100)==0)
            std::cout << ks << " : (" << i << "/" << imax << ")\n";
        //sample k elements
        mpz_set_ui(indicesW,0);
        for(int j=0;j<ks;++j) {
            mpz_urandomm(ind, rs, cs.tW);
            indices.push_back(cs.draw(ind));
            mpz_add(indicesW,indicesW,cs.getCandidate(indices.back()).w);
        }

        std::pair<int,ParamPoint> t = cs.findNonCovered(indices);
        if(t.first == -1){
            for(auto index : indices){
                std::cout << "(" << index.first << "," << index.second << ") ";
            }
            std::cout << "\nFound one\n";
            cs.showCovering(indices);
            return indices;
        }
        mpz_mul_ui(indicesW,indicesW,r);
        if(ks >= kU || mpz_cmp(cs.tW,indicesW)>=0){
            std::cout << "   reweighted\n";
            cs.reweight(t);
            i++;
            ks = kL;
        }else{
            ks ++;
        }
        indices.clear();
    }
    return {};
}

Clustering approxCoverRound(Curves& curves,CandidateSet& cs, double delta) {

    int k=16;
    int r = 1;
    //doubling dimension is 110*d + 412... as d=96, that would be unreasonably large steps
    int gamma = 4;
    std::vector<std::pair<int,int>> indices;
    while(true){
        std::cout << "k = " << k << "\n";
        int ks = 16*k*gamma*((int)log2(16*k*gamma));
        int imax = 5*k*((int)log2(cs.totalCount/k));
        indices = kApproxCover(curves,cs,delta,r,gamma*k,2*gamma*k,imax);
        if(!indices.empty())
            break;
        cs.resetWeights();
        k=k*2;
        r=k*2;
    }

    //return indices;
}
*/
Clustering approxCover(Curves& curves, double delta,int l,int max_rounds){
    //something about rounds
    //double const guarantee = 1.0+1.0+2*(7.0/3.0);
    //CandidateSet cs = CandidateSet(curves,guarantee*delta);
    //cs.computeCandidates(l);

    //return approxCoverRound(curves,cs,delta);
}

bool cmpLeftLower(std::pair<int,Subcurve> l, std::pair<int,Subcurve> r){
    return (l.first < r.first) || (l.first == r.first && (l.second.start < r.second.start));
}

bool cmpLeftUpper(std::pair<int,Subcurve> l, std::pair<int,Subcurve> r){
    return (l.first < r.first) || (l.first == r.first && (l.second.start < r.second.end));
}

bool cmpRightUpper(std::pair<int,Subcurve> l, std::pair<int,Subcurve> r){
    return (l.first < r.first) || (l.first == r.first && (l.second.end < r.second.end));
}
//TODO: something here is fucked.
void updateCandidate(Curves& curves, Candidate& c, std::vector<std::pair<int,Subcurve>> covering, std::vector<double> suffixLengths, int roundID){
    double newLength = 0;

    for(auto matching : c.matchings){

        auto l = std::upper_bound(covering.begin(),covering.end(),matching,cmpLeftUpper);
        int lIdx = l - covering.begin();
        //lIdx points to the first interval that can lie right of s

        auto r = std::upper_bound(covering.begin(),covering.end(),matching, cmpRightUpper);
        int rIdx = r - covering.begin();
        //rIdx points to the first interval that can lie right of t

        bool leftContains = (l!=covering.end()) && (l->first == matching.first) && (l->second.contains(matching.second.start));
        bool rightContains = (r!=covering.end()) && (r->first == matching.first) && (r->second.contains(matching.second.end));

        if(leftContains && rightContains && lIdx == rIdx)
            continue;

        newLength += curves[matching.first].subcurve_length(matching.second.start,matching.second.end);

        if(leftContains){
            newLength -= curves[matching.first].subcurve_length(matching.second.start,l->second.end);
            lIdx ++;
        }
        //now lIdx points to the first interval that lies strictly to the right of s

        if(rightContains){
            newLength -= curves[matching.first].subcurve_length(r->second.start,matching.second.end);
            //rIdx --;
        }
        //now rIdx points to the last inerval that lies strictly to the left of t
        if(rIdx > lIdx)
            newLength -= (suffixLengths[rIdx] - suffixLengths[lIdx]);
    }
    if(c.semiUpdatedCoverLength + EPSILON < newLength){
        std::cout << "?" << c.semiUpdatedCoverLength << " -> " << newLength<<std::endl;
        //updateCandidate(curves,c,covering,suffixLengths,roundID);
    }
    c.semiUpdatedCoverLength = newLength;
    c.roundOfUpdate = roundID;
}

void printFirst50(CandidateSetPQ& cs){
    std::vector<std::pair<int,Candidate>> temp;
    for(int i=0;i<5;++i){
        std::pair<int,Candidate> c = cs.top();
        cs.pop();
        temp.push_back(c);
        std::cout << "(" << c.second.roundOfUpdate << ";" << c.second.semiUpdatedCoverLength << ") ";
    }
    std::cout << "\n";
    for(std::pair<int,Candidate> c : temp)
        cs.push(c);
}

double lengthOfUncovered(Curves curves, std::vector<std::pair<int,Candidate>> candidateSet){
    //step 1: merge
    std::vector<std::pair<int,Subcurve>> presorted;
    for(auto c:candidateSet){
        presorted.insert(presorted.end(),c.second.matchings.begin(),c.second.matchings.end());
    }
    std::sort(presorted.begin(), presorted.end(), cmpLeftLower);


    double uncovered = 0;
    std::pair<int, ParamPoint> pcur = {0, {0, 0}};

    std::vector<std::pair<int,Subcurve>> covering;


    if(!presorted.empty()) {
        std::pair<int, Subcurve> cur = presorted[0];
        for (int covI = 1; covI < presorted.size(); covI++) {
            std::pair<int, Subcurve> next = presorted[covI];
            if ((next.first == cur.first) && (next.second.start < cur.second.end)) {
                cur.second.end = std::max(next.second.end, cur.second.end);
            } else {
                covering.push_back(cur);
                cur = next;
            }
        }
        covering.push_back(cur);

        //now covering is a disjoint set of intervals
        for (auto cov: covering) {
            while (pcur.first < cov.first) {
                uncovered += curves[pcur.first].subcurve_length(pcur.second,
                                                                {(int) (curves[pcur.first].size()) - 2, 1.0});
                pcur = {pcur.first + 1, {0, 0}};
            }
            if (cov.first == pcur.first) {
                uncovered += curves[cov.first].subcurve_length(pcur.second, cov.second.start);
                pcur = {pcur.first, cov.second.end};
            }
        }
    }
    while(pcur.first < curves.size()){
        uncovered += curves[pcur.first].subcurve_length(pcur.second,{(int)(curves[pcur.first].size())-2,1.0});
        pcur = {pcur.first+1,{0,0}};
    }
    return uncovered;
}

Curves greedyCover(Curves& curves, double delta, int l, int max_rounds, bool show){
    std::vector<std::pair<int,Candidate>> bestResultVisualizer = greedyCoverUnsanitizedOutput(curves,delta,l,max_rounds,show);
    Curves bestresult;

    for (const auto &sub: bestResultVisualizer) {
        bestresult.push_back(
                curves[sub.first].constructSubcurve(sub.second.getStart(), sub.second.getEnd()));
    }

    std::cout << "\nImportances: ";
    for(auto c : bestResultVisualizer){
        std::cout << c.second.importance << " ";
    }
    std::cout << "\n";
    for(int count = 0;count < ((20<bestResultVisualizer.size())?20:bestResultVisualizer.size());count++) {
        std::pair<int, Candidate> bestCandidate = bestResultVisualizer[count];
        Candidate bC = bestCandidate.second;
        for (int i = 0; i < bC.visualMatchings.size(); ++i) {
            auto matching = bC.visualMatchings[i];
            io::exportSubcurve(
                    "/Users/styx/data/curveclustering/results/cluster/matching" + std::to_string(count) +"/interval" + std::to_string(i) + ".txt",
                    curves[matching.first], matching.second.start, matching.second.end, 100);
        }
    }
    return bestresult;
}


std::vector<std::pair<int,Candidate>> greedyCoverUnsanitizedOutput(Curves& curves, double delta, int l, int max_rounds, bool show){
    Curves bestresult;
    stdc::SWatch swatchinsert;
    stdc::SWatch swatchupdate;
    std::vector<std::pair<int,Candidate>> bestResultVisualizer;
    double const guarantee = 1.0+1.0+2*(7.0/3.0);
    CandidateSetPQ cs = CandidateSetPQ(curves,guarantee*delta);
    //cs.computeCandidates(l);
    cs.ultrafastComputeCandidates(l,(int)(l/4));


    io::exportSubcurve("/Users/styx/data/curveclustering/results/bestcandidate.txt",curves[cs.top().first],cs.top().second.getStart(),cs.top().second.getEnd());

    for(int r = 0;r<max_rounds;++r) {

        //swatchinsert.reset();
        //swatchupdate.reset();

        std::vector<std::pair<int,Candidate>> result;
        std::vector<std::pair<int,Subcurve>> covering;
        std::vector<double> suffixLengths;

        std::cout << "\nRound " << r;

        for (int i = 0; i < 501; ++i) {
            auto ID = 0;
            auto roundID = i;
            std::cout << ".";
            //update top, until the roundID matches
            while (cs.top().second.roundOfUpdate != roundID) {
                ID++;
                std::pair<int, Candidate> c = cs.top();
                cs.pop();
                swatchupdate.start();
                updateCandidate(curves, c.second, covering, suffixLengths, roundID);
                swatchupdate.stop();
                swatchinsert.start();
                cs.push(c);
                swatchinsert.stop();
            }
            if (lengthOfUncovered(curves,result) <= EPSILON || cs.top().second.semiUpdatedCoverLength <= EPSILON) {
                std::cout << "\nTrying to refine... ";
                int deletecount = 0;
                for(int igni=result.size()-1;igni>=0;--igni){
                    //i is the index to be ignored
                    std::vector<std::pair<int,Candidate>> temp;
                    for(int j=0;j<result.size();++j){
                        if(igni != j){
                            temp.push_back(result[j]);
                        }
                    }
                    double importance = lengthOfUncovered(curves,temp);
                    if(importance <= EPSILON) {
                        std::cout << " ( " << igni << " , " << importance << " )";
                        result.erase(result.begin() + igni);
                        deletecount += 1;
                    }else{
                        result[igni].second.importance = importance;
                    }
                }
                if (deletecount == 0)
                    std::cout << "nothing";
                std::cout << " greedely deleted\n";
                for(int igni=result.size()-1;igni>=0;--igni){
                    //i is the index to be ignored
                    std::vector<std::pair<int,Candidate>> temp;
                    for(int j=0;j<result.size();++j){
                        if(igni != j){
                            temp.push_back(result[j]);
                        }
                    }
                    double importance = lengthOfUncovered(curves,temp);
                    result[igni].second.importance = importance;
                }
                //std::sort(result.begin(), result.end(), [](auto& a, auto& b){return a.second.importance > b.second.importance;});
                std::cout << " Solution of size " << result.size() << " found. ";
                if (bestresult.empty() || bestresult.size() > result.size()) {
                    if(!bestresult.empty()){
                        std::cout << "This improves on the current best by " << bestresult.size() - result.size() << "! ";
                    }
                    bestresult.clear();
                    for (const auto &sub: result) {
                        bestresult.push_back(
                                curves[sub.first].constructSubcurve(sub.second.getStart(), sub.second.getEnd()));
                        bestResultVisualizer.push_back(sub);
                    }
                }
                break;

            }

            //printFirst50(cs);

            double addedweight;
            //cs.showCovering(result);
            //std::cout << cs.top().second.semiUpdatedCoverLength << " ";
            //std::pair<int,Candidate> ctemp = cs.top();
            //std::cout << cs.top().second.matchings.size();
            //updateCandidate(curves, ctemp.second, covering, suffixLengths, roundID);


            //cs.top() will be added to c
            if (i == 0) {
                std::vector<std::pair<int, Candidate>> temp;
                for (int j = 0; j < r; ++j) {
                    temp.push_back(cs.top());
                    cs.pop();
                }
                std::pair<int, Candidate> c = cs.top();
                result.push_back(c);

                cs.pop();
                addedweight = c.second.semiUpdatedCoverLength;
                c.second.semiUpdatedCoverLength = 0;
                cs.push(c);
                for (const auto &ctemp: temp)
                    cs.push(ctemp);
            } else {
                std::pair<int, Candidate> c = cs.top();
                result.push_back(c);

                //io::exportSubcurve("/Users/styx/data/curveclustering/results/bestcandidate.txt",curves[0],cs.top().second.getStart(),cs.top().second.getEnd());
                //cs.showCovering(result);

                cs.pop();
                addedweight = c.second.semiUpdatedCoverLength;
                c.second.semiUpdatedCoverLength = 0;
                cs.push(c);
            }

            //std::cout << " added weight: " << addedweight << ". Updated " << ID << " lengths\n";


            std::pair<int, Candidate> c = result.back();

            //update covering
            std::vector<std::pair<int, Subcurve>> temp(covering);
            temp.insert(temp.end(), c.second.matchings.begin(), c.second.matchings.end());
            std::sort(temp.begin(), temp.end(), cmpLeftLower);
            covering.clear();
            suffixLengths.clear();
            suffixLengths.push_back(0);
            std::pair<int, Subcurve> cur = temp[0];
            for (int covI = 1; covI < temp.size(); covI++) {
                std::pair<int, Subcurve> next = temp[covI];
                if ((next.first == cur.first) && (next.second.start < cur.second.end)) {
                    cur.second.end = std::max(next.second.end, cur.second.end);
                } else {
                    covering.push_back(cur);
                    suffixLengths.push_back(
                            suffixLengths.back() + curves[cur.first].subcurve_length(cur.second.start, cur.second.end));
                    cur = next;
                }
            }
            covering.push_back(cur);
            suffixLengths.push_back(
                    suffixLengths.back() + curves[cur.first].subcurve_length(cur.second.start, cur.second.end));
            //if((i+1)%100==0)
            //cs.showCovering(result);
        }
        std::cout << std::endl;
        std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(swatchinsert.elapsed()).count() << std::endl;
        std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(swatchupdate.elapsed()).count() << std::endl;
        std::cout << "Cleaning up for next round ... ";
        cs.resetWeights();
        std::cout << "Done";
    }
    if(show) {
#ifdef HASVISUAL
        cs.showCovering(bestResultVisualizer);
#endif
    }
    return bestResultVisualizer;
}
