
#include "center_clustering_algs.h"

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

bool cmpLeftLower(Subcurve l, Subcurve r){
    return (l.curveIdx < r.curveIdx) || (l.curveIdx == r.curveIdx && (l.start < r.start));
}

bool cmpLeftUpper(Subcurve l, Subcurve r){
    return (l.curveIdx < r.curveIdx) || (l.curveIdx == r.curveIdx && (l.start < r.end));
}

bool cmpRightUpper(Subcurve l, Subcurve r){
    return (l.curveIdx < r.curveIdx) || (l.curveIdx == r.curveIdx && (l.end < r.end));
}
//TODO: something here is fucked.
void updateCandidate(Curves& curves, Candidate& c, std::vector<Subcurve> covering, std::vector<double> suffixLengths, int roundID){
    double newLength = 0;

    for(auto matching : c.matchings){

        auto l = std::upper_bound(covering.begin(),covering.end(),matching,cmpLeftUpper);
        int lIdx = l - covering.begin();
        //lIdx points to the first interval that can lie right of s

        auto r = std::upper_bound(covering.begin(),covering.end(),matching, cmpRightUpper);
        int rIdx = r - covering.begin();
        //rIdx points to the first interval that can lie right of t

        bool leftContains = (l!=covering.end()) && (l->curveIdx == matching.curveIdx) && (l->contains(matching.start));
        bool rightContains = (r!=covering.end()) && (r->curveIdx == matching.curveIdx) && (r->contains(matching.end));

        if(leftContains && rightContains && lIdx == rIdx)
            continue;

        newLength += curves[matching.curveIdx].subcurve_length(matching.start,matching.end);

        if(leftContains){
            newLength -= curves[matching.curveIdx].subcurve_length(matching.start,l->end);
            lIdx ++;
        }
        //now lIdx points to the first interval that lies strictly to the right of s

        if(rightContains){
            newLength -= curves[matching.curveIdx].subcurve_length(r->start,matching.end);
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
    std::vector<Candidate> temp;
    for(int i=0;i<5;++i){
        Candidate c = cs.top();
        cs.pop();
        temp.push_back(c);
        std::cout << "(" << c.roundOfUpdate << ";" << c.semiUpdatedCoverLength << ") ";
    }
    std::cout << "\n";
    for(Candidate c : temp)
        cs.push(c);
}

double lengthOfUncovered(Curves curves, std::vector<Candidate> candidateSet){
    //step 1: merge
    std::vector<Subcurve> presorted;
    for(auto c:candidateSet){
        presorted.insert(presorted.end(),c.matchings.begin(),c.matchings.end());
    }
    std::sort(presorted.begin(), presorted.end(), cmpLeftLower);


    double uncovered = 0;
    std::pair<int, ParamPoint> pcur = {0, {0, 0}};

    std::vector<Subcurve> covering;

    if(!presorted.empty()) {
        Subcurve cur = presorted[0];
        for (int covI = 1; covI < presorted.size(); covI++) {
            Subcurve next = presorted[covI];
            if ((next.curveIdx == cur.curveIdx) && (next.start < cur.end)) {
                cur.end = std::max(next.end, cur.end);
            } else {
                covering.push_back(cur);
                cur = next;
            }
        }
        covering.push_back(cur);

        //now covering is a disjoint set of intervals
        for (auto cov: covering) {
            while (pcur.first < cov.curveIdx) {
                uncovered += curves[pcur.first].subcurve_length(pcur.second,
                                                                {(int) (curves[pcur.first].size()) - 2, 1.0});
                pcur = {pcur.first + 1, {0, 0}};
            }
            if (cov.curveIdx == pcur.first) {
                uncovered += std::max(0.0,curves[cov.curveIdx].subcurve_length(pcur.second, cov.start));
                pcur = {pcur.first, cov.end};
                if((ParamPoint){(int)(curves[pcur.first].size()-2),1.0} <= pcur.second ){
                    pcur = {pcur.first + 1, {0, 0}};
                }
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
    std::vector<Candidate> bestResultVisualizer = greedyCoverUnsanitizedOutput(curves,delta,l,max_rounds,show,[=](const Candidate& a){return a.getEnd().id - a.getStart().id > l/4;});
    Curves bestresult;

    for (auto &sub: bestResultVisualizer) {
        bestresult.push_back(
                curves[sub.getIndex()].constructSubcurve(sub.getStart(), sub.getEnd()));
    }

    std::cout << "\nImportances: ";
    for(auto c : bestResultVisualizer){
        std::cout << c.importance << " ";
    }
    std::cout << "\n";
    for(int count = 0;count < ((20<bestResultVisualizer.size())?20:bestResultVisualizer.size());count++) {
        Candidate bestCandidate = bestResultVisualizer[count];
        Candidate bC = bestCandidate;
        for (int i = 0; i < bC.visualMatchings.size(); ++i) {
            auto matching = bC.visualMatchings[i];
            io::exportSubcurve(
                    "/Users/styx/data/curveclustering/results/cluster/matching" + std::to_string(count) +"/interval" + std::to_string(i) + ".txt",
                    curves[matching.curveIdx], matching.start, matching.end, 100);
        }
    }
    return bestresult;
}