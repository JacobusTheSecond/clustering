//
// Created by Jacobus Conradi on 17.04.23.
//

#include <vector>
#include "io.h"
#include <fstream>
#include <sstream>
#include "defs.h"

void io::exportSubcurve(std::string const& filename, Curve curve, CPoint s, CPoint t){

    std::ofstream file(filename);
    if (!file.is_open()) {
        ERROR("Couldn't open clustering file: " << filename);
    }
    for(auto coord:curve.eval(s)){
        file << coord << " ";
    }
    file << "\n";

    int step = (s>t)?-1:1;

    for(int i = s.getPoint()+std::max(step,0);(s>t)?i>t.getPoint():i<=t.getPoint();i+=step){
        for(auto coord:curve[i]){
            file << coord << " ";
        }
        file << "\n";
    }
    for(auto coord:curve.eval(t)){
        file << coord << " ";
    }
    file << "\n";
    file.close();
}

void io::exportSubcurve(std::string const& filename, Curve curve, CPoint s, CPoint t, int interpol){

    std::ofstream file(filename);
    if (!file.is_open()) {
        ERROR("Couldn't open clustering file: " << filename);
    }
    for(int i=0;i<=interpol;++i){
        //std::cout << "(" << s.interpol(t,i,interpol).id << "," << s.interpol(t,i,interpol).t << ") ";
        for(auto coord:curve.eval(s.interpol(t,i,interpol))){
            file << coord << " ";
        }
        file << "\n";
    }
    std::cout << std::endl;
    file.close();
}