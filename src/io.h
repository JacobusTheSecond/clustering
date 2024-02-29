//
// Created by Jacobus Conradi on 17.04.23.
//

#ifndef CLUSTERING_IO_H
#define CLUSTERING_IO_H

#include <utility>
#include "basic_types.h"
#include "geometry_types.h"
#include "Curve.h"

//enum Label {_labelstart,transition,walk,jump,punch,leg_kick,squat,run,stand,arm_up,drink,stretch,slap,turn,_labelend};


namespace io
{
    void exportSubcurve(std::string const& filename, Curve curve, CPoint s, CPoint t);
    void exportSubcurve(std::string const& filename, Curve curve, CPoint s, CPoint t, int interpol);
};


#endif //CLUSTERING_IO_H
