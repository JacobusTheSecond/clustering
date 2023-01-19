#pragma once

#include "curve.h"
#include "io.h"

Curve simplify(Curve const& curve, int l);

Curve good_simplification(const Curve &c, distance_t delta, std::vector<int>* indices = nullptr);
Curve good_simplification(const Curve &c, const distance_t delta, std::vector<std::pair<Label,int>>& groundthruth, std::vector<std::pair<Label,ParamPoint>>* simplifiedLabels);