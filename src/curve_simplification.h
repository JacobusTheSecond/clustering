#pragma once

#include "curve.h"
#include "io.h"

Curve simplify(Curve const& curve, int l);

//TODO: simplifier

Curve good_simplification(const Curve &c, distance_t delta, std::vector<int>* indices = nullptr);
Curve good_simplification(const Curve &c, const distance_t delta, int maxSkip, std::vector<int>* indices = nullptr);
Curve good_simplification(const Curve &c, const distance_t delta, FrameLabeling& groundthruth, ParamLabeling* simplifiedLabels, std::vector<int>* indices = nullptr);
Curve good_simplification(const Curve &c, const distance_t delta, FrameLabeling& groundthruth, ParamLabeling* simplifiedLabels, int maxSkip, std::vector<int>* indices = nullptr);