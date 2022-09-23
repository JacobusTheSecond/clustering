#pragma once

#include "curve.h"

Curve simplify(Curve const& curve, int l);

Curve good_simplification(const Curve &c, const distance_t delta, std::vector<int>* indices = nullptr);
