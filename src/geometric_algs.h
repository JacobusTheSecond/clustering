#pragma once

#include "geometry_basics.h"

namespace unit_tests { void testMEC(); }

// call-by-value do to random shuffling
Circle calcMinEnclosingCircle(Points points);
