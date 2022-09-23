// Make sure to always print in the unit tests
#ifdef NVERBOSE
#	undef NVERBOSE
#endif
#include "unit_tests.h"

#include "defs.h"
#include "geometric_algs.h"

//
// Define some helpers
//

#define TEST(x)                                                                \
	do {                                                                       \
		if (!(x)) {                                                            \
			std::cout << "\n";                                                 \
			std::cout << "TEST_FAILED!\n";                                     \
			std::cout << "File: " << __FILE__ << "\n";                         \
			std::cout << "Line: " << __LINE__ << "\n";                         \
			std::cout << "Function: " << __func__ << "\n";                     \
			std::cout << "Test: " << #x << "\n";                               \
			std::cout << "\n";                                                 \
			std::cout << std::flush;                                           \
			std::abort();                                                      \
		}                                                                      \
	} while (0)

//
// Unit Tests
//

void unit_tests::testAll()
{
	testMEC();
}

void unit_tests::testMEC()
{
	auto circle1 = calcMinEnclosingCircle(
		Points{{0.,0.}, {0.,2.}, {1., 1.}, {0.3, 1.3}, {.6,.6}});
	TEST(circle1.center == Point(0., 1.) && circle1.radius == 1.);

	auto circle2 = calcMinEnclosingCircle(
		Points{{0.,0.}, {4.,0.}, {2., .5}, {2., -.5}, {.1,.1}});
	TEST(circle2.center == Point(2., 0.) && circle2.radius == 2.);
}

// just in case anyone does anything stupid with this file...
#undef TEST
