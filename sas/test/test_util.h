// Minimal test helpers (no external test framework), in the style of the
// existing third_party/Test*.cc programs: a test binary returns non-zero on failure.
#ifndef SAS_TEST_UTIL_H
#define SAS_TEST_UTIL_H

#include <cmath>
#include <iostream>

namespace sastest {

	inline int &failures() {
		static int n = 0;
		return n;
	}

	inline int finish(const char *name) {
		if (failures() == 0) {
			std::cout << "[  PASSED  ] " << name << std::endl;
			return 0;
		}
		std::cerr << "[  FAILED  ] " << name << ": " << failures() << " check(s) failed" << std::endl;
		return 1;
	}

} // namespace sastest

#define CHECK(cond)                                                              \
	do {                                                                         \
		if (!(cond)) {                                                           \
			std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond \
					  << std::endl;                                              \
			++sastest::failures();                                               \
		}                                                                        \
	} while (0)

#define CHECK_EQ(a, b)                                                           \
	do {                                                                         \
		if (!((a) == (b))) {                                                     \
			std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK_EQ failed: " #a \
					  << " == " #b << " (" << (a) << " vs " << (b) << ")"        \
					  << std::endl;                                              \
			++sastest::failures();                                               \
		}                                                                        \
	} while (0)

#define CHECK_NEAR(a, b, tol)                                                   \
	do {                                                                        \
		double va_ = (a), vb_ = (b);                                            \
		if (!(std::fabs(va_ - vb_) <= (tol))) {                                 \
			std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK_NEAR failed: " \
					  << #a << " = " << va_ << ", " << #b << " = " << vb_       \
					  << " (tol " << (tol) << ")" << std::endl;                 \
			++sastest::failures();                                              \
		}                                                                       \
	} while (0)

#endif
