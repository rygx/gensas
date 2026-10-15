// Elementary single-argument functions of the Modelica subset.
//
// Native functions (sin, cos, exp, sqrt, inv, log, atan, asin, acos) are solved by their own series
// recurrences (see SasComputation.cpp). Derived functions (tan, sinh, cosh, tanh, log10) are rewritten
// into native ones before solving (see HybridModel.cpp). Both kinds are evaluated numerically here when
// they appear in parameters, start values, conditions and event expressions.
#ifndef SAS_MATH_FUNCTIONS_H
#define SAS_MATH_FUNCTIONS_H

#include <cmath>
#include <string>

namespace che {
	namespace core {

		// Evaluates a native or derived function of one argument; false if the name is not one of them.
		inline bool evalMathFunction(const std::string &n, double a, double &out) {
			if (n == "sin") {
				out = std::sin(a);
			} else if (n == "cos") {
				out = std::cos(a);
			} else if (n == "exp") {
				out = std::exp(a);
			} else if (n == "sqrt") {
				out = std::sqrt(a);
			} else if (n == "inv") {
				out = 1.0 / a;
			} else if (n == "log") {
				out = std::log(a);
			} else if (n == "atan") {
				out = std::atan(a);
			} else if (n == "asin") {
				out = std::asin(a);
			} else if (n == "acos") {
				out = std::acos(a);
			} else if (n == "tan") {
				out = std::tan(a);
			} else if (n == "sinh") {
				out = std::sinh(a);
			} else if (n == "cosh") {
				out = std::cosh(a);
			} else if (n == "tanh") {
				out = std::tanh(a);
			} else if (n == "log10") {
				out = std::log10(a);
			} else {
				return false;
			}
			return true;
		}

	} // namespace core
} // namespace che

#endif
