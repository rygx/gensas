// Elementary single-argument functions of the Modelica subset.
//
// Native functions are solved by their own series recurrences (see SasComputation.cpp).
// Derived functions are rewritten into native ones before solving (see HybridModel.cpp).
// Both kinds are evaluated numerically here when they appear in parameters, start values,
// conditions and event expressions.
#ifndef SAS_MATH_FUNCTIONS_H
#define SAS_MATH_FUNCTIONS_H

#include <cmath>
#include <cstring>
#include <string>

namespace che {
	namespace core {

		// sin, cos, exp, sqrt, inv and the inverse functions log, atan, asin and acos.
		inline bool isNativeMathFunction(const std::string &n) {
			static const char *const names[] = {"sin", "cos", "exp", "sqrt", "inv", "log", "atan", "asin", "acos", NULL};
			for (int i = 0; names[i] != NULL; i++) {
				if (n == names[i]) {
					return true;
				}
			}
			return false;
		}

		// tan, sinh, cosh, tanh and log10, expressed with the native functions.
		inline bool isDerivedMathFunction(const std::string &n) {
			static const char *const names[] = {"tan", "sinh", "cosh", "tanh", "log10", NULL};
			for (int i = 0; names[i] != NULL; i++) {
				if (n == names[i]) {
					return true;
				}
			}
			return false;
		}

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
