/* Compiled as C: checks that gensas.h is valid C and that the basic flow works. */
#include "capi/gensas.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;
#define CHECK(cond)                                                \
	do {                                                           \
		if (!(cond)) {                                             \
			fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
			failures++;                                            \
		}                                                          \
	} while (0)

static int count = 0;
static int stop_after_first(double t, double end, void *user) {
	(void)t;
	(void)end;
	(void)user;
	count++;
	return 1;
}

int main(void) {
	char err[256];
	const char *job =
		"{\"gensas\":{\"format\":\"job\",\"version\":1},\"kind\":\"modelica\","
		"\"model\":{\"source\":\"model D Real x(start=1); equation der(x) = -x; end D;\"},"
		"\"solver\":{\"endTime\":2}}";
	gensas_simulation *sim;
	size_t n;

	CHECK(gensas_api_version() == GENSAS_API_VERSION);
	CHECK(gensas_create("{", NULL, err, sizeof err) == NULL);
	CHECK(strstr(err, "invalid JSON") != NULL);
	CHECK(gensas_create(NULL, NULL, err, sizeof err) == NULL);

	sim = gensas_create(job, NULL, err, sizeof err);
	CHECK(sim != NULL);
	if (sim == NULL) {
		return 1;
	}
	CHECK(gensas_sample_count(sim) == 0);
	CHECK(gensas_time(sim) == NULL);

	CHECK(gensas_run(sim, stop_after_first, NULL) == GENSAS_CANCELLED);
	CHECK(count == 1);
	CHECK(gensas_sample_count(sim) == 0);

	CHECK(gensas_run(sim, NULL, NULL) == GENSAS_OK);
	n = gensas_sample_count(sim);
	CHECK(n > 10);
	CHECK(gensas_variable_count(sim) == 1);
	CHECK(strcmp(gensas_variable_name(sim, 0), "x") == 0);
	CHECK(gensas_variable_name(sim, 1) == NULL);
	if (n > 0) {
		CHECK(fabs(gensas_time(sim)[n - 1] - 2.0) < 1e-9);
		CHECK(fabs(gensas_variable_values(sim, 0)[n - 1] - exp(-2.0)) < 1e-4);
	}
	CHECK(gensas_event_count(sim) == 0);
	CHECK(gensas_discrete_count(sim) == 0);
	CHECK(gensas_write_result(sim, "xml", "x", err, sizeof err) != 0);
	CHECK(strstr(err, "unknown output format") != NULL);
	gensas_destroy(sim);
	gensas_destroy(NULL);

	if (failures == 0) {
		printf("[  PASSED  ] capi_c_test\n");
		return 0;
	}
	return 1;
}
