#include <stdint.h>
#include <stdlib.h>

#include "test_general.h"
#include "test_graph.h"
#include "test_entropy.h"
#include "test_equivalence.h"

#ifdef TEST_ALL
#define TEST_GRAPH_RELATIVE_PROPORTION
#define TEST_ENTROPY_SHANNON_WEAVER
#define TEST_ENTROPY_RENYI
#define TEST_ENTROPY_PATIL_TAILLIE
#define TEST_ENTROPY_Q_LOGARITHMIC
#define TEST_EQUIVALENCE_ENTROPY
#endif

static int32_t num_calls_info;
static int32_t num_calls_warning;
static int32_t num_calls_error;

struct test_case cases[] = {
	#ifdef TEST_GRAPH_RELATIVE_PROPORTION
	{test_compute_graph_relative_proportions, 0},
	#endif
	#ifdef TEST_ENTROPY_SHANNON_WEAVER
	{test_shannon_weaver_entropy, 0},
	#endif
	#ifdef TEST_ENTROPY_RENYI
	{test_renyi_entropy, 0},
	#endif
	#ifdef TEST_ENTROPY_PATIL_TAILLIE
	{test_patil_taillie_entropy, 0},
	#endif
	#ifdef TEST_ENTROPY_Q_LOGARITHMIC
	{test_q_logarithmic_entropy, 0},
	#endif
	#ifdef TEST_EQUIVALENCE_ENTROPY
	{test_equivalence_entropy, 0},
	#endif
};

int32_t main(void){
	int32_t n = sizeof(cases) / sizeof(struct test_case);
	int32_t num_failures = 0;
	for(int32_t i = 0 ; i < n ; i++){
		int32_t local_result = cases[i].function();
		if(!cases[i].can_fail){
			num_failures += (local_result != 0);
		}
	}

	printf("INFO: %i\n", num_calls_info);
	printf("WARNING: %i\n", num_calls_warning);
	printf("ERROR: %i\n", num_calls_error);
	printf("num_failures: %i\n", num_failures);

	return num_failures > 0;
}
