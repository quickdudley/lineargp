#ifndef _pareto_h
#define _pareto_h 1

#include "vmgenome.h"

typedef struct _candidate {
	genome *genome;
	int *evaluations;
	int dominated;
} candidate;

typedef struct _genepool {
	int num_candidates;
	int num_criteria;
	candidate candidates[0];
} genepool;

// Returns the pareto optimal set, puts the remainder into remainder.
// Does not preserve original.
genepool* pareto_front(genepool *original, genepool **remainder);

void delete_genepool(genepool *scrap);

#endif
