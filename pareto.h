#ifndef _pareto_h
#define _pareto_h 1

#include "vmgenome.h"
#include "vm.h"

typedef struct _candidate {
	genome *genome;
	int *evaluations; //first the crowding penalty, followed by results of evaluation functions
	double dominated;
	int age;
} candidate;

typedef struct _genepool {
	int num_candidates;
	int num_criteria;
	candidate candidates[];
} genepool;

typedef int (*eval_func)(genome *g, void *context);

typedef struct _eval_closure {
	eval_func func;
	void *context;
	double pressure;
	int stop;
	char *label;
} eval_closure;

extern int poolsize_hover;
extern int spawn_factor;

genome* selection_loop(genepool *m, eval_closure* crit, int num_criteria);
void evaluate_pool(genepool *pool, eval_closure* crit, int num_criteria);
void pool_age(genepool *pool);

// Returns the pareto optimal set, puts the remainder into remainder.
// Does not preserve original.
genepool* pareto_front(genepool *original, genepool **remainder, int resultsize, eval_closure *conditions);

genepool* initial_genepool(int size, int length);
genepool* spawn_genepool(genepool* parents, int size);

void delete_genepool(genepool *scrap);
genepool* concat_genepool(genepool *a, genepool *b);

void save_genepool(int fd, genepool *g);
genepool* load_genepool(int fd);

#endif
