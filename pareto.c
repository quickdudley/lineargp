#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "vmgenome.h"
#include "vm.h"
#include "pareto.h"

void delete_genepool(genepool *scrap)
{
	for(int i = 0; i < scrap->num_candidates; i++) {
		free(scrap->candidates[i].evaluations);
		delete_genome(scrap->candidates[i].genome);
	}
	free(scrap);
}

genepool* pareto_front(genepool *original, genepool **remainder)
{
	int dc = 0;
	int fc = 0;
	genepool *r;
	/* Mark candidates as dominated if another candidate exists which has a lower or
	 * equal cost value for every criterion. */
	for(int i = 0; i < original->num_candidates; i++) {
		original->candidates[i].dominated = 0;
		for(int j = 0; !(original->candidates[i].dominated) && j < original->num_candidates; j++) if(i != j){
			int td = 1;
			int eq = 1;
			for(int c = 0; td && c < original->num_criteria; c++) {
				int a = original->candidates[i].evaluations[c];
				int b = original->candidates[j].evaluations[c];
				/* If a < b then b does not dominate a */
				if(a < b) {
					td = 0;
					break;
				}
				/* if b < a then b will dominate a if there is no other criterion
				 * where a < b */
				if(b < a) {
					eq = 0;
				}
			}
			original->candidates[i].dominated = td && !eq;
		}
		if(original->candidates[i].dominated) {
			dc++;
		}
	}

	fc = original->num_candidates - dc;
	r = malloc(sizeof(genepool) + fc * sizeof(candidate));
	r->num_candidates = fc;
	r->num_criteria = original->num_criteria;
	*remainder = malloc(sizeof(genepool) + dc * sizeof(candidate));
	(*remainder)->num_candidates = dc;
	(*remainder)->num_criteria = original->num_criteria;
	for(int a=0, b=0; a + b < original->num_candidates;) {
		candidate c = original->candidates[a + b];
		if(c.dominated) {
			(*remainder)->candidates[b] = c;
			b++;
		}
		else {
			r->candidates[a] = c;
			a++;
		}
	}
	free(original);
	return r;
}

genepool* concat_genepool(genepool *a, genepool *b)
{
	a = realloc(a, sizeof(genepool) + (a->num_candidates + b->num_candidates) * sizeof(candidate));
	for(int i = a->num_candidates, j = 0; j < b->num_candidates; i++, j++) {
		a->candidates[i] = b->candidates[j];
	}
	a->num_candidates += b->num_candidates;
	free(b);
	return a;
}

void evaluate_pool(genepool *pool, eval_closure* crit, int num_criteria)
{
	pool->num_criteria = num_criteria;
	for(int i = 0; i < pool->num_candidates; i++) {
		if(pool->candidates[i].evaluations != NULL) {
			free(pool->candidates[i].evaluations);
		}
		// first element in evaluation array is age.
		pool->candidates[i].evaluations = malloc(sizeof(int) * (num_criteria + 2));
		pool->candidates[i].evaluations[0] = 0;
		for(int r = 0; r < num_criteria; r++) {
			pool->candidates[i].evaluations[r + 2] = crit[r].func(pool->candidates[i].genome, crit[r].context);
		}
	}
}

genepool* initial_genepool(int size)
{
	genepool *ret = malloc(sizeof(genepool) + sizeof(candidate) * size);
	ret->num_candidates = size;
	for(int i = 0; i < size; i++) {
		ret->candidates[i].evaluations = NULL;
		ret->candidates[i].genome = random_genome(16);
	}
	return ret;
}

genepool* spawn_genepool(genepool* parents, int size)
{
	genepool *ret = malloc(sizeof(genepool) + sizeof(candidate) * size);
	ret->num_candidates = size;
	for(int i = 0; i < size; i++) {
		ret->candidates[i].evaluations = NULL;
		genome *p[2];
		for(int j = 0; j < 2; j++) {
			p[j] = parents->candidates[random() % parents->num_candidates].genome;
		}
		ret->candidates[i].genome = crossover_genome(p[0], p[1]);
		mutate_genome(ret->candidates[i].genome);
	}
	return ret;
}

void pool_age(genepool *pool)
{
	for(int i = 0; i < pool->num_candidates; i++) {
		pool->candidates[i].evaluations[0]++;
	}
}

void pool_mark_similar(genepool *pool) {
	for(int i = 0; i < pool->num_candidates; i++) {
		pool->candidates[i].evaluations[1] = 0;
		for(int j = 0; j < pool->num_candidates; j++) if(i != j) {
			pool->candidates[i].evaluations[1] +=
					genome_compare(pool->candidates[i].genome,
							pool->candidates[j].genome);
		}
	}
}

static inline genome* readystop(genepool *pool, int *stop)
{
	static int minerror = INT_MAX;
	static int mincerror = INT_MAX;
	for(int i = 0; i < pool->num_candidates; i++) {
		if(pool->candidates[i].evaluations[4] < mincerror) {
			mincerror = pool->candidates[i].evaluations[4];
			printf("\nNew minimum error: %d\n", mincerror);
		}
		if(pool->candidates[i].evaluations[3] < minerror) {
			minerror = pool->candidates[i].evaluations[3];
			printf("\nNew minimum bitwise error: %d\n", minerror);
		}
		int a = 1;
		for(int j = 1; j < pool->num_criteria; j++) {
			if(pool->candidates[i].evaluations[j] > stop[j - 1]) {
				a = 0;
				break;
			}
		}
		if(a == 1) {
			return pool->candidates[i].genome;
		}
	}
	return NULL;
}

genome* selection_loop(eval_closure* crit, int num_criteria, int *stop)
{
	int gc = 0;
	genome *r = NULL;
	genepool *m = initial_genepool(500);
	evaluate_pool(m, crit, num_criteria);
	while((r = readystop(m, stop)) == NULL) {
		genepool *n = NULL, *x = NULL;
		pool_mark_similar(m);
		m = pareto_front(m, &x);
		while(m->num_candidates < 300 && x->num_candidates > 0) {
			pool_age(x);
			genepool *f = pareto_front(x, &x);
			m = concat_genepool(m, f);
		}
//		if(m->num_candidates >= 1500)
			pool_age(m);
		printf("Generation %d: %d individuals    \r", gc++, m->num_candidates);
		fflush(stdout);
		n = spawn_genepool(m, 150); //originally I deleted the old one first, but that caused problems for memoization.
		delete_genepool(x);
		evaluate_pool(n, crit, num_criteria);
		m = concat_genepool(m, n);
	}
	return r;
}
