#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include "vmgenome.h"
#include "vm.h"
#include "pareto.h"

int poolsize_hover = 200;
int spawn_factor = 1;

void delete_genepool(genepool *scrap)
{
	for(int i = 0; i < scrap->num_candidates; i++) {
		free(scrap->candidates[i].evaluations);
		delete_genome(scrap->candidates[i].genome);
	}
	free(scrap);
}

void age_genepool(genepool *pool)
{
	for(int i = 0; i < pool->num_candidates; i++) {
		pool->candidates[i].age++;
	}
}

genepool* pareto_front(genepool *original, genepool **remainder)
{
	if(original->num_candidates <= poolsize_hover) {
		*remainder = NULL;
		return original;
	}
	int fc = 0;
	int i = 0;
	genepool *r;
	for(int i = 0; i < original->num_candidates; i++) {
		original->candidates[i].dominated = 0;
		original->candidates[i].evaluations[0] = 0;
	}
	/* Bubble sort by pareto dominance */
	for(i = 0; i < original->num_candidates && (i < poolsize_hover || fc == original->candidates[i].dominated); i++) {
		for(int j = i + 1; j < original->num_candidates; j++) {
			int td = 1;
			int eq = 1;
			for(int c = 0; c < original->num_criteria; c++) {
				int a = original->candidates[i].evaluations[c];
				int b = original->candidates[j].evaluations[c];
				/* If a < b then b does not dominate a */
				if(a < b) {
					td = 0;
				}
				/* if b < a then b will dominate a if there is no other criterion
				 * where a < b */
				if(b < a) {
					eq = 0;
				}
			}
			// If two candidates appear the same then:
			// Demote one at random to the next front. This actually prioritises newer candidates
			// because they are compared fewer times than older candidates.
			if(td && eq) {
				if(original->candidates[i].age < original->candidates[j].age ||
						(original->candidates[i].age == original->candidates[j].age && (random() & 1))) {
					td = 0;
					original->candidates[j].evaluations[0]++;
				} else {
					eq = 0;
					original->candidates[i].evaluations[0]++;
				}
			}
			if(td && !eq && original->candidates[i].dominated <= original->candidates[j].dominated)
				original->candidates[i].dominated = original->candidates[j].dominated + 1;
			if(eq && !td && original->candidates[j].dominated <= original->candidates[i].dominated)
				original->candidates[j].dominated = original->candidates[i].dominated + 1;
			if(original->candidates[j].dominated < original->candidates[i].dominated) {
				candidate t = original->candidates[i];
				original->candidates[i] = original->candidates[j];
				original->candidates[j] = t;
				j = i;
			}
			fc = original->candidates[i].dominated;
		}
	}

	fc = original->num_candidates - i;

	r = malloc(sizeof(genepool) + i * sizeof(candidate));
	r->num_candidates = i;
	r->num_criteria = original->num_criteria;
	*remainder = malloc(sizeof(genepool) + fc * sizeof(candidate));
	(*remainder)->num_candidates = fc;
	(*remainder)->num_criteria = original->num_criteria;
	memcpy(r->candidates, original->candidates, sizeof(candidate) * i);
	memcpy((*remainder)->candidates, &(original->candidates[i]), sizeof(candidate) * fc);
	free(original);
	return r;
}

genepool* concat_genepool(genepool *a, genepool *b)
{
	if(a == b)
		return a;
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
		// first element in evaluation array is special.
		pool->candidates[i].evaluations = malloc(sizeof(int) * (num_criteria + 1));
		pool->candidates[i].evaluations[0] = 0;
		for(int r = 0; r < num_criteria; r++) {
			pool->candidates[i].evaluations[r + 1] = crit[r].func(pool->candidates[i].genome, crit[r].context);
		}
	}
}

genepool* initial_genepool(int size, int length)
{
	genepool *ret = malloc(sizeof(genepool) + sizeof(candidate) * size);
	ret->num_candidates = size;
	for(int i = 0; i < size; i++) {
		ret->candidates[i].evaluations = NULL;
		ret->candidates[i].genome = random_genome(length);
	}
	return ret;
}


/*
 * Just selecting individuals at random does eventually get the job done, since the best
 * are always guaranteed to survive, but biasing the selection greatly reduces the number
 * of generations required.
 */
static int biased_random(genepool *source)
{
	// Initial biased index: I don't understand this part but deong from stack overflow's
	// professor does.
	double bias = 1.5;
	double r = (double)random() / (double)RAND_MAX;
	int x = (int)(source->num_candidates * (bias - sqrt(bias*bias -4.0*(bias-1.0)* r)));
	// Now take into account areas of the genepool with equal rankings.
	int l = x;
	int h = x;
	while(l > 0 && source->candidates[l - 1].dominated == source->candidates[x].dominated)
		l--;
	while(h < source->num_candidates && source->candidates[h].dominated == source->candidates[x].dominated)
		h++;
	x = l + random() % (h - l);
	return x;
}

genepool* spawn_genepool(genepool* parents, int size)
{
	genepool *ret = malloc(sizeof(genepool) + sizeof(candidate) * size);
	ret->num_candidates = size;
	for(int i = 0; i < size; i++) {
		ret->candidates[i].evaluations = NULL;
		if(random() & 1) {
			genome *p[2];
			for(int j = 0; j < 2; j++) {
				p[j] = parents->candidates[biased_random(parents)].genome;
			}
			ret->candidates[i].genome = crossover_genome(p[0], p[1]);
		} else {
			ret->candidates[i].genome = copy_genome(parents->candidates[biased_random(parents)].genome);
			mutate_genome(ret->candidates[i].genome);
		}
	}
	return ret;
}

static inline genome* readystop(genepool *pool, int *stop)
{
	static int *minerror = NULL;
	static char *errorname[] = {
			"manhattan difference",
			"error by suffix",
			"heap use",
			"runtime",
			"illegal instructions"
	};
	if(minerror == NULL) {
		int length = pool->num_criteria * sizeof(int);
		minerror = malloc(length);
		for(int i = 0; i < pool->num_criteria; i++) {
			minerror[i] = INT_MAX;
		}
	}
	for(int i = 0; i < pool->num_candidates; i++) {
		for(int em = 0; em < pool->num_criteria; em++) {
			if(pool->candidates[i].evaluations[em + 1] < minerror[em]) {
				minerror[em] = pool->candidates[i].evaluations[em + 1];
				if(em == 0) {
					printf("\nNew minimum genome size: %d\n", minerror[em]);
				} else {
					int batch = (em - 1) / 5 + 1;
					int test = (em - 1) % 5;
					printf("\nNew minimum %s case %d: %d\n", errorname[test], batch, minerror[em]);
				}
			}
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

int save_genepool(int fd, genepool *g)
{
	for(int i = 0; i < g->num_candidates; i++) {
		int count = genome_size(g->candidates[i].genome);
		write(fd, &count, sizeof(count));
		save_genome(fd, g->candidates[i].genome);
	}
}
genepool* load_genepool(int fd)
{
	genepool *r = malloc(sizeof(genepool));
	genepool *a = malloc(sizeof(genepool) + 10 * sizeof(candidate));
	int end = 0;
	r->num_candidates = 0;
	r->num_criteria = 0;
	a->num_candidates = 0;
	a->num_criteria = 0;
	do {
		int size;
		read(fd, &size, sizeof(size));
		a->candidates[a->num_candidates].genome = load_genome(fd, size);
		if(a->candidates[a->num_candidates].genome == NULL) {
			end = 1;
		}
		else {
			a->candidates[a->num_candidates].evaluations = NULL;
			a->num_candidates++;
			if(a->num_candidates == 10) {
				r = concat_genepool(r, a);
				a = malloc(sizeof(genepool) + 10 * sizeof(candidate));
				a->num_candidates = 0;
				a->num_criteria = 0;
			}
		}
	} while(end == 0);
	return concat_genepool(r,a);
}

genome* selection_loop(genepool *m, eval_closure* crit, int num_criteria, int *stop)
{
	int gc = 0;
	genome *r = NULL;
	char pfn[] = "pool99999999.cgp";
	char ofn[] = "pool99999999.cgp";
	evaluate_pool(m, crit, num_criteria);
	while((r = readystop(m, stop)) == NULL) {
		genepool *n = NULL, *x = NULL;
		int ns;
		genepool *xn = NULL;
		ns = m->num_candidates;
		m = pareto_front(m, &xn);
		if(x != NULL) {
			x = concat_genepool(x, xn);
		} else {
			x = xn;
		}
//		while(m->num_candidates < 60 && x->num_candidates > 0) {
//		for(int i = 0; i < 2 && x->num_candidates > 0; i++) {
//			genepool *f = pareto_front(x, &x);
//			m = concat_genepool(m, f);
//		}
		{
			int fd;
			strncpy(ofn, pfn, sizeof(pfn));
			sprintf(pfn, "pool%08d.cgp", gc);
			fd = creat(pfn, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
			save_genepool(fd, m);
			close(fd);
			remove(ofn);
		}
		printf("Generation %d: %d individuals    \r", gc++, m->num_candidates);
		fflush(stdout);
		n = spawn_genepool(m, m->num_candidates * spawn_factor); //originally I deleted the old one first, but that caused problems for memoization.
		if(x != NULL)
			delete_genepool(x);
		evaluate_pool(n, crit, num_criteria);
		age_genepool(m);
		m = concat_genepool(m, n);
	}
	return r;
}
