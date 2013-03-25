#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
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
	if(original->num_candidates <= 200) {
		*remainder = NULL;
		return original;
	}
	int dc = 0;
	int fc = 0;
	int fl = 1;
	genepool *r;
	for(int i = 0; i < original->num_candidates; i++) {
		original->candidates[i].dominated = 0;
	}
	/* Mark candidates as dominated if another candidate exists which has a lower or
	 * equal cost value for every criterion. */
	for(int i = 0; i < original->num_candidates; i++) {
		for(int j = 0; j < original->num_candidates; j++) if(i != j){
			int td = 1;
			int eq = 1;
			for(int c = 0; td && c < original->num_criteria; c++) {
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
			// I hear it's bad to let 2 candidates take the same spot on the pareto front
			// This hack just demotes one to the next front. LaterI'll try removing duplicates instead
			if(td && eq && original->candidates[i].dominated == original->candidates[j].dominated) {
				if(random() & 1) {
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
			if(original->candidates[i].dominated > fl || original->candidates[j].dominated > fl) {
				i = original->candidates[i].dominated;
				j = original->candidates[j].dominated;
				fl = i > j ? i : j;
				i = -1; // this will be incremented before use
				j = -1;
				break;
			}
		}
	}

	// Heapsort
	for(int i = original->num_candidates - 1; i > 0; i--) {
		int p = (i - 1) / 2;
		if(original->candidates[i].dominated > original->candidates[p].dominated) {
			candidate t = original->candidates[i];
			original->candidates[i] = original->candidates[p];
			original->candidates[p] = t;
		}
	}
	for(int i = original->num_candidates - 1; i > 0; i--) {
		int j = 0;
		candidate t;
		if(original->candidates[0].dominated > original->candidates[i].dominated)
		{
			t = original->candidates[0];
			original->candidates[0] = original->candidates[i];
			original->candidates[i] = t;
		}
		while(j < i - 1) {
			int c = 2 * j + 1;
			if(c < i - 2 && original->candidates[c].dominated < original->candidates[c + 1].dominated) {
				c++;
			}
			if(c < i - 1 && original->candidates[c].dominated > original->candidates[j].dominated) {
				t = original->candidates[c];
				original->candidates[c] = original->candidates[j];
				original->candidates[j] = t;
			} else {
				break;
			}
			j = c;
		}
	}
	fc = 200;
	fl = original->candidates[fc].dominated;
	while(original->candidates[fc + 1].dominated == fl)
		fc++;
	dc = original->num_candidates - fc;

	r = malloc(sizeof(genepool) + fc * sizeof(candidate));
	r->num_candidates = fc;
	r->num_criteria = original->num_criteria;
	*remainder = malloc(sizeof(genepool) + dc * sizeof(candidate));
	(*remainder)->num_candidates = dc;
	(*remainder)->num_criteria = original->num_criteria;
	memcpy(r->candidates, original->candidates, sizeof(candidate) * fc);
	memcpy((*remainder)->candidates, &(original->candidates[fc]), sizeof(candidate) * dc);
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
		// first element in evaluation array is age.
		pool->candidates[i].evaluations = malloc(sizeof(int) * (num_criteria + 1));
		pool->candidates[i].evaluations[0] = 0;
		for(int r = 0; r < num_criteria; r++) {
			pool->candidates[i].evaluations[r + 1] = crit[r].func(pool->candidates[i].genome, crit[r].context);
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
		if(random() & 1) {
			genome *p[2];
			for(int j = 0; j < 2; j++) {
				p[j] = parents->candidates[random() % parents->num_candidates].genome;
			}
			ret->candidates[i].genome = crossover_genome(p[0], p[1]);
		} else {
			ret->candidates[i].genome = copy_genome(parents->candidates[random() % parents->num_candidates].genome);
		}
		mutate_genome(ret->candidates[i].genome);
	}
	return ret;
}

void pool_mark_similar(genepool *pool) {
//	int *med = malloc(sizeof(int) * pool->num_candidates);
	for(int i = 0; i < pool->num_candidates; i++) {
		pool->candidates[i].evaluations[0] = 0;
//		for(int j = 0; j < pool->num_candidates; j++) {
//			med[j] = genome_compare(pool->candidates[i].genome, pool->candidates[j].genome);
//		}
//		for(int x = 0; x < pool->num_candidates; x++) {
//			for(int y = x + 1; y < pool->num_candidates; y++) {
//				if(med[y] > med[x]) {
//					int t = med[y];
//					med[y] = med[x];
//					med[x] = t;
//				}
//			}
//		}
//		int d = pool->num_candidates / 2;
//		if(pool->num_candidates % 2 == 0) {
//			pool->candidates[i].evaluations[0] = (int)(((long long int)(med[d]) + med[d + 1]) / 2);
//		}
//		else {
//			pool->candidates[i].evaluations[0] = med[d];
//		}
	}
//	free(med);
}

static inline genome* readystop(genepool *pool, int *stop)
{
	static int minerror = INT_MAX;
	static int mincerror = INT_MAX;
	for(int i = 0; i < pool->num_candidates; i++) {
		if(pool->candidates[i].evaluations[3] < mincerror) {
			mincerror = pool->candidates[i].evaluations[3];
			printf("\nNew minimum error: %d\n", mincerror);
		}
		if(pool->candidates[i].evaluations[2] < minerror) {
			minerror = pool->candidates[i].evaluations[2];
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
		do {
			genepool *xn = NULL;
			ns = m->num_candidates;
			pool_mark_similar(m);
			m = pareto_front(m, &xn);
			if(x != NULL) {
				x = concat_genepool(x, xn);
			} else {
				x = xn;
			}
		} while(m->num_candidates > 1500 && ns != m->num_candidates);
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
		n = spawn_genepool(m, m->num_candidates >= 50 ? m->num_candidates : 50); //originally I deleted the old one first, but that caused problems for memoization.
		if(x != NULL)
			delete_genepool(x);
		evaluate_pool(n, crit, num_criteria);
		m = concat_genepool(m, n);
	}
	return r;
}
