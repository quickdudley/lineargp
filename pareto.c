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


genepool* filter_genepool(genepool *original, int *filter, genepool **remainder)
{
	int num_accept = 0;
	genepool *r;
	for(int i = 0; i < original->num_candidates; i++) {
		original->candidates[i].dominated = 0;
		for(int j = 0; j < original->num_criteria; j++) {
			if(original->candidates[i].evaluations[j] > filter[j]) {
				original->candidates[i].dominated = 1;
				break;
			}
		}
		if(original->candidates[i].dominated == 0) {
			num_accept++;
		}
	}
	if(num_accept == original->num_candidates) {
		*remainder = NULL;
		return original;
	}
	r = malloc(sizeof(genepool) + num_accept * sizeof(candidate));
	*remainder = malloc(sizeof(genepool) + (original->num_candidates - num_accept) * sizeof(candidate));
	r->num_candidates = num_accept;
	(*remainder)->num_candidates = original->num_candidates - num_accept;
	r->num_criteria = original->num_criteria;
	int ri = 0;
	int mi = 0;
	for(int i = 0; i < original->num_candidates; i++) {
		if(original->candidates[i].dominated == 0) {
			r->candidates[ri] = original->candidates[i];
			ri++;
		} else {
			(*remainder)->candidates[mi] = original->candidates[i];
			mi++;
		}
	}
	free(original);
	return r;
}

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


genepool* concat_genepool(genepool *a, genepool *b)
{
	if(a == b || b == NULL)
		return a;
	if(a == NULL)
		return b;
	a = realloc(a, sizeof(genepool) + (a->num_candidates + b->num_candidates) * sizeof(candidate));
	for(int i = a->num_candidates, j = 0; j < b->num_candidates; i++, j++) {
		a->candidates[i] = b->candidates[j];
	}
	a->num_candidates += b->num_candidates;
	free(b);
	return a;
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

genepool* spawn_genepool(genepool* parents, int size)
{
	genepool *ret = malloc(sizeof(genepool) + sizeof(candidate) * size);
	ret->num_candidates = size;
	for(int i = 0; i < size; i++) {
		ret->candidates[i].evaluations = NULL;
		if(!(random() & 3)) {
			genome *p[2];
			for(int j = 0; j < 2; j++) {
				p[j] = parents->candidates[random() % parents->num_candidates].genome;
			}
			ret->candidates[i].genome = crossover_genome(p[0], p[1]);
			if(random() & 1) {
				mutate_genome(ret->candidates[i].genome);
			}
		} else {
			ret->candidates[i].genome = copy_genome(parents->candidates[random() % parents->num_candidates].genome);
			mutate_genome(ret->candidates[i].genome);
		}
	}
	return ret;
}

void save_genepool(int fd, genepool *g)
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


