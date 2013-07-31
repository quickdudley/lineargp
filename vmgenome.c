#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include "vmgenome.h"
#include "levenshtein.h"

#define MUTATION_THRESHOLD (RAND_MAX / 20)
#define SHUFFLE_THRESHOLD (RAND_MAX / 6000)

static void sort_genome(genome * unsorted, int criterion);

void sort_crossover(genome * unsorted) {
	sort_genome(unsorted, 0);
}

void sort_execute(genome * unsorted) {
	sort_genome(unsorted, 1);
}

/* Bottom up merge sort */
static void sort_genome(genome * unsorted, int criterion)
{
	if (criterion != 1 && criterion != 0) {
		return;
	}
	while(unsorted->prev != NULL)
		unsorted = unsorted->prev;
	int len = genome_size(unsorted);
	genome **gl = malloc(sizeof(genome *) * len);
	genome *current = unsorted;
	/* Break genome into list of single gene genomes */
	for(int i = 0; i < len; i++) {
		gl[i] = current;
		current = current->next;
		gl[i]->next = NULL;
		gl[i]->prev = NULL;
	}
	/* Iterate until only one genome remains */
	while(len != 1) {
		int nl = len / 2 + (len & 1);
		genome **ngl = malloc(sizeof(genome *) * nl);
		if(len & 1) {
			ngl[len / 2] = gl[len - 1];
		}
		for(int i = 0; i < (len & ~1); i += 2) {
			genome *a = gl[i], *b = gl[i + 1];
			genome *n = NULL;
			while(a != NULL || b != NULL) {
				genome *t = a == NULL ? b : (b == NULL ? a :
						(criterion == 0 || a->first.execution_position == b->first.execution_position ?
						(a->first.crossover_position <= b->first.crossover_position ? a : b) :
						(a->first.execution_position <= b->first.execution_position ? a : b)));
				if(t == a)
					a = a->next;
				else
					b = b->next;
				if(n != NULL) {
					n->next = t;
					t->prev = n;
				}
				else {
					ngl[i / 2] = t;
				}
				n = t;
			}
		}
		free(gl);
		gl = ngl;
		len = nl;
	}
	free(gl);
}

void delete_genome(genome * scrap)
{
	while(scrap->prev != NULL) {
		scrap = scrap->prev;
	}
	while(scrap != NULL) {
		genome *n = scrap;
		scrap = scrap->next;
		free(n);
	}
}

genome * random_genome(int size)
{
	genome *r = NULL;
	while(0 < size) {
		genome *n, *i, *t;
		n = malloc(sizeof(genome));
		n->next = NULL;
		n->prev = NULL;
		n->first.crossover_position = random() / 2;
		n->first.execution_position = random() / 2;
		for(int x = 0; x < sizeof(n->first.instructions); x++) {
			n->first.instructions[x] = (char)(random() & 0xff);
		}
		i = r;
		while(i != NULL) {
			if(n->first.crossover_position >= i->first.crossover_position)
				n->first.crossover_position++;
			if(n->first.execution_position >= i->first.execution_position)
				n->first.execution_position++;
			t = i;
			i = i->next;
		}
		if(r == NULL) {
			r = n;
		}
		else {
			t->next = n;
			n->prev = t;
		}
		size--;
	}
	return r;
}

genome * copy_genome(genome *original)
{
	if(original == NULL)
		return NULL;
	genome *n, *t, *i;
	n = NULL;
	t = NULL;
	while(original->prev != NULL)
		original = original->prev;
	while(original != NULL) {
		i = malloc(sizeof(genome));
		if(n == NULL)
			n = i;
		memcpy(&(i->first), &(original->first), sizeof(i->first));
		i->next = NULL;
		i->prev = t;
		if(t != NULL)
			t->next = i;
		t = i;
		original = original->next;
	}
	return n;
}

void mutate_genome(genome *x)
{
	while(x->prev != NULL)
		x = x->prev;
	genome *t = x;
	int local_rate = random() % MUTATION_THRESHOLD;
	while(t != NULL) {
		if(random() < SHUFFLE_THRESHOLD) {
			t->first.crossover_position = random() / 2;
			genome *c = x;
			while(c != NULL) {
				if(c != t && t->first.crossover_position >= c->first.crossover_position) {
					t->first.crossover_position++;
				}
				c = c->next;
			}
		}
		if(random() < SHUFFLE_THRESHOLD) {
			t->first.execution_position = random() / 2;
			genome *c = x;
			while(c != NULL) {
				if(c != t && t->first.execution_position >= c->first.execution_position) {
					t->first.execution_position++;
				}
				c = c->next;
			}
		}
		if(random() < (local_rate / 4)) {
			char tmp[sizeof(t->first.instructions)];
			int rs = random() & 0xFF;
			for(int i = 0; i < sizeof(tmp); i++) {
				tmp[(i + rs) % sizeof(tmp)] = t->first.instructions[i];
			}
			for(int i = 0; i < sizeof(tmp); i++) {
				t->first.instructions[i] = tmp[i];
			}
		}
		for(int i = 0; i < sizeof(t->first.instructions); i++) {
			if(random() < local_rate) {
				t->first.instructions[i] = (char)(random() & 0xFF);
			}
		}
		t = t->next;
	}
	// Mutation to add one random gene: mutation to remove one doesn't seem necessary
	// because the crossover operator can remove or duplicate genes.
	if(random() < (local_rate / 4)) {
		t = random_genome(1);
		t->next = x;
		x->prev = t;
		while(x != NULL) {
			if(t->first.execution_position >= x->first.execution_position)
				t->first.execution_position++;
			if(t->first.crossover_position >= x->first.crossover_position)
				t->first.crossover_position++;
			x = x->next;
		}
	}
}

genome * crossover_genome(genome *parent1, genome *parent2)
{
	sort_crossover(parent1);
	sort_crossover(parent2);
	while(parent1->prev != NULL)
		parent1 = parent1->prev;
	while(parent2->prev != NULL)
		parent2 = parent2->prev;
	if(parent1 == parent2) {
		return copy_genome(parent1);
	}
	genome *ret = NULL;
	genome *rt = NULL;
	while(ret == NULL) {
		genome *s1 = parent1;
		genome *s2 = parent2;
		while(s1 != NULL || s2 != NULL) {
			int pp1 = s1 == NULL ? INT_MAX : s1->first.crossover_position;
			int pp2 = s2 == NULL ? INT_MAX : s2->first.crossover_position;
			if(s2 == NULL || pp1 < pp2) {
				if(random() & 1) {
					genome *n = malloc(sizeof(genome));
					n->next = NULL;
					if(rt != NULL)
						rt->next = n;
					n->prev = rt;
					n->first = s1->first;
					rt = n;
					if(ret == NULL)
						ret = n;
				}
				s1 = s1->next;
			}
			else if(s1 == NULL || pp2 < pp1) {
				if(random() & 1) {
					genome *n = malloc(sizeof(genome));
					n->next = NULL;
					if(rt != NULL)
						rt->next = n;
					n->prev = rt;
					n->first = s2->first;
					rt = n;
					if(ret == NULL)
						ret = n;
				}
				s2 = s2->next;
			} else {
				genome *n = malloc(sizeof(genome));
				n->next = NULL;
				if(rt != NULL)
					rt->next = n;
				n->prev = rt;
				n->first = random() & 1 ? s2->first : s1->first;
				if(random() & 1) {
					for(int i = 0; i < sizeof(n->first.instructions); i++) {
						n->first.instructions[i] = random() & 1 ?
								s2->first.instructions[i] :
								s1->first.instructions[i];
					}
				}
				rt = n;
				if(ret == NULL)
					ret = n;
				s1 = s1->next;
				s2 = s2->next;
			}
		}
	}
	return ret;
}

int genome_size(genome * s)
{
	int count = 0;
	while (s->prev != NULL)
		s = s->prev;
	while (s != NULL) {
		count++;
		s = s->next;
	}
	return count;
}

int eval_genome_size(genome *s, void *unused)
{
	return genome_size(s);
}

static inline int isqrt(int o)
{
	int r = 0;
	for(long long int i = 1 << 30; i != 0; i = i >> 1) {
		long long int n = r | i;
		if(n * n <= o)
			r = n;
	}
	return r;
}

int genome_compare(genome *a, genome *b)
{
	while(a->prev != NULL)
		a = a->prev;
	while(b->prev != NULL)
		b = b->prev;
	// Begin using the hungarian algorithm to match the most similar segments
	// between the two genomes.
	int w = genome_size(a);
	{
		int h = genome_size(b);
		if(w < h)
			w = h;
	}
	int ms = sizeof(int) * w * w;
	int max_match = 0;
	int *cm = malloc(ms);
	int *hm = malloc(ms);
	int *lx = calloc(w, sizeof(int));
	int *ly = calloc(w, sizeof(int));
	int *xy = calloc(w, sizeof(int));
	int *yx = calloc(w, sizeof(int));
	int *cw = calloc(w, sizeof(int));
	int *pw = calloc(w, sizeof(int));
	int *sx = calloc(w, sizeof(int));
	int *sy = calloc(w, sizeof(int));
	genome *ai = a;
	genome *bi;
	int x, y;
	for(x = 0; ai != NULL; x++, ai = ai->next) {
		bi = b;
		for(y = 0; bi != NULL; y++, bi = bi->next) {
			int idx = x * w + y;
			cm[idx] = 0;
			if(ai->first.crossover_position != bi->first.crossover_position)
				cm[idx] = 1;
			if(ai->first.execution_position != bi->first.execution_position)
				cm[idx] += 1;
			cm[idx] += bytewiseLevenshtein(ai->first.instructions, sizeof(ai->first.instructions),
					bi->first.instructions, sizeof(bi->first.instructions));
		}
		for(; y < w; y++)
			cm[x * w + y] = 0;
	}
	for(; x < w; x++)
		for(y = 0; y < w; y++)
			cm[x * w + y] = 0;
	memcpy(hm, cm, ms);

	// Cost matrix is set up and copied to hm. Now reduce hm
	for(x = 0; x < w; x++) {
		int min = INT_MAX;
		for(y = 0; y < w; y++) {
			if(hm[x * w + y] < min)
				min = hm[x * w + y];
		}
		for(y = 0; y < w; y++) {
			hm[x * w + y] -= min;
		}
	}
	for(y = 0; y < w; y++) {
		int min = INT_MAX;
		for(x = 0; x < w; x++) {
			if(hm[x * w + y] < min)
				min = hm[x * w + y];
		}
		for(x = 0; x < w; x++) {
			hm[x * w + y] -= min;
		}
	}

	// Setup initial feasible labeling
	memset(ly, 0, sizeof(int) * w);
	for(x = 0; x < w; x++) {
		lx[x] = 0;
		for(y = 0; y < w; y++) {
			if(hm[x * w + y] < lx[x])
				lx[x] = hm[x * w + y];
		}
	}

	// Perform initial greedy matching
	for(x = 0; x < w; x++) {
		xy[x] = -1;
		yx[x] = -1;
	}
	for(x = 0; x < w; x++) {
		for(y = 0; y < w; y++) {
			if(xy[x] == -1 && yx[y] == -1 && hm[x * w + y] - lx[x] - ly[y] == 0) {
				xy[x] = y;
				yx[y] = x;
			}
		}
	}

	int u;

	// Fetch unmatched x;
	for(u = 0; u < w; u++) {
		if(xy[u] == -1) {
			break;
		}
	}

	memset(sx, 0, sizeof(int) * w);
	memset(sy, 0xFF, sizeof(int) * w);

	while(u < w) {
		memset(cw, 0, sizeof(int) * w);
		memset(pw, 0xFF, sizeof(int) * w);
		cw[u] = 1;
		for(y = 0; y < w; y++) {
			sy[y] = hm[u * w + y] - lx[u] - ly[y];
			sx[y] = u;
		}
		while(1) {
			int minsx = -1, minsy = -1;
			int minsv = INT_MAX;
			for(int y = 0; y < w; y++) {
				if(pw[y] == -1 && sy[y] < minsv) {
					minsv = sy[y];
					minsx = sx[y];
					minsy = y;
				}
			}
			if(minsv > 0) {
				for(x = 0; x < w; x++) {
					if(cw[x]) {
						lx[x] += minsv;
					}
				}
				for(y = 0; y < w; y++) {
					if(pw[y] != -1) {
						ly[y] -= minsv;
					}
					else {
						sy[y] -= minsv;
					}
				}
			}
			pw[minsy] = minsx;
			if(yx[minsy] == -1) {
				int cy = minsy;
				int pwk = pw[cy];
				while(1) {
					int tmp = xy[pwk];
					xy[pwk] = cy;
					yx[cy] = pwk;
					cy = tmp;
					if(cy == -1)
						break;
					pwk = pw[cy];
				}
				goto hungarian_reiter;
			}
			else {
				int worker = yx[minsy];
				cw[worker] = 1;
				for(int y = 0; y < w; y++) {
					if(pw[y] == -1) {
						int slack = hm[worker * w + y] - lx[worker] - ly[y];
						if(sy[y] > slack) {
							sy[y] = slack;
							sx[y] = worker;
						}
					}
				}
			}
		}
hungarian_reiter:
		for(u = 0; u < w; u++) {
			if(xy[u] == -1) {
				break;
			}
		}
	}

	// Free memory used only for computing the matching
	free(hm);
	free(lx);
	free(ly);
	free(yx);
	free(cw);
	free(pw);
	free(sx);
	free(sy);

	int r = 0;
	for(x = 0; x < w; x++) {
		r += (2 + sizeof(a->first.instructions)) - cm[x * w + xy[x]];
	}

	free(cm);
	free(xy);

	return r * 100 / w;
}

int genome_equal(genome *a, genome *b) {
	sort_crossover(a);
	sort_crossover(b);
	while(a->prev != NULL)
		a = a->prev;
	while(b->prev != NULL)
		b = b->prev;
	while(a != NULL && b != NULL) {
		if(a->first.crossover_position != b->first.crossover_position)
			return 0;
		if(b->first.execution_position != b->first.execution_position)
			return 0;
		for(int i = 0; i < sizeof(a->first.instructions); i++)
			if(a->first.instructions[i] != b->first.instructions[i])
				return 0;
		a = a->next;
		b = b->next;
	}
	if(a != NULL || b != NULL)
		return 0;
	return 1;
}

int save_genome(int fd, genome *g)
{
	int count;
	while(g->prev != NULL)
		g = g->prev;
	while(g != NULL) {
		count += write(fd, &(g->first), sizeof(g->first)) >= 0;
		g = g->next;
	}
	return count;
}

genome* load_genome(int fd, int size)
{
	gene n;
	genome *r = NULL;
	genome *t = NULL;
	while(size != 0) {
		if(read(fd, &n, sizeof(n)) <= 0)
			size = 0;
		else {
			genome *x = malloc(sizeof(genome));
			memcpy(&(x->first), &n, sizeof(n));
			x->next = NULL;
			x->prev = t;
			if(t != NULL)
				t->next = x;
			t = x;
			if(r == NULL)
				r = t;
		}
		if(size > 0)
			size--;
	}
	return r;
}
