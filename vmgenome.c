#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "vmgenome.h"

#define MUTATION_THRESHOLD (RAND_MAX / 17)

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
						(criterion == 0 ?
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
	while(t != NULL) {
		if(random() < MUTATION_THRESHOLD) {
			t->first.crossover_position = random() / 2;
			genome *c = x;
			while(c != NULL) {
				if(c != t && t->first.crossover_position >= c->first.crossover_position) {
					t->first.crossover_position++;
				}
				c = c->next;
			}
		}
		if(random() < MUTATION_THRESHOLD) {
			t->first.execution_position = random() / 2;
			genome *c = x;
			while(c != NULL) {
				if(c != t && t->first.execution_position >= c->first.execution_position) {
					t->first.execution_position++;
				}
				c = c->next;
			}
		}
		if(random() < MUTATION_THRESHOLD) {
			if(random() & 1) {
				for(int i = 0; i < (sizeof(t->first.instructions) - 1); i++) {
					t->first.instructions[i] = t->first.instructions[i + 1];
				}
				t->first.instructions[sizeof(t->first.instructions) - 1] = (char)(random() & 0xFF);
			}
			else {
				for(int i = sizeof(t->first.instructions) - 1; i > 0; i--) {
					t->first.instructions[i] = t->first.instructions[i - 1];
				}
				t->first.instructions[0] = (char)(random() & 0xFF);
			}
		}
		for(int i = 0; i < sizeof(t->first.instructions); i++) {
			if(random() < MUTATION_THRESHOLD) {
				t->first.instructions[i] = (char)(random() & 0xFF);
			}
		}
		t = t->next;
	}
	// Mutation to add one random gene: mutation to remove one doesn't seem necessary
	// because the crossover operator can remove or duplicate genes.
	if(random() < MUTATION_THRESHOLD && random() < MUTATION_THRESHOLD) {
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

int genome_compare(genome *a, genome *b)
{
	int r = 0;
	int s = genome_size(a);
	sort_execute(a);
	sort_execute(b);
	while(a->prev != NULL)
		a = a->prev;
	while(b->prev != NULL)
		b = b->prev;
	while(a != NULL && b != NULL) {
		if(a->first.execution_position < b->first.execution_position) {
			a = a->next;
		}
		else if(b->first.execution_position < a->first.execution_position) {
			b = b->next;
		}
		else {
			r +=actualLevenshtein(a->first.instructions, sizeof(a->first.instructions), b->first.instructions, sizeof(b->first.instructions));
			a = a->next;
			b = b->next;
		}
	}
	return r * 100 / s;
}
