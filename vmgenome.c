#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "vmgenome.h"

#define MUTATION_THRESHOLD (RAND_MAX * 0.05)

static void sort_genome(genome * unsorted, int criterion);

typedef struct _genome_list {
	genome * data;
	struct _genome_list * next;
} genome_list;

void sort_crossover(genome * unsorted) {
	sort_genome(unsorted, 0);
}

void sort_execute(genome * unsorted) {
	sort_genome(unsorted, 1);
}

/* Bottom up merge sort */
static void sort_genome(genome * unsorted, int criterion) {
	genome * tail;
	genome_list * fragments;
	genome * fin;
	if (criterion != 1 && criterion != 0) {
		return;
	}
	fin = malloc(sizeof(genome)); // copy the first node to prevent making a circular list when returning
	*fin = *unsorted;
	tail = unsorted->next;
	tail->prev = fin;
	unsorted->next = NULL;
	tail->prev = NULL;
	fragments = malloc(sizeof(genome_list));
	fragments->data = fin;
	fragments->next = NULL;
	/* Break genome into list of single gene genomes */
	while(tail != NULL) {
		genome_list * old_fragments;
		fragments->data->next = NULL;
		fragments->data->prev = NULL;
		old_fragments = fragments;
		fragments = malloc(sizeof(genome_list));
		fragments->data = tail;
		fragments->next = old_fragments;
		tail = tail->next;
	}
	/* Iterate until only one genome remains */
	while(fragments->next != NULL) {
		genome_list * iter = fragments;
		genome_list * del;
		/* Merge every second genome */
		while(iter != NULL) {
			genome * a;
			a = iter->data;
			if(iter->next != NULL) {
				genome * b;
				genome * r;
				genome * t;
				t = NULL;
				b = iter->next->data;
				/* It makes for slightly shorter code to have a blank gene at the start and delete it when finished */
				r = malloc(sizeof(genome));
				t = r;
				while(a != NULL && b != NULL) {
					t->next = malloc(sizeof(genome));
					t->next->prev = t;
					t = t->next;
					if((criterion == 0 && (a->first.crossover_position <= b->first.crossover_position)) ||
					   (criterion == 1 && (a->first.execution_position <= b->first.execution_position))) {
						t->first = a->first;
						if(a->next != NULL) {
							a = a->next;
							free(a->prev);
						} else {
							free(a);
							a = NULL;
						}
					} else {
						t->first = b->first;
						if(b->next != NULL) {
							b = b->next;
							free(b->prev);
						} else {
							free(b);
							b = NULL;
						}
					}
				}
				if(a != NULL) {
					a->prev = t;
					t->next = a;
				} else if (b != NULL) {
					b->prev = t;
					t->next = b;
				}
				t = r->next;
				free(r);
				t->prev = NULL;
				iter->data = t;
				del = iter->next;
				iter->next = del->next;
				free(del);
				iter = iter->next;
			} else {
				iter = NULL;
			}
		}
	}
	*unsorted = *(fragments->data);
	unsorted->next->prev = unsorted;
	free(fragments->data);
	free(fragments);
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
			}
		}
		if(random() < MUTATION_THRESHOLD) {
			t->first.execution_position = random() / 2;
			genome *c = x;
			while(c != NULL) {
				if(c != t && t->first.execution_position >= c->first.execution_position) {
					t->first.execution_position++;
				}
			}
		}
		for(int i = 0; i < sizeof(t->first.instructions); i++) {
			if(random() < MUTATION_THRESHOLD) {
				t->first.instructions[i] = (char)(random() & 0xFF);
			}
		}
		t = t->next;
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
			if(s2 == NULL || s1->first.crossover_position < s2->first.crossover_position) {
				if(random() ^ 1) {
					genome *n = malloc(sizeof(genome));
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
			else if(s1 == NULL || s2->first.crossover_position < s1->first.crossover_position) {
				if(random() ^ 1) {
					genome *n = malloc(sizeof(genome));
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
				if(rt != NULL)
					rt->next = n;
				n->prev = rt;
				n->first = random() ^ 1 ? s2->first : s1->first;
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
