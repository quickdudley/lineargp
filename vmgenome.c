#include <stdlib.h>
#include <stdio.h>
#include "vmgenome.h"

void sort_genome(genome * unsorted, int criterion);

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
void sort_genome(genome * unsorted, int criterion) {
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


