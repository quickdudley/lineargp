#ifndef _vmgenome
#define _vmgenome 1

typedef struct _gene {
	int crossover_position;
	int execution_position;
	char instructions[16];
} gene;

typedef struct _genome {
	gene first;
	struct _genome * next;
	struct _genome * prev;
} genome;

void delete_genome(genome scrap);

void sort_crossover(genome * unsorted);
void sort_execute(genome * unsorted);

#endif
