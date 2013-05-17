#ifndef _vmgenome
#define _vmgenome 1

typedef struct _gene {
	int crossover_position;
	int execution_position;
	unsigned char instructions[16];
} gene;

typedef struct _genome {
	gene first;
	struct _genome * next;
	struct _genome * prev;
} genome;

void delete_genome(genome * scrap);

int genome_size(genome * s);
int eval_genome_size(genome *s, void *unused);

// genome_compare() returns a high value when comparing similar genomes,
// low value when comparing more different ones. It was hard to write it
// so that it wasn't biased by genome length
int genome_compare(genome *a, genome *b);
int genome_equal(genome *a, genome *b);

void sort_crossover(genome * unsorted);
void sort_execute(genome * unsorted);

genome * random_genome(int size);
genome * copy_genome(genome *original);
void mutate_genome(genome *x);
genome * crossover_genome(genome *parent1, genome *parent2);

int save_genome(int fd, genome *g);
genome* load_genome(int fd, int size);

#endif
