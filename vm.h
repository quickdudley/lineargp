#ifndef _vm_h
#define _vm_h 1

#include <setjmp.h>

#include "vmgenome.h"

jmp_buf jmpbuf;

typedef struct _heapPage {
	int base;
	struct _heapPage * next;
	struct _heapPage * prev;
	int data[1024];
} heapPage;

typedef int (*inputFunction)(void *);
// return value for outputFunction: see vmStep()
typedef int (*outputFunction)(void *, int);

typedef struct _output {
	outputFunction sink;
	void * context;
	int closed;
} output;

typedef struct _input {
	inputFunction source;
	void * context;
	int exhausted;
} input;

typedef struct _environment {
	input ii[16];
	output oo[16];
	heapPage *heap;
	int rgs[0x10];
	int s1;
	int s2;
	int fd;
	unsigned int pc;
	unsigned int pcn;
} environment;

void init_environment(environment *env);
int delete_heap(heapPage *heap);

//return value for vmStep: 0 for neutral,
// negative for penalty, positive to abort (greater than one
// to both abort and give penalty)
int vmStep(genome **g, environment *env);
int vmRun(genome *g, environment *env, long long int *steps);

typedef struct _evalset {
	genome *last_genome;
	int heap_pages;
	long long int steps;
	int manhattan_error;
	int suffix_error;
	int illegal;
	char *target;
	int target_len;
	char *input;
	int input_len;
} evalset;

extern void eval_run(genome *g, evalset *eval);
extern int eval_runtime(genome *g, evalset *eval);
extern int eval_memory(genome *g, evalset *eval);
extern int eval_manhattan(genome *g, evalset *eval);
extern int eval_suffix(genome *g, evalset *eval);
extern int eval_illegal(genome *g, evalset *eval);

#endif
