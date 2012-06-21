#ifndef _vm_h
#define _vm_h 1
#include "vmgenome.h"

typedef struct _heapPage {
	int base;
	struct _heapPage * next;
	struct _heapPage * prev;
	int data[1024];
} heapPage;

typedef unsigned int (*inputFunction)(void *);
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
	int cb;
	int s1;
	int s2;
	int fd;
	unsigned int pc;
	unsigned int pcn;
} environment;

typedef environment * (*setupEnvironmentFunction)(environment *);
typedef void (*takedownEnvironmentFunction)(environment *);

typedef struct _environmentTemplate {
	setupEnvironmentFunction setup;
	takedownEnvironmentFunction takedown;
} environmentTemplate;

//return value for vmStep: 0 for neutral,
// negative for penalty, positive to abort (greater than one
// to both abort and give penalty)
int vmStep(genome *g, environment *env);
int vmRun(genome *g, environment *env, long long int steps);

#endif
