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
typedef unsigned int (*outputFunction)(void *, int);

typedef struct _output {
	outputFunction sink;
	void * context;
} output;

typedef struct _input {
	inputFunction source;
	void * context;
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

typedef environment (*setupEnvironmentFunction)(void);
typedef void (*takedownEnvironmentFunction)(environment);

typedef struct _environmentTemplate {
	setupEnvironmentFunction setup;
	takedownEnvironmentFunction takedown;
} environmentTemplate;

int vmStep(genome *g, environment *env);

#endif
