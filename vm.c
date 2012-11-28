#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "vm.h"
#include "levenshtein.h"

int * getMem(heapPage ** m, int address) {
	int block = address / 1024;
	int offset = address % 1024;
	while((*m)->prev != NULL && (*m)->prev->base >= block) {
		*m = (*m)->prev;
	}
	while((*m)->next != NULL && (*m)->next->base <= block) {
		*m = (*m)->next;
	}
	if((*m)->prev == NULL && (*m)->base > block) {
		(*m)->prev = (heapPage *)malloc(sizeof(heapPage));
		(*m)->prev->next = (*m);
		(*m)->prev->base = block;
		(*m)->prev->prev = NULL;
		(*m) = (*m)->prev;
	} else if((*m)->next == NULL && (*m)->base < block) {
		(*m)->next = (heapPage *)malloc(sizeof(heapPage));
		(*m)->next->prev = (*m);
		(*m)->next->base = block;
		(*m)->next->next = NULL;
		(*m) = (*m)->next;
	} else if((*m)->base > block) {
		heapPage * n = (heapPage *)malloc(sizeof(heapPage));
		n->base = block;
		n->next = *m;
		n->prev = (*m)->prev;
		n->prev->next = n;
		(*m)->prev = n;
		(*m) = n;
	} else if((*m)->base < block) {
		heapPage * n = (heapPage *)malloc(sizeof(heapPage));
		n->base = block;
		n->prev = *m;
		n->next = (*m)->next;
		n->next->prev = n;
		(*m)->next = n;
		(*m) = n;
	}
	return (*m)->data + offset;
}

unsigned char vmFetch(genome **g, unsigned int *pc) {
	int block = *pc / 16;
	while((*g)->prev != NULL && (*g)->prev->first.execution_position >= block) {
		*g = (*g)->prev;
	}
	while((*g)->next != NULL && (*g)->next->first.execution_position <= block) {
		*g = (*g)->next;
	}
	if((*g)->first.execution_position > block) {
		*pc = block * 16;
		return (*g)->first.instructions[0];
	} else if((*g)->first.execution_position < block) {
		*pc = UINT_MAX;
		return 0x0F; // halt
	} else {
		int x = *pc % 16;
		return (*g)->first.instructions[x];
	}
}

int vmStep(genome *g, environment *env) {
	unsigned char instruction;
	unsigned char h, l;
	instruction = vmFetch(&g, &(env->pc));
	env->pcn = env->pc + 1;
	h = instruction & 0xF0;
	l = instruction & 0x0F;
	if(h == 0) {
		if(l == 0x0F) {
			return 1;
		}
		return 0;
	}
	else if(h == 0x10) {
		env->cb <<= 4;
		env->cb |= l;
		return 0;
	}
	else if(h == 0x20) {
		env->s1 = l;
		return 0;
	}
	else if(h == 0x30) {
		env->s2 = l;
		return 0;
	}
	else if(h == 0x40) {
		env->rgs[l] = env->cb;
		env->cb = 0;
		return 0;
	}
	else if(h == 0x50) {
		env->rgs[l] = *(getMem(&(env->heap), env->rgs[env->s1]));
		return 0;
	}
	else if(h == 0x60) {
		*(getMem(&(env->heap), env->rgs[env->s2])) = env->rgs[l];
		return 0;
	}
	else if(h == 0x70) {
		env->rgs[l] = env->rgs[env->s1];
		return 0;
	}
	else if(h == 0x80 || (h == 0x90 && env->rgs[env->s1] > env->rgs[env->s2])
			  || (h == 0xA0 && env->rgs[env->s1] == env->rgs[env->s2])) {
		env->rgs[0x0F] = env->pcn;
		env->pcn = env->rgs[l];
		env->cb = 0;
		return 0;
	}
	else if(h == 0xB0) {
		env->rgs[env->s1] = env->rgs[l];
		return 0;
	} 
	else if(h == 0xC0) {
		if(env->ii[env->fd].source != NULL) {
			env->rgs[l] = env->ii[env->fd].source(env->ii[env->fd].context);
			if(env->rgs[l] == -1) {
				if(env->ii[env->fd].exhausted) {
					return -1;
				} else {
					env->ii[env->fd].exhausted = 1;
				}
			}
			return 0;
		}
	}
	else if(h == 0xD0) {
		if(env->oo[env->fd].sink != NULL) {
			if(!env->oo[env->fd].closed) {
				int rv = env->oo[env->fd].sink(env->oo[l].context, env->rgs[l]);
				if(env->rgs[l] == -1) {
					env->oo[env->fd].closed = 1;
				}
				return rv;
			}
			return -2;
		}
	}
	else if(h == 0xE0) {
		env->fd = l;
		return 0;
	}
	else if(h == 0xF0) {
		switch(l) {
		case 0x0:
			env->rgs[env->s1] = env->rgs[env->s1] + env->rgs[env->s2];
			break;
		case 0x1:
			env->rgs[env->s1] = env->rgs[env->s1] - env->rgs[env->s2];
			break;
		case 0x2:
			env->rgs[env->s1] = env->rgs[env->s1] * env->rgs[env->s2];
			break;
		case 0x3:
			env->rgs[env->s1] = env->rgs[env->s1] / env->rgs[env->s2];
			break;
		case 0x4:
			env->rgs[env->s1] = env->rgs[env->s1] % env->rgs[env->s2];
			break;
		case 0x5:
			env->rgs[env->s1] = ~(env->rgs[env->s2]);
			break;
		case 0x6:
			env->rgs[env->s1] = env->rgs[env->s1] & env->rgs[env->s2];
			break;
		case 0x7:
			env->rgs[env->s1] = env->rgs[env->s1] | env->rgs[env->s2];
			break;
		case 0x8:
			env->rgs[env->s1] = env->rgs[env->s1] << env->rgs[env->s2];
			break;
		case 0x9:
			env->rgs[env->s1] = env->rgs[env->s1] >> env->rgs[env->s2];
			break;
		case 0xA:
		case 0xC:
		case 0xE:
			env->rgs[env->s1] = 0;
			break;
		case 0xB:
		case 0xD:
			env->rgs[env->s1] = 1;
			break;
		case 0xF:
			env->rgs[env->s1] = 0xFFFFFFFF;
			break;
		}
		return 0;
	}
	return -1;
}

int vmRun(genome *g, environment *env, long long int *steps){
	int penalty;
	int s = *steps;
	sort_execute(g);
	while(s > 0) {
		int result = vmStep(g, env);
		env->pc = env->pcn;
		if(result > 0) {
			s = 0;
			penalty += result - 1;
		} else {
			penalty -= result;
		}
		s--;
		*steps = s;
	}
	return penalty;
}

void init_environment(environment *env)
{
	memset(env, 0, sizeof(environment));
	env->heap = malloc(sizeof(heapPage));
	memset(env->heap, 0, sizeof(heapPage));
}

int delete_heap(heapPage *heap)
{
	while(heap->prev != NULL)
		heap = heap->prev;
	int r = 0;
	heapPage* n;
	while(heap != NULL) {
		n = heap->next;
		free(heap);
		heap = n;
		r++;
	}
}

struct buffercontext {
	char *buffer;
	int length;
	int pos;
};

static int bufferInput(void *context)
{
	if(((struct buffercontext*)context)->pos >= ((struct buffercontext*)context)->length)
		return -1;
	int r = ((struct buffercontext*)context)->buffer[((struct buffercontext*)context)->pos];
	((struct buffercontext*)context)->pos++;
	return r;
}

static int bufferOutput(void *context, int value)
{
	if(value & (~0xff))
		return 1;
	((struct buffercontext*)context)->buffer[((struct buffercontext*)context)->pos] = (char)value;
	((struct buffercontext*)context)->pos++;
	return ((struct buffercontext*)context)->pos == ((struct buffercontext*)context)->length;
}

void eval_run(genome *g, evalset *eval)
{
	if(g == eval->last_genome)
		return;
	eval->last_genome = g;
	environment env;
	struct buffercontext inputbuffer;
	struct buffercontext outputbuffer;
	init_environment(&env);
	inputbuffer.pos = 0;
	inputbuffer.length = eval->input_len;
	inputbuffer.buffer = eval->input;
	outputbuffer.pos = 0;
	outputbuffer.length = (int)(eval->target_len * 1.75) + 1;
	outputbuffer.buffer = malloc(outputbuffer.length);
	env.ii[0].source = bufferInput;
	env.ii[0].context = &inputbuffer;
	env.oo[0].sink = bufferOutput;
	env.oo[0].context = &outputbuffer;
	eval->steps = 200000;
	vmRun(g, &env, &(eval->steps));
	eval->steps = 200000 - eval->steps;
	eval->heap_pages = delete_heap(env.heap);
	eval->difference = bitwiseLevenshtein(outputbuffer.buffer, outputbuffer.pos, eval->target, eval->target_len);
	free(outputbuffer.buffer);
}

int eval_error(genome *g, evalset *eval)
{
	eval_run(g, eval);
	return eval->difference;
}

int eval_runtime(genome *g, evalset *eval)
{
	eval_run(g, eval);
	return eval->steps;
}

int eval_memory(genome *g, evalset *eval)
{
	eval_run(g, eval);
	return eval->heap_pages;
}
