#include <limits.h>
#include <stdlib.h>

#include "vm.h"

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
		return 0;
	} else {
		int x = *pc % 16;
		return (*g)->first.instructions[x];
	}
}

int vmStep(genome *g, environment *env) {
	unsigned char instruction;
	unsigned char h, l;
	instruction = vmFetch(&g, &(env->pc));
	h = instruction & 0xF0;
	l = instruction & 0x0F;
	if(h == 0) {
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
	}
	else if(h == 0x40) {
		env->rgs[l] = env->cb;
		env->cb = 0;
	}
	else if(h == 0x50) {
		env->rgs[l] = *(getMem(&(env->heap), env->rgs[env->s1]));
	}
	else if(h == 0x60) {
		*(getMem(&(env->heap), env->rgs[env->s2])) = env->rgs[l];
	}
	else if(h == 0x70) {
		env->rgs[l] = env->rgs[env->s1];
	}
	else if(h == 0x80 || (h == 0x90 && env->rgs[env->s1] > env->rgs[env->s2])) {
		env->rgs[0x0F] = env->pcn;
		env->pcn = env->rgs[l];
		env->cb = 0;
	}
	else if(h == 0xA0) {
		env->rgs[env->s1] = env->rgs[l];
	}
	return -1;
}


