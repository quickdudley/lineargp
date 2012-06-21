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
			return 0;
		}
	}
	else if(h == 0xD0) {
		if(env->oo[env->fd].sink != NULL) {
			env->oo[env->fd].sink(env->oo[l].context, env->rgs[l]);
			return 0;
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


