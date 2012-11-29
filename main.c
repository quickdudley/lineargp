#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>

#include "vmgenome.h"
#include "vm.h"
#include "pareto.h"

int main(int argc, char ** args)
{
	evalset eval;
	eval_closure crit[5];
	genome *r;
	int stop[] = {INT_MAX, 0, 0, INT_MAX, INT_MAX};
	srand(time(0));
	memset(&eval, 0, sizeof(eval));
	crit[0].func = eval_genome_size;
	crit[1].func = (eval_func)eval_error;
	crit[2].func = (eval_func)eval_cerror;
	crit[3].func = (eval_func)eval_memory;
	crit[4].func = (eval_func)eval_runtime;
	for(int i = 1; i < 5; i++) {
		crit[i].context = (void*)&eval;
	}
	eval.input = "";
	eval.input_len = 0;
	eval.target = "Hello";
	eval.target_len = 5;
	r = selection_loop(crit, 5, stop);
}


