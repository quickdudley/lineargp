#include <stdio.h>

#include "vmgenome.h"
#include "vm.h"
#include "pareto.h"

int main(int argc, char ** args)
{
	evalset eval;
	eval_closure crit[4];
	genome *r;
	int stop[] = {200, 0, 30, 200};
	crit[0].func = eval_genome_size;
	crit[1].func = (eval_func)eval_error;
	crit[2].func = (eval_func)eval_memory;
	crit[3].func = (eval_func)eval_runtime;
	for(int i = 1; i < 4; i++) {
		crit[i].context = (void*)&eval;
	}
	eval.input = "";
	eval.input_len = 0;
	eval.target = "Hello World!\n";
	eval.target_len = 13;
	r = selection_loop(crit, 4, stop);
}


