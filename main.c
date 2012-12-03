#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "vmgenome.h"
#include "vm.h"
#include "pareto.h"

int main(int argc, char ** args)
{
	evalset eval;
	eval_closure crit[5];
	genome *r;
	genepool *e = NULL;
	int stop[] = {INT_MAX, 0, 0, INT_MAX, INT_MAX};
	for(int i = 1; i < argc; i++) {
		if(!strcmp(args[i], "-C")) {
			int fd;
			i++;
			fd = open(args[i], O_RDONLY);
			e = load_genepool(fd);
			close(fd);
		}
	}
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
	eval.target = "Taumatawhakatangihangakoauauotamateapokaiwhenuakitanatahu";
	eval.target_len = strlen(eval.target);
	if(e == NULL) {
		e = initial_genepool(50);
	}
	r = selection_loop(e, crit, 5, stop);
}


