#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "vm.h"
#include "pareto.h"

int stdInput(void *file)
{
	int r = getc((FILE *)file);
	if(!feof(file))
		return r;
	else
		return -1;
}

int stdOutput(void *file, int c)
{
	if(c <= 0xFF) {
		putc(c, (FILE *)file);
		return 0;
	} else {
		return 1;
	}
}

int main(int argc, char **args)
{
	if(argc < 2) {
		fprintf(stderr, "usage: runvm file [index]\n");
		return 1;
	}
	genepool *e;
	environment env;
	int fd;
	int ri = 0;
	int si;
	long long int steps;
	fd = open(args[1], O_RDONLY);
	e = load_genepool(fd);
	close(fd);
	if(argc >= 3 && sscanf(args[2], "%d", &si) && si > 0 && si < e->num_candidates) {
		ri = si;
	}
	init_environment(&env);
	env.ii[0].source = stdInput;
	env.ii[0].context = stdin;
	env.oo[0].sink = stdOutput;
	env.oo[0].context = stdout;
	steps = 200000;
	vmRun(e->candidates[ri].genome, &env, &steps);
}
