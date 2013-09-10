#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#include "vmgenome.h"
#include "vm.h"
#include "pareto.h"

int main(int argc, char ** args)
{
	evalset *eval;
	eval_closure *crit;
	int batches = 0;
	int batchn = 0;
	genome *r;
	genepool *e = NULL;
	for(int i = 1; i < argc; i++) {
		if(!strcasecmp(args[i], "-io")) {
			batches++;
			i += 2;
		}
	}
	crit = malloc((1 + 5 * batches) * sizeof(eval_closure));
	eval = malloc(5 * batches * sizeof(evalset));
	memset(eval, 0, 5 * batches * sizeof(evalset));
	crit[0].stop = INT_MAX;
	for(int i = 0; i < 5 * batches; i++) {
		if(i % 5 < 2) {
			crit[i + 1].stop = 0;
		} else {
			crit[i + 1].stop = INT_MAX;
		}
	}
	for(int i = 1; i < argc; i++) {
		if(!strcmp(args[i], "-C")) {
			int fd;
			genepool *t;
			i++;
			fd = open(args[i], O_RDONLY);
			t = load_genepool(fd);
			if(e) {
				e = concat_genepool(e, t);
			}
			else {
				e = t;
			}
			close(fd);
		}
		else if(!strcmp(args[i], "-R")) {
			int s = 1, l = 1;
			int success;
			genepool *t;
			i++;
			success = sscanf(args[i], "%d", &s);
			i++;
			success = success && sscanf(args[i], "%d", &l);
			if(success) {
				t = initial_genepool(s, l);
				if(e) {
					e = concat_genepool(e, t);
				}
				else {
					e = t;
				}
			}
		}
		else if(!strcmp(args[i], "-f")) {
			int success, value;
			i++;
			success = sscanf(args[i], "%d", &value);
			if(success) {
				spawn_factor = value;
			}
		}
		else if(!strcasecmp(args[i], "-io")) {
			int imap = isupper(args[i][1]);
			int omap = isupper(args[i][2]);
			i++;
			if(!imap) {
				eval[batchn].input = args[i];
				eval[batchn].input_len = strlen(args[i]);
			} else {
				struct stat stats;
				int fd = open(args[i], O_RDONLY);
				void *content;
				fstat(fd, &stats);
				content = mmap(NULL, (size_t)stats.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
				eval[batchn].input = content;
				eval[batchn].input_len = (int)stats.st_size;
			}
			i++;
			if(!omap) {
				eval[batchn].target = args[i];
				eval[batchn].target_len = strlen(args[i]);
			} else {
				struct stat stats;
				int fd = open(args[i], O_RDONLY);
				void *content;
				fstat(fd, &stats);
				content = mmap(NULL, (size_t)stats.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
				eval[batchn].target = content;
				eval[batchn].target_len = (int)stats.st_size;
			}
			batchn++;
		}
	}
	srand(time(0));
	crit[0].func = eval_genome_size;
	crit[0].pressure = 1.0 / 128;
	crit[0].label = "genome size";
	for(int i = 1; i < batches * 5; i += 5) {
		int b = (i - 1) / 5 + 1;
		crit[i].func = (eval_func)eval_manhattan;
		crit[i].pressure = 1;
		crit[i].label = malloc(38);
		sprintf(crit[i].label, "manhattan distance: batch %d", b);
		crit[i + 1].func = (eval_func)eval_suffix;
		crit[i + 1].pressure = 1;
		crit[i + 1].label = malloc(38);
		sprintf(crit[i + 1].label, "error by suffix: batch %d", b);
		crit[i + 2].func = (eval_func)eval_memory;
		crit[i + 2].pressure = 0.1;
		crit[i + 2].label = malloc(38);
		sprintf(crit[i + 2].label, "memory blocks used: batch %d", b);
		crit[i + 3].func = (eval_func)eval_runtime;
		crit[i + 3].pressure = 0.1;
		crit[i + 3].label = malloc(38);
		sprintf(crit[i + 3].label, "runtime: batch %d", b);
		crit[i + 4].func = (eval_func)eval_illegal;
		crit[i + 4].pressure = 0.3;
		crit[i + 4].label = malloc(38);
		sprintf(crit[i + 4].label, "illegal instructions: batch %d", b);
	}
	for(int i = 1; i < 6 * batches; i++) {
		crit[i].context = (void*)&(eval[(i - 1) / 5]);
	}
	if(e == NULL) {
		e = initial_genepool(50, 16);
	}
	poolsize_hover = e->num_candidates;
	r = selection_loop(e, crit, 1 + 5 * batches);
	int len = genome_size(r);
	int of = creat("finished.cgp", S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
	write(of, &len, sizeof(len));
	save_genome(of, r);
}


