#include <stdlib.h>
#include "vmgenome.h"
#include "pareto.h"

void delete_genepool(genepool *scrap)
{
	for(int i = 0; i < scrap->num_candidates; i++) {
		free(scrap->candidates[i].evaluations);
		delete_genome(scrap->candidates[i].genome);
	}
	free(scrap);
}

genepool* pareto_front(genepool *original, genepool **remainder)
{
	int dc = 0;
	int fc = 0;
	genepool *r;
	/* Mark candidates as dominated if another candidate exists which has a lower or
	 * equal cost value for every criterion. */
	for(int i = 0; i < original->num_candidates; i++) {
		original->candidates[i].dominated = 0;
		for(int j = 0; !(original->candidates[i].dominated) && j < original->num_candidates; j++) if(i != j){
			int td = 1;
			int eq = 1;
			for(int c = 0; td && c < original->num_criteria; c++) {
				int a = original->candidates[i].evaluations[c];
				int b = original->candidates[j].evaluations[c];
				/* If a < b then b does not dominate a */
				if(a < b) {
					td = 0;
					break;
				}
				/* if b < a then b will dominate a if there is no other criterion
				 * where a < b */
				if(b < a) {
					eq = 0;
				}
			}
			original->candidates[i].dominated = td || !eq;
		}
		if(original->candidates[i].dominated) {
			dc++;
		}
	}

	fc = original->num_candidates - dc;
	r = malloc(sizeof(genepool) + fc * sizeof(candidate));
	r->num_candidates = fc;
	r->num_criteria = original->num_criteria;
	*remainder = malloc(sizeof(genepool) + dc * sizeof(candidate));
	(*remainder)->num_candidates = dc;
	(*remainder)->num_criteria = original->num_criteria;
	for(int a=0, b=0; a + b < original->num_candidates;) {
		candidate c = original->candidates[a + b];
		if(c.dominated) {
			(*remainder)->candidates[b] = c;
			b++;
		}
		else {
			r->candidates[a] = c;
			a++;
		}
	}
	free(original);
	return r;
}
