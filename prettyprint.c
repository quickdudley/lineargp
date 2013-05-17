#include "vmgenome.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **args)
{
	int r = open(args[1], O_RDONLY);
	int l;
	gene c;
	if(argc != 2) {
		fprintf(stderr, "usage: prettyprint file\n");
		return 1;
	}
	while(read(r, &l, sizeof(l)) == sizeof(l)) {
		printf("{\n");
		for(int g = 0; g < l; g++) {
			memset(&c, 0, sizeof(c));
			read(r, &c, sizeof(c));
			printf("\t\t{%d, %d\n", c.crossover_position, c.execution_position);
			printf("\t\t\t{\n");
			for(int i = 0; i < 16; i++) {
				if(i > 0) {
					printf(",");
				}
				if(i % 4 == 0) {
					if(i != 0) {
						printf("\n");
					}
					printf("\t\t\t\t");
				} else {
					printf(" ");
				}
				printf("0x%02X", c.instructions[i]);
			}
			printf("\n\t\t\t}\n\t\t}\n");
		}
		printf("};\n\n");
	}
	return 0;
}
