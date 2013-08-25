#include "vmgenome.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

/*
 * Register use:
 * 0: current letter
 * 1: 1
 * 2: z
 * A: loopback address*/

gene sample[] = {
		{0, 0,
			{
				0x13, 0x11, 0x21, 0xFB,
				0x22, 0x13, 0x16, 0x2A,
				0x19, 0xD0, 0x20, 0x31,
				0xF0, 0x22, 0x30, 0x9A
			}
		}
};

int main(int argc, char **args)
{
	int out = creat("sample.cgp", S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
	int len = sizeof(sample) / sizeof(gene);
	int i;
	write(out, &len, sizeof(len));
	for(i = 0; i < len; i++) {
		write(out, &(sample[i]), sizeof(gene));
	}
	return 0;
}
