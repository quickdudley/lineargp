#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "levenshtein.h"

int bitwiseLevenshtein(char *a, size_t asize, char *b, size_t bsize)
{
	if(asize > bsize) {
		char *t = a;
		size_t tsize = asize;
		a = b;
		asize = bsize;
		b = t;
		bsize = tsize;
	}
	int width = asize * 8;
	int height = bsize * 8;
	if(width == 0) {
		return height;
	}
	int *m[2];

	for(int i = 0; i < 2; i++) {
		m[i] = malloc(width * sizeof(int));
	}
	for(int y = 0; y < height; y++) {
		for(int x = 0; x < width; x++) {
			int od, ou, ol;
			int cost, p, c;
			od = y == 0 ? x : (x == 0 ? y : m[y & 1][x - 1]);
			ou = y == 0 ? x + 1 : m[y & 1][x];
			ol = x == 0 ? y + 1 : m[!(y & 1)][x - 1];
			cost = (((a[x / 8]) >> (8 - x % 8)) & 1) != (((b[y / 8]) >> (8 - y % 8)) & 1);
			p = od + cost;
			ou++;
			ol++;
			if(p > ou)
				p = ou;
			if(p > ol)
				p = ol;
			m[!(y & 1)][x] = p;
		}
	}
	width = m[height & 1][width - 1]; // reusing the variable
	for(int i = 0; i < 2; i++) {
		free(m[i]);
	}
	return width;
}

// Modified to make a better fitness function: the edit operations are now:
// increment, decrement, insert 0, remove 0
int bytewiseLevenshtein(char *a, size_t asize, char *b, size_t bsize)
{
	if(asize > bsize) {
		char *t = a;
		size_t tsize = asize;
		a = b;
		asize = bsize;
		b = t;
		bsize = tsize;
	}
	if(asize == 0) {
		return bsize * 255;
	}
	int *m[2];

	for(int i = 0; i < 2; i++) {
		m[i] = malloc(asize * sizeof(int));
	}
	for(int y = 0; y < bsize; y++) {
		for(int x = 0; x < asize; x++) {
			int od, ou, ol;
			int cost, p, c;
			od = y == 0 ? x * 255 : (x == 0 ? y * 255 : m[y & 1][x - 1]);
			ou = y == 0 ? (x + 1) * 255 : m[y & 1][x];
			ol = x == 0 ? (y + 1) * 255 : m[!(y & 1)][x - 1];
			cost = a[x] - b[y];
			cost = cost >= 0 ? cost : -cost;
			p = od + cost;
			ou += 255;
			ol += 255;
			if(p > ou)
				p = ou;
			if(p > ol)
				p = ol;
			m[!(y & 1)][x] = p;
		}
	}
	asize = m[bsize & 1][asize - 1]; // reusing the variable
	for(int i = 0; i < 2; i++) {
		free(m[i]);
	}
	if(asize <= 12) {
		int x = 0;
	}
	return asize;
}

int actualLevenshtein(char *a, size_t asize, char *b, size_t bsize)
{
	if(asize > bsize) {
		char *t = a;
		size_t tsize = asize;
		a = b;
		asize = bsize;
		b = t;
		bsize = tsize;
	}
	if(asize == 0) {
		return bsize;
	}
	int *m[2];

	for(int i = 0; i < 2; i++) {
		m[i] = malloc(asize * sizeof(int));
	}
	for(int y = 0; y < bsize; y++) {
		for(int x = 0; x < asize; x++) {
			int od, ou, ol;
			int cost, p, c;
			od = y == 0 ? x : (x == 0 ? y : m[y & 1][x - 1]);
			ou = y == 0 ? x + 1 : m[y & 1][x];
			ol = x == 0 ? y + 1 : m[!(y & 1)][x - 1];
			cost = a[x] == b[y];
			p = od + cost;
			ou += 1;
			ol += 1;
			if(p > ou)
				p = ou;
			if(p > ol)
				p = ol;
			m[!(y & 1)][x] = p;
		}
	}
	asize = m[bsize & 1][asize - 1]; // reusing the variable
	for(int i = 0; i < 2; i++) {
		free(m[i]);
	}
	if(asize <= 12) {
		int x = 0;
	}
	return asize;
}

/* Main function for testing */
//int main(int argc, char **args)
//{
//	int r;
//	if(argc != 3) {
//		fprintf(stderr, "Need 2 arguments to calculate Levenshtein distance\n");
//		return -1;
//	}
//	r = bitwiseLevenshtein(args[1], strlen(args[1]), args[2], strlen(args[2]));
//	printf("Levenshtein difference: %d\n", r);
//	return 0;
//}
