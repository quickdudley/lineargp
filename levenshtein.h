#ifndef _levenshtein_h
#define _levenshtein_h 1
#include <stddef.h>

extern int bitwiseLevenshtein(char *a, size_t asize, char *b, size_t bsize);
extern int bytewiseLevenshtein(char *a, size_t asize, char *b, size_t bsize);
extern int errorfreeprogress(char *a, size_t asize, char *b, size_t bsize);

#endif
