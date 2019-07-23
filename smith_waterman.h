#ifndef _smith_waterman
#define _smith_waterman

#include <stddef.h>

// Simple match / mismatch / gap penalties scores
// No affine gap penalties
// pointers and scores must already have been allocated and
// have dimensions of (l1+1) * (l2+1) (the sequence lengths)
// Note that that pointers table is using ints. That is terrible,
// But it makes it easier to handle within R

void smith_waterman_id(const char *seq_1, const char *seq_2, size_t l1, size_t l2,
		       int match, int mis_match, int gap,
		       int *pointers, int *scores);

void smith_waterman_cm(const char *seq_1, const char *seq_2, size_t l1, size_t l2,
		       int match, int mis_match, int gap,
		       int *scores, int *col_maxes);


#endif
