#include "smith_waterman.h"
#include <string.h>
#include <R.h>
#include <Rinternals.h>


// it would be better to make max_set a struct
// as that would be much clearer to read in the code
void which_max(int score_set[3], int max_set[2]){
  memset((void*)max_set, 0, sizeof(int) * 2);
  for(size_t i=0; i < 3; ++i){
    if(max_set[1] < score_set[i]){
      max_set[0] = i+1;
      max_set[1] = score_set[i];
    }
  }
}

// note that match, mis_match, and gap should all be positive
// but that this is left to the caller
// l1 and l2 are the lengths of the sequences. Note that the scoring
// matrix is larger than this.. 
void smith_waterman_id(const char *seq_1, const char *seq_2, size_t l1, size_t l2,
		       int match, int mis_match, int gap,
		       int *pointers, int *scores)
{
  if(!l1 || !l2)
    return;
  size_t nrow=l1 + 1;
  size_t ncol=l2 + 1;
  memset((void*)pointers, 0, sizeof(int) * nrow * ncol);
  memset((void*)scores, 0, sizeof(int) * nrow * ncol);
  
  // left, down, diagonal --> 1, 2, 3; 0 if not set.. 
  int score_set[3] = {0, 0, 0};
  // the offset and the value
  int max_set[2] = {0, 0};

  // fill in the first row
  // Note that the following would be better done in the opposite orientation.
  for(size_t row=1; row < nrow; ++row){
    for(size_t column=1; column < ncol; ++column){
      score_set[0] = scores[ row + nrow * (column-1)] - gap;
      score_set[1] = scores[ row + nrow * column -1 ] - gap;
      score_set[2] = scores[ row - 1 + nrow * (column-1) ] +
	(seq_1[row-1] == seq_2[column-1] ? match : -mis_match);
      /* Rprintf("row: %d column %d\tscores %d, %d, %d  prev.score %d  s1: %c  s2: %c\n", row, column, */
      /* 	      score_set[0], score_set[1], score_set[2], scores[ row - 1 + nrow * (column-1) ], */
      /* 	      seq_1[row-1], seq_2[column-1]); */
      which_max( score_set, max_set );
      scores[ row + nrow * column ] = max_set[1];
      pointers[ row + nrow * column] = max_set[0];
    }
  }
  // and that should be enough..
}
