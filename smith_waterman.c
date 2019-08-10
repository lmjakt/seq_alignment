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
  // keep a record of the position of the maximum score
  // score, row, column
  int max_score = {0, 0, 0};
  
  // Weirdly, I have no difference between using row or column
  // major iteration. Strangely row major seems to work marginally
  // faster. That might change if both sequences were long;
  for(size_t row=1; row < nrow; ++row){
    for(size_t column=1; column < ncol; ++column){
      score_set[0] = scores[ row + nrow * (column-1)] - gap;
      score_set[1] = scores[ row + nrow * column -1 ] - gap;
      score_set[2] = scores[ row - 1 + nrow * (column-1) ] +
	(seq_1[row-1] == seq_2[column-1] ? match : -mis_match);
      which_max( score_set, max_set );
      scores[ row + nrow * column ] = max_set[1];
      pointers[ row + nrow * column] = max_set[0];
      if(max_score[0] < max_set[1]){
	max_score[0] = max_set[1];
	max_score[1] = row;
	max_score[2] = column;
      }
    }
  }
  // and that should be enough..
}

// A function which simply obtains column maxes.. 
// col_maxes must have been allocated with l2 + 1 members
void smith_waterman_cm(const char *seq_1, const char *seq_2, size_t l1, size_t l2,
		       int match, int mis_match, int gap,
		       int *scores, int *col_maxes)
{
  if(!l1 || !l2)
    return;
  size_t nrow=l1 + 1;
  size_t ncol=l2 + 1;
  memset((void*)scores, 0, sizeof(int) * nrow * ncol);
  // memset across the whole array isn't necessary so we can try
  memset((void*)col_maxes, 0, sizeof(int) * ncol);
  size_t offset = 0;
  register int score = 0;
  register int column_max = 0;
  for(size_t column=1; column < ncol; ++column){
    column_max = 0;
    for(size_t row=1; row < nrow; ++row){
      offset = row + nrow * column;
      scores[offset] = 0;
      // get the left score
      score = scores[ offset - nrow ] - gap;
      scores[offset] = score > scores[offset] ? score : scores[offset];
      // get the right score
      score = scores[ offset -1 ] - gap;
      scores[offset] = score > scores[offset] ? score : scores[offset];
      // and the diagonal
      score = scores[ offset - (1 + nrow) ] +
	(seq_1[row-1] == seq_2[column-1] ? match : -mis_match);
      scores[offset] = score > scores[offset] ? score : scores[offset];
      //    col_maxes[column] = col_maxes[column] < scores[offset] ? scores[offset] : col_maxes[column];
      column_max = column_max < scores[offset] ? scores[offset] : column_max;
    }
    col_maxes[column] = column_max;
  }
}
