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
  int max_score[3] = {0, 0, 0};
  
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


// As above, but optimised for memory usage. Rather than remember the full table we only remember
// the last column; If nrow >> ncol
// This means that we can scan the complet genome relatively quickly
// Again we do not allocate memory here, but leave that to the calling function
// In order for the calling function to be able to optimise memory allocation..
// Note that it would be relatively easy to modify this function to handle affine gap penalties
// The idea of this function is to use it to scan for tandem repeats. If we give a tandem
// repeat sequence as the short sequence where they occur should give rise to a periodic
// maximum signal, with a periodicity of the length of the tandem repeat. This can then
// be used combined with a running sum to find regions rich for such repeats in a manner
// similar to how you can find CpG repeats, but using scores rather than simple occurences.
//
// If one sets up the running sum, such that scores below some threshold count as negative
// then one should be able to distinguish from repeats of parts of the repeat sequence.
// That could be used by following up this with a peak detection, rather than by
// providing the full set of scores to R (4 * chr_length) bytes (for L. piscatorius
// around 100 megabytes per chromosome; not so bad.. )
//
// scores should be a pointer to 2 * (l1 + 1) integers
// seq_1 should be short, and seq_2 should be long
void smith_waterman_cm_mo(const char *seq_1, const char *seq_2, size_t l1, size_t l2,
			  int match, int mis_match, int gap,
			  int *scores, int *col_maxes)
{
  size_t nrow = l1 + 1;
  size_t ncol = l2 + 1;

  memset((void*)col_maxes, 0, sizeof(int) * ncol);
  memset((void*)scores, 0, sizeof(int) * 2 * nrow);
  
  for(size_t i=1; i < ncol; ++i){
    // current and previous columns
    int *ccol = scores + nrow * (i % 2);
    int *pcol = scores + nrow * ((i + 1) % 2);
    memset((void*)ccol, 0, sizeof(int) * nrow);
    col_maxes[i] = 0;
    ccol[0] = 0;
    for(size_t row=1; row < nrow; ++row){
      ccol[row] = 0;
      ccol[row] = ccol[row] < ccol[row-1] - gap ? ccol[row-1] - gap : ccol[row];
      ccol[row] = ccol[row] < pcol[row] - gap ? pcol[row] - gap : ccol[row];
      int match_p = seq_1[row-1] == seq_2[i-1] ? match : -mis_match;
      ccol[row] = pcol[row-1] + match_p > ccol[row] ? pcol[row-1] + match_p : ccol[row];
      col_maxes[i] = col_maxes[i] < ccol[row] ? ccol[row] : col_maxes[i];
    }
  }
}
