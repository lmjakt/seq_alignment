#include <R.h>
#include <Rinternals.h>
#include "smith_waterman.h"

// Creates smith waterman score and pointer matrices
// for two sequences. Does not extract an alignment.

SEXP smith_water_matrices(SEXP seqs_r, SEXP penalties_r){
  if( TYPEOF(seqs_r) != STRSXP )
    error("first argument should be a character vector");
  if( TYPEOF(penalties_r) != INTSXP )
    error("second argument should be an integer vector");
  if(length(seqs_r) != 2 || length(penalties_r) != 3)
    error("First argument should contain two sequences, second argument three numbers");
  
  int *penalties = INTEGER(penalties_r);
  SEXP seq1_r = STRING_ELT(seqs_r, 0);
  SEXP seq2_r = STRING_ELT(seqs_r, 1);
  int l1 = length(seq1_r);
  int l2 = length(seq2_r);
  if(!l1 || !l2)
    error("Both sequences must have non-zero length");

  const char *seq1 = CHAR(seq1_r);
  const char *seq2 = CHAR(seq2_r);

  // Allocate the data for the return vectors. We will return a list containing
  // two matrices;
  SEXP ret_data = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ret_data, 0, allocMatrix(INTSXP, l1+1, l2+1));
  SET_VECTOR_ELT(ret_data, 1, allocMatrix(INTSXP, l1+1, l2+1));

  int *pointers = INTEGER(VECTOR_ELT(ret_data, 0));
  int *scores = INTEGER(VECTOR_ELT(ret_data, 1));

  smith_waterman_id(seq1, seq2, l1, l2,
		    penalties[0], penalties[1], penalties[2],
		    pointers, scores);
  UNPROTECT(1);
  return(ret_data);
}
