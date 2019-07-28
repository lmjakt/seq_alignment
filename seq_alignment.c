#include <R.h>
#include <Rinternals.h>
#include <pthread.h>
#include "smith_waterman.h"

static const int sw_matrix_size = 30 * 10000;
static const int MAX_THREADS = 16;

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

// A function that returns the column maxima for two sets of sequences
// Returns a matrix of column maximums for each sequence...

SEXP smith_water_col_max(SEXP seq_short_r, SEXP seq_long_r, SEXP penalties_r){
  if( TYPEOF(seq_short_r) != STRSXP || TYPEOF(seq_long_r) != STRSXP )
    error("The first two arguments should both be character vectors");
  if( TYPEOF(penalties_r) != INTSXP )
    error("second argument should be an integer vector");
  if( length(seq_short_r) < 1 || length(seq_long_r) < 1 || length(penalties_r) != 3)
    error("Arguments should have lengths: l1 > 0, l2 > 0, l3 == 3");
  
  int *penalties = INTEGER(penalties_r);
  int n_short = length(seq_short_r);
  int n_long = length(seq_long_r);

  // for each long sequence we wish to determine column maxes. That means
  // we need to return a list of n_long items, containing matrices with
  // n_short columns and nrow equal to the long sequence... 
  SEXP ret_data = PROTECT(allocVector(VECSXP, n_long));
  
  // we define a single matrix that we will try to reuse
  // This will take up about 2.4MB of memory and should be sufficient
  int *score_buffer = malloc(sizeof(int) * 30 * 10000);

  // Then iterate through the long sequences:
  for(int i=0; i < n_long; ++i){
    SEXP seq_l_r = STRING_ELT(seq_long_r, i);
    int seq_length = length(seq_l_r);
    const char *seq_l = CHAR(seq_l_r);
    if(seq_length < 1)
      continue;
    SET_VECTOR_ELT(ret_data, i, allocMatrix(INTSXP, seq_length + 1, n_short));
    int *c_max = INTEGER(VECTOR_ELT(ret_data, i));
    for(int j=0; j < n_short; ++j){
      SEXP seq_s_r = STRING_ELT(seq_short_r, j);
      int primer_length = length(seq_s_r);
      const char *seq_s = CHAR(seq_s_r);
      // call realloc to make sure that we have sufficient memory for the
      // buffer
      score_buffer = realloc((void*)score_buffer, sizeof(int) * (1 + seq_length) * (1 + primer_length));
      smith_waterman_cm(seq_s, seq_l, primer_length, seq_length,
			penalties[0], penalties[1], penalties[2],
			score_buffer, c_max + (seq_length + 1) * j);
    }
  }
  free( score_buffer );
  UNPROTECT(1);
  return(ret_data);
}


// long_seqs and short_seqs are the sets of sequences to be aligned
// ret_dat is a vector of matrices which has been allocated with
// nrow=length(long_seqs[i]), ncol=length(short_seqs[i])
// start is the first long sequence which should be compared.
// end is the index of the last sequence + 1
struct col_max_thread_arg {
  SEXP long_seqs;
  SEXP short_seqs;
  int *score_buffer;
  int *penalties;
  SEXP ret_data;
  int start;
  int end;
};

void* smith_water_col_max_thread(void *args_ptr){
  struct col_max_thread_arg args = *(struct col_max_thread_arg*)args_ptr;
  int n_short = length( args.short_seqs );
  for(int i=args.start; i < args.end; ++i){
    SEXP seq_l_r = STRING_ELT(args.long_seqs, i);
    int seq_l_l = length(seq_l_r);
    const char *seq_l = CHAR(seq_l_r);
    if(seq_l_l < 1)
      continue;
    int *c_max = INTEGER(VECTOR_ELT(args.ret_data, i));
    for(int j=0; j < n_short; ++j){
      SEXP seq_s_r = STRING_ELT(args.short_seqs, j);
      int seq_s_l = length(seq_s_r);
      const char *seq_s = CHAR(seq_s_r);
      args.score_buffer = realloc((void*)args.score_buffer,
				  sizeof(int) * (1 + seq_l_l) * (1 + seq_s_l));
      smith_waterman_cm(seq_s, seq_l, seq_s_l, seq_l_l,
			args.penalties[0], args.penalties[1], args.penalties[2],
			args.score_buffer, c_max + (seq_l_l + 1) * j);
    }
  }
  pthread_exit(NULL);
}
  
// A function that returns the column maxima for two sets of sequences
// Returns a matrix of column maximums for each sequence...
// multithreaded for speed.. 
SEXP smith_water_col_max_mt(SEXP seq_short_r, SEXP seq_long_r, SEXP penalties_r, SEXP thread_n_r){
  if( TYPEOF(seq_short_r) != STRSXP || TYPEOF(seq_long_r) != STRSXP )
    error("The first two arguments should both be character vectors");
  if( TYPEOF(penalties_r) != INTSXP || TYPEOF(thread_n_r) != INTSXP )
    error("third and fourth arguments should be integer vectors");
  if( length(seq_short_r) < 1 || length(seq_long_r) < 1 || length(penalties_r) != 3 || length(thread_n_r) != 1)
    error("Arguments should have lengths: l1 > 0, l2 > 0, l3 == 3, l4 == 1");
  
  int *penalties = INTEGER(penalties_r);
  int n_short = length(seq_short_r);
  int n_long = length(seq_long_r);
  int thread_n = INTEGER(thread_n_r)[0];
  if(thread_n < 1 || thread_n >= MAX_THREADS)
    error("bad thread number %d", thread_n);

  // Since we do not wish to allocate mamory within the worker threads we will
  // first go through all of the sequences and allocate memory for te
  // return data structure. Then we will create the threads and run them
  // in parallel.
  
  // for each long sequence we wish to determine column maxes. That means
  // we need to return a list of n_long items, containing matrices with
  // n_short columns and nrow equal to the long sequence... 
  SEXP ret_data = PROTECT(allocVector(VECSXP, n_long));
  

  // Then iterate through the long sequences:
  for(int i=0; i < n_long; ++i){
    SEXP seq_l_r = STRING_ELT(seq_long_r, i);
    int seq_length = length(seq_l_r);
    if(seq_length < 1)
      continue;
    SET_VECTOR_ELT(ret_data, i, allocMatrix(INTSXP, seq_length + 1, n_short));
  }

  // we will define a separate score buffer for each worker thread
  // This will take up about 2.4MB of memory and should be sufficient
  int **score_buffers = malloc( sizeof(int*) * thread_n);
  pthread_t *threads = malloc(sizeof(pthread_t) * thread_n);
  struct col_max_thread_arg *thread_arg = malloc(sizeof(struct col_max_thread_arg)
						 * thread_n);
  for(int i=0; i < thread_n; ++i){
    score_buffers[i] = malloc(sizeof(int) * sw_matrix_size);
    thread_arg[i].long_seqs = seq_long_r;
    thread_arg[i].short_seqs = seq_short_r;
    thread_arg[i].score_buffer = score_buffers[i];
    thread_arg[i].penalties = penalties;
    thread_arg[i].ret_data = ret_data;
    thread_arg[i].start = i * n_long / thread_n;
    thread_arg[i].end = (i == thread_n - 1) ? n_long : (i + 1) * n_long / thread_n;
    // and then we create the thread and start it..
    threads[i] = pthread_create(&threads[i], NULL,
				smith_water_col_max_thread,
				(void *)&(thread_arg[i]) );
  }
  // wait for the threads to finish using join
  // we do not consider any errors here. Probably we cannot use
  // error() in the worker thread, so we have to be careful..
  void *status;
  for(int i=0; i < thread_n; ++i){
    pthread_join(threads[i], &status );
    free(score_buffers[i]);
  }
  free(score_buffers);
  free(thread_arg);
  free(threads);

  UNPROTECT(1);
  return( ret_data );
}


// A function to return some peaks. This is simple enough.
// only handles integers
// do two passes to avoid having to do funky memory arrangements
// This would be better with std::vector<> and push_back()
SEXP get_peaks(SEXP v_r, SEXP min_v_r, SEXP min_dist_r){
  if(TYPEOF(v_r) != INTSXP || TYPEOF(min_v_r) != INTSXP || TYPEOF(min_dist_r) != INTSXP)
    error("Bad argument types");
  if(length(v_r) < 2 || length(min_v_r) != 1 || length(min_dist_r) != 1)
    error("Bad argument lengths");
  int n = length(v_r);
  int *v = INTEGER(v_r);
  int min_v = *INTEGER(min_v_r);
  int min_dist = *INTEGER(min_dist_r);
  
  unsigned char *peaks = malloc(sizeof(char) * n);
  size_t peak_n = 0;
  memset((void*)peaks, 0, sizeof(char) * n);
  int last_peak = -1;
  for(int i=1; i < n; ++i){
    if(v[i] > min_v && v[i] > v[i-1] && (i == n-1 || v[i] > v[i+1])){
      if(last_peak < 0 || i - last_peak > min_dist){
	peaks[i] = 1;
	peak_n++;
	last_peak = i;
      }else{
	if(v[i] > v[last_peak]){
	  peaks[last_peak] = 0;
	  peaks[i] = 1;
	  last_peak = i;
	}
      }
    }
  }
  SEXP peak_pos_r = PROTECT(allocVector(INTSXP, peak_n));
  int *peak_pos = INTEGER(peak_pos_r);
  int peak_i = 0;
  for(int i=0; i < n && peak_i < peak_n; ++i){
    if(peaks[i]){
      peak_pos[peak_i] = i+1;
      ++peak_i;
    }
  }
  UNPROTECT(1);
  return( peak_pos_r );
}
