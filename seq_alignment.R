## Wrapper and other functions for seq_alignment

dyn.load("seq_alignment.so")


## return matrices for Smith-Waterman matrices
## optionally include the actual alignement
SW.matrices  <- function(seqs, penalties, do.align=FALSE){
    if(length(seqs) < 2)
        stop("Supply at least two sequences")
    if(!is.character(seqs))
        stop("sequences must be character vectors")
    if(!is.numeric(penalties))
        stop("penalties must be numeric")
    if(length(penalties) != 3)
        stop("Supply three scores / penalties: match, mismatch, gap")
    swm  <- .Call("smith_water_matrices", seqs[1:2], as.integer(penalties))
    if(length(swm) == 2)
        names(swm)  <- c('scores', 'pointers')
    if(length(swm) == 3)
        names(swm)  <- c('scores', 'pointers', 'seqs')
    swm
}

## Returns a list of matrices of dimensions:
## nchar(seqs.long[i]) x length(seqs.short)
## containing the maximum for each column in a SW score
## matrix where the seqs.short make up the rows.
SW.colMax  <- function(seqs.short, seqs.long, penalties){
    if(!is.integer(penalties))
        warning("Penalties will be converted to integers")
    cmax  <- .Call("smith_water_col_max", seqs.short, seqs.long,
                   as.integer(penalties))
    cmax
}

## As SW.colMax, but splits the long sequences into groups
## which are processed in parallel
SW.colMax  <- function(seqs.short, seqs.long, penalties, nthreads){
    if(!is.integer(penalties))
        warning("Penalties will be converted to integers")
    cmax  <- .Call("smith_water_col_max", seqs.short, seqs.long,
                   as.integer(penalties), as.integer(nthreads))
    cmax
}


## a very simple peak identification algorithm
## no smoothing, but does use a minimum peak
## distance parameter
## values: the values
## min.value: the minimum peak value
## min.distance: the minimum distance between peaks
int.peaks  <- function(values, min.value, min.distance){
    if(!is.integer(values) || !is.integer(min.value) || !is.integer(min.distance))
        warning("All arguments will be converted to integer values")
    .Call("get_peaks", as.integer(values), as.integer(min.value),
          as.integer(min.distance) )
}
