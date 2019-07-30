## try to use smith waterman to split sequences
## this done in R as we can use R to visualise
## statistics. These can then be passed down
## to other languages...

source("~/R/experimental_R/general_functions.R")

p.ff <- "TGGATTGATATGTAATACGACTCACTATAG"
p.fr <- "TCTCAGGCGTTTTTTTTTTTTTTTTTT"
p.bf <- "AAAAAAAAAAAAAAAAAACGCCTGAGA"
p.br <- "CTATAGTGAGTCGTATTACATATCAATCCA"

fastq.file <- "CCSreads.fastq";

## read in sequences from CCSreads.fastq
fastq <- readLines( fastq.file )
fastq.id <- fastq[ seq(1, length(fastq), 4) ]
fastq.seq <- fastq[ seq(2, length(fastq), 4) ]

## try the smith-waterman alignment..
dl <- dyn.load("seq_alignment.so")

p.seqs <- c(p.ff, p.fr, p.bf, p.br)

## check the time required... 
system.time(
for(i in 1:length(fastq.seq)){
    tmp <- lapply( p.seqs, function(p){
        .Call("smith_water_matrices", c(p, fastq.seq[i]), as.integer(c(4, 4, 5)))
    })
})
## 0.67, 0.728, 0.723, 0.723

par(mfrow=c(5,1))
par(mar=c(1, 4.1, 1, 1))
for(i in 1:length(fastq.seq)){
    tmp <- lapply( p.seqs, function(p){
        .Call("smith_water_matrices", c(p, fastq.seq[i]), as.integer(c(4, 4, 5)))
    })
    cmaxes <- sapply( tmp, function(x){
        apply(x[[2]], 2, max) })
    ##
    invisible( lapply(tmp, function(x){
        image( t(x[[2]][ nrow(x[[2]]):1, ]), col=hsvScale( 0:1000 ) )
    }))
    plot(1:nrow(cmaxes), ylim=range(cmaxes), type='n', xlab='position', ylab='max score', xaxs='i')
    for(j in 1:ncol(cmaxes))
        lines(1:nrow(cmaxes), cmaxes[,j], col=j, lwd=2)
    abline(h=c(100,90,80))
    inpt <- readline(paste("next:", i))
    if(inpt == 'q')
        break
}

## a nice example of a chimeric sequence at 16
i <- 16
png("chimeric_sequence.png", width=1200, height=1000)
par(mfrow=c(5,1))
par(mar=c(3.1, 4.1, 1.1, 1.1))
tmp <- lapply( p.seqs, function(p){
    .Call("smith_water_matrices", c(p, fastq.seq[i]), as.integer(c(4, 4, 5)))
})

cmaxes <- sapply( tmp, function(x){
    apply(x[[2]], 2, max) })
##
invisible( lapply(tmp, function(x){
    image( t(x[[2]][ nrow(x[[2]]):1, ]), col=hsvScale( 0:1000 ) )
}))
plot(1:nrow(cmaxes), ylim=range(cmaxes), type='n', xlab='position', ylab='max score', xaxs='i')
for(j in 1:ncol(cmaxes))
    lines(1:nrow(cmaxes), cmaxes[,j], col=j, lwd=2)
abline(h=c(100,90,80))
dev.off()

### It takes about 0.7 seconds for the above loop
### to go through 250 sequences (with four alignments
### per sequence. However, we don't really need the full matrices
### as we are more interested in the column max values
### and we should be able to vectorise the operations reasonably
### easily. Note that we still have to return a list

## smith_water_col_max
## simply returns the column maximum values.

system.time( tmp2 <- .Call("smith_water_col_max", p.seqs, fastq.seq, as.integer(c(4, 4, 5))) )
## compiled with -O3 this now gives 0.2 seconds or less.. Much better!
## 0.212, 0.19, 0.217, 0.179, 0.214
## note that the speed is about three times faster with -O3, then with
## a default non-optimised compile.
##
## on my laptop I am getting 0.101 seconds... 

get.peaks <- function(v, min.v=80L, min.sep=27L){
    .Call("get_peaks", v, min.v, min.sep)
}

## Make sure that the vectorised column max function gives the same result
## as the one based on the full tables.
par(mfrow=c(6,1))
par(mar=c(3.1, 4.1, 1.1, 1.1))
for(i in 1:length(fastq.seq)){
    tmp <- lapply( p.seqs, function(p){
        .Call("smith_water_matrices", c(p, fastq.seq[i]), as.integer(c(4, 4, 5)))
    })
    cmaxes <- sapply( tmp, function(x){
        apply(x[[2]], 2, max) })
    ##
    invisible( lapply(tmp, function(x){
        image( t(x[[2]][ nrow(x[[2]]):1, ]), col=hsvScale( 0:1000 ) )
    }))
    plot(1:nrow(cmaxes), ylim=range(cmaxes), type='n', xlab='position', ylab='max score', xaxs='i')
    for(j in 1:ncol(cmaxes))
        lines(1:nrow(cmaxes), cmaxes[,j], col=j, lwd=2)
    abline(h=c(100,90,80))
    plot(1:nrow(tmp2[[i]]), tmp2[[i]][,1], type='l', xaxs='i', lwd=2, ylim=range(cmaxes) )
    lines(1:nrow(tmp2[[i]]), tmp2[[i]][,2], type='l', col=2, lwd=2)
    lines(1:nrow(tmp2[[i]]), tmp2[[i]][,3], type='l', col=3, lwd=2)
    lines(1:nrow(tmp2[[i]]), tmp2[[i]][,4], type='l', col=4, lwd=2)
    peaks <- apply( tmp2[[i]], 2, get.peaks, min.v=75L )
    abline(h=75)
    for(j in 1:length(peaks))
        abline(v=peaks[[j]], col=j, lwd=4)
    inpt <- readline(paste("next:", i))
    if(inpt == 'q')
        break
}


system.time(
tmp3 <- lapply(tmp2, function(x){
    apply(x, 2, get.peaks, min.v=75L) })
)
## 0.029.. it is sufficiently fast even using lapply to

## I have also a multithreaded function.

## try with 2 threads first
#system.time(

tmp4 <- .Call("smith_water_col_max_mt", p.seqs, fastq.seq, as.integer(c(4,4,5)), 1L)
## segmentation error!

