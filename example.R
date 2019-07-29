## try to use smith waterman to split sequences
## this done in R as we can use R to visualise
## statistics. These can then be passed down
## to other languages...

source("~/R/general_functions.R")

p.ff <- "TGGATTGATATGTAATACGACTCACTATAG"
p.fr <- "TCTCAGGCGTTTTTTTTTTTTTTTTTT"
p.bf <- "AAAAAAAAAAAAAAAAAACGCCTGAGA"
p.br <- "CTATAGTGAGTCGTATTACATATCAATCCA"

fastq.file <- "../star_mapping_R/CCSreads.fastq";

## read in a 250 sequences
fastq <- readLines( fastq.file, n=1000 )
fastq.id <- fastq[ seq(1, length(fastq), 4) ]
fastq.seq <- fastq[ seq(2, length(fastq), 4) ]

## try the smith-waterman alignment..
dl <- dyn.load("~/R/seq_alignment/seq_alignment.so")

p.seqs <- c(p.ff, p.fr, p.bf, p.br)


par(mfrow=c(5,1))
par(mar=c(1, 4.1, 1, 1))
system.time(
for(i in 1:length(fastq.seq)){
    tmp <- lapply( p.seqs, function(p){
        .Call("smith_water_matrices", c(p, fastq.seq[i]), as.integer(c(4, 4, 5)))
    })
    ## zlim <- range(unlist(sapply(tmp, function(x){ x[[2]] })))
    ## cmaxes <- sapply( tmp, function(x){
    ##     apply(x[[2]], 2, max) })
    ## ##
    ## invisible( lapply(tmp, function(x){
    ##     image( t(x[[2]][ nrow(x[[2]]):1, ]), col=hsvScale( 0:1000 ) )
    ## }))
    ## plot(1:nrow(cmaxes), ylim=range(cmaxes), type='n', xlab='position', ylab='max score', xaxs='i')
    ## for(j in 1:ncol(cmaxes))
    ##     lines(1:nrow(cmaxes), cmaxes[,j], col=j, lwd=2)
    ## abline(h=c(100,90,80))
    ## inpt <- readline("next:")
}
)


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

### It takes about 1.4 seconds for the above loop
### to go through 250 sequences (with four alignments
### per sequence. However, we don't really need the full matrices
### as we are more interested in the column max values
### and we should be able to vectorise the operations reasonably
### easily. Note that we still have to return a list

### Still, it is not too bad a situation... As it would take about
### 1.5 hours to do the full set of a million at that speed.

## I now have another function: smith_water_col_max
## which simply returns the column maximum values.
dl <- dyn.load("~/R/seq_alignment/seq_alignment.so")

system.time( tmp2 <- .Call("smith_water_col_max", p.seqs, fastq.seq, as.integer(c(4, 4, 5))) )
## compiled with -O3 this now gives 0.2 seconds or less.. Much better!
## 0.

## with default compilation the below is what you get
## .678 seconds
## 0.678 / 250
##   0.002712
## 2.7 milliseconds per set of four
## mean seq lengths: 28.5 * 1572.14
## 28.5 * 1572.14 * 4 * 250 / 1e6
## [1] 44.80599
## which means we have about
## 66 million score entries per second..
## we have four possible scores for each position so we can
## say we have to do 66 * 4 comparisons
## that is 264 million operations per second, or
## 0.26e9 ops / second. Well, I guess that is about 5 operations
## per 

## OK, so we now need to work out how to extract the peaks from this.
## I have some peak finding algorithms, but it might work out alright to do:

get.peaks <- function(v, min.v=80L, min.sep=27L){
    .Call("get_peaks", v, min.v, min.sep)
}

dl <- dyn.load("~/R/seq_alignment/seq_alignment.so")

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
## 0.029.. its ok.. 

## OK, it seems that we are getting OK performances from these. So we can go ahead
## and handle the full data set and then exporting some numbers for cutting off..

## read in the full data set in one go.. Terrible..
## read in a 250 sequences
fastq <- readLines( fastq.file )
fastq.id <- fastq[ seq(1, length(fastq), 4) ]
fastq.seq <- fastq[ seq(2, length(fastq), 4) ]

system.time( fastq.cm <- .Call("smith_water_col_max", p.seqs, fastq.seq, as.integer(c(4, 4, 5))) )
##    user  system elapsed 
## 675.996  14.500 690.548 
## 11 minutes. Not so bad..

system.time( fastq.cm.peaks <- lapply(fastq.cm, function(x){
    apply(x, 2, function(y){
        paste(get.peaks(y, min.v=75L), collapse=",")})
}))
##    user  system elapsed 
## 145.477   5.527 151.022 

## so that got the regions in less than 15 minutes, involving
## determining 4 million SW.. things.. Not so bad..

## lets output to text..
## identifier: (from fastq.id)
## peak positions:  p1,p2,p3,p4

fastq.cm.peak.table <- t(sapply( fastq.cm.peaks, eval))
## because terminal primer sequences are more likely, and because it
## is not a big issue to remove sequences from termini even if incorrect
## we want to also include the terminal scores for each sequence.

p.lengths <- nchar(p.seqs)

fastq.cm.term.scores <- t(sapply(fastq.cm, function(x){
    sapply(1:length(p.lengths), function(i){ ## get begin and end scores
        paste( x[p.lengths[i],i], x[nrow(x),i], sep="," )})
}))

## that is enough data to write out:
write.table(cbind(fastq.id, fastq.cm.peak.table, fastq.cm.term.scores),
            file='pacbio_primer_pos.txt', quote=FALSE, sep="\t")

