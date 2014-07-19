winIdx = function(distm, i, winsize, mafidx=NULL) {
    disti = distm[,i]

    # index of SNPs within the winsize range of the ith SNP (downstream)
    winidx = i+which(disti > 0 & disti < winsize)
    # if a mafidx provided, then filter by that also
    if(!is.null(mafidx)) {
        winidx = winidx[which(winidx %in% mafidx)]
    } 

    res = list(winidx = winidx,
               startsnp = i)
    class(res) = 'winIdx'
    return(res)
}
