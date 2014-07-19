hapinfo = function(machphasef) {
    hapinfo = read.table(machphasef, colClasses = c('NULL', 'NULL', 'character'), stringsAsFactors = FALSE)
    # convert 'a t c g' to 'A T C G'
    hapmat  = matrix(apply(hapinfo, 1, toupper), ncol=1)
    # split a genotype string into raw hex code
    hapmat  = apply(hapmat, 1, charToRaw)
    # convert them to 0 or 1
    hapmat  = apply(hapmat, 1, basecode)
    # calculate maf
    # maf = colMeans(hapmat)

    nr_hapmat = nrow(hapmat)
    # since there are 2 haplotypes for each person, the number of rows must be even
    stopifnot(nr_hapmat %% 2 == 0)
    # odd and even rows
    idx = as.logical((1:nr_hapmat) %% 2)
    idx1 = which(idx)
    idx2 = which(!idx)

    hapmat1 = t(hapmat[idx1, ])
    hapmat2 = t(hapmat[idx2, ])

    res = list(h1=hapmat1, h2=hapmat2)
    class(res) = 'hapinfo'

    res
}
