collapseGen = function(...) UseMethod('collapseGen')

collapseGen.default = function(hapm1, hapm2=NULL, win, collmethod=2) {

    stopifnot(all(dim(hapm1) == dim(hapm2)))

    if(inherits(win, 'winIdx')) {
        startsnp = win$startsnp
        winidx = win$winidx
    } else {
        nsnps = nrow(hapm1)
        win = as.integer(win)
        stopifnot(length(win) == 1 & win > 0 & win <= nsnps)
        startsnp = win
        if(win < nsnps) winidx = (win+1):nsnps
        else winidx = integer(0)
    }
    winwidth = length(winidx)

    if(is.null(hapm2)) {
        cat('Only one matrix given, assuming unphased data.\n')
        snpfix = snp_no_collapse = hapm1[startsnp, , drop=FALSE]

        if(winwidth > 0) {
            snpvar = hapm1[winidx, , drop=FALSE]
            snpfix = snpfix[rep(1, winwidth), ]

            cgen = collapseGen.geno(snpfix, snpvar, collmethod, snp_no_collapse)

            re = list(snpidx=c(startsnp, winidx),
                      cgen=cgen)
        } else {
            re = list(snpidx = startsnp,
                      cgen = t(snp_no_collapse))
        }

    } else {
        cat('Two matrices given, assuming phased data.\n')

        snpfix1 = hapm1[startsnp, , drop=FALSE]
        snpfix2 = hapm2[startsnp, , drop=FALSE]
        snp_no_collapse = snpfix1 + snpfix2

        if(winwidth > 0) {
            snpvar1 = hapm1[winidx, , drop=FALSE]
            snpvar2 = hapm2[winidx, , drop=FALSE]
            snpfix1 = snpfix1[rep(1, winwidth), ]
            snpfix2 = snpfix2[rep(1, winwidth), ]

            cgen = collapseGen.haplo(snpfix1, snpfix2, snpvar1, snpvar2, collmethod, snp_no_collapse)

            re = list(snpidx=c(startsnp, winidx),
                      cgen=cgen)
        } else {
            re = list(snpidx = startsnp,
                      cgen = t(snp_no_collapse))
        }
    }

    class(re) = c('collapsegen', class(re))
    re
}


collapseGen.geno = function(snpfix, snpvar, collmethod, snp_no_collapse) {
    if(is.null(dim(snpfix))) snpfix = matrix(snpfix, nrow=1)
    if(is.null(dim(snpvar))) snpvar = matrix(snpvar, nrow=1)

    gsum = round(snpfix + snpvar)

    cgen = (gsum > 1) + 0
    idx = which(gsum > 2)

    if(collmethod == 2) {

    } else if(collmethod == 1) {
        cgen[idx] = 0
    } else if(collmethod == 3) {
        cgen[idx] = NA
    } else {
        stop('Unknown collmethod option!')
    }

    if(!missing(snp_no_collapse)) 
        cgen = rbind(snp_no_collapse, cgen)
    dimnames(cgen) = NULL
    cgen = t(cgen)
    cgen
}


collapseGen.haplo = function(snpfix1, snpfix2, snpvar1, snpvar2, collmethod, snp_no_collapse) {
    if(is.null(dim(snpfix1))) snpfix1 = matrix(snpfix1, nrow=1)
    if(is.null(dim(snpfix2))) snpfix2 = matrix(snpfix2, nrow=1)
    if(is.null(dim(snpvar1))) snpvar1 = matrix(snpvar1, nrow=1)
    if(is.null(dim(snpvar2))) snpvar2 = matrix(snpvar2, nrow=1)

    gsum = round((snpfix1 + snpvar1) * (snpfix2 + snpvar2))
    cgen = (gsum > 0) + 0
    idx = which(gsum > 1)

    if(collmethod == 2) {
    } else if(collmethod == 1) {
        cgen[idx] = 0
    } else if(collmethod == 3) {
        cgen[idx] = NA
    } else {
        stop('Unknown collmethod option!')
    }

    if(!missing(snp_no_collapse)) 
        cgen = rbind(snp_no_collapse, cgen)
    dimnames(cgen) = NULL

    cgen = t(cgen)
    cgen
}
