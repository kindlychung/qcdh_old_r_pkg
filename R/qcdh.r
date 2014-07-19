qcdh = function(...) UseMethod('qcdh')

# qcdh for single files
qcdh.default = function(
                        phe, family, 
                        hapm1f, hapm2f=NULL, 
                        distmf, infof, 
                        maflim=NULL, 
                        winsize=1e5, 
                        collmethod=2,
                        inputtype
                        ) {
    cat('Haplotype 1 file: ', hapm1f, '\n')
    cat('Haplotype 2 file: ', hapm2f, '\n')
    cat('Distance file: ', distmf, '\n')
    cat('Info file: ', infof, '\n')

    if (inputtype == "bigmatrixfile") {
        hapm1 = abigm(hapm1f)

        if(is.null(hapm2f)) hapm2=NULL
        else hapm2 = abigm(hapm2f)

        distm = abigm(distmf)
        info = read.table(infof, head=TRUE)
    } else if (inputtype == "matrix") {
        hapm1 = hapm1f
        hapm2 = hapm2f
        distm = distmf
        info = infof
    }

    mafv = info$MAF

    if (is.null(maflim)) mafidx = 1:nrow(hapm1)
    else mafidx = which(mafv > maflim[1] & mafv <= maflim[2])

    if(length(mafidx) == 0) {
        cat('No SNP with MAF between ', maflim[1], ' and ', maflim[2],  
            'see the info file: ', infof, '\n')
        return(NULL)
    }

    # get chr number
    chr = chrnumber(infof)

    wins_stats = list()
    for(i in mafidx) {
        win = winIdx(distm, i, winsize, mafidx)
        collgen = collapseGen(hapm1, hapm2, win, collmethod=collmethod)
        winstats = qcdhWin(phe, collgen, family=family)
        wins_stats[[length(wins_stats) + 1]] = winstats
    }
    wins_stats = do.call('rbind', wins_stats)

    snp1idx =  wins_stats$firstsnp
    snp2idx =  wins_stats$secondsnp
    info1 = extractInfo(info, snp1idx, 1)
    info2 = extractInfo(info, snp2idx, 2)

    wins_stats = wins_stats[, c('Estimate', 'Pr(>|t|)', 'b1', 'p1', 'Ntests')]
    wins_stats = cbind(info1, info2,  wins_stats, chr=chr)

    # fix the awkward "Pr(>|t|)" colname 
    pidx = grep('Pr.>.t.', colnames(wins_stats))
    colnames(wins_stats)[pidx] = 'p'
    pidx = grep('Estimate', colnames(wins_stats))
    colnames(wins_stats)[pidx] = 'b'

    rownames(wins_stats) = NULL

    wins_stats
}

qcdh.batch = function(
                      phe, family, 
                      hapm1fs, hapm2fs=NULL, 
                      distmfs, infofs, 
                      ...
                      ) {
    h1list = Sys.glob(hapm1fs)

    if(is.null(hapm2fs)) h2list = rep(NULL, length(h1list))
    else h2list = Sys.glob(hapm2fs)

    distlist = Sys.glob(distmfs)
    infolist = Sys.glob(infofs)

    n_files = sapply(list(h1list, h2list, distlist, infolist), length)
    stopifnot(all(n_files[1] == n_files))

    # file_stats = list()
    # for(i in 1:(n_files[1])) {
    #     filestat = qcdh.default(phe, family, h1list[i], h2list[i], distlist[i], infolist[i], ...)
    #     file_stats[[length(file_stats) + 1]] = filestat
    # }
    # file_stats = do.call('rbind', file_stats)

    registerDoMC()
    options(cores=6)
    file_stats = foreach(i=1:(n_files[1]), .combine='rbind') %dopar% {
        qcdh.default(phe, family, h1list[i], h2list[i], distlist[i], infolist[i], ...)
    }

    rownames(file_stats) = NULL
    file_stats
}

qcdh.formula = function(
                        formula, covm, family=gaussian(), 
                        batch=FALSE, 
                        ...
                        ) {
    yresid = resid(glm(formula=formula, family=family, data=covm))
    if(batch == TRUE) {
        return(qcdh.batch(yresid, family, ...))
    } else {
        return(qcdh.default(yresid, family, ...))
    }
}

