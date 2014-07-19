phewins = function(h1f, h2f, distf, infof, maflim, b=.5) {
    info = read.table(infof, head=TRUE)
    mafidx = which(maflim[1] < info$MAF & info$MAF < maflim[2])

    mydist = abigm(distf)
    h1 = abigm(h1f)
    h2 = abigm(h2f)

    if(length(mafidx) > 60) winspots = sample(mafidx, 60)
    else winspots = mafidx

    wins = list()
    for (i in mafidx) {
        # get the window
        winI = winIdx(mydist, i, 3e5)
        snpidx = c(winI$startsnp, winI$winidx)

        # sample causal_snps
        causal_snps = sample(snpidx, 2)
        xidx = causal_snps[1]
        yidx = causal_snps[2]
        x1 = h1[xidx, ]
        x2 = h2[xidx, ]
        y1 = h1[yidx, ]
        y2 = h2[yidx, ]

        # get the causal collapsed genotypes
        causal_cgI = collapseGen.haplo(x1, x2, y1, y2, collmethod=2)

        # simulate the phenotypes
        nindiv = length(causal_cgI)
        phe1 = b * causal_cgI + rnorm(nindiv)
        phe0 = rnorm(nindiv)

        # remove the causal_snps
        snpidx = snpidx[!(snpidx %in% causal_snps)]

        h1winI = h1[snpidx, , drop=FALSE]
        h2winI = h1[snpidx, , drop=FALSE]

        res = list(phe1=phe1, phe0=phe0, h1win=h1winI, h2win=h2winI)
        wins[[length(wins) + 1]] = res
    }

    wins
}

qcdhp = function(phe, cg) qcdhWin(phe, cg, verbose=F)$"Pr(>|t|)"

pwin = function(phe, collmethod, h1win, h2win) {
    nsnps = nrow(h1win)

    res = numeric()
    for(i in 1:nsnps) {
        cgmI = collapseGen(h1win, h2win, win=i, collmethod=collmethod)
        pminI = qcdhp(phe, cgmI)
        res[length(res) + 1] = pminI
    }

    min(res)
}

# todo: two for loops over phe0,1 and collmethod1,2,3
pswin = function(phe0, phe1, h1win, h2win) {
    dat = expand.grid(method=1:3, phe=c('phe0', 'phe1'))
    pval = numeric()
    for (phe in list(phe0, phe1)) {
        for (collmethod in 1:3) {
            pval[length(pval) + 1] = pwin(phe, collmethod, h1win, h2win)
        }
    }
    dat = cbind(dat, pval)
    dat
}

psfile = function(h1f, h2f, distf, infof, b=.5) {
    filewins = phewins(h1f, h2f, distf, infof, b)
    fileps = lapply(filewins, function(w) pswin(w$phe0, w$phe1, w$h1win, w$h2win))
    fileps = do.call('rbind', fileps)
    fileps
}

psfilelist = function(h1flist, h2flist, distflist, infoflist, b=.5) {
    stopifnot(length(h1flist) == length(h2flist) &&
              length(h2flist) == length(distflist) &&
              length(distflist) == length(infoflist))

    # listps = lapply(1:length(h1flist) function(i) {
    #                 h1f = h1flist[i]
    #                 h2f = h2flist[i]
    #                 distf = distflist[i]
    #                 infof = infoflist[i]
    #                 fileps = psfile(h1f, h2f, distf, infof, b)
    #                 fileps
    #           })

    registerDoMC()
    options(cores=6)

    listps = foreach(i=1:length(h1flist), .combine='rbind') %dopar% {
        h1f = h1flist[i]
        h2f = h2flist[i]
        distf = distflist[i]
        infof = infoflist[i]
        fileps = psfile(h1f, h2f, distf, infof, b)
        fileps
    }

    listps
}

