# for simulation purposes
collgen1 = function(snp11, snp21, snp12, snp22) {
    cgen = (xor(snp11, snp21) & xor(snp12, snp22)) + 0
    cgen
}


# for simulation purposes
collgen2 = function(snp11, snp21, snp12, snp22) {
    cgen = ((snp11 | snp21) & (snp12 | snp22)) + 0
    cgen
}

rsnp12 = function(snp11, snp21, snp12, snp22) {
    snp1 = snp11 + snp12
    snp2 = snp21 + snp22
    cor(snp1, snp2)
}

simSnps = function(infofile, maflim=.05, winsizelim=5e5) {
    info = read.table(infofile, head=TRUE)
    s = which(info$MAF < maflim)
    s1 = s[-length(s)]
    s2 = s[-1]

    maf1 = info$MAF[s1]
    maf2 = info$MAF[s2]
    pos1 = info$position[s1]
    pos2 = info$position[s2]
    winsize = pos2 - pos1

    simsnps = data.frame(ss1=info$snp[s1], 
               ss2=info$snp[s2], 
               s1, s2, 
               winsize)

    simsnps = na.omit(simsnps)
    simsnps[which(simsnps$winsize < winsizelim), ]
}

psimPmin = function(s1, s2, h1, h2,                     
                    h1files, h2files, distfiles, infofiles, 
                    collmethod=1, 
                    maflim=.05, 
                    b=c(.1, .5, 2.5), 
                    winsizelim=5e5) {
    snp11 = h1[s1, ]
    snp12 = h2[s1, ]
    snp21 = h1[s2, ]
    snp22 = h2[s2, ]

    collfunctions = list(collgen1, collgen2)
    collfunc = collfunctions[[collmethod]]
    cg = collfunc(snp11, snp21, snp12, snp22)

    phe1 = cg * b[1] + rnorm(length(cg))
    phe2 = cg * b[2] + rnorm(length(cg))
    phe3 = cg * b[3] + rnorm(length(cg))
    phelist = list(phe1, phe2, phe3)

    psim = sapply(phelist, function(phe) {
                        tryCatch({
                            summary(glm(phe~cg))$coefficients[2, 4]
                        }, error = function(e) NA
                        )
                    })

    qcdhlist = list()
    for(phe in list(phe1, phe2, phe3)) {
        qcdhres = qcdh.batch(phe, 
                             gaussian(), 
                             h1files, 
                             h2files,
                             distfiles,
                             infofiles,
                             maflim=c(0, maflim), 
                             collmethod=collmethod, 
                             winsize=winsizelim)    

        qcdhres = simpleQcdhres(qcdhres)
        qcdhlist[[length(qcdhlist) + 1]] = qcdhres
    }             

    pickups = do.call('rbind', qcdhlist)
    cbind(pickups, psim, b)
}


simpleQcdhres = function(qcdhres) {
    qcdhres = qcdhres[, c('snp1', 'snp2', 'p', 'p1', 'Ntests')]
    minidx = which.min(qcdhres$p)
    minidx1 = which.min(qcdhres$p1)
    pickup = qcdhres[minidx, ]
    pickup1 = qcdhres[minidx1, ]
    names(pickup1)[1:2] = c('snp3', 'snp4')

    pickup = cbind(pickup, pickup1)
    return(pickup)
}
