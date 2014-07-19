infolist = Sys.glob('~/Phased/info/chr22_0?.info')
h1list = Sys.glob('~/Phased/bigmat/chr22_0?_mach1.hap1.desc')
h2list = Sys.glob('~/Phased/bigmat/chr22_0?_mach1.hap2.desc')
distlist = Sys.glob('~/Phased/distm/chr22_0?.dist.desc')
stopifnot(length(infolist) == length(h1list) & 
          length(h1list) == length(h2list))

filedat = data.frame(info=infolist, h1=h1list, h2=h2list)

scanOneInfo = function(infofile, h1file, h2file,  
                       maflim=.05, winsizelim=3e5,
                       collmethod=1, 
                       b=c(.1, .5, 2.5)
                       # b=c(4, 5, 6)
                       ) {
    simsnps = simSnps(infofile, maflim=maflim, winsizelim=winsizelim)
    
    h1 = abigm(h1file)
    h2 = abigm(h2file)

    pminlist = list()
    for(i in 1:nrow(simsnps)) {
        snpdat = simsnps[i, ]
        s1 = snpdat$s1
        s2 = snpdat$s2

        pminlist[[length(pminlist) + 1]] = psimPmin(s1, s2, h1, h2, 
                                                    h1list, h2list, distlist, infolist,
                                                    collmethod=collmethod,
                                                    maflim=maflim,
                                                    winsizelim=winsizelim, 
                                                    b=b
                                                    )
    }
    pmins = do.call(rbind, pminlist)

    simsnps = simsnps[rep(1:nrow(simsnps), each=3), ]
    stopifnot(nrow(simsnps) == nrow(pmins))
    cbind(simsnps, pmins)
}

dosimu = function() {
    allpmins = list()
    for(i in 1:nrow(filedat)) {
        infofile = filedat[i, 'info']
        h1file = filedat[i, 'h1']
        h2file = filedat[i, 'h2']
        cat('======================\n')
        print(infofile)
        print(h1file)
        print(h2file)
        cat('======================\n')
        pmins = scanOneInfo(infofile, h1file, h2file)
        allpmins[[length(allpmins) + 1]] = pmins
    }
    allpmins = do.call(rbind, allpmins)
    row.names(allpmins) = NULL
    names(allpmins)[8:9] = c('qmp', 'qmp1')
    save(allpmins, file='allpmins.rda')
    allpmins
}

anaSimu = function(allpmins) {
    plist = split(allpmins, as.factor(allpmins$b))

    pickratelist = list()
    for(i in 1:length(plist)) {
        dat = plist[[i]]
        b = dat$b[1]
        ss1 = pmin(dat$ss1, dat$ss2)
        ss2 = pmax(dat$ss1, dat$ss2)
        snp1 = pmin(dat$snp1, dat$snp2)
        snp2 = pmax(dat$snp1, dat$snp2)
        snp3 = dat$snp3

        qp1 = mean(ss1 == snp1 | ss2 == snp2)
        qp2 = mean(ss1 == snp1 & ss2 == snp2)
        sp1 = mean(ss1 == snp3 | ss2 == snp3)

        pickratelist[[length(pickratelist) + 1]] = data.frame(b, qp1, qp2, sp1, prr=qp1/sp1)
    }
    pickrates = do.call(rbind, pickratelist)
    pickrates
}
