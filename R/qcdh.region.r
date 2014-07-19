qcdh.region = function(...) UseMethod('qcdh.region')

# This does not require a distance matrix
# it iterate from the first snp until the last
# thus it's very time consuming, and only suitable 
# for a small region, hence the name
qcdh.region.default = function(phe, family=gaussian(), h1win, h2win=NULL) {
    nsnps = nrow(h1win)

    res = lapply(1:nsnps, function(i) {
        cgmI = collapseGen(h1win, h2win, win=i, collmethod=2)
        qcdhWin(phe=phe, collgen=cgmI, family=family, verbose=F)
    })

    res = do.call('rbind', res)
    res
}

qcdh.region.formula = function(formula, covm, family=gaussian(), h1win, h2win=NULL) {
    yresid = resid(glm(formula=formula, family=family, data=covm))
    qcdh.region(yresid, family=family, h1win, h2win)
}
