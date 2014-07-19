singleCollTest = function(b, SD, maf, nobser, nsnp) {
    g1 = sample(0:2, nobser*nsnp, replace=T, prob=c((1-maf)^2, 2*(1-maf)*maf, maf^2))
    g2 = sample(0:2, nobser*nsnp, replace=T, prob=c((1-maf)^2, 2*(1-maf)*maf, maf^2))
    cg = fg2(g1, g2)

    phe = cg * b + rnorm(length(cg))

    gp1 = g1 + rnorm(length(g1), 0, SD)
    gp2 = g2 + rnorm(length(g2), 0, SD)
    cgp = fg2(gp1, gp2)

    cgpm = matrix(cgp, nobser, nsnp)
    phem = matrix(phe, nobser, nsnp)

    getp = function(x, y) summary(lm(y~x))$coefficients[2, 4]

    ps = numeric()
    for(i in 1:nsnp) {
        g = cgpm[, i]
        phe = phem[, i]
        ps[length(ps) + 1] = getp(g, phe)
    }

    thresh = 5e-2 / nsnp
    mean(ps < thresh)
}


simuCollTest = function() {
    b = c(.1, .5, 2.5)
    SD = c(.01, .05, .15)
    maf = c(.01, .05, .1)
    conds = expand.grid(b=b, SD=SD, maf=maf)

    require(foreach)
    require(doMC)
    registerDoMC()
    options(cores=6)

    # powers = numeric()
    powers = foreach(i=1:nrow(conds), .combine='c') %dopar% {
        cond = conds[i, ]
        # powers[length(powers) + 1] = singleCollTest(cond$b, cond$SD, cond$maf, 5000, 200)
        singleCollTest(cond$b, cond$SD, cond$maf, 5000, 500)
    }

    res = cbind(conds, powers)
    res
}
