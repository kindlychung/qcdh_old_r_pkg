qcdhWin = function(...) UseMethod('qcdhWin')

# qcdhWin.formula = function(formula, covm, collgen, ...) {
#     yresid = resid(glm(formula=formula, data=covm, ...))
#     qcdhWin.default(yresid, collgen, ...)
# }

qcdhWin.default = function(phe, collgen, family=gaussian(), verbose=F) {
    stopifnot(inherits(collgen, 'collapsegen'))

    winstats = list()
    for(i in 1:ncol(collgen$cgen)) {
        g = collgen$cgen[, i]
        dat = data.frame(g=g, phe=phe)
        dat = dat[complete.cases(dat), ]
        coefs = tryCatch({
            res = summary(glm(dat$phe~dat$g, family=family))$coefficients[2, c(1, 2, 4)]
            # if(res[2] < 1e-10) save(phe, g, file=paste('/tmp/log', i, '.rda', sep=''))
            res
        }, error = function(e) {
            if(sd(g) != 0) {
                warning('Error in glm model,
                        but SD of collapse genotype not equal to zero,
                        you need to look into this.
                        ')
            }
            # c(Estimate=NA, 'Pr(>|t|)'=NA)
            c(rep(NA, 3))
        })
        winstats[[length(winstats) + 1]] = coefs
    }
    winstats = do.call('rbind', winstats)

    secondsnp = collgen$snpidx
    firstsnp = rep(secondsnp[1], length(secondsnp))
    stopifnot(length(secondsnp) == nrow(winstats))
    winstats = cbind(firstsnp, secondsnp, winstats)

    if(verbose) {
        winstats
    } else {
        winstats1 = na.omit(as.data.frame(winstats))
        Ntests = nrow(winstats1)

        if(nrow(winstats1) > 0) {
            # index of minimum p
            pminidx = which.min(winstats1[, 'Pr(>|t|)'])[1]
            p1idx = which(winstats1$firstsnp == winstats1$secondsnp)
            if(length(p1idx) > 0) {
                stopifnot(length(p1idx) == 1 & p1idx == 1)
                p1dat = winstats1[p1idx, c("Std. Error", "Pr(>|t|)")]
                names(p1dat) = c('b1', 'p1')
            } else {
                p1dat = data.frame(b1=NA, p1=NA)
            }
            winstats1 = cbind(winstats1[pminidx, ], p1dat)
            rownames(winstats1) = NULL
        } else {
            p1dat = data.frame(b1=NA, p1=NA)
            winstats1 = cbind(winstats[1, ], p1dat)
        }

        winstats1$Ntests = Ntests
        winstats1
    }
}
