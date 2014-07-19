compP = function(...) UseMethod('compP')

compP.default = function(qcdhres) {
    require(foreach)
    splitres = split(qcdhres, qcdhres$firstsnp)
    reduced = foreach(i=1:length(splitres), .combine='rbind') %do% {
        res_i = splitres[[i]]
        Psingle = res_i[1, 4]
        Pbest = min(res_i[, 4])
        c(Psingle=Psingle, Pbest=Pbest, 
          logPsingle=-log10(Psingle), logPbest=-log10(Pbest))
    }
    class(reduced) = 'compP'
    reduced
}

plot.compP = function(Ps, pch=20, cex=.6, ...) {
    Ps = Ps[, 3:4]
    Ps = Ps[complete.cases(Ps), ]
    ymin = min(Ps)
    ymax = max(Ps)
    plot(Ps[, 1], ylim=c(ymin, ymax), ylab='-log P', xlab='Position', col='darkcyan', pch=pch, cex=cex, ...)
    points(Ps[, 2], col='brown4', pch=pch, cex=cex)
}


