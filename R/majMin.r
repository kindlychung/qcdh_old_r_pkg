# calculate major and minor alleles from info (data frame)
majMin = function(info) {
    alleles = info[, c('Al1', 'Al2')]

    majoridx = (info$Freq1 == info$MAF) + 1
    minoridx = 3 - majoridx
    majoridx = cbind(1:length(majoridx), majoridx)
    minoridx = cbind(1:length(minoridx), minoridx)
    majA = alleles[majoridx]
    minA = alleles[minoridx]

    reorder_alleles = cbind(majA, minA)
    res = cbind(subset(info, select=-c(Al1, Al2)), reorder_alleles)
    res
}
