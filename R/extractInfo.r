# extract relevant info into result
extractInfo = function(info, snpidx, suffix) {
    info = info[snpidx,]
    info = majMin(info)
    info = subset(info, select=-c(Freq1, Quality, Rsq))
    names(info) = paste(names(info), suffix, sep='')
    info
}
