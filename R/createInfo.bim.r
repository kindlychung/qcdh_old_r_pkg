createInfo.bim = function(bimfile) {
    bimdat = read.csv(bimfile)
}
createInfo.file = function(machdat, infofile, legendfile, destfile) {
    infosmall = read.table(machdat, head=F)
    infobig = read.table(infofile, head=T)
    legendbig = read.table(legendfile, head=T)

    # eliminate multple headers caused by concatenating files
    idx = which(legendbig$rs == 'rs')
    if(length(idx) > 0) {
        infobig = infobig[-idx, ]
        legendbig = legendbig[-idx, ]
    }

    # first merge
    infomerge = merge_sameord(infosmall, infobig, by.x='V2', by.y='SNP')
    names(infomerge)[1] = 'snp'
    # the V1 column is all "M" from mach output, of no use to us.
    infomerge$V1 = NULL

    # second merge
    allmerge = merge_sameord(infomerge, legendbig, by.x='snp', by.y='rs')
    # check consistency between info and legend files
    if(!(all(allmerge$Al1 == allmerge$X0 & allmerge$Al2 == allmerge$X1)))
        stop('Inconsistency in alleles between the info file and the legend file')
    # remove redundant allele info
    allmerge$X0 = allmerge$X1 = NULL

    write.table(allmerge, file=destfile, quote=F, row.names=F, col.names=T)

    NULL
}
