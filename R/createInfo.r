# Given the data.dat file (provided by mach prephasing),
# the .info and .legend file (result of imputation),
# create a single info file, containing, among other things,
# MAF and bp position
createInfo = function(arg1, ...) {
    if(file.info(arg1)$isdir) {
        createInfo.dir(arg1, ...)
    } else if(file.exists(arg1)) {
        createInfo.file(arg1, ...)
    } else {
        stop('File does not exist!')
    }

    NULL
}

createInfo.dir = function(datdir, infofiledir, legendfiledir, dest) {
    # get lists of files
    datlist = Sys.glob(paste(datdir, '/*.dat', sep=''))
    infolist = Sys.glob(paste(infofiledir, '/*.info', sep=''))
    legendlist = Sys.glob(paste(legendfiledir, '/*.legend', sep=''))

    if(!file.exists(dest)) dir.create(dest)
    
    # get number of segments for each chromosome
    getnAppear = function(key) sum(grepl(key, datlist))

    # data.dat files start with these prefixes
    datkeys1 = paste('chr0', 1:9, '_', sep='')
    datkeys2 = paste('chr', 10:22, '_', sep='')
    datkeys = c(datkeys1, datkeys2)

    nAppear = sapply(datkeys, getnAppear)
    require(foreach)
    # create an index for the info and legend files
    idx = foreach(i=1:length(nAppear), .combine='c') %do% {
        rep(i, nAppear[i])
    }
    infoExpand = infolist[idx]
    legendExpand = legendlist[idx]
    stopifnot(length(infoExpand) == length(datlist))
    stopifnot(length(legendExpand) == length(datlist))

    # names for the new info files
    m = regexpr('chr\\d+_\\d+', datlist, perl=TRUE)
    destfiles = regmatches(datlist, m)
    destfiles = paste(dest, '/', destfiles, '.info', sep='')

    require(doMC)
    registerDoMC()
    options(cores=6)
    foreach(i=1:length(datlist)) %dopar% {
        createInfo.file(datlist[i], infoExpand[i], legendExpand[i], destfiles[i])
    }

    NULL
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
