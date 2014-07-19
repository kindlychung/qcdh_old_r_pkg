createDist = function(...) UseMethod('createDist')

createDist.dir = function(infodir, dirn, depth) {
    infolist = Sys.glob(file.path(infodir, '[a-zA-Z]*.info'))
    registerDoMC()
    options(cores=6)
    foreach(i=1:length(infolist), .combine='c') %dopar% {
        postable = read.table(infolist[i], head=T, colClasses=c(rep('NULL', 7), 'integer'))
        posvec = postable$position
        distfn = basename(file_path_sans_ext(infolist[i]))
        createDist.default(posvec, depth, dirn, distfn)
    }
}

createDist.file = function(infofile, ...) {
    postable = read.table(infofile, head=T, colClasses=c(rep('NULL', 7), 'integer'))
    posvec = postable$position
    createDist.default(posvec, ...)
}

# To create a distance matrix from a given vector
# of bp positions, save it in a big.matrix,
# returns the big.matrix.
# The result should be used by columns
createDist.default = function(posvec, depth, dirn, distfn) {
    distbin = paste(distfn, '.dist.bin', sep='')
    distdesc = paste(distfn, '.dist.desc', sep='')

    owd = getwd()
    setwd(dirn)
    if(file.exists(distbin)) file.remove(distbin)
    if(file.exists(distdesc)) file.remove(distdesc)
    setwd(owd)

    require(bigmemory)
    require(biganalytics)
    posvec = as.integer(posvec)
    posvec1 = c(posvec, rep(NA, depth))
    posm = foreach(i=1:length(posvec), .combine='cbind') %do% {
        posvec1[(i+1):(i+depth)]
    }
    posvecm = t(as.matrix(posvec))
    posvecm = posvecm[rep(1, depth), ]
    distm = posm - posvecm
    distm[is.na(distm)] = 0

    n = length(posvec)
    stopifnot(all(dim(distm) == c(depth, n)))

    options(bigmemory.typecast.warning=FALSE)
    distbigm = filebacked.big.matrix(nrow=depth, ncol=n,
                                  init=0,
                                  backingpath=dirn,
                                  backingfile=distbin,
                                  descriptorfile=distdesc,
                                  type='integer')
    distbigm[1:depth, 1:n] = distm
    file.path(dirn, distdesc)
}
