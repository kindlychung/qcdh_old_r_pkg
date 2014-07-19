# attach.big.matrix wrapper
abigm = function(path) {
    # getm = function(fn) {
    #     withCallingHandlers(attach.big.matrix(dget(fn)), warning=roHandler)
    # }
    # chroute(getm, path)
    withCallingHandlers(attach.big.matrix(path), warning=roHandler)
}

# Change working dir, then execute function, then switch
# to old working dir.
chroute = function(fun, path) {
    fn = basename(path)
    fdir = dirname(path)

    owd = getwd()
    setwd(fdir)

    res = fun(fn)
    setwd(owd)
    res
}

# no read-only warnings from attach.big.matrix
roHandler = function(w) {
    if(any(grepl('read-only', w))) {
        invokeRestart('muffleWarning')
    }
}
