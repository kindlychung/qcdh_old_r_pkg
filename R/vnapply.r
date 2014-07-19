vnapply = function(...) UseMethod('vnapply')

vnapply.character = function(xn, yn, func, expand, mcomb, ...) {
    stopifnot(is.character(yn))
    stopifnot(is.logical(expand))

    if(expand) {
        dat = expand.grid(xn=xn, yn=yn, stringsAsFactors=FALSE)
    } else {
        stopifnot(length(xn) == length(yn))
        dat = data.frame(xn=xn, yn=yn)
    }

    res = vnapply(dat, func, ...)
    res = do.call(mcomb, res)
    if(mcomb == 'c' & length(res) == nrow(dat)) {
        res = cbind(dat, res)
    }

    res
}

vnapply.matrix = vnapply.data.frame = function(dat, func, ...)  {
    stopifnot(ncol(dat) == 2)

    res = list()
    for(i in 1:nrow(dat)) {
        xni = dat[i, 1]
        yni = dat[i, 2]

        # get variables from names
        xi = get(xni)
        yi = get(yni)
        
        resi = func(xi, yi, ...)
        res[[length(res) + 1]] = resi
    }

    res
}

