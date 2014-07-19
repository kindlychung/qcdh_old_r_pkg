fh2 = function(...) UseMethod('fh2')

# polynomial formula for collapsing haplotypes, method 2
fh2.default = function(x1, x2, y1, y2) {
    x = (x1 + y1)*(x2 + y2)
    0.125*x^3 - 0.875*x^2 + 1.750*x
}

fh2.data.frame = function(d) {
    stopifnot(ncol(d) == 4)
    d[, ncol(d) + 1] = fh2.default(d[, 1], d[, 2], d[, 3], d[, 4])
    d
}

fh2.matrix = function(d) {
    stopifnot(ncol(d) == 4)
    d[, ncol(d) + 1] = fh2.default(d[, 1], d[, 2], d[, 3], d[, 4])
    d
}
