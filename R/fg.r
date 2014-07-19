fg1 = function(...) UseMethod('fg1')
fg2 = function(...) UseMethod('fg2')

# x1 = 1:4
# A1 = cbind(x1^4, x1^3, x1^2, x1)
# b1 = c(0, 1, 0, 0)
# coefs1 = solve(A1, b1)

# polynomial formula for collapsing genotypes, method 1
fg1.default = function(s1, s2) .25*(s1 + s2)^4 - 2*(s1 + s2)^3 + 4.75*(s1 + s2)^2 - 3*(s1 + s2)

fg1.data.frame = function(d) {
    stopifnot(ncol(d) == 2)
    d[, ncol(d) + 1] = fg1.default(d[, 1], d[, 2])
    d
}

fg1.matrix = function(d) {
    stopifnot(ncol(d) == 2)
    d[, ncol(d) + 1] = fg1.default(d[, 1], d[, 2])
    d
}

# x2 = 1:4
# A2 = cbind(x2^4, x2^3, x2^2, x2)
# b2 = c(0, 1, 1, 1)
# coefs2 = solve(A2, b2)

# polynomial formula for collapsing genotypes, method 2
fg2.default = function(s1, s2) 0.1250*(s1 + s2)^4 - 1.0833*(s1 + s2)^3 + 2.8750*(s1 + s2)^2 - 1.9167*(s1 + s2)

fg2.data.frame = function(d) {
    stopifnot(ncol(d) == 2)
    d[, ncol(d) + 1] = fg2.default(d[, 1], d[, 2])
    d
}

fg2.matrix = function(d) {
    stopifnot(ncol(d) == 2)
    d[, ncol(d) + 1] = fg2.default(d[, 1], d[, 2])
    d
}
