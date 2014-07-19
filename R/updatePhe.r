# Update the given phenotype file with the given genotype vector,
# returns the phenotype table at the same time.
updatePhe = function(fn, geno) {
    geno = as.vector(geno)
    phe_dat = read.table(fn, T)
    stopifnot(nrow(phe_dat) == length(geno))
    phe_dat$gen = geno
    phe = sapply(1:10, function(i) phe_dat$gen * i / 10 + rnorm(length(phe_dat$gen)))
    phe_dat[, 1:10] = phe
    write.table(phe_dat, file=fn, quote=F, row.names=F, col.names=T)
    phe_dat
}
