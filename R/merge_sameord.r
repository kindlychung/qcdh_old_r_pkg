# merge two data frames while conserving the order of the first one
merge_sameord = function(x, y, ...) {
    UseMethod('merge_sameord')
}

merge_sameord.data.frame = function(x, y, ...) {
    # a random string for the name of the new column
    rstr = paste(sample(c(0:9, letters, LETTERS), 12, replace=TRUE), collapse='')
    # the new column recording the row order
    x[, rstr] = 1:nrow(x)
    res = merge(x, y, all.x=TRUE, sort=FALSE, ...)
    res = res[order(res[, rstr]), ]
    res[, rstr] = NULL
    res
}


# shorttable = read.table('data.dat', colClasses=c('NULL', 'character'), stringsAsFactors=FALSE)
# longtable = read.table('GoNL.chr22_1.legend', colClasses=c('character', 'integer', 'NULL', 'NULL'), header=TRUE)
# colnames(shorttable)[1] = 'rs'
# colnames(longtable)[1] = 'rs'
# snptable = merge_sameord(shorttable, longtable)
# write.table(snptable, file='snptable.csv', quote=FALSE, row.names=FALSE, col.names=FALSE)
