# get the chromosome number from the info file name
chrnumber = function(infof) {
    chr_string = regmatches(infof, regexpr('chr\\d+.*\\.info', infof, perl=TRUE))
    chr = regmatches(chr_string, regexpr('\\d+', chr_string))
    chr = as.integer(chr)
    chr
}

