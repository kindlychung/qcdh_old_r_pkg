haplom.bm.dir = function(machdir, bmdir) {
    machlist = Sys.glob(file.path(machdir, 'chr*mach1.out'))

    registerDoMC()
    options(cores=6)

    foreach(i=1:length(machlist)) %dopar% {
        haplom.bm(machlist[i], bmdir)
    }
}

haplom.bm = function(machphasef, bmdir) {
    if(!file.exists(bmdir)) {
        dir.create(bmdir)
    } else if(!file.info(bmdir)$isdir) {
        stop('Target is not a directory!')
    }

    # big.matrix backing filenames
    require(tools)
    bn = basename(file_path_sans_ext(machphasef))
    binf1  = paste(bn, '.hap1.bin', sep  = '')
    binf2  = paste(bn, '.hap2.bin', sep  = '')
    descf1 = paste(bn, '.hap1.desc', sep = '')
    descf2 = paste(bn, '.hap2.desc', sep = '')

    # remove old files, if they exist
    targetfiles = file.path(bmdir, c(binf1, binf2, descf1, descf2))
    for(f in targetfiles) {
        if(file.exists(f)) {
            print(targetfiles)
            stop('Target file exists!')
        }
    }


    hapmats = hapinfo(machphasef)
    hapmat1 = hapmats$h1
    hapmat2 = hapmats$h2

    require(bigmemory)
    require(biganalytics)
    as.big.matrix(hapmat1, type='short',
                  backingpath=bmdir,
                  backingfile=binf1,
                  descriptorfile=descf1)
    as.big.matrix(hapmat2, type='short',
                  backingpath=bmdir,
                  backingfile=binf2,
                  descriptorfile=descf2)

    # Sys.chmod(targetfiles, mode='440')

    res = list(hapm1=paste(bmdir, descf1, sep='/'),
               hapm2=paste(bmdir, descf2, sep='/'))
    class(res) = 'haplom.bm'
    res
}
