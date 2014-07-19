# the haplom function for databel format
haplom.da.dir = function(machdir, dadir) {
    machlist = Sys.glob(file.path(machdir, 'chr*mach1.out'))

    registerDoMC()
    options(cores=6)

    foreach(i=1:length(machlist)) %dopar% {
        haplom.da(machlist[i], dadir)
    }
}

haplom.da = function(machphasef, dadir) {
    dadir = Sys.glob(dadir)
    if(!file.exists(dadir)) {
        dir.create(dadir)
    } else if(!file.info(dadir)$isdir) {
        stop('Target is not a directory!')
    }

    # big.matrix backing filenames
    require(tools)
    bn = basename(file_path_sans_ext(machphasef))
    bn1 = paste(bn, '.hap1', sep='')
    bn2 = paste(bn, '.hap2', sep='')
    datafile1  = paste(bn1, '.fvd', sep  = '')
    indexfile1 = paste(bn1, '.fvi', sep = '')
    datafile2  = paste(bn2, '.fvd', sep  = '')
    indexfile2 = paste(bn2, '.fvi', sep = '')
    n1 = file.path(dadir, bn1)
    n2 = file.path(dadir, bn2)

    # # remove old files, if they exist
    # targetfiles = file.path(dadir, c(datafile1, indexfile1, datafile2, indexfile2))
    # for(f in targetfiles) {
    #     if(file.exists(f)) {
    #         print(targetfiles)
    #         stop('Target file exists!')
    #     }
    # }


    hapmats = hapinfo(machphasef)
    hapmat1 = hapmats$h1
    hapmat2 = hapmats$h2

    require(DatABEL)
    hapm1 = matrix2databel(hapmat1, filename=n1, 
                           type='SHORT_INT', readonly=TRUE)
    hapm2 = matrix2databel(hapmat2, filename=n2, 
                           type='SHORT_INT', readonly=TRUE)


    # Sys.chmod(targetfiles, mode='440')

    res = list(hapm1=hapm1, hapm2=hapm2)
    class(res) = 'haplom.da'
    res
}
