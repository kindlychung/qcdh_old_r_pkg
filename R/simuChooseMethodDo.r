doChooseMethod = function() {
    distflist = Sys.glob("~/Phased/distm/chr*.dist.desc")
    h1flist = Sys.glob("~/Phased/bigmat/chr*mach1.hap1.desc")
    h2flist = Sys.glob("~/Phased/bigmat/chr*mach1.hap2.desc")
    infoflist = Sys.glob("~/Phased/info/chr*.info")
    # distf = "~/Phased/distm/chr22_01.dist.desc"
    # h1f = "~/Phased/bigmat/chr22_01_mach1.hap1.desc"
    # h2f = "~/Phased/bigmat/chr22_01_mach1.hap2.desc"
    # infof = "~/Phased/info/chr22_01.info"

    res = psfilelist(h1flist, h2flist, distflist, infoflist, b=.5)
    res
}
