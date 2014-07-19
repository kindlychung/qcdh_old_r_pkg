doskin20_692 = function() {
    g = read.csv('~/data/skin20_692/skn.gen.csv')
    g = g[, -(1:2)]
    g = t(g)

    phecovm = read.csv('~/data/skin20_692/skn.phe.csv')
    formulas = c("hskn~age+sex", "sskn~age+sex")

    registerDoMC()
    options(cores=2)

    res = foreach(fml=formulas) %do% {
        resfml = qcdh.region(as.formula(fml), covm=phecovm, family=gaussian(), h1win=g)
        resfml$firstsnp = row.names(g)[resfml$firstsnp]
        resfml$secondsnp = row.names(g)[resfml$secondsnp]
        resfml
    }

    names(res) = c('hskn', 'sskn')
    res
}
