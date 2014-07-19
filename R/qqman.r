# plotqcdh = function(qcdhres, pch=20, cex=.6, ...) {
#     ptype = qcdhres$ptype
#     stopifnot(sum(ptype=='single') == sum(ptype=='qcdhmin'))

#     logp  = -log10(qcdhres$'Pr(>|t|)')
#     ymin = min(logp)
#     ymax = max(logp)

#     idx = (ptype == 'single')
#     idx1 = which(idx)
#     idx2 = which(!idx)

#     plot(logp[idx1], ylim=c(ymin, ymax), ylab='-log P', xlab='Position', col='darkcyan', pch=pch, cex=cex, ...)
#     points(logp[idx2], col='brown4', pch=pch, cex=cex)
# }
manhattan <- function(dataframe, limitchromosomes=NULL,
                      pt.col=c('gray10','gray50'), pt.bg=c('gray10','gray50'), pt.cex=0.7,
                      pch=20, cex=.2, cex.axis=0.6,
                      gridlines=F, gridlines.col='gray83', gridlines.lty=1, gridlines.lwd=1,
                      ymax=NULL, ymax.soft=T,
                      annotate=TRUE, annotate.cex=0.7, annotate.font=3,
                      threshlines=FALSE,
                      suggestiveline=-log10(1e-5), suggestiveline.col='blue',
                      suggestiveline.lwd=1.5, suggestiveline.lty=1,
                      genomewideline=-log10(5e-8), genomewideline.col='red',
                      genomewideline.lwd=1.5, genomewideline.lty=1,
                      highlight=NULL, highlight.col=c('green3','magenta'),
                      highlight.bg=c('green3','magenta'),
                      filename=NULL, ...) {
	#============================================================================================
	######## Check data and arguments
    d = na.omit(dataframe) # omit NAs
    d$p = as.numeric(d$p)
    print(min(d$p))
    print(max(d$p))

    if (!("chr" %in% names(d) & "position" %in% names(d) & "p" %in% names(d))) {
        stop("Make sure your data frame contains columns chr, position, and p")
    }

	if(is.logical(annotate)){
		if(annotate == TRUE){
            idx = which(d$p < 5e-8)
			annotate <- as.character(d$snp[idx])
            if(length(annotate) == 0) annotate = NULL
			highlight <- annotate
            cat('These SNPs are annotated: \n')
            print(d[which(d$snp %in% annotate), ])
		} else {
			annotate <- NULL
		}
	} else if(is.vector(annotate)) {
		if ('snp' %in% names(d)){
			if (!all(annotate %in% d$snp)) {
                stop ("Annotate vector must be a subset of the snp column.")
            }
		} else {
			stop("Dataframe must have a column $snp with rs_ids to use annotate feature.")
		}
    }

    if(!is.null(limitchromosomes)){
        limitchromosomes = suppressWarnings(as.numeric(limitchromosomes))
    	if (TRUE %in% is.na(limitchromosomes)){
    		stop('limitchromosomes argument is not numeric')
    	} else {
    		d = d[d$chr %in% limitchromosomes, ]
    	}
    }


    ######################

    # Set positions, ticks, and labels for plotting
    d=subset(d[order(d$chr, d$position), ], (p>0 & p<=1)) # sort, and keep only 0<p<=1
    d$logp = -log10(d$p)
    d$pos=NA


    # Ymax
    if(is.null(ymax)){  # not numeric
    	ymax = ceiling(max(-log10(d$p)))
    }

    if (ymax.soft==T){ #if soft, ymax is just the lower limit for ymax
    	ymax = max(ymax, ceiling(max(-log10(d$p))))

    	# make ymax larger if top annotate snp is very high
    	if (!is.null(annotate)){
    		annotate.max = max(d[which(d$snp %in% annotate),]$logp)
    		if ((ymax - annotate.max) < 0.18*ymax){
    			ymax = annotate.max + 0.18*ymax
    		}
    	}
    } #else, ymax = ymax

	## Fix for the bug where one chromosome is missing. Adds index column #####
	d$index=NA
	ind = 0
	for (i in unique(d$chr)){
		ind = ind + 1
		d[d$chr==i,]$index = ind
	}
	########

    nchr=length(unique(d$chr))
    if (nchr==1) {
        d$pos=d$position
        ticks=floor(length(d$pos))/2+1
        xlabel = paste('Chromosome',unique(d$chr),'position')
        labs = ticks
    } else {
    	ticks = rep(NA,length(unique(d$chr))+1)
    	ticks[1] = 0
        for (i in 1:max(d$index)) {
          	d[d$index==i, ]$pos   =    (d[d$index==i, ]$position - d[d$index==i,]$position[1]) +1 +ticks[i]
    		ticks[i+1] = max(d[d$index==i,]$pos)
    	}
    	xlabel = 'Chromosome'
    	labs = append(unique(d$chr),'')
	}

    # Initialize plot
    xmax = max(d$pos) * 1.03
    xmin = max(d$pos) * -0.03
    #ymax = ceiling(ymax * 1.03)
    ymin = -ymax*0.03

    if(!is.null(filename)) {
        filename = as.character(filename)
        png(filename)

        plot(0,col=F,xaxt='n',bty='n',xaxs='i',yaxs='i',xlim=c(xmin,xmax), ylim=c(ymin,ymax),
                xlab=xlabel,ylab=expression(-log[10](italic(p))),las=1,cex.axis=cex.axis, cex=cex)

        # stagger labels
        blank = rep('',length(labs))
        lowerlabs = rep('',length(labs))
        upperlabs = rep('',length(labs))

        for (i in 1:length(labs)){
            if (i %% 2 == 0){
                lowerlabs[i] = labs[i]
            } else{
                upperlabs[i] = labs[i]
            }
        }

        axis(1,at=ticks,labels=blank,lwd=0,lwd.ticks=0,cex.axis=cex.axis, col='red')
        axis(1,at=ticks,labels=upperlabs,lwd=1,lwd.ticks=1,cex.axis=cex.axis,line=-0.25)
        axis(1,at=ticks,labels=lowerlabs,lwd=0,lwd.ticks=0,cex.axis=cex.axis,line=0.25, col='blue')

        yvals = par('yaxp')
        yinterval = par('yaxp')[2] / par('yaxp')[3]
        axis(2,at= (seq(0,(ymax+yinterval/2),yinterval) - yinterval/2),labels=F,lwd=0,lwd.ticks=1,cex.axis=cex.axis)

        # Gridlines
        if (isTRUE(gridlines)){

            abline(v=ticks,col=gridlines.col[1],lwd=gridlines.lwd,lty=gridlines.lty) #at ticks
            abline(h=seq(0,ymax,yinterval),col=gridlines.col[1],lwd=gridlines.lwd,lty=gridlines.lty) # at labeled ticks
            #abline(h=(seq(0,ymax,yinterval) - yinterval/2),col=gridlines.col[1],lwd=1.0) # at unlabeled ticks
        }

        # Points, with optional highlighting
        pt.col = rep(pt.col,max(d$chr))[1:max(d$chr)]
        pt.bg = rep(pt.bg,max(d$chr))[1:max(d$chr)]
        d.plain = d
        if (!is.null(highlight)) {
            if(class(highlight)!='character' & class(highlight)!='list'){
                stop('"highlight" must be a char vector (for 1 color) or list (for multi color).')
            }

            if (class(highlight)=='character'){ #if char vector, make list for consistency in plotting below
                highlight = list(highlight)
            }

            if ('snp' %in% names(d)){
                for (i in 1:length(highlight)){
                    if (FALSE %in% (highlight[[i]] %in% d$snp)) stop ("D'oh! Highlight vector/list must be a subset of the snp column.")
                }
            } else {
                stop("D'oh! Dataframe must have a column $snp with rs_ids to use highlight feature.")
            }

            highlight.col = rep(highlight.col,length(highlight))[1:length(highlight)]
            highlight.bg = rep(highlight.bg,length(highlight))[1:length(highlight)]

            for (i in 1:length(highlight)){
                d.plain = d.plain[which(!(d.plain$snp %in% highlight[[i]])), ]
            }
        }

        icol=1
        for (i in unique(d.plain$chr)) {
            with(d.plain[d.plain$chr==i, ],points(pos, logp, col=pt.col[icol],bg=pt.bg[icol],cex=pt.cex,pch=pch,...))
            icol=icol+1
        }

        if (!is.null(highlight)){
            for (i in 1:length(highlight)){
                d.highlight=d[which(d$snp %in% highlight[[i]]), ]
                with(d.highlight, points(pos, logp, col=highlight.col[i],bg=highlight.bg[i],cex=pt.cex,pch=pch,...))
            }
        }

        # Significance lines
        if(threshlines == TRUE) {
            if (is.numeric(suggestiveline)) abline(h=suggestiveline, col=suggestiveline.col[1],lwd=suggestiveline.lwd,lty=suggestiveline.lty)
            if (is.numeric(genomewideline)) abline(h=genomewideline, col=genomewideline.col[1],lwd=genomewideline.lwd,lty=genomewideline.lty)
        }

        # Annotate
        if (!is.null(annotate)){
            d.annotate = d[which(d$snp %in% annotate),]
            text(d.annotate$pos,
                (d.annotate$logp + 0.019*ymax),
                labels=d.annotate$snp,
                srt=45,
                cex=annotate.cex,adj=c(0,0.48),
                font=annotate.font
            )
        }


        dev.off()
    } else {
        stop('You should give me a filename to save the plot to!')
    }


	# Box
	# box()
}







## Make a pretty QQ plot of p-values ### Add ymax.soft
qq = function(d,
              gridlines=FALSE, gridlines.col='gray83',
              gridlines.lwd=1, gridlines.lty=1,
              confidence=T, confidence.col='gray81', 
              pt.cex=0.5, pt.col='black', pt.bg='black',
              pch=20, abline.col='red', abline.lwd=1.8, abline.lty=1,
              ymax=8, ymax.soft=T, 
              highlight=NULL, highlight.col=c('green3','magenta'),
              highlight.bg=c('green3','magenta'), 
              annotate=TRUE, annotate.cex=0.7, 
              annotate.font=3, cex.axis=0.95,
              filename, ...) {
	#======================================================================================================
	######## Check data and arguments; create observed and expected distributions
    d = na.omit(d) # remove NA, and non-numeric [which were converted to NA during as.numeric()]
    d = d[which(d$p>0 & d$p<1), ] # only Ps between 0 and 1

	if(is.logical(annotate)){
		if(annotate == TRUE){
            idx = which(d$p < 5e-8)
			annotate <- as.character(d$snp[idx])
            if(length(annotate) == 0) annotate = NULL
			highlight <- annotate
            cat('These SNPs are annotated: \n')
            print(d[which(d$snp %in% annotate), ])
		} else {
			annotate <- NULL
		}
	} else if(is.vector(annotate)) {
		if ('snp' %in% names(d)){
			if (!all(annotate %in% d$snp)) {
                stop ("Annotate vector must be a subset of the snp column.")
            }
		} else {
			stop("Dataframe must have a column $snp with rs_ids to use annotate feature.")
		}
    }


	d = d[order(d$p,decreasing=F), ] # sort
    d$logp = -log10(d$p)
	# o = -log10(d)
    d$expectedP = -log10(ppoints(length(d)))
    # e = -log10( ppoints(length(d) ))

	# Ymax
	if (!is.numeric(ymax) | ymax<max(d$logp)) ymax <- max(d$logp)

	################################

	# Initialize plot
	#print('Setting up plot.')
	#print(ymax)
	xspace = 0.078
	xmax = max(d$expectedP) * 1.019
    xmin = max(d$expectedP) * -0.035
    #ymax = ceiling(ymax * 1.03)
    ymin = -ymax*0.03
	plot(0,xlab=expression(Expected~~-log[10](italic(p))),ylab=expression(Observed~~-log[10](italic(p))),
			col=F,las=1,xaxt='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty='n',xaxs='i',yaxs='i',cex.axis=cex.axis)
	axis(side=1,labels=seq(0,max(e),1),at=seq(0,max(e),1),cex.axis=cex.axis,lwd=0,lwd.ticks=1)

	# Grid lines
	if (isTRUE(gridlines)){
		yvals = par('yaxp')
		yticks = seq(yvals[1],yvals[2],yvals[2]/yvals[3])
		abline(v=seq(0,max(e),1),col=gridlines.col[1],lwd=gridlines.lwd,lty=gridlines.lty)
		abline(h=yticks,col=gridlines.col[1],lwd=gridlines.lwd,lty=gridlines.lty)
	}

	 #Confidence intervals
	 find_conf_intervals = function(row){
	 	i = row[1]
	 	len = row[2]
	 	if (i < 10000 | i %% 100 == 0){
	 		return(c(-log10(qbeta(0.95,i,len-i+1)), -log10(qbeta(0.05,i,len-i+1))))
	 	} else { # Speed up
	 		return(c(NA,NA))
	 	}
	 }

	 # Find approximate confidence intervals
	if (isTRUE(confidence)){
		#print('Plotting confidence intervals.')
		ci = apply(cbind( 1:length(e), rep(length(e),length(e))), MARGIN=1, FUN=find_conf_intervals)
	 	bks = append(seq(10000,length(e),100),length(e)+1)
		for (i in 1:(length(bks)-1)){
	 		ci[1, bks[i]:(bks[i+1]-1)] = ci[1, bks[i]]
	 		ci[2, bks[i]:(bks[i+1]-1)] = ci[2, bks[i]]
		}
		colnames(ci) = names(e)
		# Extrapolate to make plotting prettier (doesn't affect intepretation at data points)
		slopes = c((ci[1,1] - ci[1,2]) / (e[1] - e[2]), (ci[2,1] - ci[2,2]) / (e[1] - e[2]))
		extrap_x = append(e[1]+xspace,e) #extrapolate slightly for plotting purposes only
		extrap_y = cbind( c(ci[1,1] + slopes[1]*xspace, ci[2,1] + slopes[2]*xspace), ci)

		polygon(c(extrap_x, rev(extrap_x)), c(extrap_y[1,], rev(extrap_y[2,])),col = confidence.col[1], border = confidence.col[1])
	}

	# Points (with optional highlighting)
	#print('Plotting data points.')
	fills = rep(pt.bg,length(o))
	borders = rep(pt.col,length(o))
	names(fills) = names(borders) = names(o)
	if (!is.null(highlight)){
		borders[highlight] = rep(NA,length(highlight))
		fills[highlight] = rep(NA,length(highlight))
	}
	points(e,o,pch=pch,cex=pt.cex,col=borders,bg=fills)

	if (!is.null(highlight)){
		points(e[highlight],o[highlight],pch=pch,cex=pt.cex,col=highlight.col,bg=highlight.bg)
	}

	#Abline
	abline(0,1,col=abline.col,lwd=abline.lwd,lty=abline.lty)

	# Annotate snps
	if (!is.null(annotate)){
		x = e[annotate] # x will definitely be the same
		y = -0.1 + apply(rbind(o[annotate],ci[1,annotate]),2,min)
		text(x,y,labels=annotate,srt=90,cex=annotate.cex,adj=c(1,0.48),font=annotate.font)
	}
	# Box
	box()
}
