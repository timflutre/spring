suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(plyr))
suppressMessages(library(reshape2))
suppressMessages(library(Hmisc))
suppressMessages(library(scales))

cov.plot <- function(cov.res, corr=FALSE) {

  if (corr) {
    cov.res <- cov2cor(cov.res)
  }

  cov.plot <- ggplot(melt(as.matrix(cov.res)), aes(Var1, Var2, fill = value)) + geom_tile(colour="gray80") +
      scale_fill_gradient2(name="correlation",limits=c(-max(cov.res),max(cov.res)),low = muted("blue"), mid = "white", high = muted("red"),
                           midpoint = 0.04, space = "Lab", na.value = "grey50", guide = "colourbar") +
    theme_bw() +  labs(x="", y="") +
        theme(legend.title  = element_text(size = 12), axis.text  = element_text(size = 12),
              plot.background = element_blank()
              ,panel.grid.major = element_blank()
              ,panel.grid.minor = element_blank()
              ,panel.border = element_blank()
              ,panel.background = element_blank()
              )
  print(cov.plot)
}

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

multiplot <- function(..., legend=FALSE, plotlist=NULL, cols) {
    require(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # Make the panel
    plotCols = cols                       # Number of columns of plots
    plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
        curCol = (i-1) %% plotCols + 1
        print(plots[[i]], vp = vplayout(curRow, curCol ))
    }
    if (legend) {
      grid.draw(g_legend(plots[[i]]))
    }

}

plot.map <- function(par.hat, gen.map, title="Selected features along the outcomes", plot=TRUE, show.names=TRUE,
                     chr=1:max(gen.map$chrom)) {

  ## restricting to some specified chromosoms
  par.hat <- par.hat[gen.map$chrom %in% chr, ]
  gen.map <- gen.map[gen.map$chrom %in% chr, ]

  pheno  <- colnames(par.hat)
  npheno <- length(pheno)

  ## retrieve information related to the selected gen.map
  dx <- do.call(rbind, apply(par.hat, 2, function(sel) gen.map[sel!=0, ]))
  uchrom <- sort(unique(dx$chrom))
  nchrom <- length(uchrom)
  ochrom <- rep(0, max(uchrom))
  ochrom[uchrom] <- seq(nchrom)
  dx$outcomes  <- factor(rep(pheno, sapply(apply(par.hat, 2, function(sel) gen.map[sel!=0, ]), nrow)), levels=pheno)
  dx$position  <-  ochrom[dx$chrom] + (as.numeric(dx$outcomes) * 2 - npheno - 1) *.35 /npheno
  dx$intensity <-  par.hat[par.hat !=0]
  dx$sign      <-  rep("+", length(dx$intensity))
  dx$sign[sign(dx$intensity) < 0] <- "-"
  dx$sign <- as.factor(dx$sign)
  dx$intensity <-  abs(dx$intensity)

  ## Constructing ther map for chromosome representation with annotation
  map.x <-rep(seq(nchrom), tapply(gen.map$loci,
                                  factor(gen.map$chrom, ordered=TRUE, levels=1:max(uchrom)), length)[uchrom])
  map.y <- gen.map$loci[gen.map$chrom %in% uchrom]
  map.n <- gen.map$names[gen.map$chrom %in% uchrom]
  map <- data.frame(xs.map=map.x-0.475, xe.map=map.x+0.475, y=map.y,
                    xs.names=map.x-0.5, label=map.n)

  ## Build the genetic Map and put annotation of the gen.map
  g <- ggplot(map) + geom_segment(aes(x=xs.map,xend=xe.map,y=y,yend=y), alpha=0.6, colour="#00cc33")
  if (show.names)
    g <- g + geom_text(data=map, aes(x=xs.names,y=y,label=label,hjust=0,vjust=0), size=4, alpha=0.6, colour="#00cc33")
  ## Set the points where features have been selected by
##  g <- g + geom_point(data=dx, aes(x=position,y=loci,group=outcomes,colour=outcomes, size=intensity,shape=sign))
  g <- g + geom_point(data=dx, aes(x=position,y=loci,group=outcomes,colour=outcomes, size=intensity))
  ## Various adjustements (axis, legends, labels...)
  g <- g + scale_x_discrete(breaks=1:nchrom, labels=uchrom, limits=1:nchrom)
  g <- g + guides(alpha=FALSE) + labs(x="chromosomes",y="loci",title=title)

  if (plot)
    print(g)

  return(g)
}
