.defaults <- list(
    "PopGenome" = list(
        "summary" = list(which = c("n.sites", "trans.transv.ratio")),
        "neutrality" = list(which = c("Tajima.D", "Fu.Li.F", "Fu.Li.D")),
        "fixed" = list(wrap.formula = "~ population"),
        "shared" = list(wrap.formula = "~ population"),
        "diversity" = list(which = c("nuc.diversity.within", "nuc.F_ST.vs.all", "Pi")),
        "diversity.between" = list(wrap.formula = "~ population"),
        "F_ST" = list(wrap.formula = "~ statistic", which = c ("nucleotide.F_ST", "Nei.G_ST")),
        "F_ST.pairwise" = list(which = c("nuc.F_ST.pairwise", "hap.F_ST.pairwise")),
        "segregating.sites" = list(wrap.formula = "~ population")
    )
)

.getOption <- function(application, statistics, key, default) {
    if (is.na(application))
        return (default)
    if (!(statistics %in% names(.defaults[[application]])))
        return (default)
    if (!(key %in% names(.defaults[[application]][[statistics]])))
        return (default)
    .defaults[[application]][[statistics]][[key]]
}



##' Plot GStats objects
##'
##' Plot GStats objects
##'
##' @param data GStats object
##' @param x variable to map to x aestethic
##' @param y variable to map to y aestethic
##' @param type plot type
##' @param xlim x limit
##' @param ylim y limit
##' @param main plot title
##' @param xlab x label
##' @param ylab y label
##' @param size plot size
##' @param colour colours to use
##' @param colour.var variable to map to colour aestethic
##' @param wrap wrap plots
##' @param wrap.formula wrap formula
##' @param wrap.ncol number of columns in facet wrap
##' @param compact.facet compact facet representation
##' @param strip.position strip position
##' @param scales scales
##' @param hide.legend whether or not to hide legend
##' @param hide.xaxis hide x axis tick marks and labels
##' @param grid include grid lines
##' @param text.size text size
##' @param text.x.angle x text angle
##' @param text.x.hjust x text horizontal justification
##' @param which statistic to plot
##' @param per.site plot statistics normalized by window length
##' @param zscore Plot Z-scores
##' @param ... extra arguments
##' @return ggplot
##' @author Per Unneberg
##' @export
##'
##' @import ggplot2
##' @importFrom RColorBrewer brewer.pal
##'
##' @describeIn gplot
##'
setMethod("gplot", c(data="GStats"),
          function(data, x="feature_id", y="value",
                   type="point",
                   xlim=NULL, ylim=NULL, main=paste(data@statistics, "statistics"),
                   xlab="bp", ylab=NULL, size=1,
                   colour=brewer.pal(3, "Dark2"), colour.var="seqnames",
                   wrap=TRUE, wrap.formula="",
                   pos="coordinate",
                   wrap.ncol=1, compact.facet=TRUE,
                   strip.position="right", scales="free_y",
                   hide.legend=TRUE, hide.xaxis=TRUE, grid=FALSE,
                   text.size=14, text.x.angle=45, text.x.hjust=1,
                   which=NULL, per.site=FALSE, zscore=FALSE,
                   ...) {
    if (wrap.formula == "")
        wrap.formula <- .getOption(data@application, data@statistics, "wrap.formula", "statistic ~ population")
    if (is.null(which))
        which <- .getOption(data@application, data@statistics, "which", which)
    df <- as.data.frame(asGRanges(data, per.site = per.site))
    which <- match.arg(which, levels(factor(df$statistic)), several.ok = TRUE)
    df <- subset(df, statistic %in% which)
    type <- match.arg(type, c("point", "line"))
    ## Make sure factor is ordered according to order of occurrence if not numeric
    if (!is.numeric(df[[x]])) {
        if (!is.factor(df[[x]])) df[[x]] <- factor(df[[x]], levels = unique(df[[x]]))
    }
    ## Convert to superscaffold coordinates
    pos <- match.arg(pos, c("coordinate", "index", "chr"))
    if (pos == "coordinate") {
        w0 <- width(rowRanges(data))[1]
        df$pos <- cumsum(width(rowRanges(data))) - w0 + 1 + width(rowRanges(data)) / 2
        x <- "pos"
    } else if (pos == "index") {
        xlab <- "window"
    } else if (pos == "chr") {
        xlab <- "seqnames"
    }
    ## By default scale by statistic
    if (zscore) df[[y]] <- unlist(tapply(df[[y]], df$statistic, scale))
    p <- ggplot(df, aes_string(x = x, y = y, colour = colour.var))
    if (wrap) p <- p + facet_wrap(as.formula(wrap.formula), ncol = wrap.ncol, strip.position = strip.position, scales = scales, ...)
    if (type == "point") {
        p <- p + geom_point(size = size)
    } else if (type == "line") {
        p <- p + geom_line(size = size)
    }
    if (!is.null(xlab)) p <- p + xlab(xlab)
    if (!is.null(ylab)) p <- p + ylab(ylab)
    if (!is.null(main)) {
        p <- p + ggtitle(main)
    }
    if (compact.facet) {
        p <- p + theme(panel.spacing = unit(0, "lines"))
    }
    if (hide.legend) p <- p + theme(legend.position = "none")
    if (!grid) p <- p + theme(panel.grid = element_blank())
    if (hide.xaxis) {
        p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    } else {
        p <- p + theme(axis.text.x = element_text(angle = text.x.angle, hjust = text.x.hjust))
    }
    p <- p + theme(text = element_text(size = text.size))
    # Colouring setup
    nlevels <- length(levels(df[[colour.var]]))
    nc <- length(colour)
    n <- floor(nlevels / nc)
    nmod <- nlevels %% nc
    if (nmod == 0) {
        p <- p + scale_colour_manual(values = rep(colour, n))
    } else {
        p <- p + scale_colour_manual(values = c(rep(colour, n), colour[1:nmod]))
    }
    p
})
