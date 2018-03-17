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
    if (!(statistics %in% names(.defaults[[application]])))
        return (default)
    if (!(key %in% names(.defaults[[application]][[statistics]])))
        return (default)
    .defaults[[application]][[statistics]][[key]]
}

##' Make box/violin plot of GStats object
##'
##' Make a box/violin plot of data
##'
##' @param data long format from GStats function
##' @param formula a formula, such as y ~ grp
##' @param type plot type, either point or line
##' @param xlim x limit
##' @param ylim y limit
##' @param main plot title
##' @param xlab x label
##' @param ylab y label
##' @param colour colours to use
##' @param colour.var variable to map to colour aestethic
##' @param wrap wrap plots
##' @param wrap.formula wrap formula
##' @param wrap.ncol number of columns in facet wrap
##' @param compact.facet compact facet representation
##' @param strip.position strip position
##' @param scales scales
##' @param hide.legend whether or not to hide legend
##' @param grid include grid lines
##' @param text.size text size
##' @param text.x.angle angle of x tick labels
##' @param text.x.hjust horizontal adjustment of x tick labels
##' @param which statistic to plot
##' @param per.site normalize statistic by window length
##' @param ... arguments passed on to facet_wrap
##'
##' @return ggplot object
##' @author Per Unneberg
##'
##' @export
##' @rdname gboxplot
##'
##' @import ggplot2
##' @importFrom stats as.formula
##' @importFrom RColorBrewer brewer.pal
##'
setMethod("gboxplot", "GStats",
          function(data, formula = "value ~ population",
                   type="box",
                   xlim=NULL, ylim=NULL, main=paste(data@statistics, "statistics"),
                   xlab="group", ylab=NULL,
                   colour=brewer.pal(3, "Dark2"), colour.var="seqnames",
                   wrap=FALSE, wrap.formula="",
                   wrap.ncol=1,
                   compact.facet=TRUE, strip.position="right",
                   scales="free_y", hide.legend=TRUE, grid=FALSE,
                   text.size=14, text.x.angle=45, text.x.hjust=1,
                   which=NULL, per.site=FALSE,
                   ...) {
    stopifnot(data@application %in% names(.defaults))
    if (wrap.formula == "")
        wrap.formula <- .getOption(data@application, data@statistics, "wrap.formula", "~ statistic")
    if (is.null(which))
        which <- .getOption(data@application, data@statistics, "which", which)
    df <- as.data.frame(asGRanges(data, per.site = per.site))
    which <- match.arg(which, levels(factor(df$statistic)), several.ok = TRUE)
    df <- subset(df, statistic %in% which)
    type <- match.arg(type, c("box", "violin"))
    y.var <- as.character(as.list(as.formula(formula))[[2]])
    x.var <- as.character(as.list(as.formula(formula))[[3]])
    colour.var <- colour.var %||% x.var
    df[[colour.var]] <- factor(df[[colour.var]], levels = unique(df[[colour.var]]))
    p <- ggplot(df, aes_string(x = x.var, y = y.var, colour = colour.var))
    if (wrap) p <- p + facet_wrap(as.formula(wrap.formula), ncol = wrap.ncol, strip.position = strip.position, scales = scales, ...)
    if (type == "box") {
        p <- p + geom_boxplot()
    } else if (type == "violin") {
        p <- p + geom_violin()
    }
    if (!is.null(xlab)) p <- p + xlab(xlab)
    if (!is.null(ylab)) p <- p + ylab(ylab)
    if (!is.null(main)) p <- p + ggtitle(main)
    if (compact.facet) {
        p <- p + theme(panel.spacing = unit(0, "lines"))
    }
    if (hide.legend) {
         p <- p + theme(legend.position = "none")
    }
    if (!grid) p <- p + theme(panel.grid = element_blank())
    p <- p + theme(text = element_text(size = text.size), axis.text.x = element_text(angle=text.x.angle, hjust=text.x.hjust))
    # Colouring setup
    nlevels <- length(levels(df[[colour.var]]))
    nc <- length(colour)
    n <-  floor(nlevels / nc)
    nmod <- nlevels %% nc
    if (nmod == 0) {
        p <- p + scale_colour_manual(values = rep(colour, n))
    } else {
        p <- p + scale_colour_manual(values = c(rep(colour, n), colour[1:nmod]))
    }
    p
})
