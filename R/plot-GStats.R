##' Plot GStats objects
##'
##' Plot GStats objects
##'
##' @param data GStats object
##' @param x variable to map to x aestethic
##' @param y variable to map to y aestethic
##' @param colour colours to use
##' @param colour.var variable to map to colour aestethic
##' @param wrap wrap plots
##' @param wrap.formula wrap formula
##' @param wrap.ncol number of columns in facet wrap
##' @param plot.type plot type, either point or line
##' @param x.lab x label
##' @param y.lab y label
##' @param main plot title
##' @param compact.facet compact facet representation
##' @param strip.position strip position
##' @param scales scales
##' @param size plot size
##' @param hide.legend whether or not to hide legend
##' @param hide.xaxis hide x axis tick marks and labels
##' @param grid include grid lines
##' @param text.size text size
##' @param text.x.angle x text angle
##' @param text.x.hjust x text horizontal justification
##' @param ... extra arguments
##' @param which statistic to plot
##' @return ggplot
##' @author Per Unneberg
##' @export
setMethod("plot", "GStats",
          function(data, x="ranges", y="value",
                   colour=brewer.pal(3, "Dark2"), colour.var="seqnames",
                   wrap=TRUE, wrap.formula="key ~ population",
                   wrap.ncol=1, plot.type="point", x.lab="window",
                   y.lab=NULL, main=NULL,
                   compact.facet=TRUE, strip.position="right",
                   scales="free_y", size=1, hide.legend=TRUE,
                   hide.xaxis=TRUE, grid=FALSE,
                   text.size=14, text.x.angle=45, text.x.hjust=1,
                   ...) {
    plot.type <- match.arg(plot.type, c("point", "line"))
    data[[colour.var]] <- factor(data[[colour.var]], levels=unique(data[[colour.var]]))
    ## Make sure factor is ordered according to order of occurrence if not numeric
    if (!is.numeric(data[[x]])) {
        if (!is.factor(data[[x]])) data[[x]] <- factor(data[[x]], levels=unique(data[[x]]))
    }

    p <- ggplot(data, aes_string(x = x, y = y, colour = colour.var))
    if (wrap) p <- p + facet_wrap(as.formula(wrap.formula), ncol = wrap.ncol, strip.position = strip.position, scales = scales, ...)
    if (plot.type == "point") {
        p <- p + geom_point(size = size)
    } else if (plot.type == "line") {
        p <- p + geom_line(size = size)
    }
    if (!is.null(x.lab)) p <- p + xlab(x.lab)
    if (!is.null(y.lab)) p <- p + ylab(y.lab)
    if (!is.null(main)) p <- p + ggtitle(main)
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
    nlevels <- length(levels(data[[colour.var]]))
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
##' @title plot.GRanges
##' @describeIn plot.pg Make plot of GRanges object
##' @param per.site normalize statistics to per-site values
##' @export
plot.GRanges <- function(data, which=levels(factor(data$key)), per.site=TRUE, y="value", ...) {
    which <- match.arg(which, levels(factor(data$key)), several.ok = TRUE)
    df <- cbind(seqnames = as.character(seqnames(data)), as.data.frame(values(data)))
    df$ranges <- paste(as.character(seqnames(data)), start(data), "-", end(data), ":", sep = " ")
    if (per.site) df[[y]] <- df[[y]] / width(data)
    plot.pg(subset(df, key %in% which), y = y, ...)
}
##' @title plot.pg.summary
##' @describeIn plot.pg Make plot of summary
##' @export
plot.pg.summary <- function(data, x="seqnames", y="value", wrap.formula="~ key", which=c("n.sites", "trans.transv.ratio"), x.lab="scaffold", hide.xaxis=FALSE, ...) {
    which <- match.arg(which, c("n.sites", "n.biallelic.sites", "n.gaps", "n.unknowns", "n.valid.sites", "n.polyallelic.sites", "trans.transv.ratio"), several.ok = TRUE)
    plot.pg(subset(data, key %in% which), x = x, y = y, wrap.formula = wrap.formula, x.lab = x.lab, hide.xaxis = hide.xaxis, ...)
}
##' @title plot.pg.detail
##' @describeIn plot.pg Make plot of details
##' @export
plot.pg.detail <- function(data, which=c("MDG1", "MDG2", "MDSD"), ...) {
    which = match.arg(which, c("MDG1", "MDG2", "MDSD"), several.ok = TRUE)
    plot.pg(subset(data, key %in% which), ...)
}
##' @title plot.pg.neutrality
##' @describeIn plot.pg Make plot of neutrality statistics
##' @export
plot.pg.neutrality <- function(data, which=c("Tajima.D", "Fu.Li.F", "Fu.Li.F"), ...) {
    which = match.arg(which, c("Tajima.D", "n.segregating.sites", "Rozas.R_2", "Fu.Li.F", "Fu.Li.D", "Fu.F_S", "Fay.Wu.H", "Zeng.E", "Strobeck.S"), several.ok = TRUE)
    plot.pg(subset(data, key %in% which), ...)
}
##' @title plot.pg.fixed
##' @describeIn plot.pg Make plot of fixed sites
##' @export
plot.pg.fixed <- function(data, wrap.formula="~ population", ...) {
    plot.pg(data, wrap.formula = wrap.formula, ...)
}
##' @title plot.pg.shared
##' @describeIn plot.pg Make plot of shared sites
##' @export
plot.pg.shared <- function(data, wrap.formula="~ population", ...) {
    plot.pg(data, wrap.formula = wrap.formula, ...)
}
##' @title plot.pg.diversity
##' @describeIn plot.pg Make plot of diversity
##' @export
plot.pg.diversity <- function(data, which=c("nuc.diversity.within", "nuc.F_ST.vs.all", "Pi"),  ...) {
    which <- match.arg(which, c("hap.diversity.within", "hap.F_ST.vs.all", "nuc.diversity.within", "nuc.F_ST.vs.all", "Pi"), several.ok = TRUE)
    plot.pg(subset(data, key %in% which), ...)
}
##' @title plot.pg.diversity.between
##' @describeIn plot.pg Make plot of between population diversity
##' @export
plot.pg.diversity.between <- function(data, wrap.formula="~ population", colour="name", ...) {
    plot.pg(data, wrap.formula = wrap.formula, colour = colour, ...)
}
##' @title plot.pg.F_ST
##' @describeIn plot.pg Make plot of F_ST
##' @export
plot.pg.F_ST <- function(data, wrap.formula="~ key", which=c("nucleotide.F_ST", "Nei.G_ST"), ...) {
    which <- match.arg(which, c("haplotype.F_ST", "nucleotide.F_ST", "Nei.G_ST", "Hudson.G_ST", "Hudson.H_ST", "Hudson.K_ST"), several.ok = TRUE)
    plot.pg(subset(data, key %in% which), wrap.formula = wrap.formula, ...)
}
##' @title plot.pg.F_ST.pairwise
##' @describeIn plot.pg Make plot of pairwise F_ST
##' @export
plot.pg.F_ST.pairwise <- function(data, ...) {
    plot.pg(data, ...)
}
##' @title plot.pg.segregating.sites
##' @describeIn plot.pg Make plot of segregating sites
##' @export
plot.pg.segregating.sites <- function(data, wrap.formula="~ population", ...) {
    plot.pg(data, wrap.formula = wrap.formula, ...)
}

