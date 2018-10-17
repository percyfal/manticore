#!/usr/bin/env Rscript
#
## Copyright (C) 2018 by Per Unneberg
##
## Author: Per Unneberg
##
## Description:
## Calculate genome statistics for a vcf file with PopGenome
##

library("getopt")
library("PopGenome")
library("nonmodelr")
library("tidyr")
library("dplyr")

spec <- matrix(
    c("help"         , "h", 0, "logical",
      "scaffold"     , "s", 1, "character",
      "vcf"          , "v", 1, "character",
      "gff"          , "a", 1, "character",
      "metadata"     , "m", 1, "character",
      "windowsize"   , "w", 2, "integer",
      "wdir"         , "d", 2, "character",
      "outfile"      , "o", 2, "character",
      "cpus"         , "c", 2, "integer"),
    ncol = 4, byrow = TRUE)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

## Setup analysis
curdir <- getwd()


if (is.null(opt$windowsize)){
    opt$windowsize <- 100000
}
if (is.null(opt$cpus)){
    opt$cpus <- 1
}
if (!is.null(opt$wdir)){
    wdir <- file.path(curdir, opt$wdir)
    message("Setting working directory to ", wdir)
    if (!dir.exists(wdir)) {
        dir.create(wdir, recursive = TRUE)
    }
    setwd(wdir)
}

## Load parallel package if cpus > 1
if (opt$cpus > 1) {
    library("parallel")
}

# Get length of scaffold
ln <- as.numeric(gsub(">", "", gsub(".*length=", "", system(paste0("bcftools view -h ", opt$vcf, " | grep \"ID=", opt$scaffold, ",\""), intern=TRUE))))

# Get samples present in vcf
samples <- system(paste0("bcftools query -l ", opt$vcf), intern=TRUE)

## Read sample metadata file and group samples by population
message("Reading sample metadata file ", opt$metadata)
metadata.df <- read.csv(opt$metadata, stringsAsFactors = FALSE)
metadata.df <- subset(metadata.df, SM %in% samples)
if (!("POOL" %in% colnames(metadata.df)))
    metadata.df$POOL <- FALSE

metadata.list <- as.list(subset(metadata.df, POOL == "False")[, c("SM", "POP")] %>% mutate(i=unlist(lapply(as.list(table(POP)), seq))) %>% spread(POP, SM) %>% select(-i))
populations.list <- lapply(names(metadata.list), function(x){as.vector(na.omit(metadata.list[[x]]))})
names(populations.list) <- names(metadata.list)


load_data <- function(scaffold) {
    message(paste("\nloading scaffold ", scaffold))
    ## Check whether scaffold exists
    gff_file <- opt$gff
    if (!scaffold %in% levels(gff.scaffolds$V1)) {{
        message(paste("Scaffold ", scaffold, " not in gff scaffolds; setting to FALSE"))
        gff_file <- FALSE
    }}
    ## Get length from vcf
    message("Opening vcf file ", opt$vcf)
    sink("/dev/null")
    vcf <- tryCatch({
        vcf <-readVCF(opt$vcf, 1000, frompos = 1, topos = ln,
                      tid = as.character(scaffold), gffpath = gff_file)
    }, error = function(err) {
        message("Scaffold not in vcf; returning NA")
        return (NA)
    })
    sink()
    if (inherits(vcf, "GENOME")) {
        vcf@region.names <- as.character(scaffold)
        vcf <- set.populations(vcf, populations.list, diploid = TRUE)
    }
    vcf
}

## Read gff file
message("Reading annotation file ", opt$gff)
gff.scaffolds <- read.table(opt$gff, sep = "\t",
                            colClasses = c("factor", rep("NULL", 8)))

## Load data
if (opt$cpus > 1) {
    message("loading data in parallel...")
    GENOME.classes <- parallel::mclapply(as.list(opt$scaffold),
                                         load_data,
                                         mc.cores = opt$cpus-1,
                                         mc.silent = TRUE,
                                         mc.preschedule = TRUE)
} else {
    message("loading data...")
    GENOME.classes <- lapply(as.list(opt$scaffold),
                             load_data)
}

## Concatenate data if more than one
message("Setting GENOME.class from GENOME.classes, length ",
        length(GENOME.classes))
if (length(GENOME.classes) > 1) {
    message("Concatenate data")
    GENOME.class <- concatenate.classes(GENOME.classes)
} else {
    GENOME.class <- GENOME.classes[[1]]
}

if (inherits(GENOME.class, "GENOME")) {
    message("Setting populations on big class: ", toString(populations.list))
    GENOME.class <- set.populations(GENOME.class, populations.list, diploid = TRUE)
    message("\nCalculating stats for GENOME.class")
    GENOME.class <- genomewide.stats(GENOME.class)
}

## Provide sliding window
message("Preparing windows")
if (ln >= opt$windowsize & !is.na(GENOME.class)) {
    GENOME.class.slide <- tryCatch({
        GENOME.class.slide <- sliding.window.transform(GENOME.class,
                                                       opt$windowsize,
                                                       opt$windowsize,
                                                       type = 2, whole.data = FALSE)
        ## Calculate statistics
        message("\nCalculating stats for GENOME.class.slide")
        GENOME.class.slide <- genomewide.stats(GENOME.class.slide)
    }, error = function(err) {
        message(paste0("Failed to create window size ", opt$windowsize, " for scaffold ", opt$scaffold, ", length ", ln))
        return(NA)
    })
} else {
    message("Scaffold too short for window analysis")
    GENOME.class.slide <- NA
}

## Finally save the output
message("Saving GENOME.class and GENOME.class.slide")
save(GENOME.class, GENOME.class.slide, file = opt$outfile)

## Cd back to original directory
message("cd back to ", curdir)
setwd(curdir)
