library(parallel)
library(mefa4)
library(gbm)
library(dismo)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
PROJ <- "run3"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2019-10-29.RData"))

## CAWA, OSFL, RCKI, RUBL, MAWA
#SPP <- colnames(yy)
#SPP <- c("ALFL", "AMRO", "BOCH", "BTNW", "CAWA",
#    "OSFL", "OVEN", "RUBL", "WCSP", "YRWA")
#SPP <- "OSFL"
SPP <- "CAWA"
SPPBCRss <- SPPBCR[grep(spp, SPPBCR)]

DONE <- sapply(strsplit(list.files(file.path(ROOT, "out", PROJ)), ".", fixed=TRUE), function(z) z[1L])
if (length(DONE) > 0) {
    cat("OK\n* Summary so far:\n")
    table(sapply(strsplit(gsub(".RData", "", DONE), "-"), "[[", 2))
}
TOGO <- setdiff(SPPBCRss, DONE)

spp <- "CAWA"
TOGO <- paste0(spp, "-BCR_", u)

#RUN <- "OSFL-BCR_60"
run_brt1 <- function(RUN, SUB=NULL, RATE=0.001) {
    ## parse input
    tmp <- strsplit(RUN, "-")[[1L]]
    spp <- tmp[1L]
    BCR <- tmp[2L]
    ## create data subset
    ss <- dd[,BCR] == 1L
    if (sum(yy[ss, spp]) < 1)
        return(structure(
            sprintf("0 detections for %s in %s", spp, BCR),
            class="try-error"))
    DAT <- data.frame(
        count=as.numeric(yy[ss, spp]),
        offset=off[ss, spp],
        weights=dd$wi[ss],
        ARU=dd$ARU[ss], # ARU added here, but not as layer
        dd2[ss, CN[[BCR]]])
    if (!is.null(SUB)) {
        DAT$weights <- 1
        SUB <- min(SUB, nrow(DAT))
        ## this is not bootstrap resampling
        ## just subset to reduce memory footprint
        DAT <- DAT[sample(nrow(DAT), SUB, prob=DAT$weights),]
    }
    cat("\nFitting gbm.step for", spp, "in", BCR, "\n")
    cat("    Sample size =", nrow(DAT), "\n\n")
    out <- try(gbm.step(DAT,
        gbm.y = 1,
        gbm.x = 4:ncol(DAT),
        offset = DAT$offset, site.weights = DAT$weights,
        family = "poisson", tree.complexity = 3, learning.rate = RATE, bag.fraction = 0.5))
    out
}
run_brt2 <- function(RUN, SUB=NULL, RATE=0.001, ntree=1000) {
    ## parse input
    tmp <- strsplit(RUN, "-")[[1L]]
    spp <- tmp[1L]
    BCR <- tmp[2L]
    ## create data subset
    ss <- dd[,BCR] == 1L
    if (sum(yy[ss, spp]) < 1)
        return(structure(
            sprintf("0 detections for %s in %s", spp, BCR),
            class="try-error"))
    DAT <- data.frame(
        count=as.numeric(yy[ss, spp]),
        offset=off[ss, spp],
        weights=dd$wi[ss],
        ARU=dd$ARU[ss], # ARU added here, but not as layer
        dd2[ss, CN[[BCR]]])
    if (!is.null(SUB)) {
        DAT$weights <- 1
        SUB <- min(SUB, nrow(DAT))
        ## this is not bootstrap resampling
        ## just subset to reduce memory footprint
        DAT <- DAT[sample(nrow(DAT), SUB, prob=DAT$weights),]
    }
    cat("\nFitting gbm for", spp, "in", BCR, "\n")
    cat("    Sample size =", nrow(DAT), "\n")
    cat("    N trees     =", ntree, "\n\n")
    out <- try(gbm::gbm(DAT$count ~ . + offset(DAT$offset),
        data=DAT[,-(1:3)],
        n.trees = ntree,
        interaction.depth = 3,
        shrinkage = RATE,
        bag.fraction = 0.5,
        weights = DAT$weights,
        distribution = "poisson",
        var.monotone = NULL,#rep(0, length(4:ncol(DAT))),
        keep.data = FALSE,
        verbose = TRUE,
        n.cores = 1))
    out
}

set.seed(as.integer(Sys.time()))
for (i in TOGO) {
    gc()
    cat("\n  -", length(TOGO), "more pieces to go -", date(), "... ")
    flush.console()
    out <- run_brt1(i, RATE=0.001)
    if (is.null(out) || inherits(out, "try-error"))
        out <- run_brt1(i, RATE=0.0001)
    cat("\n    > Saving:", i, "... ")
    saveRDS(out, file=file.path(ROOT, "out", PROJ, paste0(i, ".rds")))
    cat("OK")
    DONE <- as.character(
        sapply(strsplit(list.files(file.path(ROOT, "out", PROJ)), ".", fixed=TRUE), function(z) z[1L]))
    TOGO <- setdiff(SPPBCRss, DONE)
}

