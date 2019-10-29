## --- settings ---
## file name for data bundle, need to be in /data/ dir
fn <- "BAMdb-GNMsubset-2019-10-29.RData"
## project name for storing the output
PROJ <- "run3"
## cli arguments
bcr <- u
SUB <- NULL
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    ## if 1 arg provided, it is BCR
    if (length(args) >= 1) {
        bcr <- args[1L]
        if (as.integer(bcr) == 0L)
            bcr <- u
    }
    ## if 2nd arg provided, it is subset size
    if (length(args) >= 2) {
        SUB <- as.integer(commandArgs(trailingOnly = TRUE)[2L])
    }
}
cat("* Using BCR:\n")
print(bcr)
if (!is.null(SUB))
    cat("* Sample size set to", SUB, "\n")

## test suite uses limited sets
TEST <- FALSE
if (TEST) {
    cat("* Note: this is a test run!\n")
}

cat("* Loading packages and sourcing functions:")
library(parallel)
library(mefa4)
library(gbm)
library(dismo)
#source("~/repos/abmianalytics/birds/00-functions.R")

## Create an array from the NODESLIST environnement variable
if (interactive()) {
    nodeslist <- 2
    setwd("d:/bam/BAM_data_v2019/gnm")
} else {
    cat("OK\n* Getting nodes list ... ")
    nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))
    cat("OK\n  Nodes list:\n")
    print(nodeslist)
}

cat("OK\n* Loading data on master ... ")
load(file.path("data", fn))
#load(file.path("data", "NTREE.RData"))

## Create the cluster with the nodes name.
## One process per count of node name.
## nodeslist = node1 node1 node2 node2, means we are starting 2 processes
## on node1, likewise on node2.
cat("* Spawning workers...")
cl <- makePSOCKcluster(nodeslist, type = "PSOCK")

cat("OK\n* Loading packages on workers ... ")
tmpcl <- clusterEvalQ(cl, library(mefa4))
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(dismo))

cat("OK\n* Exporting and data loading on workers ... ")
tmpcl <- clusterExport(cl, "fn")
if (interactive())
    tmpcl <- clusterEvalQ(cl, setwd("d:/bam/BAM_data_v2019/gnm"))
#tmpcl <- clusterEvalQ(cl, load(file.path("data", fn)))
clusterExport(cl, c("dd", "dd2", "off", "yy", "CN", "PROJ"))

cat("OK\n* Establishing checkpoint ... ")
SPP <- colnames(yy)
SPPBCRss <- SPPBCR[grep(paste0("BCR_", bcr), SPPBCR)]


#DONE <- character(0)
#if (interactive() || TEST)
#    SPPBCR <- SPPBCR[1:2]

DONE <- as.character(
    sapply(strsplit(list.files(paste0("out/", PROJ)), ".", fixed=TRUE), function(z) z[1L]))
#DONE <- unique(c(DONE0, DONE))
if (length(DONE) > 0) {
    cat("OK\n* Summary so far:\n")
    table(sapply(strsplit(gsub(".RData", "", DONE), "-"), "[[", 2))
}
TOGO <- setdiff(SPPBCRss, DONE)

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
        family = "poisson",
        tree.complexity = 3,
        learning.rate = RATE,
        bag.fraction = 0.5))
    if (is.null(out)) {
        out <- try(gbm.step(DAT,
            gbm.y = 1,
            gbm.x = 4:ncol(DAT),
            offset = DAT$offset, site.weights = DAT$weights,
            family = "poisson",
            tree.complexity = 3,
            learning.rate = RATE/10,
            bag.fraction = 0.5))
    }
    out
}

cat("OK\n* Start running models:")
set.seed(as.integer(Sys.time()))
while (length(TOGO) > 0) {
    #SET <- sample(TOGO)[seq_len(min(length(TOGO), length(cl)))]
    SET <- TOGO[seq_len(min(length(TOGO), length(cl)))]
    cat("\n  -", length(TOGO), "more pieces to go -", date(), "... ")
    if (interactive())
        flush.console()
    #system.time(hhh <- run_brt2("CAWA-BCR_6", ntree=100))
    #res <- lapply(X=SET, fun=run_brt2, SUB=SUB, RATE=0.001, ntree = 100)
    res <- parLapply(cl=cl, X=SET, fun=run_brt1, SUB=SUB, RATE=0.001)
    #res <- parLapply(cl=cl, X=SET, fun=run_brt2, SUB=SUB, RATE=0.001,
    #    ntree = if (interactive() || TEST) 100 else NULL)
    names(res) <- SET
    cat("OK")
    for (i in SET) {
        cat("\n    > Saving:", i, "... ")
        out <- res[[i]]
        #saveRDS(out, file=paste0("out/", PROJ, "/", if (TEST) "00test_" else "", i, ".RData"))
        save(out, file=paste0("out/", PROJ, "/", if (TEST) "00test_" else "", i, ".RData"))
        cat("OK")
    }
    DONE <- as.character(
        sapply(strsplit(list.files(paste0("out/", PROJ)), ".", fixed=TRUE), function(z) z[1L]))
    TOGO <- setdiff(SPPBCRss, DONE)
}

## Releaseing resources.
cat("\n* Shutting down ... ")
stopCluster(cl)
cat("OK\nDONE!\n")
q("no")
