## --- settings ---
## file name for data bundle, need to be in /data/ dir
fn <- "BAMdb-GNMsubset-2019-03-01.RData"
## project name for storing the output
PROJ <- "gnm"
## cli arguments
bcr <- 4:14
SUB <- NULL
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    ## if 1 arg provided, it is BCR
    if (length(args) >= 1) {
        bcr <- args[1L]
        if (as.integer(bcr) == 0L)
            bcr <- 4:14
    }
    ## if 2nd arg provided, it is subset size
    if (length(args) >= 2) {
        SUB <- as.integer(commandArgs(trailingOnly = TRUE)[2L])
    }
}
cat("* Using BCR:\n")
print(bcr)

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

cat("OK\n* Loading data on master ... ")
load(file.path("data", fn))

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
clusterExport(cl, c("dd", "dd2", "off", "yy", "cnf", "CN", "PROJ"))



cat("OK\n* Establishing checkpoint ... ")
SPP <- colnames(yy)
#SPP <- c("ALFL", "AMRO", "BOCH", "BTNW", "CAWA",
#    "OSFL", "OVEN", "RUBL", "WCSP", "YRWA")
BCRlist <- paste0("BCR_", bcr)
tmp <- expand.grid(SPP=SPP, BCR=BCRlist)
SPPBCR <- as.character(interaction(tmp$SPP, tmp$BCR, sep="-"))



DONE <- character(0)
#if (interactive() || TEST)
#    SPPBCR <- SPPBCR[1:2]

DONE <- sapply(strsplit(list.files(paste0("out/", PROJ)), ".", fixed=TRUE), function(z) z[1L])
TOGO <- setdiff(SPPBCR, DONE)

run_brt <- function(RUN, TRY=1, TEST=FALSE, SUB=NULL) {
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
    if (TEST) {
        y <- as.numeric(yy[ss, spp])
        off <- off[ss, spp]
        w <- dd$wi[ss]
        DAT <- data.frame(dd2[ss, c(cnf, CN[[BCR]][1:5])])
        out <- try(glm(y ~ ., DAT, offset=off, weights=w, family=poisson))
    } else {
        DAT <- data.frame(
            count=as.numeric(yy[ss, spp]),
            offset=off[ss, spp],
            weights=dd$wi[ss],
            dd2[ss, c(cnf, CN[[BCR]])])
        if (!is.null(SUB)) {
            DAT$weights <- 1
            ## this is not bootstrap resampling
            ## just subset to reduce memory footprint
            DAT <- DAT[sample(nrow(DAT), SUB, prob=DAT$weights),]
            cat("Sample size =", nrow(DAT))
        }
        RATE <- 0.001
        ## fit BRT
        out <- try(gbm.step(DAT,
            gbm.y = 1,
            gbm.x = 4:ncol(DAT),
            offset = DAT$offset, site.weights = DAT$weights,
            family = "poisson", tree.complexity = 3, learning.rate = RATE, bag.fraction = 0.5))
        if (TRY > 1 && inherits(out, "try-error"))
            out <- try(gbm.step(DAT,
                gbm.y = 1,
                gbm.x = 4:ncol(DAT),
                offset = DAT$offset, site.weights = DAT$weights,
                family = "poisson", tree.complexity = 3, learning.rate = RATE/10, bag.fraction = 0.5))
        if (TRY > 2 && inherits(out, "try-error"))
            out <- try(gbm.step(DAT,
                gbm.y = 1,
                gbm.x = 4:ncol(DAT),
                offset = DAT$offset, site.weights = DAT$weights,
                family = "poisson", tree.complexity = 3, learning.rate = RATE/100, bag.fraction = 0.5))
    }
    out
}

cat("OK\n* Start running models:")
set.seed(as.integer(Sys.time()))
while (length(TOGO) > 0) {
    SET <- sample(TOGO)[seq_len(min(length(TOGO), length(cl)))]
    cat("\n  -", length(TOGO), "more species to go -", date(), "... ")
    if (interactive())
        flush.console()
    t0 <- proc.time()
    #hhh <- run_brt(SET[1], SUB=1000)
    res <- parLapply(cl=cl, X=SET, fun=run_brt, TEST=interactive() || TEST, SUB=SUB)
    names(res) <- SET
    cat("OK")
    for (i in SET) {
        cat("\n    > Saving:", i, "... ")
        out <- res[[i]]
        save(out, file=paste0("out/", PROJ, "/", if (TEST) "00test_" else "", i, ".RData"))
        cat("OK")
    }
    DONE <- sapply(strsplit(list.files(paste0("out/", PROJ)), ".", fixed=TRUE), function(z) z[1L])
    TOGO <- setdiff(SPPBCR, DONE)
}

## Releaseing resources.
cat("\n* Shutting down ... ")
stopCluster(cl)
cat("OK\nDONE!\n")
q("no")
