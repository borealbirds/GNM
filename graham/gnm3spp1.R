## --- settings ---
## file name for data bundle, need to be in /data/ dir
fn <- "BAMdb-GNMsubset-2019-10-29.RData"
## project name for storing the output
PROJ <- "run3spp1"
## no cli args
SUB <- NULL
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
SPP <- c("CAWA", "OSFL")
bcr <- u
SPPBCRss <- c(SPPBCR[grep(SPP[1], SPPBCR)], SPPBCR[grep(SPP[2], SPPBCR)])
TOGO <- SPPBCRss

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
SET <- TOGO
res <- parLapply(cl=cl, X=SET, fun=run_brt1, SUB=SUB, RATE=0.001)
names(res) <- SET
cat("OK")
for (i in SET) {
    cat("\n    > Saving:", i, "... ")
    out <- res[[i]]
    save(out, file=paste0("out/", PROJ, "/", if (TEST) "00test_" else "", i, ".RData"))
    cat("OK")
}

## Releaseing resources.
cat("\n* Shutting down ... ")
stopCluster(cl)
cat("OK\nDONE!\n")
q("no")
