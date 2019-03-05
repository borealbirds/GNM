## --- settings ---
## file name for data bundle, need to be in /data/ dir
fn <- "BAMdb-GNMsubset-2019-03-01.RData"
## project name for storing the output
PROJ <- "gnm"

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

## Create the cluster with the nodes name.
## One process per count of node name.
## nodeslist = node1 node1 node2 node2, means we are starting 2 processes
## on node1, likewise on node2.
cat("* Spawning workers...")
cl <- makePSOCKcluster(nodeslist, type = "PSOCK")

cat("OK\n* Loading data on master ... ")
load(file.path("data", fn))

cat("OK\nload packages on workers .. .")
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
#SPP <- colnames(yy)
SPP <- c("ALFL", "AMRO", "BOCH", "BTNW", "CAWA",
    "OSFL", "OVEN", "RUBL", "WCSP", "YRWA")
BCRlist <- paste0("BCR_", 4:14)
tmp <- expand.grid(BCR=BCRlist, SPP=SPP)
SPPBCR <- as.character(interaction(tmp$SPP, tmp$BCR, sep="-"))

DONE <- character(0)
if (interactive() || TEST)
    SPPBCR <- SPPBCR[1:2]

DONE <- sapply(strsplit(list.files(paste0("out/", PROJ)), ".", fixed=TRUE), function(z) z[1L])
TOGO <- setdiff(SPPBCR, DONE)

run_brt <- function(RUN, SAVE=TRUE, TEST=FALSE) {
    ## parse input
    tmp <- strsplit(RUN, "-")[[1L]]
    spp <- tmp[1L]
    BCR <- tmp[2L]
    ## create data subset
    ss <- dd[,BCR] == 1L
    DAT <- data.frame(
        count=as.numeric(yy[ss, spp]),
        offset=off[ss, spp],
        weights=dd$wi[ss],
        dd2[ss, c(cnf, CN[[BCR]])])
    if (TEST)
        DAT <- DAT[sample(nrow(DAT), 5000),]
    RATE <- 0.001
    ## fit BRT
    out <- try(gbm.step(DAT,
        gbm.y = 1,
        gbm.x = 4:ncol(DAT),
        offset = DAT$offset, site.weights = DAT$weights,
        family = "poisson", tree.complexity = 3, learning.rate = RATE, bag.fraction = 0.5))
    if (inherits(out, "try-error"))
        out <- try(gbm.step(DAT,
            gbm.y = 1,
            gbm.x = 4:ncol(DAT),
            offset = DAT$offset, site.weights = DAT$weights,
            family = "poisson", tree.complexity = 3, learning.rate = RATE/10, bag.fraction = 0.5))
    if (inherits(out, "try-error"))
        out <- try(gbm.step(DAT,
            gbm.y = 1,
            gbm.x = 4:ncol(DAT),
            offset = DAT$offset, site.weights = DAT$weights,
            family = "poisson", tree.complexity = 3, learning.rate = RATE/100, bag.fraction = 0.5))
    if (SAVE) {
        save(out, file=paste0("out/", PROJ, "/", RUN, ".RData"))
        return(!inherits(out, "try-error"))
    }
    out
}

cat("OK\n* Start running models:")
set.seed(as.integer(Sys.time()))
TOGO <- sample(TOGO)
cat("\n  -", length(DONE), "done,", length(TOGO), "more to go -", date(), "... ")

res <- parLapply(cl=cl, X=TOGO, fun=run_brt, SAVE=TRUE, TEST=interactive() || TEST)

DONE <- sapply(strsplit(list.files(paste0("out/", PROJ)), ".", fixed=TRUE), function(z) z[1L])
TOGO <- setdiff(SPPBCR, DONE)
cat("OK")

cat("\n* Result summary: ")
print(table(res))

## Releaseing resources.
cat("\n* Shutting down ... ")
stopCluster(cl)
cat("OK\nDONE!\n")
q("no")
