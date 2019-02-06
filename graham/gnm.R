## --- settings ---
## file name for data bundle, need to be in /data/ dir
fn <- "ab-birds-north-2019-01-30.RData"
## project name for storing the output
PROJ <- "north"

## CAIC = alpha * AIC + (1 - alpha) * BIC, 1: AIC, 0: BIC
CAICalpha <- 1
## Number of bootstrap runs, 100 or 240
MaxB <- 8*32 # 256

## test suite uses limited sets
TEST <- FALSE

if (TEST) {
    MaxB <- 2
    cat("* Note: this is a test run!\n")
}

cat("* Loading packages and sourcing functions:")
library(parallel)
library(mefa4)
library(opticut)
source("~/repos/abmianalytics/birds/00-functions.R")

## Create an array from the NODESLIST environnement variable
if (interactive()) {
    nodeslist <- 2
    BBB <- 2
    setwd("d:/abmi/AB_data_v2018/data/analysis/birds")
} else {
    cat("OK\n* Getting nodes list ... ")
    nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))
    BBB <- MaxB
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
tmpcl <- clusterEvalQ(cl, library(opticut))
tmpcl <- clusterEvalQ(cl, source("~/repos/abmianalytics/birds/00-functions.R"))

cat("OK\n* Exporting and data loading on workers ... ")
tmpcl <- clusterExport(cl, "fn")
if (interactive())
    tmpcl <- clusterEvalQ(cl, setwd("d:/abmi/AB_data_v2018/data/analysis/birds"))
#tmpcl <- clusterEvalQ(cl, load(file.path("data", fn)))
clusterExport(cl, c("DAT", "YY", "OFF", "BB", "SSH", "OFFmean"))

cat("OK\n* Establishing checkpoint ... ")
SPP <- colnames(YY)
DONE <- character(0)
if (interactive() | TEST)
    SPP <- SPP[1:2]

DONE <- substr(list.files(paste0("out/", PROJ)), 1, 4)
TOGO <- setdiff(SPP, DONE)

cat("OK\n* Start running models:")
set.seed(as.integer(Sys.time()))
while (length(TOGO) > 0) {
    SPP1 <- sample(TOGO, 1)
    cat("\n  -", length(DONE), "done,", length(TOGO), "more to go, doing", SPP1, "on", date(), "... ")
    if (interactive())
        flush.console()
    t0 <- proc.time()
    #z <- run_path1(1, "AMRO", mods, CAICalpha=1, wcol="vegw", ssh_class="vegc", ssh_fit="Space")
    if (interactive()) {
        res <- pblapply(cl=cl, X=1:BBB, FUN=run_path1,
            i=SPP1, mods=mods, CAICalpha=CAICalpha,
            wcol="vegw", ssh_class="vegca", ssh_fit="Space")
    } else {
        res <- parLapply(cl, 1:BBB, run_path1,
            i=SPP1, mods=mods, CAICalpha=CAICalpha,
            wcol="vegw", ssh_class="vegca", ssh_fit="Space")
    }
    attr(res, "timing") <- proc.time() - t0
    attr(res, "proj") <- PROJ
    attr(res, "spp") <- SPP1
    attr(res, "CAICalpha") <- CAICalpha
    attr(res, "date") <- as.character(Sys.Date())
    attr(res, "ncl") <- length(cl)
    save(res,
        file=paste0("out/", PROJ, "/", SPP1, ".RData"))
    DONE <- substr(list.files(paste0("out/", PROJ)), 1, 4)
    TOGO <- setdiff(SPP, DONE)
    cat("OK")
}

## Releaseing resources.
cat("\n* Shutting down ... ")
stopCluster(cl)
cat("OK\nDONE!\n")
q("no")
