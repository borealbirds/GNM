## --- settings ---
## file name for data bundle, need to be in /data/ dir
fn <- "BAMdb-GNMsubset-2020-01-08.RData"
## project name for storing the output
PROJ <- "boot"
## output directory
OUTDIR <- if (interactive())
    paste0("d:/bam/BAM_data_v2019/gnm/out/", PROJ) else paste0("/scratch/psolymos/out/", PROJ)

cat("* Loading packages and sourcing functions:")
library(parallel)
library(mefa4)
library(gbm)

run_brt_boot <- function(b, spp) {
    sppbcr <- SPPBCR[grep(spp, SPPBCR)]
    bcr <- strsplit(RUN, "-")[[1]][2]
    if (!dir.exists(paste0(OUTDIR, "/", spp)))
        dir.create(paste0(OUTDIR, "/", spp))
    if (!dir.exists(paste0(OUTDIR, "/", spp, "/",bcr)))
        dir.create(paste0(OUTDIR, "/", spp, "/",bcr))
    for (RUN in sppbcr) {
        out <- .run_brt_boot(b, RUN)
        save(out, file=paste0(OUTDIR, "/", spp, "/",bcr, "/gnmboot-",
            spp, "-", bcr, "-", b, ".RData"))
    }
    invisible(TRUE)
}

## b: bootstrap id (>= 1)
## RUN: SPP-BCR_NO tag
.run_brt_boot <- function(b, RUN, verbose=interactive()) {
    t0 <- proc.time()["elapsed"]
    ## parse input
    tmp <- strsplit(RUN, "-")[[1L]]
    spp <- tmp[1L]
    BCR <- tmp[2L]
    bcr <- as.integer(strsplit(BCR, "_")[[1L]][2L])
    if (bcr %% 100 == 0) {
        ## US all clim+topo
        cn <- CN[[BCR]]
        nt <- ntmax[spp]
    } else {
        ## Canada: based on XV
        cn <- names(xvinfo[[RUN]]$varimp)[xvinfo[[RUN]]$varimp > 0]
        nt <- xvinfo[[RUN]]$ntree
    }

    ## create data subset for BCR unit
    ss <- dd[,BCR] == 1L
    DAT <- data.frame(
        count=as.numeric(yy[ss, spp]),
        offset=off[ss, spp],
        #weights=dd$wi[ss],
        cyid=dd$cyid[ss],
        YEAR=dd$YEAR[ss],
        ARU=dd$ARU[ss], # ARU added here, but not as layer
        dd2[ss, cn])
    ## subsample based on 2.5x2.5km^2 cell x year units
    DAT <- DAT[sample.int(nrow(DAT)),]
    DAT <- DAT[!duplicated(DAT$cyid),]
    if (b > 1)
        DAT <- DAT[sample.int(nrow(DAT), replace=TRUE),]
    ## 0 detection output
    if (sum(DAT$count) < 1) {
        out <- structure(
            sprintf("0 detections for %s in %s", spp, BCR),
            class="try-error")
    } else {
        if (verbose)
            cat("\nFitting gbm::gbm for", spp, "in", BCR, "/ b =", b, "/ n =",nrow(DAT), "\n")
        out <- try(gbm::gbm(DAT$count ~ . + offset(DAT$offset),
            data=DAT[,-(1:3)],
            n.trees = nt,
            interaction.depth = 3,
            shrinkage = 1/nt,
            bag.fraction = 0.5,
            #weights = DAT$weights,
            distribution = "poisson",
            var.monotone = NULL,#rep(0, length(4:ncol(DAT))),
            keep.data = FALSE,
            verbose = verbose,
            n.cores = 1))
    }
    attr(out, "__settings__") <- list(
        species=spp, region=BCR, iteration=b, elapsed=proc.time()["elapsed"]-t0)
    out
}

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
load(file.path("data", "xvinfo.RData"))
## get max tree size for each species to inform US models
tmp <- strsplit(names(xvinfo), "-")
x <- data.frame(sppbcr=names(xvinfo),
    spp=sapply(tmp, "[[", 1),
    bcr=sapply(tmp, "[[", 2),
    nt=sapply(xvinfo, "[[", "ntree"),
    vi=sapply(xvinfo, function(z)
        length(z$varimp[z$varimp > 0])))
x <- x[x$vi>0 & x$nt>0,]
x$pocc <- colMeans(yy>0)[as.character(x$spp)]
x$pbcr <- 0
for (i in levels(x$bcr)) {
    cm <- colMeans(yy[dd[[i]] > 0,]>0)
    x$pbcr[x$bcr == i] <- cm[as.character(x$spp[x$bcr == i])]
}
ntmax <- aggregate(x$nt, list(spp=x$spp), max)
ntmax <- structure(ntmax$x, names=as.character(ntmax$spp))

## Create the cluster with the nodes name.
## One process per count of node name.
## nodeslist = node1 node1 node2 node2, means we are starting 2 processes
## on node1, likewise on node2.
cat("* Spawning workers...")
cl <- makePSOCKcluster(nodeslist, type = "PSOCK")

cat("OK\n* Loading packages on workers ... ")
tmpcl <- clusterEvalQ(cl, library(mefa4))
tmpcl <- clusterEvalQ(cl, library(gbm))

cat("OK\n* Exporting and data loading on workers ... ")
if (interactive())
    tmpcl <- clusterEvalQ(cl, setwd("d:/bam/BAM_data_v2019/gnm"))
clusterExport(cl, c("dd", "dd2", "off", "yy", "CN",
    "PROJ", "xvinfo", "ntmax", "OUTDIR", "SPPBCR", ".run_brt_boot"))

cat("OK\n* Establishing checkpoint ... ")
SPP <- colnames(yy)
DONE <- list.files(OUTDIR)
TOGO <- setdiff(SPP, DONE)

cat("OK\n* Start running models:")
set.seed(as.integer(Sys.time()))
ncl <- if (interactive())
    nodeslist else length(nodeslist)
#while (length(TOGO) > 0) {
for (counter in 1:5) { # run only 5 species (~10hrs)
    spp <- sample(TOGO, 1)
    #res <- lapply(X=seq_len(ncl), fun=run_brt_boot, spp=spp)
    parLapply(cl=cl, X=seq_len(ncl), fun=run_brt_boot, spp=spp)
    DONE <- list.files(OUTDIR)
    TOGO <- setdiff(SPP, DONE)
}

## Releaseing resources.
cat("\n* Shutting down ... ")
stopCluster(cl)
cat("OK\nDONE!\n")
q("no")
