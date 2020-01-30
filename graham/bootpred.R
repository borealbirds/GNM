cat("\n\n--- BOOTPRED started ---\n\n* Preliminary stuff:")
library(parallel)
library(gbm)
library(raster)
library(rgdal)
## variables
PROJ <- "boot"
if (interactive()) {
    ROOT1 <- "d:/bam/BAM_data_v2019/gnm"
    ROOT2 <- "d:/bam/BAM_data_v2019/gnm"
    SPP <- list.files(file.path(ROOT1, "out", PROJ))
} else {
    ROOT1 <- "/home/psolymos/bam"
    ROOT2 <- "/scratch/psolymos"
    SPP <- list.files(file.path(ROOT2, "out", PROJ))
}
B <- 32
u <- c(4, 5, 60, 61, 70, 71, 80, 81, 82, 83, 9, 10, 11, 12, 13, 14,
    200, 400, 500, 1000, 1100, 1200, 1300, 1400)
u <- u[u < 200]
if (interactive()) {
    cat("OK\n* Interactive testing mode only, using small subsets\n")
    SPP <- SPP[1:2]
    u <- u[1:2]
    B <- 2
}
## functions
predict_gbm <- function(brt, ND, r1, impute=0) {
    if (inherits(brt, "try-error") || is.null(brt)) {
        rp <- r1[[1]]
        values(rp)[!is.na(values(rp))] <- impute
    } else {
        if (!("ROAD" %in% colnames(ND$data)))
            ND$data$ROAD <- 0
        if (!("ARU" %in% colnames(ND$data)))
            ND$data$ARU <- 0
        z0 <- suppressWarnings(predict.gbm(brt, ND$data, type="response", n.trees=brt$n.trees))
        ## dealing with snow/ice
        z0[ND$data$nalc %in% c(0, 18, 19)] <- impute
        ## expand to extent
        z <- rep(NA, ND$n)
        z[ND$subset] <- z0
        zz <- matrix(t(z), ND$dim[1], ND$dim[2], byrow=TRUE)
        rp <- raster(x=zz, template=r1)
    }
    rp
}
pred_fun <- function (b, spp, BCR, PROJ) {
    fin <- file.path(ROOT2, "out", PROJ, spp, paste0("BCR_", BCR),
        paste0("gnmboot-", spp, "-BCR_", BCR, "-", b, ".RData"))
    e <- new.env()
    load(fin, envir=e)
    rrr <- predict_gbm(e$out, ND, r1, 0)
    if (!dir.exists(file.path(ROOT2, "out", "parts", spp)))
        dir.create(file.path(ROOT2, "out", "parts", spp))
    fout <- file.path(ROOT2, "out", "parts", spp,
        paste0("pred-", spp, "-BCR_", BCR, "-", PROJ, "-", b, ".tif"))
    writeRaster(rrr, fout, overwrite=TRUE)
    invisible(TRUE)
}
load_fun <- function () {
    r1 <- raster(file.path(ROOT1, "data", "templates", paste0("bcr-template-", BCR, ".grd")))
    ND <- readRDS(file.path(ROOT1, paste0("STACK-ND-BCR_", BCR, ".rds")))
    ND$data$YEAR <- 2011
    assign("r1", r1, envir=.GlobalEnv)
    assign("ND", ND, envir=.GlobalEnv)
    invisible(TRUE)
}

## Create an array from the NODESLIST environnement variable
if (interactive()) {
    nodeslist <- 2
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
tmpcl <- clusterEvalQ(cl, library(raster))
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(rgdal))

cat("OK\n* Establishing checkpoint ... ")
DONE <- list.files(file.path(ROOT2, "out", "parts"))
TOGO <- setdiff(SPP, DONE)

cat("OK\n* Start running predictions:\n")
set.seed(as.integer(Sys.time()))
ncl <- if (interactive())
    nodeslist else length(nodeslist)
#while (length(TOGO) > 0) {
for (counter in 1:10) {
    spp <- sample(TOGO, 1)
    cat("* Doing species", spp, "\n")
    for (BCR in u) {
        cat("\t- species", spp, "in BCR", BCR, "...")
        clusterExport(cl, c("load_fun", "predict_gbm", "ROOT1", "BCR", "ROOT2"))
        clusterEvalQ(cl, load_fun())
        ## predict for runs: use patLapply here
        parLapply(cl=cl, X=seq_len(ncl), fun=pred_fun,
            spp=spp, BCR=BCR, PROJ=PROJ)
        #for (b in 1:B) {
        #    pred_fun(b=b, spp=spp, BCR=BCR, PROJ=PROJ)
        #}
        cat(" OK\n")
    }
}

## Releaseing resources.
cat("\n* Shutting down ... ")
stopCluster(cl)
cat("OK\nBye!\n\n")
q("no")
