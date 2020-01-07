## extracting XV info
fl <- list.files("d:/bam/BAM_data_v2019/gnm/out/run3")
names(fl) <- gsub(".RData", "", fl)

xvinfo <- list()
for (i in names(fl)) {

    cat(i, which(names(fl) == i), "/", length(fl), "\n")
    flush.console()

    e <- new.env()
    load(paste0("d:/bam/BAM_data_v2019/gnm/out/run3/", fl[i]), envir=e)

    if (is.null(e$out) || inherits(e$out, "try-error")) {
        xvinfo[[i]] <- list(
            ntree=0,
            varimp=character(0))
    } else {
        xvinfo[[i]] <- list(
            ntree=e$out$n.trees,
            varimp=structure(e$out$contributions$rel.inf, names=rownames(e$out$contributions)))
    }

}
save(xvinfo, file="d:/bam/BAM_data_v2019/gnm/out/xvinfo.RData")

nt <- sapply(xvinfo, "[[", "ntree")
nc <- sapply(xvinfo, function(z) length(z$varimp))
summary(nt[nt>0])
sum(nt==0)
sum(nt==10000)
summary(nc[nc>0])
hist(nt[nt>0])
hist(nc[nc>0])

library(parallel)
library(mefa4)
library(gbm)
library(dismo)

fn <- "BAMdb-GNMsubset-2019-10-29.RData"
PROJ <- "run3"
load(file.path("d:/bam/BAM_data_v2019/gnm", "data", fn))
load("d:/bam/BAM_data_v2019/gnm/out/xvinfo.RData")
load("d:/bam/BAM_data_v2019/gnm/data/cid.RData")

dd$cid <- cid[rownames(dd)]
dd$cyid <- interaction(dd$cid, dd$YEAR, sep="_", drop=TRUE)

z <- data.frame(reg=colnames(dd)[grepl("BCR_", colnames(dd))],
    n=sapply(colnames(dd)[grepl("BCR_", colnames(dd))], function(i) {
        sum(dd[,i] == 1L)
    }),
    n_cy=sapply(colnames(dd)[grepl("BCR_", colnames(dd))], function(i) {
        nlevels(droplevels(dd$cyid[dd[,i] == 1L]))
    }))

RUN <- "ALFL-BCR_70"

## b: bootstrap id (>= 1)
## RUN: SPP-BCR_NO tag
run_brt_boot <- function(b, RUN, verbose=interactive()) {
    t0 <- proc.time()["elapsed"]
    ## parse input
    tmp <- strsplit(RUN, "-")[[1L]]
    spp <- tmp[1L]
    BCR <- tmp[2L]
    nt <- xvinfo[[RUN]]$ntree
    cn <- names(xvinfo[[RUN]]$varimp)[xvinfo[[RUN]]$varimp > 0]
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

x <- run_brt_boot(1, RUN)
attr(x, "__settings__")

