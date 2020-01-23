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

## bootstrap based estimation (no xv)

library(parallel)
library(mefa4)
library(gbm)

fn <- "BAMdb-GNMsubset-2020-01-08.RData"
PROJ <- "run3"
load(file.path("d:/bam/BAM_data_v2019/gnm", "data", fn))
load("d:/bam/BAM_data_v2019/gnm/out/xvinfo.RData")

## number of surveys vs cid/year blocks
z <- data.frame(reg=colnames(dd)[grepl("BCR_", colnames(dd))],
    n=sapply(colnames(dd)[grepl("BCR_", colnames(dd))], function(i) {
        sum(dd[,i] == 1L)
    }),
    n_cy=sapply(colnames(dd)[grepl("BCR_", colnames(dd))], function(i) {
        nlevels(droplevels(dd$cyid[dd[,i] == 1L]))
    }))

## explore XV results
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
summary(x)

if (FALSE) {
hist(x$nt)
hist(x$vi)
plot(nt ~ pocc,x)
plot(nt ~ pbcr,x)

m10 <- lm(nt ~ 1, x)
m11 <- lm(nt ~ bcr, x)
m12 <- lm(nt ~ spp, x)
m13 <- lm(nt ~ spp + bcr, x)
m14 <- lm(nt ~ pocc, x)

m20 <- lm(vi ~ 1, x)
m21 <- lm(vi ~ bcr, x)
m22 <- lm(vi ~ spp, x)
m23 <- lm(vi ~ spp + bcr, x)

AIC(m10,m11,m12, m13, m14) # ntree depends on the species, region is negligible
c(summary(m10)$r.squared,
summary(m11)$r.squared,
summary(m12)$r.squared,
summary(m13)$r.squared)

AIC(m20,m21,m22, m23) # no. of vars depends on the region more than on spp
c(summary(m20)$r.squared,
summary(m21)$r.squared,
summary(m22)$r.squared,
summary(m23)$r.squared)
}

ntmax <- aggregate(x$nt, list(spp=x$spp), max)
ntmax <- structure(ntmax$x, names=as.character(ntmax$spp))

#RUN <- "ALFL-BCR_70"
#b=1

## b: bootstrap id (>= 1)
## RUN: SPP-BCR_NO tag
run_brt_boot <- function(b, RUN, verbose=interactive()) {
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

x <- run_brt_boot(1, RUN)
attr(x, "__settings__")

spp <- "OVEN"
res <- list()
for (i in u) {
    cat("Fitting gbm::gbm for", spp, "in BCR", i, "\n")
    RUN <- paste0(spp, "-BCR_", i)
    res[[paste0("BCR_", i)]] <- run_brt_boot(1, RUN, verbose=FALSE)
}
save(res, file="d:/bam/BAM_data_v2019/gnm/out/test/oven.RData")

spp <- "OSFL"
res <- list()
for (i in u) {
    cat("Fitting gbm::gbm for", spp, "in BCR", i, "\n")
    RUN <- paste0(spp, "-BCR_", i)
    res[[paste0("BCR_", i)]] <- run_brt_boot(1, RUN, verbose=FALSE)
}
save(res, file="d:/bam/BAM_data_v2019/gnm/out/test/osfl.RData")

spp <- "CAWA"
res <- list()
for (i in u) {
    cat("Fitting gbm::gbm for", spp, "in BCR", i, "\n")
    RUN <- paste0(spp, "-BCR_", i)
    res[[paste0("BCR_", i)]] <- run_brt_boot(1, RUN, verbose=FALSE)
}
save(res, file="d:/bam/BAM_data_v2019/gnm/out/test/cawa.RData")


load("d:/bam/BAM_data_v2019/gnm/out/test/oven.RData")
oven <- data.frame(time=sapply(res, function(z) attr(z, "__settings__")$elapsed),
    size=sapply(res, object.size))

load("d:/bam/BAM_data_v2019/gnm/out/test/osfl.RData")
osfl <- data.frame(time=sapply(res, function(z) attr(z, "__settings__")$elapsed),
    size=sapply(res, object.size))

load("d:/bam/BAM_data_v2019/gnm/out/test/cawa.RData")
cawa <- data.frame(time=sapply(res, function(z) attr(z, "__settings__")$elapsed),
    size=sapply(res, object.size))

## checking completeness and running missing pieces
library(mefa4)
library(gbm)
## load objects
fn <- "BAMdb-GNMsubset-2020-01-08.RData"
PROJ <- "boot"
load(file.path("d:/bam/BAM_data_v2019/gnm", "data", fn))
load("d:/bam/BAM_data_v2019/gnm/out/xvinfo.RData")
## explore XV results
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

SPP <- colnames(yy)
DIRS <- list.files(file.path("d:/bam/BAM_data_v2019/gnm", "out", PROJ))
BCR <- paste0("BCR_", u)
B <- 32

CHUNK <- list()
for (spp in DIRS) {
    d1 <- list.files(file.path("d:/bam/BAM_data_v2019/gnm", "out", PROJ, spp))
    for (bcr in BCR) {
        i <- paste0(spp, "-", bcr)
        if (dir.exists(file.path("d:/bam/BAM_data_v2019/gnm", "out", PROJ, spp, bcr))) {
            d2 <- list.files(file.path("d:/bam/BAM_data_v2019/gnm", "out", PROJ, spp, bcr))
            if (length(d2)) {
                b <- sort(as.integer(sapply(strsplit(sapply(strsplit(d2, "\\."), "[[", 1), "-"), "[[", 4)))
                MISSING <- setdiff(1:B, b)
            } else {
                MISSING <- 1:B
            }
            if (length(MISSING))
                CHUNK[[i]] <- MISSING
        } else {
            CHUNK[[i]] <- 1:B
        }
    }
}

.run_brt_boot_test <- function(b, RUN, verbose=interactive()) {
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
        out <- TRUE
    }
    attr(out, "__settings__") <- list(
        species=spp, region=BCR, iteration=b, elapsed=proc.time()["elapsed"]-t0)
    out
}

## now only saving 0 detection output
for (i in names(CHUNK)) {
    cat("\n", i)
    flush.console()
    tmp <- strsplit(i, "-")[[1]]
    spp <- tmp[1]
    bcr <- tmp[2]
    dr <- file.path("d:/bam/BAM_data_v2019/gnm", "out", PROJ, spp, bcr)
    if (!dir.exists(dr))
        dir.create(dr)
    for (b in CHUNK[[i]]) {
        cat(".")
        flush.console()
        out <- .run_brt_boot_test(b, i)
        if (inherits(out, "try-error"))
            save(out, file=file.path(dr, paste0("gnmboot-", spp, "-", bcr, "-", b, ".RData")))
    }
}

## scan again
CHUNK <- list()
for (spp in DIRS) {
    d1 <- list.files(file.path("d:/bam/BAM_data_v2019/gnm", "out", PROJ, spp))
    for (bcr in BCR) {
        i <- paste0(spp, "-", bcr)
        if (dir.exists(file.path("d:/bam/BAM_data_v2019/gnm", "out", PROJ, spp, bcr))) {
            d2 <- list.files(file.path("d:/bam/BAM_data_v2019/gnm", "out", PROJ, spp, bcr))
            if (length(d2)) {
                b <- sort(as.integer(sapply(strsplit(sapply(strsplit(d2, "\\."), "[[", 1), "-"), "[[", 4)))
                MISSING <- setdiff(1:B, b)
            } else {
                MISSING <- 1:B
            }
            if (length(MISSING))
                CHUNK[[i]] <- MISSING
        } else {
            CHUNK[[i]] <- 1:B
        }
    }
}

## these are the missing ones (Timeout on Graham)
for (i in names(CHUNK)) {
    tmp <- strsplit(i, "-")[[1]]
    spp <- tmp[1]
    bcr <- tmp[2]
    dr <- file.path("d:/bam/BAM_data_v2019/gnm", "out", PROJ, spp, bcr)
    if (!dir.exists(dr))
        dir.create(dr)
    for (b in CHUNK[[i]]) {
        cat(i, b, "\n")
        flush.console()
        out <- .run_brt_boot(b, i)
        save(out, file=file.path(dr, paste0("gnmboot-", spp, "-", bcr, "-", b, ".RData")))
    }
}



