## take BRT outputs and stacks to make a density raster
## and then stitch them together

## libraries
library(mefa4)
library(gbm)
library(raster)
library(dismo)

## change this as needed
ROOT <- "d:/bam/BAM_data_v2019/gnm"

## check if we have all the pieces
DONE <- sapply(strsplit(list.files(paste0(ROOT, "/out/gnm")), ".", fixed=TRUE), function(z) z[1L])
table(sapply(strsplit(gsub(".RData", "", DONE), "-"), "[[", 2))

## load training data
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2019-03-01.RData"))

## prepare stacks
STACK <- list()
for (i in 4:14) {
    cat("\nLoading stack for BCR", i, "\n")

    pr <- stack(file.path(ROOT, "data", "stacks", paste0("bcr", i, "_1km.grd")))
    prc <- stack(file.path(ROOT, "data", "stacks", paste0("bcr", i, "cat_1km.grd")))

    names(prc) <- c("bcr", "nalc", "lf")
    prc[["nalc"]] <- as.factor(prc[["nalc"]])
    prc[["lf"]] <- as.factor(prc[["lf"]])
    pr <- pr[[which(names(pr) != "bcr")]]
    prc <- prc[[which(names(prc) != "bcr")]]
    pr <- addLayer(pr, prc)

    print(compare_sets(c(cnf, CN[[paste0("BCR_", i)]]), names(pr)))

    STACK[[paste0("BCR_", i)]] <- pr
}
rm(prc, pr)

## extracting model stats

SPP <- sort(colnames(yy))
OK <- matrix(1, length(SPP), 11)
dimnames(OK) <- list(SPP, 4:14)
NULLS <- character(0)

for (spp in SPP) {

    fl <- paste0(spp, "-BCR_", 4:14, ".RData")
    names(fl) <- 4:14
    VARIMP <- CVSTATS <- list()

    if (!dir.exists(file.path(ROOT, "artifacts", spp)))
        dir.create(file.path(ROOT, "artifacts", spp))

    ## extract stats
    cat(spp, ": ", sep="")
    pdf(file.path(ROOT, "artifacts", spp, paste0("brtplot-", spp, ".pdf")),
        height=9, width=12, onefile=TRUE)
    for (i in 4:14) {
        cat(i, ", ", sep="")
        e <- new.env()
        tmp <- try(load(file.path(ROOT, "out", "gnm", fl[as.character(i)]), envir=e))
        if (inherits(tmp, "try-error")) {
            OK[spp, as.character(i)] <- -2 # 0-byte input file
        } else {
            brt <- e$out
            rm(e)
            if (is.null(brt)) {
                NULLS <- c(NULLS, unname(fl[as.character(i)]))
                OK[spp, as.character(i)] <- -1 # NULL object
            } else {
                if (inherits(brt, "try-error")) {
                    VARIMP[[paste0("BCR_", i)]] <- NULL
                    CVSTATS[[paste0("BCR_", i)]] <- NULL
                    plot.new()
                    text(0.5,0.5,paste("BRT for", spp, "in BCR", i, "failed."))
                    OK[spp, as.character(i)] <- 0 # try-error model fit
                } else {
                    varimp <- as.data.frame(brt$contributions)
                    VARIMP[[paste0("BCR_", i)]] <- data.frame(BCR=i, varimp)
                    cvstats <- t(as.data.frame(brt$cv.statistics))
                    CVSTATS[[paste0("BCR_", i)]] <- data.frame(BCR=i, cvstats)

                    vi <- varimp[,2]
                    names(vi) <- varimp[,1]
                    vi <- vi[vi > 0]
                    op <- par(mfrow=c(1,1), mar=c(5,15,2,2), las=1)
                    barplot(rev(vi), horiz=TRUE, cex.names=0.6, xlab="Relative influence",
                        main=paste(spp, " in BCR", i))
                    par(op)
                    gbm.plot(brt, n.plots=12, smooth=TRUE, common.scale=FALSE)
                }
            }
        }
    }
    cat("\n")
    dev.off()

    VI <- do.call(rbind, VARIMP)
    CS <- do.call(rbind, CVSTATS)

    write.csv(VI, row.names = FALSE,
        file=file.path(ROOT, "artifacts", spp, paste0("varimp-", spp, ".csv")))
    write.csv(CS, row.names = FALSE,
        file=file.path(ROOT, "artifacts", spp, paste0("cvstats-", spp, ".csv")))
}

save(OK, file=file.path(ROOT, "OK.RData"))

failed <- sapply(which(OK == 0), function(i)
    paste0(rownames(OK)[row(OK)[i]], "-BCR_", colnames(OK)[col(OK)[i]]))
msg <- list()
for (h in failed) {
    e <- new.env()
    tmp <- try(load(file.path(ROOT, "out", "gnm", paste0(h, ".RData")), envir=e))
    msg[[h]] <- as.character(e$out)
}
msg <- unlist(msg)
msg <- msg[!startsWith(msg, "0 detections")]

## picking up the leftover
leftover <- sapply(which(OK < 0), function(i)
    paste0(rownames(OK)[row(OK)[i]], "-BCR_", colnames(OK)[col(OK)[i]]))

leftover <- sort(c(leftover, names(msg))) # 158

run_brt2 <- function(RUN, SUB=NULL, RATE=0.001) {
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
        dd2[ss, c(cnf, CN[[BCR]])])
    if (!is.null(SUB)) {
        DAT$weights <- 1
        SUB <- min(SUB, nrow(DAT))
        ## this is not bootstrap resampling
        ## just subset to reduce memory footprint
        DAT <- DAT[sample(nrow(DAT), SUB, prob=DAT$weights),]
        cat("Sample size =", nrow(DAT), "\n")
    }
    out <- try(gbm::gbm(count ~ . + offset(offset),
        data=DAT,
        n.trees = 1000,
        interaction.depth = 3,
        shrinkage = RATE,
        bag.fraction = 0.5,
        weights = DAT$weights,
        distribution = "poisson",
        var.monotone = NULL,#rep(0, length(4:ncol(DAT))),
        verbose = TRUE)
    )
    out
}

for (h in leftover) {
    cat("\n", h, which(leftover == h), "/", length(leftover), "\n")
    out <- run_brt2(h, SUB=15000, RATE=0.001)
    save(out, file=file.path(ROOT, "out", "leftover", paste0(h, ".RData")))
}

## collecting ntree info

SPP <- sort(colnames(yy))
NTREE <- matrix(0, length(SPP), 11)
dimnames(NTREE) <- list(SPP, 4:14)

for (spp in SPP) {

    fl <- paste0(spp, "-BCR_", 4:14, ".RData")
    names(fl) <- 4:14

    cat("\n", spp, ": ", sep="")
    for (i in 4:14) {
        cat(i, ", ", sep="")
        flush.console()
        e <- new.env()
        tmp <- try(load(file.path(ROOT, "out", "gnm", fl[as.character(i)]), envir=e))
        if (!inherits(tmp, "try-error")) {
            brt <- e$out
            if (!is.null(brt)) {
                if (!inherits(brt, "try-error")) {
                    NTREE[spp, as.character(i)] <- brt$n.trees
                }
            }
        }
    }
}
colnames(NTREE) <- paste0("BCR_",4:14)
save(NTREE, file=file.path(ROOT, "NTREE.RData"))


## now predicting

predict_gbm <- function(brt, ppp, impute=0) {
    if (inherits(brt, "try-error")) {
        rp <- ppp[[1]]
        values(rp)[!is.na(values(rp))] <- impute
    } else {
        ctb <- brt$contributions
        pset <- names(ppp)
        n <- length(values(ppp[[1]]))
        nd <- matrix(0, n, length(pset))
        colnames(nd) <- pset
        for (i in pset) {
            nd[,i] <- values(ppp[[i]])
        }

        notNA <- rowSums(is.na(nd)) == 0
        nd <- as.data.frame(nd[notNA,,drop=FALSE])
        nd$nalc <- as.factor(nd$nalc)
        nd$lf <- as.factor(nd$lf)
        nd$offset <- 0
        nd$weights <- 1
        z0 <- suppressWarnings(predict.gbm(brt, nd, type="response", n.trees=brt$n.trees))
        z <- rep(NA, n)
        z[notNA] <- z0
        zz <- matrix(t(z), dim(ppp)[1], dim(ppp)[2], byrow=TRUE)
        rp <- raster(x=zz, template=ppp)
    }
    rp
}

load(file.path(ROOT, "OK.RData"))

SPP <- rownames(OK[rowSums(OK)==ncol(OK)])


spp <- "CAWA"
PIECES <- list()
t0 <- proc.time()
for (i in 4:14) {
    gc()
    BCR <- paste0("BCR_", i)
    cat("\n", spp, BCR)
    flush.console()
    ## load spp/bcr model object
    e <- new.env()
    tmp <- try(load(file.path(ROOT, "out", "gnm", paste0(spp, "-", BCR, ".RData")), envir=e))
    #tmp <- try(load(file.path(ROOT, "out", "leftover", "ATSP-BCR_11.RData"), envir=e))
    brt <- e$out
    rm(e)

    cat(" -", if (inherits(brt, "try-error")) NA else brt$n.trees)

    ## predict
    PIECES[[BCR]] <- predict_gbm(brt, STACK[[BCR]], 0)
    cat(" @", round((proc.time() - t0)[3]/60, 2), "min")
}


rast <- mosaic(
    PIECES[['BCR_4']],
    PIECES[['BCR_5']],
    PIECES[['BCR_6']],
    PIECES[['BCR_7']],
    PIECES[['BCR_8']],
    PIECES[['BCR_9']],
    PIECES[['BCR_10']],
    PIECES[['BCR_11']],
    PIECES[['BCR_12']],
    PIECES[['BCR_13']],
    PIECES[['BCR_14']],
    fun=mean)

writeRaster(rast, file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, ".tif")))

MAX <- 3 * cellStats(rast, 'mean')
bluegreen.colors <- colorRampPalette(c("#FFFACD", "lemonchiffon","#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.5)
png(file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, ".png")),
    height=2000, width=3000)
op <- par(cex.main=3, mfcol=c(1,1), oma=c(0,0,0,0), mar=c(0,0,5,0))
plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(spp, "- BCR", i), box=FALSE)
plot(rast, col=bluegreen.colors(15), zlim=c(0,MAX), axes=FALSE,
    main=spp, add=TRUE, legend.width=1.5, horizontal = TRUE,
    smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=2))
par(op)
dev.off()




t0 <- proc.time()
rast <- raster::predict(pr, brt, type="response", n.trees=brt$n.trees)
proc.time() - t0


pp <- trim(crop(pr, extent(c(-1798000/10, 13000, 6580000, 8336000))))
ppp <- pp[[rownames(brt$contributions)]]
rast <- raster::predict(pp, brt, type="response", n.trees=10, factors=c("nalc","lf"))


t0 <- proc.time()
rast <- predict_gbm(brt, pr)
proc.time() - t0

q99 <- quantile(rast, probs=c(0.99))
prev <- cellStats(rast, 'mean')
max <- 3*prev
bluegreen.colors <- colorRampPalette(c("#FFFACD", "lemonchiffon","#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.5)

op <- par(cex.main=1.8, mfcol=c(1,1), oma=c(0,0,0,0), mar=c(0,0,5,0))
plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(spp, "- BCR", i))
plot(rast, col=bluegreen.colors(15), zlim=c(0,max), axes=FALSE,
    main=spp, add=TRUE, legend.width=1.5, horizontal = TRUE,
    smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=1.5))
par(op)



value <- as.numeric(value)
r <- as.matrix(Xtab(value ~ Row + Col, rc))
r[is.na(as.matrix(rt))] <- NA
raster(x=r, template=rt)



writeRaster(rast, filename=paste(w,speclist[j],"_pred1km3",sep=""), format="GTiff",overwrite=TRUE)

## checking prediction layers

e <- new.env()
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2019-03-01.RData"), envir=e)
dd2 <- e$dd2
i <- 6
tmp <- list()
for (i in 4:14) {
    cat(i, "\n")
    CN <- e$CN[[paste0("BCR_", i)]]
    pr <- stack(file.path(ROOT, "data", "stacks", paste0("bcr", i, "_2011_1km.grd")))
    prc <- stack(file.path(ROOT, "data", "stacks", paste0("bcr", i, "cat_1km.grd")))
    names(prc) <- c("bcr", "nalc", "lf")
    prc[["nalc"]] <- as.factor(prc[["nalc"]])
    prc[["lf"]] <- as.factor(prc[["lf"]])
    pr <- pr[[which(names(pr) != "bcr")]]
    prc <- prc[[which(names(prc) != "bcr")]]
    pr <- addLayer(pr, prc)

    tmp[[i]] <- list(comp=compare_sets(CN, names(pr)),
        diff1=setdiff(CN, names(pr)),
        diff2=setdiff(names(pr), CN),
        NM=sapply(1:107,function(i) sum(!is.na(values(pr[[i]])))))
}
