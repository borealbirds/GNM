library(mefa4)
library(raster)
library(gbm)
library(dismo)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
PROJ <- "run1"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2019-06-05.RData"))

#SPP <- colnames(yy)
#SPP <- c("ALFL", "AMRO", "BOCH", "BTNW", "CAWA",
#    "OSFL", "OVEN", "RUBL", "WCSP", "YRWA")
SPP <- "OSFL"
SPPBCRss <- NULL
for (spp in SPP) {
    SPPBCRss <- c(SPPBCRss, SPPBCR[grep(spp, SPPBCR)])
}

DONE <- sapply(strsplit(list.files(file.path(ROOT, "out", PROJ)), ".", fixed=TRUE), function(z) z[1L])
if (length(DONE) > 0) {
    cat("OK\n* Summary so far:\n")
    data.frame(table(sapply(strsplit(gsub(".RData", "", DONE), "-"), "[[", 2)))
}

for (i in u) {
    cat("\nLoading stack for BCR", i, "\n")
    flush.console()
    gc()

    pr <- stack(file.path(ROOT, "data", "subunits", paste0("bcr", i, "all_1km.grd")))
    rr <- raster(file.path(ROOT, "data", "subunits", "road", paste0("bcr", i, "road_1km.grd")))
    #pr[["nalc"]] <- as.factor(pr[["nalc"]])
    #pr[["lf"]] <- as.factor(pr[["lf"]])
    try(pr <- addLayer(pr, rr))
}

## prepare stacks
for (i in u) {
    cat("\nLoading stack for BCR", i, "\n")
    flush.console()
    gc()

    pr <- stack(file.path(ROOT, "data", "subunits", paste0("bcr", i, "all_1km.grd")))
    rr <- raster(file.path(ROOT, "data", "subunits", "road", paste0("bcr", i, "road_1km.grd")))
    #pr[["nalc"]] <- as.factor(pr[["nalc"]])
    #pr[["lf"]] <- as.factor(pr[["lf"]])
    if (!compareRaster(pr, rr, stopiffalse=FALSE)) {
        rr <- resample(rr, pr, "ngb")
    }
    pr <- addLayer(pr, rr)
    names(pr)[length(names(pr))] <- "ROAD"

    print(compare_sets(CN[[paste0("BCR_", i)]], names(pr))) # diff is ROAD

    #pr <- trim(pr)

    #STACK[[paste0("BCR_", i)]] <- pr

    #pset <- names(pr)
    pset <- CN[[paste0("BCR_", i)]]

#    pset <- pset[pset != "ROAD"]

    n <- length(values(pr[[1]]))
    nd <- matrix(0, n, length(pset))
    colnames(nd) <- pset
    for (j in pset) {
        nd[,j] <- values(pr[[j]])
    }
    notNA <- rowSums(is.na(nd)) == 0
    #notNA <- !is.na(nd[,1])
    nd <- as.data.frame(nd[notNA,,drop=FALSE])
    nd$nalc <- as.factor(nd$nalc)
    nd$lf <- as.factor(nd$lf)

#    nd$ROAD <- 0

    nd$offset <- 0L
    nd$weights <- 1L
    ND <- list(data=nd, n=n, subset=which(notNA), dim=dim(pr))
    saveRDS(ND, file=file.path(ROOT, paste0("STACK-ND-BCR_", i, ".rds")))

}

## need to set NALC land cover types not used in model (snow/ice) to 0

predict_gbm_data <- function(ppp) {
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

    nd$ROAD <- 0

    nd$offset <- 0
    nd$weights <- 1
    list(data=nd, n=n, subset=which(notNA), dim=dim(ppp))
}

t0 <- proc.time()
for (i in 4:14) {
    gc()
    BCR <- paste0("BCR_", i)
    cat("\n", BCR)
    flush.console()
    ND <- predict_gbm_data(STACK[[BCR]])
    save(ND, file=file.path(ROOT, paste0("STACK-ND-BCR_", i, ".RData")))
    cat(" @", round((proc.time() - t0)[3]/60, 2), "min")
}

ppp <- predict_gbm_data(STACK[[BCR]])
pr <- predict_gbm(brt, ppp, STACK[[BCR]][[1]], 0)

## predict ---------------------------------

## Done: 4, 5, 6, 7, 8,  9,  10,  11
i <- 12 # this is BCR

library(mefa4)
library(gbm)
library(raster)
ROOT <- "d:/bam/BAM_data_v2019/gnm"
#ROOT <- "c:/p/tmp/gnm"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2019-03-01.RData"))

predict_gbm <- function(brt, ppp, r, impute=0) {
    if (inherits(brt, "try-error")) {
        rp <- r[[1]]
        values(rp)[!is.na(values(rp))] <- impute
    } else {
        if (!("ROAD" %in% colnames(ppp$data)))
            ppp$data$ROAD <- 0
        z0 <- suppressWarnings(predict.gbm(brt, ppp$data, type="response", n.trees=brt$n.trees))
        z <- rep(NA, ppp$n)
        z[ppp$subset] <- z0
        zz <- matrix(t(z), ppp$dim[1], ppp$dim[2], byrow=TRUE)
        rp <- raster(x=zz, template=r)
    }
    rp
}

#load(file.path(ROOT, "OK.RData"))
#SPP <- rownames(OK)[rowSums(OK >= 0)==ncol(OK)]

SPP <- colnames(yy)
#SPP <- rev(colnames(yy))

#PROJ <- "gnm"
PROJ <- "roadfix"

#spp <- "CAWA"
#spp <- "AMRO"
spp <- "OSFL"

r1 <- raster(file.path(ROOT, "data", "stacks", paste0("bcr", i, "_1km.grd")))
load(file.path(ROOT, paste0("STACK-ND-BCR_", i, ".RData")))
table(ND$data$nalc, useNA="a")
for (spp in SPP) {
    gc()
    fout0 <- file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-", PROJ, ".tif"))
    fout <- file.path(ROOT, "artifacts", spp,
        paste0("mosaic-", spp, "-BCR_", i, "-", PROJ, "-nalcfix.tif"))
    if (!file.exists(fout0) && !file.exists(fout)) {
        BCR <- paste0("BCR_", i)
        cat("\n", spp, BCR)
        flush.console()
        ## load spp/bcr model object
        e <- new.env()
        tmp <- try(load(file.path(ROOT, "out", PROJ, paste0(spp, "-", BCR, ".RData")), envir=e))
        brt <- e$out
        rm(e)
        ## preprocesses stack
        cat(" -", if (inherits(brt, "try-error")) NA else brt$n.trees)
        ## predict
        rrr <- predict_gbm(brt, ND, r1, 0)
        writeRaster(rrr, fout, overwrite=TRUE)
    }
}

spp <- "AMRO"
for (spp in SPP) {
    fout <- file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-", PROJ, "-nalcfix.tif"))
    fin <- file.path(ROOT, "artifacts", spp,
        paste0("mosaic-", spp, "-BCR_", 4:14, "-", PROJ, "-nalcfix.tif"))
    if (!file.exists(fout)) {
		cat(spp, "\n")
		flush.console()
		r4 <- raster(fin[1])
		r5 <- raster(fin[2])
		r6 <- raster(fin[3])
		r7 <- raster(fin[4])
		r8 <- raster(fin[5])
		r9 <- raster(fin[6])
		r10 <- raster(fin[7])
		r11 <- raster(fin[8])
		r12 <- raster(fin[9])
		r13 <- raster(fin[10])
		r14 <- raster(fin[11])
		rast <- mosaic(r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, fun=mean)
		writeRaster(rast, fout, overwrite=TRUE)
    }
}

## making png maps
library(rgdal)
library(raster)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
#ROOT <- "c:/p/tmp/gnm"

bluegreen.colors <- colorRampPalette(c("#FFFACD", "lemonchiffon","#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.5)

BCR <- readOGR(dsn=file.path(ROOT, "data", "bcr"), "bcrfinallcc")
PROV <- readOGR(dsn=file.path(ROOT, "data", "prov"), "province_state_line")
LAKES <- readOGR(dsn=file.path(ROOT, "data", "lakes"), "lakes_lcc")
LAKES <- spTransform(LAKES, proj4string(BCR))

#spp <- "AMGO"
#library(opticut)
for (spp in SPP) {
    fout0 <- file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-", PROJ, "-nalcfix.tif"))
    if (file.exists(fout0)) {
        cat(spp, "\n")
        flush.console()

        rast <- raster(file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-", PROJ, ".tif")))
#        rast <- raster(file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-", PROJ, "-nalcfix.tif")))

        ## play with Lc
        #lc <- lorenz(values(rast)[!is.na(values(rast))])
        #q <- quantile(lc, probs=c(0.05, 0.1, 0.2, 0.5, 0.8, 0.99), type="L")
        #MAX <- q["80%"]
        MAX <- cellStats(rast, 'mean')
        png(file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-", PROJ, ".png")),
            height=2000, width=3000)
#        png(file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-", PROJ, "-nalcfix.png")),
#            height=2000, width=3000)
        op <- par(cex.main=3, mfcol=c(1,1), oma=c(0,0,0,0), mar=c(0,0,5,0))
        plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(spp), box=FALSE)
        plot(rast, col=bluegreen.colors(15), zlim=c(0,MAX), axes=FALSE,
            main=spp, add=TRUE, legend.width=1.5, horizontal = TRUE,
            smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=2))
        plot(PROV, col="grey", add=TRUE)
        plot(LAKES,col="#aaaaff", border=NA, add=TRUE)
        plot(BCR, add=TRUE)
        par(op)
        dev.off()

        col1 <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8","#91BFDB","#4575B4")))(100)
        op <- par(cex.main=3, mfcol=c(1,1), oma=c(0,0,0,0), mar=c(0,0,5,0))
        plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(spp), box=FALSE)
        plot(rast, col=col1, zlim=c(0,MAX), axes=FALSE,
            main=spp, add=TRUE, legend.width=1.5, horizontal = TRUE,
            smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=2))
        plot(PROV, col="grey", add=TRUE)
        plot(LAKES,col="#aaaaff", border=NA, add=TRUE)
        plot(BCR, add=TRUE)
        par(op)


    }
}

## copying mosaice`d files
for (spp in SPP) {
    fi <- file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-", PROJ, ".tif"))
    fo <- file.path("d:/bam/BAM_data_v2019/gnm", "roadfix-mosaic", paste0("mosaic-", spp, "-", PROJ, ".tif"))
    file.copy(fi, fo)
}

## bcr level maps
i <- 6
BCRi <- BCR[BCR@data$BCR == i,]

for (spp in SPP) {

    rast <- try(raster(file.path(ROOT, "artifacts", spp,
        paste0("mosaic-", spp, "-BCR_", i, "-", PROJ, ".tif"))))
    if (!inherits(rast, "try-error")) {

        MAX <- 3 * cellStats(rast, 'mean')
        png(file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-BCR_", i, "-", PROJ, ".png")),
            height=2000, width=3000)
        op <- par(cex.main=3, mfcol=c(1,1), oma=c(0,0,0,0), mar=c(0,0,5,0))
        plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(spp), box=FALSE)
        plot(rast, col=bluegreen.colors(15), zlim=c(0,MAX), axes=FALSE,
            main=spp, add=TRUE, legend.width=1.5, horizontal = TRUE,
            smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=2))
        plot(PROV, col="grey", add=TRUE)
        plot(LAKES,col="#aaaaff", border=NA, add=TRUE)
        plot(BCRi, add=TRUE)
        par(op)
        dev.off()
    } else {
        cat(spp, "\n")
    }
}


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
