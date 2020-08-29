library(mefa4)
library(gbm)
library(raster)
library(rgdal)

## adjust these paths as needed
ROOT <- "d:/bam/BAM_data_v2019/gnm"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2020-01-08.RData"))
SUB <- readOGR(dsn=file.path(ROOT, "data", "sub-bound"), "BCRSubunits")
BCR <- readOGR(dsn=file.path(ROOT, "data", "bcr"), "bcrfinallcc")
PROV <- readOGR(dsn=file.path(ROOT, "data", "prov"), "province_state_line")
LAKES <- readOGR(dsn=file.path(ROOT, "data", "lakes"), "lakes_lcc")
LAKES <- spTransform(LAKES, proj4string(BCR))

PROJ <- "boot"
B <- 32 # max is 32
BCRs <- u[u < 200] # Canada only
SPP <- colnames(yy)

bluegreen.colors <- function(n) {
    x <- hcl.colors(n, "Lajolla")
    c(x, rep(x[length(x)], ceiling(n*0.1)))
}

## -------- calculating Canada-wide mosaiced layers for each boot run -----------

for (spp in SPP) {
    for (i in 1:B) {
        BCRr <- sample(BCRs)
        for (bcr in BCRr) {
            cat(spp, "run", i, "in BCR", bcr, "\n")
            flush.console()
            bi <- if (B > 32)
                sample(B, 1) else i
            p <- SUB[SUB@data[,1] == bcr,]
            ## linear (uniform): plot(ecdf(runif(10^3)))
            #w <- runif(1, 0, 1)
            ## sigmoid cdf: plot(ecdf(rbeta(10^3, 2, 2)))
            w <- rbeta(1, 2, 2)
            pb <- rgeos::gBuffer(p, width=w * 100*10^3)
            ri <- raster(file.path(ROOT, "out", "parts", spp,
                paste0("pred-", spp, "-BCR_", bcr, "-boot-", bi, ".tif")))
            pb <- rgeos::gBuffer(p, width=w * 100 * 10^3)
            if (bcr == BCRr[1]) {
                rout <- mask(ri, pb)
                #pb0 <- rgeos::gBuffer(p, width=100*10^3)
                #plot(ri2)
                #plot(p, add=TRUE)
                #plot(pb0, add=TRUE, border=2)
            } else {
                rout <- mosaic(rout, mask(ri, pb), fun=mean)
            }
        }
        if (!dir.exists(file.path(ROOT, "artifacts", spp)))
            dir.create(file.path(ROOT, "artifacts", spp))
        fout <- file.path(ROOT, "artifacts", spp, paste0("pred-", spp, "-CAN-boot-", bi, ".tif"))
        writeRaster(rout, overwrite=TRUE, filename=fout)
    }
}


## -------- making Canada-wide mean & SD layers -----------

for (spp in SPP) {

    cat(spp, "\n")
    fi <- file.path(ROOT, "artifacts", spp, paste0("pred-", spp, "-CAN-boot-", 1:B, ".tif"))
    st <- list()
    for (b in 1:B) {
        st[[b]] <- raster(fi[b])
    }

    # mean
    rast <- mosaic(st[[1]], st[[2]], st[[3]], st[[4]], st[[5]], st[[6]], st[[7]], st[[8]],
        st[[9]], st[[10]], st[[11]], st[[12]], st[[13]], st[[14]], st[[15]], st[[16]],
        st[[17]], st[[18]], st[[19]], st[[20]], st[[21]], st[[22]], st[[23]], st[[24]],
        st[[25]], st[[26]], st[[27]], st[[28]], st[[29]], st[[30]], st[[31]], st[[32]],
        fun=mean)
    fo <- file.path(ROOT, "artifacts", spp, paste0("pred-", spp, "-CAN-Mean.tif"))
    writeRaster(rast, fo, overwrite=TRUE)

    # SD
    rast <- mosaic(st[[1]], st[[2]], st[[3]], st[[4]], st[[5]], st[[6]], st[[7]], st[[8]],
        st[[9]], st[[10]], st[[11]], st[[12]], st[[13]], st[[14]], st[[15]], st[[16]],
        st[[17]], st[[18]], st[[19]], st[[20]], st[[21]], st[[22]], st[[23]], st[[24]],
        st[[25]], st[[26]], st[[27]], st[[28]], st[[29]], st[[30]], st[[31]], st[[32]],
        fun=sd)
    fo <- file.path(ROOT, "artifacts", spp, paste0("pred-", spp, "-CAN-SD.tif"))
    writeRaster(rast, fo, overwrite=TRUE)
}

## -------- png's of Canada-wide mean & SD layers -----------

for (spp in SPP) {
    cat(spp, "\n")
    flush.console()

    rast <- raster(file.path(ROOT, "artifacts", spp, paste0("pred-", spp, "-CAN-Mean.tif")))
    q <- quantile(values(rast), 0.999, na.rm=TRUE)
    values(rast)[!is.na(values(rast)) & values(rast)>q] <- q
    jpeg(file.path(ROOT, "artifacts", spp, paste0("pred-", spp, "-CAN-Mean.jpg")),
        height=4000, width=6000, res=300)
    op <- par(cex.main=3, mfcol=c(1,1), oma=c(0,0,0,0), mar=c(0,0,5,0))
    plot(rast, col=bluegreen.colors(15), axes=FALSE, legend=TRUE, main=paste(spp, "Mean"), box=FALSE)
    plot(PROV, col="#aaaaaa", add=TRUE)
    plot(LAKES,col="#aaaaff", border=NA, add=TRUE)
    plot(BCR, add=TRUE)
    par(op)
    dev.off()

    # SD
    rast <- raster(file.path(ROOT, "artifacts", spp, paste0("pred-", spp, "-CAN-SD.tif")))
    q <- quantile(values(rast), 0.999, na.rm=TRUE)
    values(rast)[!is.na(values(rast)) & values(rast)>q] <- q
    jpeg(file.path(ROOT, "artifacts", spp, paste0("pred-", spp, "-CAN-SD.jpg")),
        height=4000, width=6000, res=300)
    op <- par(cex.main=3, mfcol=c(1,1), oma=c(0,0,0,0), mar=c(0,0,5,0))
    plot(rast, col=bluegreen.colors(15), axes=FALSE, legend=TRUE, main=paste(spp, "SD"), box=FALSE)
    plot(PROV, col="#aaaaaa", add=TRUE)
    plot(LAKES,col="#aaaaff", border=NA, add=TRUE)
    plot(BCR, add=TRUE)
    par(op)
    dev.off()

}

