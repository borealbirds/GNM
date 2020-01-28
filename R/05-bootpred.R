library(mefa4)
library(gbm)
library(raster)
ROOT <- "d:/bam/BAM_data_v2019/gnm"

fn <- "BAMdb-GNMsubset-2020-01-08.RData"
PROJ <- "boot"
load(file.path(ROOT, "data", fn))


## prepare stacks for US
if (FALSE) {
for (i in u[u>=200]) {
    cat("\nLoading stack for BCR", i, "\n")
    flush.console()
    gc()

    pr <- stack(file.path(ROOT, "data", "subunits", paste0("bcr", i, "all_1km.grd")))
    plot(pr[["ROAD"]], main=paste("ROAD", i))

    pset <- CN[[paste0("BCR_", i)]]
    print(compare_sets(pset, names(pr)))

    n <- length(values(pr[[1]]))
    nd <- matrix(0, n, length(pset))
    colnames(nd) <- pset
    for (j in pset) {
        #cat(j, "\n");flush.console()
        nd[,j] <- values(pr[[j]])
    }
    notNA <- rowSums(is.na(nd[,colnames(nd) != "ROAD"])) == 0
    #notNA <- !is.na(nd[,1])
    nd <- as.data.frame(nd[notNA,,drop=FALSE])
    nd$nalc <- as.factor(nd$nalc)
#    nd$lf <- as.factor(nd$lf)
    nd$ROAD[is.na(nd$ROAD)] <- 0

    nd$offset <- 0L
    nd$weights <- 1L
    ND <- list(data=nd, n=n, subset=which(notNA), dim=dim(pr))
    saveRDS(ND, file=file.path(ROOT, paste0("STACK-ND-BCR_", i, ".rds")))

}
}

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

SPP <- colnames(yy)

B <- 32
spp <- "CAWA"

for (BCR in u) {
    cat("loading stack: BCR", BCR, "\n")
    flush.console()
    #spp <- "OSFL"
    r1 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr", BCR, "all_1km.grd")))
    ND <- readRDS(file.path(ROOT, paste0("STACK-ND-BCR_", BCR, ".rds")))

    ND$data$YEAR <- 2011

    ## predict for runs
    for (b in 1:B) {
        cat("\t", spp, "@ BCR", BCR, ">", b, "/", max(B))
        flush.console()
        fin <- file.path(ROOT, "out", PROJ, spp, paste0("BCR_", BCR), paste0("gnmboot-", spp, "-BCR_", BCR, "-", b, ".RData"))
        if (file.exists(fin)) {
            e <- new.env()
            aa <- try(load(fin, envir=e))
            if (inherits(aa, "try-error")) {
                writeLines(as.character(aa),
                    file.path(ROOT, "out", "error",
                        paste0(PROJ, "-", spp, "-BCR_", BCR, ".txt")))
            } else {
                brt <- e$out
            }
            rm(e)
        } else {
            brt <- NULL
        }
        rrr <- predict_gbm(brt, ND, r1, 0)
        if (!dir.exists(file.path(ROOT, "artifacts2", spp)))
            dir.create(file.path(ROOT, "artifacts2", spp))
        if (!dir.exists(file.path(ROOT, "artifacts2", spp, "parts")))
            dir.create(file.path(ROOT, "artifacts2", spp, "parts"))
        fout <- file.path(ROOT, "artifacts2", spp, "parts",
            paste0("mosaic-", spp, "-BCR_", BCR, "-", PROJ, "-", b, ".tif"))
        writeRaster(rrr, fout, overwrite=TRUE)
        cat(" > OK\n")
    }
    ## combine runs: mean & SD
    st <- list()
    for (b in 1:B) {
        fin <- file.path(ROOT, "artifacts2", spp, "parts",
            paste0("mosaic-", spp, "-BCR_", BCR, "-", PROJ, "-", b, ".tif"))
        st[[b]] <- raster(fin)
    }

    Mean <- mosaic(st[[1]], st[[2]], st[[3]], st[[4]], st[[5]], st[[6]], st[[7]], st[[8]],
        st[[9]], st[[10]], st[[11]], st[[12]], st[[13]], st[[14]], st[[15]], st[[16]],
        st[[17]], st[[18]], st[[19]], st[[20]], st[[21]], st[[22]], st[[23]], st[[24]],
        st[[25]], st[[26]], st[[27]], st[[28]], st[[29]], st[[30]], st[[31]], st[[32]],
        fun=mean)
    fout <- file.path(ROOT, "artifacts2", spp,
        paste0(spp, "-BCR_", BCR, "-", PROJ, "-Mean.tif"))
    writeRaster(Mean, fout, overwrite=TRUE)
    Sd <- mosaic(st[[1]], st[[2]], st[[3]], st[[4]], st[[5]], st[[6]], st[[7]], st[[8]],
        st[[9]], st[[10]], st[[11]], st[[12]], st[[13]], st[[14]], st[[15]], st[[16]],
        st[[17]], st[[18]], st[[19]], st[[20]], st[[21]], st[[22]], st[[23]], st[[24]],
        st[[25]], st[[26]], st[[27]], st[[28]], st[[29]], st[[30]], st[[31]], st[[32]],
        fun=sd)
    fout <- file.path(ROOT, "artifacts2", spp,
        paste0(spp, "-BCR_", BCR, "-", PROJ, "-SD.tif"))
    writeRaster(Sd, fout, overwrite=TRUE)

}

## mosaic the mean and SD maps across regions

mosaic_fun <- function(fin, fout) {
    r4 <- raster(fin[1])
    r5 <- raster(fin[2])
    r60 <- raster(fin[3])
    r61 <- raster(fin[4])
    r70 <- raster(fin[5])
    r71 <- raster(fin[6])
    r80 <- raster(fin[7])
    r81 <- raster(fin[8])
    r82 <- raster(fin[9])
    r83 <- raster(fin[10])
    r9 <- raster(fin[11])
    r10 <- raster(fin[12])
    r11 <- raster(fin[13])
    r12 <- raster(fin[14])
    r13 <- raster(fin[15])
    r14 <- raster(fin[16])
    r200 <- raster(fin[17])
    r400 <- raster(fin[18])
    r500 <- raster(fin[19])
    r1000 <- raster(fin[20])
    r1100 <- raster(fin[21])
    r1200 <- raster(fin[22])
    r1300 <- raster(fin[23])
    r1400 <- raster(fin[24])
    rast <- mosaic(r4, r5, r60, r61, r70, r71,
        r80, r81, r82, r83, r9, r10, r11, r12, r13, r14,
        r200, r400, r500, r1000, r1100, r1200, r1300, r1400,
        fun=mean)
    writeRaster(rast, fout, overwrite=TRUE)
    invisible()
}

SPP <- c("CAWA", "OSFL", "OVEN", "BBWA")
#spp <- "OSFL"
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    fi1 <- file.path(ROOT, "artifacts2", spp, paste0(spp, "-BCR_", u, "-boot-Mean.tif"))
    fo1 <- file.path(ROOT, "artifacts2", spp, paste0(spp, "-BCR_ALL-boot-Mean.tif"))
    mosaic_fun(fi1, fo1)
    fi2 <- file.path(ROOT, "artifacts2", spp, paste0(spp, "-BCR_", u, "-boot-SD.tif"))
    fo2 <- file.path(ROOT, "artifacts2", spp, paste0(spp, "-BCR_ALL-boot-SD.tif"))
    mosaic_fun(fi2, fo2)
}

## now let's make B mosaiced layers 1st, and then average those
B <- 32
for (spp in SPP) {
    #b <- 1
    cat(spp, "\n")
    for (b in 1:B) {
        cat("\t", b, "\n")
        flush.console()
        fi <- file.path(ROOT, "artifacts2", spp, "parts", paste0("mosaic-", spp, "-BCR_", u, "-boot-", b, ".tif"))
        fo <- file.path(ROOT, "artifacts2", spp, "b", paste0("mosaic-", spp, "-BCR_ALL-boot-", b, ".tif"))
        mosaic_fun(fi, fo)
    }
}

for (spp in SPP) {
    #b <- 1
    cat(spp, "\n")
    fi <- file.path(ROOT, "artifacts2", spp, "b", paste0("mosaic-", spp, "-BCR_ALL-boot-", 1:B, ".tif"))
    st <- list()
    for (b in 1:B) {
        st[[b]] <- raster(fi[b])
    }

    rast <- mosaic(st[[1]], st[[2]], st[[3]], st[[4]], st[[5]], st[[6]], st[[7]], st[[8]],
        st[[9]], st[[10]], st[[11]], st[[12]], st[[13]], st[[14]], st[[15]], st[[16]],
        st[[17]], st[[18]], st[[19]], st[[20]], st[[21]], st[[22]], st[[23]], st[[24]],
        st[[25]], st[[26]], st[[27]], st[[28]], st[[29]], st[[30]], st[[31]], st[[32]],
        fun=mean)
    fo <- file.path(ROOT, "artifacts2", spp, paste0(spp, "-BCR_ALL-boot-Mean2.tif"))
    writeRaster(rast, fo, overwrite=TRUE)

    rast <- mosaic(st[[1]], st[[2]], st[[3]], st[[4]], st[[5]], st[[6]], st[[7]], st[[8]],
        st[[9]], st[[10]], st[[11]], st[[12]], st[[13]], st[[14]], st[[15]], st[[16]],
        st[[17]], st[[18]], st[[19]], st[[20]], st[[21]], st[[22]], st[[23]], st[[24]],
        st[[25]], st[[26]], st[[27]], st[[28]], st[[29]], st[[30]], st[[31]], st[[32]],
        fun=sd)
    fo <- file.path(ROOT, "artifacts2", spp, paste0(spp, "-BCR_ALL-boot-SD2.tif"))
    writeRaster(rast, fo, overwrite=TRUE)
}

## making nice maps

library(rgdal)
library(raster)

ROOT <- "d:/bam/BAM_data_v2019/gnm"

#bluegreen.colors <- colorRampPalette(c("#FFFACD", "lemonchiffon","#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=0.5)
bluegreen.colors <- function(n) {
    x <- hcl.colors(n, "Lajolla")
    c(x, rep(x[length(x)], ceiling(n*0.1)))
}

SUB <- readOGR(dsn=file.path(ROOT, "data", "sub-bound"), "BCRSubunits")
BCR <- readOGR(dsn=file.path(ROOT, "data", "bcr"), "bcrfinallcc")
PROV <- readOGR(dsn=file.path(ROOT, "data", "prov"), "province_state_line")
LAKES <- readOGR(dsn=file.path(ROOT, "data", "lakes"), "lakes_lcc")
LAKES <- spTransform(LAKES, proj4string(BCR))

for (spp in SPP) {
    cat(spp, "\n")
    flush.console()

    f <- file.path(ROOT, "artifacts2", spp, paste0(spp, "-BCR_ALL-boot-Mean.tif"))
    #f <- file.path(ROOT, "artifacts2", spp, paste0(spp, "-BCR_ALL-boot-Mean2.tif"))
    rast <- raster(f)
    MAX <- 2*cellStats(rast, 'mean')
    png(gsub("\\.tif", "\\.png", f), height=2000, width=3000)
    op <- par(cex.main=3, mfcol=c(1,1), oma=c(0,0,0,0), mar=c(0,0,5,0))
    plot(rast, col=bluegreen.colors(15), axes=FALSE, legend=TRUE, main=paste(spp, "mean"), box=FALSE)
    plot(PROV, col="grey", add=TRUE)
    plot(LAKES,col="#aaaaff", border=NA, add=TRUE)
    plot(BCR, add=TRUE)
    par(op)
    dev.off()

    f <- file.path(ROOT, "artifacts2", spp, paste0(spp, "-BCR_ALL-boot-SD.tif"))
    #f <- file.path(ROOT, "artifacts2", spp, paste0(spp, "-BCR_ALL-boot-SD2.tif"))
    rast <- raster(f)
    MAX <- 2*cellStats(rast, 'mean')
    png(gsub("\\.tif", "\\.png", f), height=2000, width=3000)
    op <- par(cex.main=3, mfcol=c(1,1), oma=c(0,0,0,0), mar=c(0,0,5,0))
    plot(rast, col=bluegreen.colors(15), axes=FALSE, legend=TRUE, main=paste(spp, "mean"), box=FALSE)
    plot(PROV, col="grey", add=TRUE)
    plot(LAKES,col="#aaaaff", border=NA, add=TRUE)
    plot(BCR, add=TRUE)
    par(op)
    par(op)
    dev.off()
}


## experiment with random buffering
for (spp in c("OSFL", "OVEN", "BBWA")) {
#BCRs <- c(71, 80, 81)
BCRs <- u[u < 200]

B <- 32
st <- list()
for (i in 1:B) {
    cat(spp, i, "\n");flush.console()
    for (bcr in BCRs) {
        bi <- if (B > 32)
            sample(B, 1) else i
        p <- SUB[SUB@data[,1] == bcr,]
        ## linear (uniform): plot(ecdf(runif(10^3)))
        #w <- runif(1, 0, 1)
        ## sigmoid cdf: plot(ecdf(rbeta(10^3, 2, 2)))
        w <- rbeta(1, 2, 2)
        pb <- rgeos::gBuffer(p, width=w * 100*10^3)
        ri <- raster(file.path(ROOT, "artifacts2", spp, "parts",
            paste0("mosaic-", spp, "-BCR_", bcr, "-boot-", bi, ".tif")))
        pb <- rgeos::gBuffer(p, width=w * 100 * 10^3)
        if (bcr == BCRs[1]) {
            rout <- mask(ri, pb)
            #pb0 <- rgeos::gBuffer(p, width=100*10^3)
            #plot(ri2)
            #plot(p, add=TRUE)
            #plot(pb0, add=TRUE, border=2)
        } else {
            rout <- mosaic(rout, mask(ri, pb), fun=mean)
        }
    }
    st[[i]] <- rout
    if (i == 1) {
        rast <- rout
    } else {
        #rast <-  ((i-1)/i) * rast + (1-((i-1)/i)) * rout
        #rast <- overlay(rast, rout, fun=function(x,y) {((i-1)/i) * x + (1-((i-1)/i)) * y})
        #isNA <- is.na(values(rast))
        rast <- mosaic(((i-1)/i) * rast, (1-((i-1)/i)) * rout, fun=sum)
        #values(rast)[isNA & is.na(rout)] <- NA
    }
}
#rast2 <- mosaic(st[[1]], st[[2]], st[[3]], st[[4]], st[[5]], st[[6]], st[[7]], st[[8]],
#    st[[9]], st[[10]], fun=mean)

#plot(rast, main=1)
#plot(rast2, main=2)

writeRaster(rast, overwrite=TRUE,
    filename=file.path(ROOT, "artifacts2", spp, paste0(spp, "-BCR_ALL-boot-MeanRND.tif")))

q <- quantile(values(rast), 0.999, na.rm=TRUE)
values(rast)[!is.na(values(rast)) & values(rast)>q] <- q
    png(file.path(ROOT, "artifacts2", spp, paste0(spp, "-BCR_ALL-boot-MeanRND.png")),
        height=2000, width=3000)
    op <- par(cex.main=3, mfcol=c(1,1), oma=c(0,0,0,0), mar=c(0,0,5,0))
    plot(rast, col=bluegreen.colors(15), axes=FALSE, legend=TRUE, main=paste(spp, "feathered"), box=FALSE)
    plot(PROV, col="grey", add=TRUE)
    plot(LAKES,col="#aaaaff", border=NA, add=TRUE)
    plot(BCR, add=TRUE)
    par(op)
    dev.off()
}
