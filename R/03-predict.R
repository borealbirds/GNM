library(mefa4)
library(raster)
library(gbm)
library(dismo)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
PROJ <- "run3"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2019-10-29.RData"))

#SPP <- colnames(yy)
#SPP <- c("ALFL", "AMRO", "BOCH", "BTNW", "CAWA",
#    "OSFL", "OVEN", "RUBL", "WCSP", "YRWA")
SPP <- "OSFL"
SPPBCRss <- NULL
for (spp in SPP) {
    SPPBCRss <- c(SPPBCRss, SPPBCR[grep(spp, SPPBCR)])
}

DONE <- sapply(strsplit(list.files(file.path(ROOT, "out", PROJ)), ".", fixed=TRUE), function(z) z[1L])
done <- data.frame(table(sapply(strsplit(gsub(".RData", "", DONE), "-"), "[[", 2)))
tab <- data.frame(table(sapply(strsplit(SPPBCR, "-"), "[[", 2)))
tab$done <- done$Freq[match(tab$Var1, done$Var1)]
tab$done[is.na(tab$done)] <- 0
tab$perc <- 100 * tab$done / tab$Freq
tab


## prepare stacks
for (i in u) {
    cat("\nLoading stack for BCR", i, "\n")
    flush.console()
    gc()

    pr <- stack(file.path(ROOT, "data", "subunits", paste0("bcr", i, "all_1km.grd")))

    pset <- CN[[paste0("BCR_", i)]]
    print(compare_sets(pset, names(pr)))
#    if (i == 4)
#        pset <- c(pset, "MAP", "FFP", "Landsc750_Cham_Noo_v1")
#    if (i == 11)
#        pset <- c(pset, "Landsc750_Ostr_Vir_v1", "Landsc750_Fagu_Gra_v1",
#            "Landsc750_Tsug_Can_v1", "Landsc750_Quer_Rub_v1",
#           "Species_Acer_Sah_v1")
#    if (i == 83)
#        pset <- unique(unlist(CN))

    n <- length(values(pr[[1]]))
    nd <- matrix(0, n, length(pset))
    colnames(nd) <- pset
    for (j in pset) {
        #cat(j, "\n");flush.console()
        nd[,j] <- values(pr[[j]])
    }
    notNA <- rowSums(is.na(nd)) == 0
    #notNA <- !is.na(nd[,1])
    nd <- as.data.frame(nd[notNA,,drop=FALSE])
    nd$nalc <- as.factor(nd$nalc)
#    nd$lf <- as.factor(nd$lf)

    nd$offset <- 0L
    nd$weights <- 1L
    ND <- list(data=nd, n=n, subset=which(notNA), dim=dim(pr))
    saveRDS(ND, file=file.path(ROOT, paste0("STACK-ND-BCR_", i, ".rds")))

}

## need to set NALC land cover types not used in model (snow/ice) to 0

## predict ---------------------------------

library(mefa4)
library(gbm)
library(raster)
ROOT <- "d:/bam/BAM_data_v2019/gnm"
#ROOT <- "c:/p/tmp/gnm"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2019-10-29.RData"))

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


PROJ <- "run3"
#BCR <- 83 # this is BCR
for (BCR in u) {
    cat("loading stack:", BCR, "\n")
    flush.console()
    #spp <- "OSFL"
    r1 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr", BCR, "all_1km.grd")))
    ND <- readRDS(file.path(ROOT, paste0("STACK-ND-BCR_", BCR, ".rds")))
    SPPss <- substr(SPPBCR[grep(paste0("BCR_", BCR), SPPBCR)], 1, 4)
    #SPPss <- c("OSFL","CAWA", "BBWA", "BLPW")
    #spp <- SPPss[1]
    #SPPss <- SPPss[which(SPPss==spp):length(SPPss)]
    for (spp in SPPss) {
        gc()
        cat("\t", spp, "@ BCR", BCR, ">", which(SPPss == spp), "/", length(SPPss))
        flush.console()
        fin <- file.path(ROOT, "out", PROJ, paste0(spp, "-BCR_", BCR, ".RData"))
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
        if (!dir.exists(file.path(ROOT, "artifacts", spp)))
            dir.create(file.path(ROOT, "artifacts", spp))
        fout <- file.path(ROOT, "artifacts", spp,
            paste0("mosaic-", spp, "-BCR_", BCR, "-", PROJ, ".tif"))
        writeRaster(rrr, fout, overwrite=TRUE)
        cat(" > OK\n")
    }
}

## filling in the missing pieces
SPP <- colnames(yy)
SPPBCRall <- levels(interaction(as.factor(SPP)[1], paste0("BCR_", u), sep="-"))
str(SPPBCRall)
length(SPP) * length(u)
for (BCR in u) {
    cat("loading stack:", BCR, "\n")
    flush.console()
    r1 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr", BCR, "all_1km.grd")))
    ND <- readRDS(file.path(ROOT, paste0("STACK-ND-BCR_", BCR, ".rds")))
    SPPss1 <- substr(SPPBCR[grep(paste0("BCR_", BCR), SPPBCR)], 1, 4)
    SPPss2 <- substr(SPPBCRall[grep(paste0("BCR_", BCR), SPPBCRall)], 1, 4)
    SPPss <- setdiff(SPPss2, SPPss1)
    for (spp in SPPss) {
        gc()
        cat("\t", spp, "@ BCR", BCR, ">", which(SPPss == spp), "/", length(SPPss))
        flush.console()
        rrr <- predict_gbm(NULL, ND, r1, 0)
        fout <- file.path(ROOT, "artifacts", spp,
            paste0("mosaic-", spp, "-BCR_", BCR, "-", PROJ, ".tif"))
        writeRaster(rrr, fout, overwrite=FALSE)
        cat(" > OK\n")
    }
}


SPPss <- colnames(yy)
mosaic_fun <- function(spp) {
    fout <- file.path(ROOT, "artifacts", "00mosaic", paste0("mosaic-", spp, "-", PROJ, ".tif"))
    fin <- file.path(ROOT, "artifacts", spp,
        paste0("mosaic-", spp, "-BCR_", u, "-", PROJ, ".tif"))
    if (!file.exists(fout)) {
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
        rast <- mosaic(r4, r5, r60, r61, r70, r71,
            r80, r81, r82, r83, r9, r10, r11, r12, r13, r14, fun=mean)
        writeRaster(rast, fout, overwrite=TRUE)
    }
    invisible()
}
#SPPss <- SPPss[32:150]
for (spp in SPPss) {
    cat(spp)
    flush.console()
    aa <- try(mosaic_fun(spp))
    if (!inherits(aa, "try-error"))
        cat(" - OK\n") else cat(" -------- ERROR\n")
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
for (spp in c("CAWA", "MAWA", "OSFL", "RCKI", "RUBL")) {
    #fout0 <- file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-", PROJ, ".tif"))
    if (file.exists(fout0)) {
        cat(spp, "\n")
        flush.console()

        rast <- raster(file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-", PROJ, ".tif")))

        ## play with Lc
        #lc <- lorenz(values(rast)[!is.na(values(rast))])
        #q <- quantile(lc, probs=c(0.05, 0.1, 0.2, 0.5, 0.8, 0.99), type="L")
        #MAX <- q["80%"]
        MAX <- 2*cellStats(rast, 'mean')
        png(file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-", PROJ, ".png")),
            height=2000, width=3000)
#        png(file.path(ROOT, "artifacts", spp, paste0("mosaic-", spp, "-", PROJ, "-nalcfix.png")),
#            height=2000, width=3000)
        op <- par(cex.main=3, mfcol=c(1,1), oma=c(0,0,0,0), mar=c(0,0,5,0))
        plot(rast, col="blue", axes=FALSE, legend=FALSE, main=paste(spp), box=FALSE)
        plot(rast, col=bluegreen.colors(15), zlim=c(0,MAX), axes=FALSE,
            main=spp, add=TRUE, legend.width=1.5, horizontal = TRUE)#,
            #smallplot = c(0.60,0.85,0.82,0.87), axis.args=list(cex.axis=2))
        plot(PROV, col="grey", add=TRUE)
        plot(LAKES,col="#aaaaff", border=NA, add=TRUE)
        plot(BCR, add=TRUE)
        par(op)
        dev.off()
if (FALSE) {
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
