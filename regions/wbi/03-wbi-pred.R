library(mefa4)
library(dismo)
library(gbm)
library(qs)
library(raster)
library(sf)
library(sp)
library(rgdal)
library(segmented)
library(ggplot2)

source("~/repos/GNM/regions/wbi/functions.R")

lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
wb <- readOGR("d:/bam/2021/wbi/study-area")
wb <- st_as_sf(wb)
wb <- st_transform(wb, lcc_crs)
levels(wb$HASC_1) <- gsub("CA\\.", "WB", levels(wb$HASC_1))

#qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-11.qRData")
qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-25.qRData")
qload("d:/bam/2021/wbi/pred.qRData")
#load("~/repos/GNM/regions/wbi/subsets.RData")


## !!! Important !!!
## - make sure NALC is factor
## - make sure ARU is all 0
## - make sure YEAR is set to something that represent the training data well: 2011
prd$nalc <- as.factor(prd$nalc)
prd <- prd[rowSums(is.na(prd))==0,]

r <- raster(paste0("d:/bam/2021/wbi/mosaics/AHM.tif"))
values(r)[!is.na(values(r))] <- 0



#spp <- "ALFL"
i <- 1
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    gc()

    fo <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-pred-", i, ".png")
    if (!file.exists(fo)) {

        ## load overall models and make predictions
        fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
        qload(fn)
        ## predict
        X <- model.matrix(~.,prd[,RES$vars])
        pr_gbm <- predict(RES$gbm, prd, type="response")
        pr_glm <- exp(drop(X[,names(RES$glm)] %*% RES$glm))
        pr_ssa <- exp(drop(X[,names(RES$ssa)] %*% RES$ssa))
        pr <- data.frame(prd[,c("index", "region", "nalc")],
            gbm=pr_gbm,
            glm=pr_glm,
            ssa=pr_ssa)
        ## smoothing based on pred
        m_ssa2 <- lm(log(pr$gbm) ~ ., data=prd[,RES$vars])
        pr$ssa2 <- exp(fitted(m_ssa2))
        ## clip GNM to extent
        rspp <- raster(paste0("d:/bam/BAM_data_v2019/gnm/artifacts/",
            spp, "/pred-", spp, "-CAN-Mean.tif"))
        rspp <- crop(rspp, r)
        rspp <- mask(rspp, r)
        vrspp <- values(rspp)
        pr$gnm <- vrspp[pr$index]

        W <- c("gnm", "gbm", "glm", "ssa", "ssa2")
#        fo <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-pred-", i, ".png")
        png(fo, width=1500, height=1000)
        op <- par(mfrow=c(2,3))
        for (WHAT in W) {
            rs <- recast_pred(pr[[WHAT]])
            plot(rs, axes=FALSE, box=FALSE, main=WHAT,
                col=hcl.colors(20, "lajolla"))
            plot(wb$geometry,add=TRUE,border="#00000044")
        }
        par(op)
        dev.off()

        fo2 <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-resp-", i, ".png")
        png(fo2, width=1500, height=1000)
        plot_fun(RES$gbm)
        dev.off()

        fo3 <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-resp-nalc-", i, ".png")
        png(fo3)
        plot(RES$gbm, "nalc")
        dev.off()

    }
}

#qsavem(pr, file=paste0("d:/bam/2021/wbi/res/", spp, "/", "WB-", spp, "-preds.qRData"))


## TODO:
## - explain how to assemble df for prediction
## - code to do GBM prediction + super smoothing


## check bootstrap OCCC

library(epiR)
library(mefa4)
library(dismo)
library(gbm)
library(qs)
library(raster)

qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-25.qRData")
#load("~/repos/GNM/regions/wbi/subsets4.RData")

B <- 10

pfun <- function(...) {
    suppressMessages(suppressWarnings(predict(...)))
}

#spp <- "ALFL"
oo <- list()
for (spp in SPP) {

    i <- 1
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qload(fn)
    prd <- cbind(dd, ddvar)[RES$pk,RES$vars]

    prd$YEAR <- 2011
    prd$ARU <- 0
    PR <- matrix(0, nrow(prd), B)

    cat(spp, "\t.")
    pr_gbm <- pfun(RES$gbm, prd, type="response")
    PR[,i] <- pr_gbm

    for (i in 2:B) {
        cat(".")
        flush.console()
        fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
        qload(fn)
        pr_gbm <- pfun(RES$gbm, prd, type="response")
        PR[,i] <- pr_gbm
    }
    cat("\tOK\n")

    # http://web1.sph.emory.edu/observeragreement/OCCC.pdf
    # precision (degree of variation)
    # accuracy  (degree of location or scale shift)
    # overall: prec * acc
    o <- epi.occc(PR)

    j <- sample(nrow(prd), 10*10^3)
    z <- PR[j,]
    rm <- rowMeans(z)
    z <- z[order(rm),]
    fo <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-occc-", i, ".png")
    png(fo, width=1500, height=1000)
    matplot(1:nrow(z), z, type="l", lty=1, col="#00000022",
        ylim=quantile(z, c(0.001, 0.999)),
        axes=FALSE, ann=FALSE)
#    box()
    axis(2)
    title(main=paste(spp, "GBM in WB"), ylab="Expected density")
    lines(1:nrow(z), rowMeans(z), col="white", lwd=3)
    lines(1:nrow(z), rowMeans(z), col=2)
    legend("topleft", bty="n", col=c(NA, NA, NA),
        legend=paste(c("OCCC", "Prec", "Accu"), "=",
            round(unlist(o[1:3]), 3)))
    dev.off()
    oo[[spp]] <- unlist(o[1:3])
}

save(oo, file="d:/bam/2021/wbi/occc.RData")

ooo <- do.call(rbind, oo)
ooo <- ooo[order(ooo[,1]),]


