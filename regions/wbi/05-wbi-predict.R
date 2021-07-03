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

qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-30.qRData")
qload("d:/bam/2021/wbi/pred.qRData")


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
for (spp in SPPx) {
    cat(spp, "\n")
    flush.console()
    gc()


    ## load overall models and make predictions
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qload(fn)
    ## predict
    X <- model.matrix(~.,prd[,RES$vars])
    pr_gbm <- predict(RES$gbm, prd, type="response")
    pr <- data.frame(
        prd[,c("index", "region", "nalc")],
        gbm=pr_gbm)
    ## smoothing based on pred
    ## clip GNM to extent
    rspp <- raster(paste0("d:/bam/BAM_data_v2019/gnm/artifacts/",
        spp, "/pred-", spp, "-CAN-Mean.tif"))
    rspp <- crop(rspp, r)
    rspp <- mask(rspp, r)
    vrspp <- values(rspp)
    pr$gnm <- vrspp[pr$index]

    if (!dir.exists(paste0("d:/bam/2021/wbi/pred/", spp)))
        dir.create(paste0("d:/bam/2021/wbi/pred/", spp))
    fn <- paste0("d:/bam/2021/wbi/pred/", spp, "/", "WB-", spp, "-pred.qRData")
    qsavem(pr, file=fn)

    if (FALSE) {
        fo <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-pred-", i, ".png")
        png(fo, width=1200, height=500)
        op <- par(mfrow=c(1,2), mai=c(0.1,0.1,0.5,1.5))
        for (WHAT in c("gnm", "gbm")) {
            rs <- recast_pred(pr[[WHAT]])
            plot(rs, axes=FALSE, box=FALSE, main=toupper(WHAT),
                col=hcl.colors(20, "lajolla"))
            plot(wb$geometry,add=TRUE, border="#00000044")
        }
        par(op)
        dev.off()
    }

}

for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qload(fn)

    fo2 <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-resp-", i, ".png")
    png(fo2, width=1500, height=1000)
    p <- plot_fun(RES$gbm)
    print(p)
    dev.off()

    fo3 <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-resp-nalc-", i, ".png")
    png(fo3)
    print(a)
    dev.off()
}

## NALC w proper labels

## regional lc
LCCs <- c(
    "Conifer"=1,
    "Taiga Conifer"=2,
    "Deciduous"=5,
    "Mixedwood"=6,
    "Shrub"=8,
    "Grass"=10,
    "Arctic Shrub"=11,
    "Arctic Grass"=12,
    "Arctic Barren"=13,
    "Wetland"=14,
    "Cropland"=15,
    ## exclude ---
    "Barren Lands"=16, # 0
    "Urban and Built-up"=17, # ?
    "Water"=18, # 0
    "Snow and Ice"=19) # 0
i <- 1
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    fn <- paste0("d:/bam/2021/wbi/pred/", spp, "/", "WB-", spp, "-pred.qRData")
    qload(fn)
    pr <- pr[as.numeric(as.character(pr$nalc)) <= 17,]

    s <- as.data.frame(sum_by(pr$gbm, pr$nalc))
    s$Density <- s$x/s$by
    s$nalc <- as.numeric(rownames(s))
    s <- s[order(s$nalc),]
    s$NALC <- factor(names(LCCs)[match(s$nalc, LCCs)], names(LCCs))
    s$Source <- factor("GBM", c("GNM", "GBM"))
    s1 <- s

    s <- as.data.frame(sum_by(pr$gnm, pr$nalc))
    s$Density <- s$x/s$by
    s$nalc <- as.numeric(rownames(s))
    s <- s[order(s$nalc),]
    s$NALC <- factor(names(LCCs)[match(s$nalc, LCCs)], names(LCCs))
    s$Source <- factor("GNM", c("GNM", "GBM"))
    s <- rbind(s, s1)
    rm(s1)


    p <- ggplot(s, aes(x=NALC, y=Density)) +
    #    geom_bar(stat="identity", fill="#95B6C1") +
        geom_bar(stat="identity") +
        coord_flip() +
        facet_grid(cols=vars(Source)) +
        theme_minimal()
    fo <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-avg-nalc-", i, ".png")
    png(fo, width=800, height=400)
    print(p)
    dev.off()
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
SPPx <- SPP[1:20]
SPPx <- SPP[21:40]
SPPx <- SPP[41:60]
SPPx <- SPP[61:80]
SPPx <- SPP[81:100]
SPPx <- SPP[101:117]

for (spp in SPPx) {

    i <- 1
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qload(fn)

    prd <- cbind(dd, ddvar)[RES$pk,]
    prd <- nonDuplicated(prd, SS)
    prd$YEAR <- 2011
    prd$ARU <- 0
    prd <- prd[,RES$vars]

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
    o <- epi.occc(PR)
    save(o, PR, file=paste0("d:/bam/2021/wbi/occc/", spp, ".RData"))

}

for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    load(paste0("d:/bam/2021/wbi/occc/", spp, ".RData"))

    j <- sample(nrow(PR), 10*10^3)
    z <- PR[j,]
    rm <- rowMeans(z)
    z <- z[order(rm),]
    fo <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-occc-10.png")
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
}

OCCC <- list()
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    load(file=paste0("d:/bam/2021/wbi/occc/", spp, ".RData"))
    OCCC[[spp]] <- o
}
save(OCCC, file="d:/bam/2021/wbi/OCCC.RData")


AUC <- list()
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    i <- 1
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qload(fn)
    AUC[[spp]] <- RES$AUC
}
save(AUC, file="d:/bam/2021/wbi/AUC.RData")


i <- 1
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qload(fn)
    prd <- cbind(dd, ddvar)[RES$pk,RES$vars]
    pr_gbm <- pfun(RES$gbm, prd, type="response")
    y <- as.numeric(Y[RES$pk,spp])
    o <- as.numeric(O[RES$pk,spp])
    Mx <- max(y[y<100])
    y[y>Mx] <- Mx
    Predicted_GBM_w_Offset <- exp(o)*pr_gbm
    Observed_Count <- y
    fo <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-box-", i, ".png")
    png(fo, width=1500, height=1000)
    boxplot(Predicted_GBM_w_Offset ~ Observed_Count, range=0, col="gold", main=spp)
    dev.off()
}


sf <- sf::st_as_sf(dd, coords = c("X","Y"))
sf <- st_set_crs(sf, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sf <- st_transform(sf, lcc_crs)

for (spp in SPP) {
    cat(spp, "\n")
    flush.console()

    y <- ifelse(as.numeric(Y[,spp])>0, 1, 0)
    i0 <- y==0 & !duplicated(dd$SS)
    i1 <- y==1 & !duplicated(dd$SS)

    fo <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-det.png")
    png(fo, width=1500, height=1000)
    plot(sf[i0 | i1,"geometry"], pch=19, cex=0.5, col="grey", main=spp)
    plot(wb$geometry, border=1, add=TRUE)
    plot(sf[i1,"geometry"], pch=19, cex=0.5, col=2, add=TRUE)
    dev.off()
}


## pull all the info

library(mefa4)
library(jsonlite)
qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-30.qRData")
load("d:/bam/2021/wbi/occc.RData")
load("d:/bam/2021/wbi/AUC.RData")
tab <- fromJSON("https://borealbirds.github.io/api/v4/species/")
rownames(tab) <- tab$id
tab <- tab[SPP,c("scientific", "english")]

LCCs <- c(
    "Conifer"=1,
    "Taiga Conifer"=2,
    "Deciduous"=5,
    "Mixedwood"=6,
    "Shrub"=8,
    "Grass"=10,
    "Arctic Shrub"=11,
    "Arctic Grass"=12,
    "Arctic Barren"=13,
    "Wetland"=14,
    "Cropland"=15,
    ## exclude ---
    "Barren Lands"=16, # 0
    "Urban and Built-up"=17, # ?
    "Water"=18, # 0
    "Snow and Ice"=19) # 0

#spp <- "ALFL"
load("d:/bam/2021/wbi/OCCC.RData")
load("d:/bam/2021/wbi/AUC.RData")
ALLRES <- NULL
for (spp in SPP) {
    i <- 1
    cat(spp, "\n")
    flush.console()

    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qload(fn)

    fn <- paste0("d:/bam/2021/wbi/pred/", spp, "/", "WB-", spp, "-pred.qRData")
    qload(fn)
    pr <- pr[as.numeric(as.character(pr$nalc)) <= 17,]

    Ngnm <- 2*100*sum(pr$gnm)/10^6
    Ngbm <- 2*100*sum(pr$gbm)/10^6
    AUC <- unname(RES$AUC["gbm"])
    ALLRES <- rbind(ALLRES, data.frame(
        ID=spp, Species=tab[spp,"english"],
        N_GNM=Ngnm, N_GBM=Ngbm, AUC=AUC,
        P_Occ=mean(ifelse(Y[,spp]>0,1,0)),
        OCCC=OCCC[[spp]]["occc"],
        Precision=OCCC[[spp]]["oprec"],
        Accuracy=OCCC[[spp]]["oaccu"]))
}

write.csv(ALLRES, row.names = FALSE, file="d:/bam/2021/wbi/whi-stats.csv")


x <- ALLRES

plot(x[,c("oaccu", "oprec")], xlim=c(0,1), ylim=c(0,1))
abline(h=0.5, v=0.5, lty=2)

plot(x[,c("P_Occ", "AUC")], xlim=c(0,1), ylim=c(0,1))
abline(h=0.5, lty=2)

plot(x[,c("N_GNM", "N_GBM")])
abline(0, 1, lty=2)

x[order(x$AUC),][1:10,c(1:2,5:9)]
x[order(x$oprec),][1:10,c(1:2,5:9)]
x[order(x$P_Occ),][1:10,c(1:2,5:9)]
