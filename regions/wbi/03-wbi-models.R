library(mefa4)
library(dismo)
library(gbm)
library(qs)
source("~/repos/GNM/regions/wbi/functions.R")
qs::qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-30.qRData")
load("d:/bam/2021/wbi/daic.RData")

SU <- list("WBAB"=60,
    "WBBC"=c(4,60),
    "WBYT"=4,
    "WBMB"=c(60,70,80),
    "WBSK"=c(60,80),
    "WBNT"=c(61, 70))


if (FALSE) {
## this is the part where we screen the important predictors

dim(ddvar)
M <- as.matrix(ddvar[,colnames(ddvar) != "nalc"])
CM <- cor(M)
SD <- apply(M, 2, sd)
MN <- colMeans(M)
CoV <- SD/MN
M2 <- M[,SD>0.001]
Uv <- apply(M2, 2, function(z) {
    q <- quantile(z, c(0.01, 0.99))
    z <- z[z>=q[1] & z<=q[2]]
    length(unique(z))
})
table(Uv)
M3 <- M2[,Uv > 99]

get_cn <- function(z, rmax=0.9) {
    SD <- apply(z, 2, sd)
    COR <- cor(z[,SD > 0])
    cr <- mefa:::stack.dist(as.dist(COR), dim.names = TRUE)
    cr <- cr[order(abs(cr$dist), decreasing = TRUE),]
    cr[,1] <- as.character(cr[,1])
    cr[,2] <- as.character(cr[,2])
    cr$Landsc1 <- startsWith(cr[,1], "Landsc750_")
    cr$Landsc2 <- startsWith(cr[,2], "Landsc750_")
    cr1 <- cr[cr$Landsc1 | cr$Landsc2,]
    cr2 <- cr[!(cr$Landsc1 | cr$Landsc2),]
    while(any(abs(cr1$dist) > rmax)) {
        i <- if (cr1$Landsc1[1])
            cr1[1,1] else cr1[1,2]
        j <- cr1[,1] == i | cr1[,2] == i
        cr1 <- cr1[!j,]
    }
    cr3 <- rbind(cr1, cr2)
    cr3 <- cr3[order(abs(cr3$dist), decreasing = TRUE),]
    while(any(abs(cr3$dist) > rmax)) {
        i <- if (cr3$Landsc1[1])
            cr3[1,1] else cr3[1,2]
        j <- cr3[,1] == i | cr3[,2] == i
        cr3 <- cr3[!j,]
    }
    union(as.character(cr3[,1]), as.character(cr3[,2]))
}

CN <- sort(get_cn(M3))

CN <- c("AHM", "CMD", "dev750", "EMT", "LandCover_Veg_v1", "LandCover_Veg_v1.1",
    "LandCover_VegNonTreed_v1", "LandCover_VegNonTreed_v1.1", "LandCover_VegTreed_v1.1",
    "Landsc750_Abie_Bal_v1", "Landsc750_Abie_Las_v1", "Landsc750_Abie_Spp_v1",
    "Landsc750_Acer_Neg_v1", "Landsc750_Acer_Rub_v1", "Landsc750_Betu_Pap_v1",
    "Landsc750_Betu_Spp_v1", "Landsc750_Frax_Ame_v1", "Landsc750_Frax_Nig_v1",
    "Landsc750_Frax_Pen_v1", "Landsc750_Genc_Spp_v1", "Landsc750_Lari_Lar_v1",
    "Landsc750_Lari_Lya_v1", "Landsc750_Needleleaf_Spp_v1", "Landsc750_Pice_Eng_v1",
    "Landsc750_Pice_Gla_v1", "Landsc750_Pice_Mar_v1", "Landsc750_Pice_Spp_v1",
    "Landsc750_Pinu_Ban_v1", "Landsc750_Pinu_Res_v1", "Landsc750_Pinu_Spp_v1",
    "Landsc750_Pinu_Str_v1", "Landsc750_Popu_Bal_v1", "Landsc750_Popu_Spp_v1",
    "Landsc750_Popu_Tri_v1", "Landsc750_Prun_Pen_v1", "Landsc750_Pseu_Men_v1",
    "Landsc750_Quer_Mac_v1", "Landsc750_Sali_Spp_v1", "Landsc750_Stand_Age_v1",
    "Landsc750_Thuj_Occ_v1", "Landsc750_Ulmu_Ame_v1", "led750", "MAP",
    "PPT_sm", "SHM", "slope", "Species_Abie_Bal_v1", "Species_Abie_Las_v1",
    "Species_Acer_Neg_v1", "Species_Betu_Pap_v1", "Species_Betu_Spp_v1",
    "Species_Lari_Lar_v1", "Species_Pice_Gla_v1", "Species_Pice_Mar_v1",
    "Species_Pice_Spp_v1", "Species_Pinu_Ban_v1", "Species_Pinu_Con_v1",
    "Species_Pinu_Spp_v1", "Species_Popu_Bal_v1", "Species_Popu_Spp_v1",
    "Species_Quer_Mac_v1", "Species_Sali_Spp_v1", "SpeciesGroups_Broadleaf_Spp_v1",
    "SpeciesGroups_Needleleaf_Spp_v1", "Structure_Stand_Age_v1",
    "TD", "TPI")

#spp <- "ALFL"
#DAIC <- list()
for (spp in SPP) {
    d <- get_data_all(spp, replace=FALSE, cn0=CN)

    m0 <- mgcv::gam(count ~ 1, data=d, family=poisson, offset=d$offset)
    aic <- c(Null=AIC(m0))

#    for (i in seq_len(length(CN))) {
    for (i in 47:67) {
        V <-  CN[i]
        cat(spp, i, V, "\n")
        flush.console()
        fm <- as.formula(sprintf("count ~ s(%s)", V))
        m <- try(mgcv::gam(formula=fm, data=d, family=poisson, offset=d$offset))
        newaic <- if (inherits(m, "try-error"))
            Inf else AIC(m)
        #plot(m)
        aic <- c(aic, newaic)
        names(aic)[i+1] <- V
    }
    daic <- sort(aic-min(aic))
    save(daic, aic, file=paste0("d:/bam/2021/wbi/daic/", spp, ".RData"))
    #DAIC[[spp]] <- daic
}

for (spp in names(DAIC)) {
    daic <- DAIC[[spp]]
    save(daic, file=paste0("d:/bam/2021/wbi/daic/", spp, ".RData"))

}

DAIC <- list()
for (spp in SPP) {
    load(paste0("d:/bam/2021/wbi/daic/", spp, ".RData"))
    DAIC[[spp]] <- daic
}
save(CN, DAIC, file="d:/bam/2021/wbi/daic.RData")


}




dim(get_data_by_reg("OVEN", "WBAB"))
dim(get_data_by_reg("OVEN", "WBBC"))
dim(get_data_by_reg("OVEN", "WBMB"))
dim(get_data_by_reg("OVEN", "WBNT"))
dim(get_data_by_reg("OVEN", "WBSK"))
dim(get_data_by_reg("OVEN", "WBYT"))
dim(get_data_all("OVEN", cn0=CN))
dim(get_data_all("OVEN", cn0=CN, replace=TRUE))

## one run for each spp
## use all data
SPPx <- SPP[1:30]
#SPPx <- SPP[31:60]
#SPPx <- SPP[61:90]
#SPPx <- SPP[91:117]

i <- 1
for (spp in SPPx) {
    gc()
    cat(i, spp, "\n")
    flush.console()
    tmp <- try(fit_fun(i, spp, reg=NULL))
    if (inherits(tmp, "try-error"))
        tmp <- structure(as.character(tmp), class="try-error")
    RES <- tmp
    dir.create(paste0("d:/bam/2021/wbi/out/", spp))
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qsavem(RES, file=fn)
}


## some boot
SPPx <- SPP[1:30]
SPPx <- SPP[31:60]
SPPx <- SPP[61:90]
SPPx <- SPP[91:117]

Bv <- 2:10
for (i in Bv) {
    for (spp in SPPx) {
        gc()
        cat(i, spp, "\n")
        flush.console()
        fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
        tmp <- try(fit_fun(i, spp, reg=NULL))
        if (inherits(tmp, "try-error"))
            tmp <- structure(as.character(tmp), class="try-error")
        RES <- tmp
        if (!dir.exists(paste0("d:/bam/2021/wbi/out/", spp)))
            dir.create(paste0("d:/bam/2021/wbi/out/", spp))
        #fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
        qsavem(RES, file=fn)
    }
}

# 100 run for 1 spp
spp <- "ALFL"
for (i in 2:20) {
    gc()
    cat(i, spp, "\n")
    flush.console()
    tmp <- try(fit_fun(i, spp, reg=NULL))
    if (inherits(tmp, "try-error"))
        tmp <- structure(as.character(tmp), class="try-error")
    RES <- tmp
    dir.create(paste0("d:/bam/2021/wbi/out/", spp))
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qsavem(RES, file=fn)
}

## use all the variables for certain species to see how to improve models

SPP2 <- c("RECR", "PUFI", "LCSP", "LISP", "MAWR", "SWSP", "TOSO",
    "DOWO", "PIGR", "ATTW", "BEKI", "EAKI", "EUST", "GCSP", "HAWO", "HETH", "LALO")

i <- 1
for (spp in SPPx) {
    gc()
    cat(i, spp, "\n")
    flush.console()
    if (!(spp %in% SPP2)) {
        tmp <- try(fit_fun(i, spp, reg=NULL, cn=CN, gbm_only=TRUE))
        if (inherits(tmp, "try-error"))
            tmp <- structure(as.character(tmp), class="try-error")
        RES <- tmp
        dir.create(paste0("d:/bam/2021/wbi/out/", spp))
        fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-fullCN-", i, ".qRData")
        qsavem(RES, file=fn)
    }
}

## find top variables
spp <- "RECR"
fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-fullCN-", i, ".qRData")
qload(fn)

r <- rel_inf(RES$gbm)
sort(names(DAIC[[spp]][1:10]))
sort(rownames(r)[1:10])
mefa4::compare_sets(names(DAIC[[spp]][1:10]), rownames(r)[1:10])

## compare AUC from the 3 options

#spp <- "ALFL"
auc <- list()
for (spp in SPP2) {
    cat(spp, "\n")
    flush.console()

    fn0 <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    fn1 <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-fullCN-", i, ".qRData")
    fn2 <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-refit10-", i, ".qRData")
    qload(fn0)
    RES0 <- RES
    qload(fn1)
    RES1 <- RES
    qload(fn2)
    RES2 <- RES

    auc[[spp]] <- c(RES0$AUC[c("y", "off", "gbm")],
        gbm1=unname(RES1$AUC["gbm"]),
        gbm2=unname(RES2$AUC["gbm"]))

}


## refit


i <- 1
for (spp in SPP2) {
    gc()
    cat(i, spp, "\n")
    flush.console()
    fn1 <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-fullCN-", i, ".qRData")
    qload(fn1)
    RES0 <- RES
    rm(RES)
    if (!inherits(RES0, "try-error")) {
        r <- rel_inf(RES0$gbm)
        cn <- rownames(r)[1:10]
        tmp <- try(fit_fun(i, spp, reg=NULL, cn=cn))
        if (inherits(tmp, "try-error"))
            tmp <- structure(as.character(tmp), class="try-error")
        RES <- tmp
        fn2 <- paste0("d:/bam/2021/wbi/out/", spp, "/",
            "WB-", spp, "-refit10-", i, ".qRData")
        qsavem(RES, file=fn2)
    }
}

