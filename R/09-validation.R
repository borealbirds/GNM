library(mefa4)
library(raster)
library(gbm)
library(sf)
library(rgdal)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2020-01-08.RData"))
SUB <- readOGR(dsn=file.path(ROOT, "data", "sub-bound"), "BCRSubunits")
SUB <- st_as_sf(SUB)
st_crs(SUB)$units

PROJ <- "boot"
B <- 32 # max is 32
BCRs <- u[u < 200] # Canada only
SPP <- colnames(yy)

#lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
sf <- sf::st_as_sf(dd[,c("PKEY","X","Y")], coords = c("X","Y"))
sf <- st_set_crs(sf, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sf <- st_transform(sf, st_crs(SUB))

## calculate 600 km bands
dd3 <- dd[,"PKEY",drop=FALSE]
for (bcr in BCRs) {
    cat("\n\n", bcr, "\n")
    flush.console()
    tmp <- rep(999, nrow(dd))
    names(tmp) <- rownames(dd)
    p <- SUB[SUB$BCR == bcr,]
    p1 <- st_buffer(p, 1*100*10^3)
    p2 <- st_buffer(p, 2*100*10^3)
    p3 <- st_buffer(p, 3*100*10^3)
    p4 <- st_buffer(p, 4*100*10^3)
    p5 <- st_buffer(p, 5*100*10^3)
    p6 <- st_buffer(p, 6*100*10^3)
    o <- st_intersection(sf, p6)
    tmp[rownames(o)] <- 6
    o <- st_intersection(o, p5)
    tmp[rownames(o)] <- 5
    o <- st_intersection(o, p4)
    tmp[rownames(o)] <- 4
    o <- st_intersection(o, p3)
    tmp[rownames(o)] <- 3
    o <- st_intersection(o, p2)
    tmp[rownames(o)] <- 2
    o <- st_intersection(o, p1)
    tmp[rownames(o)] <- 1
    o <- st_intersection(o, p)
    tmp[rownames(o)] <- 0
    dd3[[paste0("buf_", bcr)]] <- tmp
    print(table(tmp))
}

save(dd3, file=file.path(ROOT, "data", "bcr-buffers-0-600m.RData"))

## now make predictions

library(mefa4)
library(gbm)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2020-01-08.RData"))
#load(file.path(ROOT, "data", "bcr-buffers-0-600m.RData"))

interp <- function(x, value, v1="${", v2="}") {
    x <- as.character(x)
    s <- strsplit(x, v1, fixed=TRUE)[[1L]]
    k <- unique(sapply(strsplit(s[grep(v2, s, fixed=TRUE)], v2, fixed=TRUE), "[[", 1L))
    km <- paste0(v1, k, v2)
    if (!all(k %in% names(l)))
        stop("missing values")
    for (i in seq_along(k)) {
        x <- gsub(km[i], as.character(l[[k[i]]]), x, fixed=TRUE)
    }
    x
}

pa <- "s:/Peter/bam/BAM_data_v2019/gnm/out/boot/${spp}/${bcr}/gnmboot-${spp}-${bcr}-${run}.RData"
#ss <- rep(TRUE, nrow(dd))
ss <- !duplicated(dd$cyid)

SPP <- as.character(
    jsonlite::fromJSON(
        "https://borealbirds.github.io/api/v4/species/")$id)
#SPP <- colnames(yy)

SPP <- SPP[1:35]
SPP <- SPP[36:70]
SPP <- SPP[71:105]
SPP <- SPP[106:143]


#spp <- "OVEN"
for (spp in SPP) {
    ## get bcrs where spp occurs
    sppbcr <- SPPBCR[startsWith(SPPBCR, spp)]
    bcrs <- sapply(strsplit(sppbcr, "-"), "[[", 2L)
    uu <- as.integer(sapply(strsplit(bcrs, "_"), "[[", 2L))
    bcrs <- sort(bcrs[uu < 200]) # no US, Canada only
    ## assemble data
    DAT <- data.frame(
        #count=as.numeric(yy[ss, spp]),
        #offset=off[ss, spp],
        YEAR=dd$YEAR[ss],
        ARU=dd$ARU[ss],
        dd2[ss,])
    ## this stores E[Y]=lambda=DC values
    OUT <- matrix(NA, nrow(DAT), length(bcrs))
    rownames(OUT) <- rownames(DAT)
    colnames(OUT) <- bcrs
    attr(OUT, "spp") <- spp
    ## this stores intercept=log(D) values
    INIT <- matrix(NA, length(bcrs), 32)
    rownames(INIT) <- bcrs
    attr(INIT, "spp") <- spp


    #bcr <- "BCR_60"
    for (bcr in bcrs) {
        PR <- matrix(0, nrow(DAT), 32)
        #run <- 1
        for (run in 1:32) {
            l <- list(spp=spp, bcr=bcr, run=run)
            cat(paste(names(l), l, sep=": ", collapse=", "), "\n")
            flush.console()
            fn <- interp(pa, l)
            load(fn) # object called out
            ## PR stays 0 if error or null
            if (!inherits(out, "try-error") && !is.null(out)) {
                INIT[bcr, run] <- out$initF
                pr <- suppressWarnings(
                    predict.gbm(out, DAT, type="link", n.trees=out$n.trees))
                PR[,run] <- exp(pr + off[ss,spp])
            }
        }
        OUT[,bcr] <- rowMeans(PR)
        gc()
    }
    save(OUT, INIT, file=file.path(ROOT, "out", "validation", paste0("validation-", spp, ".RData")))
}

## get OCCC for internal points

library(mefa4)
library(gbm)
library(epiR)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2020-01-08.RData"))
load(file.path(ROOT, "data", "bcr-buffers-0-600m.RData"))

interp <- function(x, value, v1="${", v2="}") {
    x <- as.character(x)
    s <- strsplit(x, v1, fixed=TRUE)[[1L]]
    k <- unique(sapply(strsplit(s[grep(v2, s, fixed=TRUE)], v2, fixed=TRUE), "[[", 1L))
    km <- paste0(v1, k, v2)
    if (!all(k %in% names(l)))
        stop("missing values")
    for (i in seq_along(k)) {
        x <- gsub(km[i], as.character(l[[k[i]]]), x, fixed=TRUE)
    }
    x
}

pa <- "s:/Peter/bam/BAM_data_v2019/gnm/out/boot/${spp}/${bcr}/gnmboot-${spp}-${bcr}-${run}.RData"

SPP <- as.character(
    jsonlite::fromJSON(
        "https://borealbirds.github.io/api/v4/species/")$id)

#MAX <- 10^3
#B <- 5
MAX <- Inf
B <- 32
OCCC <- NULL
#spp <- "OVEN"
for (spp in SPP) {
    ## get bcrs where spp occurs
    sppbcr <- SPPBCR[startsWith(SPPBCR, spp)]
    bcrs <- sapply(strsplit(sppbcr, "-"), "[[", 2L)
    uu <- as.integer(sapply(strsplit(bcrs, "_"), "[[", 2L))
    bcrs <- sort(bcrs[uu < 200]) # no US, Canada only

    #bcr <- "BCR_60"
    for (bcr in bcrs) {
        ss <- which(!duplicated(dd$cyid) & dd3[,sub("BCR", "buf", bcr)] <= 1)
        if (length(ss) > MAX)
            ss <- sample(ss, MAX)
        DAT <- data.frame(
            #count=as.numeric(yy[ss, spp]),
            #offset=off[ss, spp],
            YEAR=dd$YEAR[ss],
            ARU=dd$ARU[ss],
            dd2[ss,])
        PR <- NULL
        #run <- 1
        for (run in 1:B) {
            l <- list(spp=spp, bcr=bcr, run=run)
            cat(paste(names(l), l, sep=": ", collapse=", "), "\n")
            flush.console()
            fn <- interp(pa, l)
            load(fn) # object called out
            ## PR stays 0 if error or null
            if (!inherits(out, "try-error") && !is.null(out)) {
                pr <- suppressWarnings(
                    predict.gbm(out, DAT, type="link", n.trees=out$n.trees))
                PR <- cbind(PR, exp(pr + off[ss,spp]))
            }
        }
        if (!is.null(PR)) {
            PR <- data.matrix(PR)
            oc <- epi.occc(PR)
            oc <-  data.frame(spp=spp,
                bcr=bcr,
                nrun=ncol(PR),
                as.data.frame(oc[1:3]))
            OCCC <- rbind(OCCC, oc)
        }
        gc()
    }
}
write.csv(OCCC, row.names = FALSE, file="www/gnm-validation-occc.csv")


## compare 250m pred in BCR 60 to c4i

library(raster)
library(cure4insect)

load_common_data()
tab <- read.csv("~/repos/abmispecies/_data/birds.csv")
rownames(tab) <- tab$AOU
fl <- list.files("~/tmp/250m/")
tmp <- strsplit(fl, "-")
SPP <- sapply(tmp, "[[", 2)
SPP <- intersect(SPP, as.character(tab$AOU))
tab <- tab[SPP,]
spc4i <- get_all_species("birds")
tab <- tab[tab$sppid %in% spc4i,]
SPP <- rownames(tab)

res <- list()
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()

    y <- load_species_data(as.character(tab[spp, "sppid"]))
    r <- raster(paste0("~/tmp/250m/pred250-", spp, "-BCR_60-boot-1.tif"))

    rc <- rasterize_results(y)
    rc <- rc[["NC"]]
    r <- try(projectRaster(r, rc))
    if (!inherits(r, "try-error")) {
        rc <- mask(rc, r)
        r <- mask(r, rc)

        vc <- data.frame(c4i=values(rc), gnm=values(r))
        vc <- vc[rowSums(is.na(vc))==0,]

        if (nrow(vc)) {
            res[[spp]] <- c(N = colSums(vc) * 100 * 2/10^6, # M inds total
                LM = coef(lm(gnm ~ c4i, vc)),
                r = cor.test(vc$c4i, vc$gnm)$estimate)
        }
    }
}

res <- do.call(rbind, res)
res <- data.frame(Species=rownames(res), res)
write.csv(res, row.names = FALSE, file="www/gnm-vs-c4i-bcr6ab.csv")


hist(res$r.cor)

op <- par(mfrow=c(1,2), mar=c(4,4,3,3))
plot(N.gnm ~ N.c4i, res)
abline(0,1)
plot(log(N.gnm) ~ log(N.c4i), res)
abline(0,1)
par(op)


## calculate validation metrics in distance bands

library(mefa4)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2020-01-08.RData"))
load(file.path(ROOT, "data", "bcr-buffers-0-600m.RData"))
BCRs <- u[u < 200] # Canada only
SPP <- as.character(
    jsonlite::fromJSON(
        "https://borealbirds.github.io/api/v4/species/")$id)

simple_roc <- function(labels, scores){
    Labels <- labels[order(scores, decreasing=TRUE)]
    data.frame(
        TPR=cumsum(Labels)/sum(Labels),
        FPR=cumsum(!Labels)/sum(!Labels),
        Labels=Labels)
}
simple_auc <- function(ROC) {
    ROC$inv_spec <- 1-ROC$FPR
    dx <- diff(ROC$inv_spec)
    sum(dx * ROC$TPR[-1]) / sum(dx)
}
## provide null when offsets are used, p is number of params if known
pseudo_r2 <- function(observed, fitted, null=NULL, p=0) {
    if (is.null(null))
        null <- mean(observed)
    ll0 <- sum(dpois(observed, null, log=TRUE))
    lls <- sum(dpois(observed, observed, log=TRUE))
    llf <- sum(dpois(observed, fitted, log=TRUE))
    n <- length(observed)
    R2 <- 1 - (lls - llf) / (lls - ll0)
    R2adj <- 1 - (1 - R2) * ((n-1) / (n-(p+1)))
    D0 <- -2 * (ll0 - lls)
    DR <- -2 * (llf - lls)
    p_value <- 1 - pchisq(DR, length(observed)-(p+1))
    c(R2=R2, R2adj=R2adj, Deviance=D0 - DR, Dev0=D0, DevR=DR, p_value=p_value)
}

#spp <- "OVEN"
VAL0 <- matrix(0, 7, 4)
dimnames(VAL0) <- list(
    c("inside", paste0(1:5*100, "-", 2:6*100, "m"), "outside"),
    c("prevalence", "AUC_init", "AUC_final", "pseudo_R2"))

RES <- NULL
for (spp in SPP) {

    load(file.path(ROOT, "out", "validation", paste0("validation-", spp, ".RData"))) # OUT, INIT

    rn <- rownames(OUT)
    y <- as.numeric(yy[rn, spp]) # counts
    o <- off[rn,spp]             # offsets
    uu <- as.integer(sapply(strsplit(colnames(OUT), "_"), "[[", 2L))

    #bcr <- 60
    for (bcr in uu) {
        cat(spp, bcr, "\n")
        flush.console()

        Bcr <- paste0("BCR_", bcr)
        v <- dd3[rn, paste0("buf_", bcr)] # buffer distances from subunit (100m's)
        v[v==0] <- 1 # inside is treated as 1
        v[v>6] <- 7 # far away

        ## bootstrap averaged prediction
        yhat <- OUT[rn,Bcr]
        ## null model (constant D & offsets)
        yhat0 <- exp(matrix(INIT[Bcr,], nrow(OUT), 32, byrow=TRUE) + o)
        yhat0[is.na(yhat0)] <- 0 # error & NULL --> set to 0 density
        yhat0 <- rowMeans(yhat0)

        VAL <- VAL0
        for (i in seq_len(nrow(VAL))) {
            s <- v==i
            VAL[i, "prevalence"] <- mean(ifelse(y[s]>0, 1, 0))
            VAL[i, "AUC_final"] <- simple_auc(simple_roc(ifelse(y[s]>0, 1, 0), yhat[s]))
            VAL[i, "AUC_init"] <- simple_auc(simple_roc(ifelse(y[s]>0, 1, 0), yhat0[s]))
            VAL[i, "pseudo_R2"] <- pseudo_r2(y[s], yhat[s], yhat0[s])[1]
        }

        RES <- rbind(RES, data.frame(spp=spp, bcr=bcr, band=rownames(VAL), dist=1:7, VAL))
    }
}

write.csv(RES, row.names = FALSE, file="www/gnm-validation.csv")


summary(RES)
RES$band <- factor(as.character(RES$band), c("inside", "100-200m", "200-300m", "300-400m",
    "400-500m", "500-600m", "outside" ))
RES$dAUC <- RES$AUC_final - RES$AUC_init

RES$pocc <- 0
for (spp in SPP) {
    b <- unique()
}
RES$sppbcr <- paste0(RES$spp, "_", RES$bcr)

vm <- RES[RES$prevalence > 0,]
table(vm$band, vm$AUC_final <= 0.5)
table(vm$band, vm$pseudo_R2 <= 0)

library(ggplot2)

p <- ggplot(vm, aes(x=band, y=AUC_final)) +
    geom_boxplot() +
    theme_minimal()

p <- ggplot(vm, aes(x=band, y=pseudo_R2)) +
    geom_boxplot() +
    theme_minimal()

vi <- RES[RES$band=="inside" & RES$prevalence > 0, ]
vi$issue <- ifelse(vi$AUC_final <= 0.5 | vi$dAUC <= 0 | vi$pseudo_R2 <= 0, 1, 0)

plot(vi$AUC_final, vi$pseudo_R2, xlim=c(0.5, 1), ylim=c(0,1))

vi2 <- droplevels(vi[vi$AUC_final <= 0.5 | vi$dAUC <= 0 | vi$pseudo_R2 <= 0,])
with(vi2, table(spp, bcr))

## problems are really associated with low sample size
m <- glm(issue ~ prevalence, data=vi, family=binomial)
plot(fitted(m) ~ prevalence, vi)

table(cut(vi$prevalence, 10), vi$issue)

## Canada wide validation metrics

library(mefa4)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2020-01-08.RData"))
load(file.path(ROOT, "data", "bcr-buffers-0-600m.RData"))
BCRs <- u[u < 200] # Canada only
SPP <- as.character(
    jsonlite::fromJSON(
        "https://borealbirds.github.io/api/v4/species/")$id)

simple_roc <- function(labels, scores){
    Labels <- labels[order(scores, decreasing=TRUE)]
    data.frame(
        TPR=cumsum(Labels)/sum(Labels),
        FPR=cumsum(!Labels)/sum(!Labels),
        Labels=Labels)
}
simple_auc <- function(ROC) {
    ROC$inv_spec <- 1-ROC$FPR
    dx <- diff(ROC$inv_spec)
    sum(dx * ROC$TPR[-1]) / sum(dx)
}
## provide null when offsets are used, p is number of params if known
pseudo_r2 <- function(observed, fitted, null=NULL, p=0) {
    if (is.null(null))
        null <- mean(observed)
    ll0 <- sum(dpois(observed, null, log=TRUE))
    lls <- sum(dpois(observed, observed, log=TRUE))
    llf <- sum(dpois(observed, fitted, log=TRUE))
    n <- length(observed)
    R2 <- 1 - (lls - llf) / (lls - ll0)
    R2adj <- 1 - (1 - R2) * ((n-1) / (n-(p+1)))
    D0 <- -2 * (ll0 - lls)
    DR <- -2 * (llf - lls)
    p_value <- 1 - pchisq(DR, length(observed)-(p+1))
    c(R2=R2, R2adj=R2adj, Deviance=D0 - DR, Dev0=D0, DevR=DR, p_value=p_value)
}


#spp <- "OVEN"

RES <- NULL
dd4 <- as.matrix(dd3[,-1])
cu <- as.integer(sapply(strsplit(colnames(dd4), "_"), "[[", 2L))
dd4 <- ifelse(dd4==0, 1, 0)
un <- colSums(t(dd4) * cu)
un <- un[un != 0]

for (spp in SPP) {

    cat(spp, "\n")
    flush.console()
    load(file.path(ROOT, "out", "validation", paste0("validation-", spp, ".RData"))) # OUT, INIT

    uu <- as.integer(sapply(strsplit(colnames(OUT), "_"), "[[", 2L))
    unn <- un[un %in% uu]
    unn <- unn[names(unn) %in% rownames(OUT)]
    y <- as.numeric(yy[names(unn), spp]) # counts
    o <- off[names(unn),spp]             # offsets
    yhat <- yhat0 <- unn * 0
    for (i in uu) {
        nn <- names(unn)[unn ==i]
        yhat[nn] <- OUT[nn, paste0("BCR_", i)]
        m <- exp(matrix(INIT[paste0("BCR_", i),], length(nn), 32, byrow=TRUE) + o[nn])
        m[is.na(m)] <- 0
        yhat0[nn] <- rowMeans(m)
    }

    RES <- rbind(RES, data.frame(spp=spp, region="Canada",
        prevalence=mean(ifelse(y>0, 1, 0)),
        AUC_final=simple_auc(simple_roc(ifelse(y>0, 1, 0), yhat)),
        AUC_init=simple_auc(simple_roc(ifelse(y>0, 1, 0), yhat0)),
        pseudo_R2=pseudo_r2(y, yhat, yhat0)[1]))
}
summary(RES)

write.csv(RES, row.names = FALSE, file="www/gnm-validation-canada.csv")


## calculate validation metrics in across subregions

library(mefa4)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2020-01-08.RData"))
load(file.path(ROOT, "data", "bcr-buffers-0-600m.RData"))
BCRs <- u[u < 200] # Canada only
SPP <- as.character(
    jsonlite::fromJSON(
        "https://borealbirds.github.io/api/v4/species/")$id)

simple_roc <- function(labels, scores){
    Labels <- labels[order(scores, decreasing=TRUE)]
    data.frame(
        TPR=cumsum(Labels)/sum(Labels),
        FPR=cumsum(!Labels)/sum(!Labels),
        Labels=Labels)
}
simple_auc <- function(ROC) {
    ROC$inv_spec <- 1-ROC$FPR
    dx <- diff(ROC$inv_spec)
    sum(dx * ROC$TPR[-1]) / sum(dx)
}
## provide null when offsets are used, p is number of params if known
pseudo_r2 <- function(observed, fitted, null=NULL, p=0) {
    if (is.null(null))
        null <- mean(observed)
    ll0 <- sum(dpois(observed, null, log=TRUE))
    lls <- sum(dpois(observed, observed, log=TRUE))
    llf <- sum(dpois(observed, fitted, log=TRUE))
    n <- length(observed)
    R2 <- 1 - (lls - llf) / (lls - ll0)
    R2adj <- 1 - (1 - R2) * ((n-1) / (n-(p+1)))
    D0 <- -2 * (ll0 - lls)
    DR <- -2 * (llf - lls)
    p_value <- 1 - pchisq(DR, length(observed)-(p+1))
    c(R2=R2, R2adj=R2adj, Deviance=D0 - DR, Dev0=D0, DevR=DR, p_value=p_value)
}


#spp <- "OVEN"
VAL0 <- matrix(0, 7, 4)
dimnames(VAL0) <- list(
    c("inside", paste0(1:5*100, "-", 2:6*100, "m"), "outside"),
    c("prevalence", "AUC_init", "AUC_final", "pseudo_R2"))

RES <- NULL
for (spp in SPP) {

    load(file.path(ROOT, "out", "validation", paste0("validation-", spp, ".RData"))) # OUT, INIT

    dd4 <- ifelse(as.matrix(dd3[,-1]) == 0, 1, 0)
    dd4 <- dd4[rowSums(dd4)>0,]
    cv <- as.integer(gsub("buf_", "", colnames(dd4)))
    dd4 <- rowSums(t(t(dd4) * cv))
    rn <- intersect(names(dd4), rownames(OUT))
    dd4 <- dd4[rn]
    uu <- as.integer(sapply(strsplit(colnames(OUT), "_"), "[[", 2L))

    #bcr <- 60
    for (bcr in uu) {
        cat(spp, bcr, "\n")
        flush.console()

        Bcr <- paste0("BCR_", bcr)

        for (bcr2 in unique(dd4)) {
            v <- names(dd4)[dd4 == bcr2]
            y <- as.numeric(yy[v, spp]) # counts
            o <- off[v,spp]             # offsets

            ## bootstrap averaged prediction
            yhat <- OUT[v,Bcr]
            ## null model (constant D & offsets)
            yhat0 <- exp(matrix(INIT[Bcr,], length(y), 32, byrow=TRUE) + o)
            yhat0[is.na(yhat0)] <- 0 # error & NULL --> set to 0 density
            yhat0 <- rowMeans(yhat0)
            RES <- rbind(RES, data.frame(spp=spp, bcr_train=bcr, bcr_test=bcr2,
                prevalence=mean(ifelse(y>0, 1, 0)),
                AUC_final=simple_auc(simple_roc(ifelse(y>0, 1, 0), yhat)),
                AUC_init=simple_auc(simple_roc(ifelse(y>0, 1, 0), yhat0)),
                pseudo_R2=pseudo_r2(y, yhat, yhat0)[1]))
        }
    }
}

write.csv(RES, row.names = FALSE, file="www/gnm-validation-cross.csv")

## find the internal & best values for each metric
RES$spp_bcr <- paste0(RES$spp, "_", RES$bcr_train)
RES$spp_bcr_test <- paste0(RES$spp, "_", RES$bcr_test)
## this has the internal validation
res <- RES[RES$bcr_train == RES$bcr_test,]
rownames(res) <- res$spp_bcr
res$bcr_test <- NULL
res$spp_bcr_test <- NULL
## now we find the best external
a <- Xtab(AUC_final ~ spp_bcr_test + bcr_test, RES)
a <- find_max(a)[rownames(res),]
r <- Xtab(pseudo_R2 ~ spp_bcr_test + bcr_test, RES)
r <- find_max(r)[rownames(res),]
res <- data.frame(res, AUC=a, R2=r)

plot(res$AUC.value, res$AUC_final)
abline(0,1)
plot(res$R2.value, res$pseudo_R2)
abline(0,1)

table(res$bcr_train == res$AUC.index)
table(res$bcr_train == res$R2.index)

## get n of trees and errors

interp <- function(x, value, v1="${", v2="}") {
    x <- as.character(x)
    s <- strsplit(x, v1, fixed=TRUE)[[1L]]
    k <- unique(sapply(strsplit(s[grep(v2, s, fixed=TRUE)], v2, fixed=TRUE), "[[", 1L))
    km <- paste0(v1, k, v2)
    if (!all(k %in% names(l)))
        stop("missing values")
    for (i in seq_along(k)) {
        x <- gsub(km[i], as.character(l[[k[i]]]), x, fixed=TRUE)
    }
    x
}
pa <- "s:/Peter/bam/BAM_data_v2019/gnm/out/boot/${spp}/${bcr}/gnmboot-${spp}-${bcr}-${run}.RData"
TREES <- NULL # -1=NULL, 0=ERROR, >0=GBM tree size
#spp <- "OVEN"
for (spp in SPP) {
    ## get bcrs where spp occurs
    sppbcr <- SPPBCR[startsWith(SPPBCR, spp)]
    bcrs <- sapply(strsplit(sppbcr, "-"), "[[", 2L)
    uu <- as.integer(sapply(strsplit(bcrs, "_"), "[[", 2L))
    bcrs <- sort(bcrs[uu < 200]) # no US, Canada only
    #bcr <- "BCR_60"
    for (bcr in bcrs) {
        #run <- 1
        for (run in 1:32) {
            l <- list(spp=spp, bcr=bcr, run=run)
            cat(paste(names(l), l, sep=": ", collapse=", "), "\n")
            flush.console()
            fn <- interp(pa, l)
            load(fn) # object called out
            ## PR stays 0 if error or null
            if (is.null(out)) {
                v <- data.frame(spp=spp, bcr=bcr, run=run, n=-1)
            } else {
                if (inherits(out, "try-error")) {
                    v <- data.frame(spp=spp, bcr=bcr, run=run, n=0)
                } else {
                    v <- data.frame(spp=spp, bcr=bcr, run=run, n=out$n.trees)
                }
            }
            TREES <- rbind(TREES, v)
        }
    }
}
write.csv(TREES, row.names = FALSE, file="www/gnm-validation-trees.csv")

## variable importance

library(mefa4)
library(gbm)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2020-01-08.RData"))
SPP <- as.character(
    jsonlite::fromJSON(
        "https://borealbirds.github.io/api/v4/species/")$id)

interp <- function(x, value, v1="${", v2="}") {
    x <- as.character(x)
    s <- strsplit(x, v1, fixed=TRUE)[[1L]]
    k <- unique(sapply(strsplit(s[grep(v2, s, fixed=TRUE)], v2, fixed=TRUE), "[[", 1L))
    km <- paste0(v1, k, v2)
    if (!all(k %in% names(l)))
        stop("missing values")
    for (i in seq_along(k)) {
        x <- gsub(km[i], as.character(l[[k[i]]]), x, fixed=TRUE)
    }
    x
}
vl <- as.character(read.csv("~/repos/api/docs/v4/BAMv4-variables-2020-02-20.csv")$variable)
pa <- "s:/Peter/bam/BAM_data_v2019/gnm/out/boot/${spp}/${bcr}/gnmboot-${spp}-${bcr}-${run}.RData"

SPP <- SPP[1:35]
SPP <- SPP[36:70]
SPP <- SPP[71:105]
SPP <- SPP[106:143]

#spp <- "OVEN"
for (spp in SPP) {
    ## get bcrs where spp occurs
    sppbcr <- SPPBCR[startsWith(SPPBCR, spp)]
    bcrs <- sapply(strsplit(sppbcr, "-"), "[[", 2L)
    uu <- as.integer(sapply(strsplit(bcrs, "_"), "[[", 2L))
    bcrs <- sort(bcrs[uu < 200]) # no US, Canada only
    #bcr <- "BCR_60"
    BIM <- matrix(NA, length(vl), length(bcrs))
    dimnames(BIM) <- list(vl, bcrs)
    for (bcr in bcrs) {
        #run <- 1
        IM <- NULL
#        IML <- list()
        for (run in 1:32) {
            l <- list(spp=spp, bcr=bcr, run=run)
            cat(paste(names(l), l, sep=": ", collapse=", "), "\n")
            flush.console()
            fn <- interp(pa, l)
            load(fn) # object called out
            if (!is.null(out) && !inherits(out, "try-error") && out$n.trees > 0) {
                im <- summary(out, plotit=FALSE)
                im$rel.inf <- im$rel.inf/sum(im$rel.inf)
                m <- im$rel.inf[match(vl, as.character(im$var))]
                names(m) <- vl
                m[is.na(m)] <- 0
                IM <- cbind(IM, m)
#                IML[[run]] <- structure(im$rel.inf, names=as.character(im$var))
            }
        }
        if (!is.null(IM)) {
            IM <- rowMeans(IM)
            BIM[names(IM), bcr] <- IM
        }
    }
    save(BIM, file=file.path(ROOT, "out", "importance", paste0("importance-", spp, ".RData")))
}


## putting all together and updating API


TAB <- jsonlite::fromJSON("https://borealbirds.github.io/api/v4/species/")
TREES <- read.csv("www/gnm-validation-trees.csv")
CROSS <- read.csv("www/gnm-validation-cross.csv")
OCCC <- read.csv("www/gnm-validation-occc.csv")
CAN <- read.csv("www/gnm-validation-canada.csv")
RES <- read.csv("www/gnm-validation.csv")


fr <- TAB$french
fr <- gsub("&eacute;", "é", fr)
fr <- gsub("&Eacute;", "É", fr)
fr <- gsub("&egrave;", "è", fr)
fr <- gsub("&agrave;", "à", fr)
fr <- gsub("&ecirc;", "ê", fr)
fr <- gsub("&acirc;", "â", fr)
fr <- gsub("&ocirc;", "ô", fr)
TAB$french <- fr
TAB$idnext <- TAB$idprevious <- TAB$show <- NULL

write.csv(TAB, row.names=FALSE, file="~/repos/api/docs/v4/BAMv4-species-2020-02-20.csv")

V <- RES[RES$dist == 1,]
V$dist <- V$band <- NULL
rownames(V) <- paste0(V$spp, "-BCR_", V$bcr)
BCRs0 <- c(
    "3 Arctic Plains & Mountains"=3,
    "4 Northwestern Interior Forest"=4,
    "5 Northern Pacific Rainforest"=5,
    "6-0 Boreal Taiga Plains, South"=60,
    "6-1 Boreal Taiga Plains, North"=61,
    "7-0 Taiga Shield & Hudson Plains, West"=70,
    "7-1 Taiga Shield & Hudson Plains, East"=71,
    "8-0 Boreal Softwood Shield, West"=80,
    "8-1 Boreal Softwood Shield, Ontario"=81,
    "8-2 Boreal Softwood Shield, East"=82,
    "8-3 Boreal Softwood Shield, Newfoundland"=83,
    "9 Great Basin"=9,
    "10 Northern Rockies"=10,
    "11 Prairie Potholes"=11,
    "12 Boreal Hardwood Transition"=12,
    "13 Lower Great Lakes/St. Lawrence Plain"=13,
    "14 Atlantic Northern Forest"=14)
V$region <- names(BCRs0)[match(V$bcr, BCRs0)]
CAN$bcr <- 0
CAN <- CAN[,colnames(V)]
rownames(CAN) <- CAN$spp
V <- rbind(CAN, V)

rownames(OCCC) <- paste0(OCCC$spp, "-", OCCC$bcr)
oa <- aggregate(OCCC[,c("occc", "oprec", "oaccu")], list(spp=OCCC$spp), mean)
rownames(oa) <- oa$spp
OCCC <- rbind(oa[,c("occc", "oprec", "oaccu")], OCCC[,c("occc", "oprec", "oaccu")])
V <- data.frame(V, OCCC[match(rownames(V), rownames(OCCC)),])

TREES$OK <- ifelse(TREES$n > 0, 1, 0)
TREES$sppbcr <- paste0(TREES$spp, "-", TREES$bcr)
x <- as.matrix(mefa4::Xtab(OK ~ sppbcr + run, TREES))
x <- rowSums(x)
V$run_complete <- NA
V[names(x), "run_complete"] <- x

V$english <- TAB$english[match(V$spp, TAB$id)]
V$scientific <- TAB$scientific[match(V$spp, TAB$id)]
V$id <- V$spp


boxplot(run_complete ~ is.na(AUC_final), V)
boxplot(prevalence ~ run_complete, V)
boxplot(AUC_final ~ run_complete, V)
boxplot(pseudo_R2 ~ run_complete, V)

cn <- c("id", "scientific", "english", "region",
    "prevalence", "run_complete",
    "AUC_init", "AUC_final", "pseudo_R2",
    "occc", "oprec", "oaccu")
V <- V[,cn]
write.csv(V, row.names=FALSE, file="~/repos/api/docs/v4/BAMv4-validation-2020-02-20.csv")


