library(mefa4)
library(gbm)
library(raster)


PROJ <- "boot"
ROOT <- "s:/Peter/bam/BAM_data_v2019/gnm"
SPP <- jsonlite::fromJSON("https://borealbirds.github.io/api/v4/species")$id

load(file.path(ROOT, "data", "BAMdb-GNMsubset-2020-01-08.RData"))

d <- data.frame(ARU=dd$ARU, YEAR=dd$YEAR, dd2)
bcrs <- c("BCR_4", "BCR_5", "BCR_60", "BCR_61", "BCR_70", "BCR_71", "BCR_80",
    "BCR_81", "BCR_82", "BCR_83", "BCR_9", "BCR_10", "BCR_11", "BCR_12",
    "BCR_13", "BCR_14")
blu <- as.matrix(dd[,bcrs])
keep <- rowSums(blu)>0
d <- d[keep,]
blu <- blu[keep,]

#b=1
#spp="OVEN"
#BCR="BCR_60"
#SPP=SPP[1:50]
#SPP=SPP[51:100]
#SPP=SPP[101:143]


for (spp in SPP) {
    PPRES <- list()
    for (BCR in bcrs) {
        d01 <- d[blu[,BCR]==1,]
        d0 <- d1 <- d01[d01$ROAD==1,]
        d0$ROAD <- 0

        pp <- matrix(NA, nrow(d01), 32)
        rownames(pp) <- rownames(d01)
        pp0 <- pp1 <- NULL
#        pp0 <- pp[rownames(d0),]
#        pp1 <- pp[rownames(d1),]

        for (b in 1:32) {
            cat(spp, BCR, b, "\n")
            flush.console()
            fin <- file.path(ROOT, "out", PROJ, spp, BCR,
                paste0("gnmboot-", spp, "-", BCR, "-", b, ".RData"))
            e <- new.env()
            load(fin, envir=e)
            if (inherits(e$out, "try-error") || is.null(e$out)) {
                pp[,b] <- 0
#                pp0[,b] <- 0
#                pp1[,b] <- 0
            } else {
                pp[,b] <- suppressWarnings(predict.gbm(e$out,
                    d01, type="response", n.trees=e$out$n.trees))
#                pp0[,b] <- suppressWarnings(predict.gbm(e$out,
#                    d0, type="response", n.trees=e$out$n.trees))
#                pp1[,b] <- suppressWarnings(predict.gbm(e$out,
#                    d1, type="response", n.trees=e$out$n.trees))
            }
        }
        PPRES[[BCR]] <- list(p=pp, p0=pp0, p1=pp1)
    }
    save(PPRES, file=paste0("d:/bam/BAM_data_v2019/gnm/predpoints/", spp, "-predpoints.RData"))
}


## calculate smooth prediction
bcrs <- c("BCR_4", "BCR_5", "BCR_60", "BCR_61", "BCR_70", "BCR_71", "BCR_80",
    "BCR_81", "BCR_82", "BCR_83", "BCR_9", "BCR_10", "BCR_11", "BCR_12",
    "BCR_13", "BCR_14")
blu <- as.matrix(dd[,bcrs])
#' Simple and fast ROC and AUC
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

spp="RCKI"
load(paste0("d:/bam/BAM_data_v2019/gnm/predpoints/", spp, "-predpoints.RData"))

y <- yy[,spp]
table(y)
o <- off[,spp]
ppres <- matrix(NA, nrow(PPRES[[1]]$p), length(bcrs))
dimnames(ppres) <- list(rownames(PPRES[[1]]$p), bcrs)
best <- ppres
for (i in 1:length(bcrs)) {
    val <- rowMeans(PPRES[[bcrs[i]]]$p)
    ppres[,i] <- val
    ss <- blu[,bcrs[i]]==1
    best[ss,i] <- val[ss]
}
best <- rowMeans(best, na.rm=TRUE)


lam <- exp(log(ppres) + o)
lamb <- exp(log(best) + o)



pred <- lamb
keep <- !is.na(pred)
ll0 <- sum(dpois(y[keep], mean(y[keep]), log=TRUE))
lls <- sum(dpois(y[keep], y[keep], log=TRUE))
llf <- sum(dpois(y[keep], pred[keep], log=TRUE))
R2 <- 1 - (lls - llf) / (lls - ll0)
R2

roc <- simple_roc(ifelse(y[keep]>0, 1, 0), pred[keep])
auc <- simple_auc(roc)
auc

boxplot(pred ~ y)

R2mat <- matrix(NA, length(bcrs), length(bcrs))
dimnames(R2mat) <- list(paste0("from_", bcrs), paste0("to_", bcrs))
AUCmat <- R2mat

for (j in 1:length(bcrs)) {      # to
    for (i in 1:length(bcrs)) {  # from
        keep <- blu[,bcrs[j]]==1 # to
        pred <- lam[, bcrs[i]]   # from
        ll0 <- sum(dpois(y[keep], mean(y[keep]), log=TRUE))
        lls <- sum(dpois(y[keep], y[keep], log=TRUE))
        llf <- sum(dpois(y[keep], pred[keep], log=TRUE))
        R2 <- 1 - (lls - llf) / (lls - ll0)
        roc <- simple_roc(ifelse(y[keep]>0, 1, 0), pred[keep])
        auc <- simple_auc(roc)
        R2mat[i,j] <- max(0, R2)
        AUCmat[i,j] <- auc
    }
}

R2mat
diag(R2mat)
AUCmat
diag(AUCmat)

## Add here AUC

## calculate road contrast

pp0 <- matrix(NA, nrow(PPRES[[1]]$p), length(bcrs))
dimnames(pp0) <- list(rownames(PPRES[[1]]$p), bcrs)
pp1 <- pp0
for (i in 1:length(bcrs)) {
    v0 <- rowMeans(PPRES[[bcrs[i]]]$p0)
    v1 <- rowMeans(PPRES[[bcrs[i]]]$p1)
    pp0[names(v0), bcrs[i]] <- v0
    pp1[names(v1), bcrs[i]] <- v1
}
pp0 <- rowMeans(pp0, na.rm=TRUE)
pp1 <- rowMeans(pp1, na.rm=TRUE)
table(is.na(pp0), is.na(pp1))
Mean0 <- mean(pp0[!is.na(pp0) & !is.na(pp1)])
Mean1 <- mean(pp1[!is.na(pp0) & !is.na(pp1)])

Mean0 # not huge diffs...
Mean1

mean(y[rowSums(blu)>0 & dd2$ROAD == 1])
mean(y[rowSums(blu)>0 & dd2$ROAD == 0])

str(dd)
library(sf)
load("d:/bam/2021/gnm/regions/BCRPROV.RData") # xy, xy1
# Nad83
ddxy <- st_as_sf(dd, coords = c("X", "Y"), crs = 4269)
ddxy <- st_transform(ddxy, st_crs(xy1))
coord <- st_coordinates(ddxy)

## how to summarize centroids:
## - bcr/jurs: this is not how models were fit
## - subunits: this is not how dete is summarized
##




