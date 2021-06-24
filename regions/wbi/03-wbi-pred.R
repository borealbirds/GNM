library(mefa4)
library(dismo)
library(gbm)
library(qs)
library(raster)

qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-11.qRData")
qload("d:/bam/2021/wbi/pred.qRData")
load("~/repos/GNM/regions/wbi/subsets4.RData")

SU <- list("WBAB"=60,
    "WBBC"=c(4,60),
    "WBYT"=4,
    "WBMB"=c(60,70,80),
    "WBSK"=c(60,80),
    "WBNT"=c(61, 70))
RID <- c(WBAB=1, WBBC=2, WBMB=3, WBNT=4, WBSK=5, WBYT=6)

## !!! Important !!!
## - make sure NALC is factor
## - make sure ARU is all 0
## - make sure YEAR is set to something that represent the training data well: 2011
prd$nalc <- as.factor(prd$nalc)
prd <- prd[rowSums(is.na(prd))==0,]

r <- raster(paste0("d:/bam/2021/wbi/mosaics/AHM.tif"))
values(r)[!is.na(values(r))] <- 0

print.wb_fit <- function(x, ...) {
    cat("WB fit:", x$spp, "/", x$reg, "/", x$i, "\n")
    print(x$AUC)
    cat("\n")
    invisible(x)
}

#spp <- "CAWA"
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    gc()

    fo <- paste0("d:/bam/2021/wbi/res/", spp, "/", "WB-", spp, "-preds.png")
    if (!file.exists(fo)) {

        ## load overall models and make predictions
        i <- 1
        fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALLRES-", i, ".qRData")
        qload(fn)

        if (!dir.exists(paste0("d:/bam/2021/wbi/res/", spp)))
            dir.create(paste0("d:/bam/2021/wbi/res/", spp))

        m <- RES
        cn <- m$vars
        X <- model.matrix(~.,prd[,m$vars])
        system.time(pr_gbm <- predict(m$gbm, prd, type="response"))
        pr_glm <- exp(drop(X[,names(m$glm)] %*% m$glm))
        pr_ssn <- exp(drop(X[,names(m$ssn)] %*% m$ssn))
        pr_ssa <- exp(drop(X[,names(m$ssa)] %*% m$ssa))
        pr <- data.frame(prd[,c("index", "region", "nalc")],
            gbm=pr_gbm,
            glm=pr_glm,
            ssn1=pr_ssn,
            ssa1=pr_ssa)

        ## clip GNM to extent
        rspp <- raster(paste0("d:/bam/BAM_data_v2019/gnm/artifacts/",
            spp, "/pred-", spp, "-CAN-Mean.tif"))
        rspp <- crop(rspp, r)
        rspp <- mask(rspp, r)
        vrspp <- values(rspp)
        pr$gnm <- vrspp[pr$index]

        ## smoothing based on pred
        m_ssn0 <- lm(pr$gnm ~ nalc, data=prd)
        m_ssa0 <- lm(pr$gnm ~ ., data=prd[,m$vars])
        m_ssn2 <- lm(pr$gbm ~ nalc, data=prd)
        m_ssa2 <- lm(log(pr$gbm) ~ ., data=prd[,m$vars])
        pr$ssn0 <- fitted(m_ssn0)
        pr$ssa0 <- fitted(m_ssa0)
        pr$ssn2 <- fitted(m_ssn2)
        pr$ssa2 <- exp(fitted(m_ssa2))

        W <- c("gnm", "gbm", "glm", "ssn0", "ssn1", "ssn2", "ssa0","ssa1", "ssa2")

        #WHAT <- "gnm"
        png(fo, width=1500, height=1500)
        op <- par(mfrow=c(3,3))
        for (WHAT in W) {
            rs <- r
            u <- pr[[WHAT]]
            q <- quantile(u, 0.9999, na.rm=TRUE)
            u[!is.na(u) & u > q] <- q
            u <- u[match(1:length(values(rs)),pr$index)]
            values(rs) <- u
            summary(u)
            summary(values(rs))
            plot(rs, main=WHAT)
        }
        par(op)
        dev.off()
    }
}

pr$ssn0 <- pr$ssn1 <- pr$ssn2 <- NULL
qsavem(pr, file=paste0("d:/bam/2021/wbi/res/", spp, "/", "WB-", spp, "-preds.qRData"))
cor(pr[,4:9])


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

qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-11.qRData")
#qload("d:/bam/2021/wbi/pred.qRData")
load("~/repos/GNM/regions/wbi/subsets4.RData")

SU <- list("WBAB"=60,
    "WBBC"=c(4,60),
    "WBYT"=4,
    "WBMB"=c(60,70,80),
    "WBSK"=c(60,80),
    "WBNT"=c(61, 70))
RID <- c(WBAB=1, WBBC=2, WBMB=3, WBNT=4, WBSK=5, WBYT=6)
prd$nalc <- as.factor(prd$nalc)
prd <- prd[rowSums(is.na(prd))==0,]

r <- raster(paste0("d:/bam/2021/wbi/mosaics/AHM.tif"))
values(r)[!is.na(values(r))] <- 0

spp <- "OVEN"
## load overall models and make predictions

B <- 20
i <- 1
fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALLRES-", i, ".qRData")
qload(fn)
prd <- cbind(dd, ddvar)[RES$pk,RES$vars]

prd$YEAR <- 2011
prd$ARU <- 0
PR <- matrix(0, nrow(prd), B)
pr_gbm <- predict(RES$gbm, prd, type="response")
PR[,i] <- pr_gbm

for (i in 2:B) {
    cat(i, "\n")
    flush.console()
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALLRES-", i, ".qRData")
    qload(fn)
    pr_gbm <- predict(RES$gbm, prd, type="response")
    PR[,i] <- pr_gbm
}

o <- epi.occc(PR)
o
# http://web1.sph.emory.edu/observeragreement/OCCC.pdf
# precision (degree of variation)
# accuracy  (degree of location or scale shift)
# overall: prec * acc

j <- sample(nrow(prd), 2*10^3)
z <- PR[j,]
rm <- rowMeans(z)
z <- z[order(rm),]
matplot(1:nrow(z), z, type="l", lty=1, col="#00000022")

