library(mefa4)
library(gbm)
library(raster)


PROJ <- "boot"
ROOT <- "s:/Peter/bam/BAM_data_v2019/gnm"
SPP <- jsonlite::fromJSON("https://borealbirds.github.io/api/v4/species")$id

load(file.path(ROOT, "data", "BAMdb-GNMsubset-2020-01-08.RData"))

d <- data.frame(ARU=dd$ARU, YEAR=dd$YEAR, dd2)

#b=1
#spp="ALFL"
#BCR=83
for (spp in SPP) {
    PPRES <- list()
    for (BCR in u) {
        d0 <- d1 <- d[dd[,paste0("BCR_", BCR)]==1,]
        d0$ROAD <- 0
        d1$ROAD <- 1

        pp <- matrix(NA, nrow(dd), 32)
        rownames(pp) <- rownames(dd)
        pp0 <- pp[rownames(d0),]
        pp1 <- pp[rownames(d1),]

        for (b in 1:32) {
            cat(spp, BCR, b, "\n")
            flush.console()
            fin <- file.path(ROOT, "out", PROJ, spp, paste0("BCR_", BCR),
                paste0("gnmboot-", spp, "-BCR_", BCR, "-", b, ".RData"))
            e <- new.env()
            load(fin, envir=e)
            if (inherits(e$out, "try-error") || is.null(e$out)) {
                pp[,b] <- 0
                pp0[,b] <- 0
                pp1[,b] <- 0
            } else {
                pp[,b] <- suppressWarnings(predict.gbm(e$out,   d, type="response", n.trees=e$out$n.trees))
                pp0[,b] <- suppressWarnings(predict.gbm(e$out, d0, type="response", n.trees=e$out$n.trees))
                pp1[,b] <- suppressWarnings(predict.gbm(e$out, d1, type="response", n.trees=e$out$n.trees))
            }
        }
        PPRES[[ paste0("BCR_", BCR)]] <- list(p=pp, p0=pp0, p1=pp1)
    }
    save(PPRES, file=paste0("d:/bam/BAM_data_v2019/gnm/predpoints/", spp, "-predpoints.RData"))
}


