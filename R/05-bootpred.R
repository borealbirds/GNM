library(mefa4)
library(gbm)
library(raster)
ROOT <- "d:/bam/BAM_data_v2019/gnm"

fn <- "BAMdb-GNMsubset-2020-01-08.RData"
PROJ <- "boot"
load(file.path(ROOT, "data", fn))


## prepare stacks for US
if (FALSE) {
for (i in u[u>2100]) {
    cat("\nLoading stack for BCR", i, "\n")
    flush.console()
    gc()

    pr <- stack(file.path(ROOT, "data", "subunits", paste0("bcr", i, "all_1km.grd")))

    pset <- CN[[paste0("BCR_", i)]]
    print(compare_sets(pset, names(pr)))

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

B <- 1:32
spp <- "CAWA"

for (BCR in u) {
    cat("loading stack: BCR", BCR, "\n")
    flush.console()
    #spp <- "OSFL"
    r1 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr", BCR, "all_1km.grd")))
    ND <- readRDS(file.path(ROOT, paste0("STACK-ND-BCR_", BCR, ".rds")))

    ND$data$YEAR <- 2011

    ## predict for runs
    for (b in B) {
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
    ## combine runs: median & IQR
    st <- list()
    for (b in B) {
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

