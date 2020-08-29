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
library(raster)
library(gbm)
library(sf)
library(rgdal)

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
#ss <- rep(TRUE, nrow(dd))
ss <- !duplicated(dd$cyid)
SPP <- colnames(yy)

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
    OUT <- matrix(0, nrow(DAT), length(bcrs))
    rownames(OUT) <- rownames(DAT)
    colnames(OUT) <- bcrs
    attr(OUT, "spp") <- spp

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
                pr <- suppressWarnings(
                    predict.gbm(out, DAT, type="link", n.trees=out$n.trees))
                PR[,run] <- exp(pr + off[ss,spp])
            }
        }
        OUT[,bcr] <- rowMeans(PR)
        gc()
    }
    save(OUT, file=file.path(ROOT, "out", "validation", paste0("validation-", spp, ".RData")))
}


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



