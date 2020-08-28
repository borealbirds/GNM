#' ---
#' title: "Predictors for BAM Generalized National Models"
#' author: "Peter Solymos <solymos@ualberta.ca>"
#' date: "`r as.Date(Sys.time())`"
#' output: pdf_document
#' ---
#+ echo=FALSE
knitr::opts_chunk$set(eval=FALSE)
ROOT <- "d:/bam/BAM_data_v2019/gnm"
#par(las=1, scipen=999)
#'
#' Making bundle with predictors and buffered BCR indicators
#'
#' # Preamble
library(mefa)
library(mefa4)
library(intrval)
library(sf)
library(sp)
library(raster)
#' Load data and set up coordinates
load(file.path(ROOT, "data", "BAMdb-patched-2019-06-04.RData"))
lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
sf <- sf::st_as_sf(dd, coords = c("X","Y"))
sf <- st_set_crs(sf, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sf <- st_transform(sf, lcc_crs)
#' Coordinates for intersections
if (FALSE) {
    ss <- nonDuplicated(dd, SS, TRUE)[,c("SS", "PCODE", "X", "Y")]
    ss$lon <- ss$X
    ss$lat <- ss$Y

    ss <- sf::st_as_sf(ss, coords = c("X","Y"))
    ss <- st_set_crs(ss, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    ss <- st_transform(ss, lcc_crs)
    str(ss)
    save(ss,file=file.path(ROOT, "data", "BAMdb-patched-xy.RData"))
}
#' Two time periods
dd$YR2 <- ifelse(dd$YEAR < 2006, 1, 2)
dd$SSYR2 <- interaction(dd$SS, dd$YR2, sep="_", drop=TRUE)
#' Randomly pick one survey from each time interval by SS
dd$o <- seq_len(nrow(dd))
set.seed(1)
dd <- dd[sample(dd$o),]
dd <- dd[!duplicated(dd$SSYR2),]
dd <- dd[order(dd$o),]
dd$o <- NULL
#' Visualize subunits
#library(rgdal)
#suu <- readOGR(file.path(ROOT, "data", "predictors", "subunits"))
#plot(suu)
#' Predictors for each time period
e <- new.env()
load(file.path(ROOT, "data", "predictors", "ss_2001attributes.RData"), envir=e)
dd2001 <- as.data.frame(e$ss)[match(dd$SS, e$ss$SS),]
rownames(dd2001) <- rownames(dd)
rm(e)
dd2001$SS <- dd2001$PCODE <- dd2001$lat <- dd2001$lon <- NULL
dd2001$geometry <- NULL
#' Second period
e <- new.env()
load(file.path(ROOT, "data", "predictors", "ss_2011attributes.RData"), envir=e)
dd2011 <- as.data.frame(e$ss)[match(dd$SS, e$ss$SS),]
rownames(dd2011) <- rownames(dd)
rm(e)
dd2011$SS <- dd2011$PCODE <- dd2011$lat <- dd2011$lon <- NULL
dd2011$geometry <- NULL
#' check consistency
stopifnot(all(colnames(dd2001)==colnames(dd2011)))
#' Combine the two tables depending on survey year
dd2 <- dd2001 # use 2001 version for -2005 data
dd2[dd$YEAR >= 2006,] <- dd2011[dd$YEAR >= 2006,] # use 2011 version for 2006- data
#' Sanity checks
#' Land facets: https://adaptwest.databasin.org/pages/adaptwest-landfacets
#dd2$lf[dd2$lf < 1] <- NA
#dd2$lf <- as.factor(as.integer(dd2$lf))
#table(dd2$lf, useNA="a")
#dd2 <- as.matrix(dd2)
#'
#' Assign BCR subunits to surveys
#bcrsu <- st_read(file.path(ROOT, "data", "predictors", "subunits", "BCRSubunits.shp"))
bcrsu <- st_read(file.path(ROOT, "data", "predictors", "subunits", "BCRunits2.shp"))
st_crs(bcrsu) <- st_crs(sf)
sf <- sf[rownames(dd),]
o <- st_join(sf, bcrsu, join = st_intersects)
o <- o[rownames(dd),]
dd$bcr <- o$BCR
table(dd$bcr, useNA="a")
#' BCR units: 6, 7, 8 are split into multiple parts in Canada (6 -> 60+61 etc)
#' BCRs in US are x100 (4 -> 400)
u <- c(4, 5, 60, 61, 70, 71, 80, 81, 82, 83, 9, 10, 11, 12, 13, 14, # Canada
    200, 400, 500, 1000, 1100, 1200, 1300, 1400)                    # US
#dd$bcrsu <- factor(dd$bcrsu, u)
#table(dd$bcrsu, useNA="a")
#' Identify points within the 100km buffers
#' Intersect US portion (no need to have 2001/2011 stacks b/c there is no 2x veg info)
for (i in u) {
    cat("BCR", i, "\n")
    ri <- raster(file.path(ROOT, "data", "subunits", paste0("bcr", i, "all_1km.gri")), 1)
    oi <- extract(ri, sf)
    oi <- !is.na(oi)
    dd[[paste0("BCR_", i)]] <- ifelse(oi, 1L, 0L)
    # read in a US BCR stack to get the list of variables there
    if (i %% 100 == 0) {
        si <- stack(file.path(ROOT, "data", "subunits", paste0("bcr", i, "all_1km.gri")))
        vi <- extract(si, sf[oi,])
        cnu <- colnames(vi)
        cnu <- cnu[cnu != "bcr"]
        dd2[oi,cnu] <- vi[,cnu]
    }
}
#'
#' Land facets dropped due to issues at places (big blobs show up on certain maps)
dd2$lf <- NULL
#' Check some variables: need to drop RH, PAS, Eref
dd2$RH <- NULL
dd2$PAS <- NULL
dd2$Eref <- NULL
#' Use ROAD layer from `dd2` not BBS or not from `dd`
dd$ROAD <- NULL
#' Dealing with ARU and BBS data
dd$isBBS <- startsWith(rownames(dd), "BBSAB")
table(dd2[,"ROAD"])
dd2[dd$isBBS,"ROAD"] <- 1
dd2$ROAD[is.na(dd2$ROAD)] <- 0
table(dd2[,"ROAD"], useNA="a")
#' ARUs
dd$ARU <- ifelse(startsWith(as.character(dd$PCODE), "BU_"), 1, 0) # BU and WildTrax
table(dd$ARU, useNA="a")
#' 0=no data, 18=water, 19=snow & ice
dd2$nalc[dd2$nalc %)(% c(1, 17)] <- NA
dd2$nalc <- as.factor(dd2$nalc)
table(dd2$nalc, useNA="a")
#' Keep points that are within the buffers and terrain/climate is not NA
m <- dd[,colnames(dd)[startsWith(colnames(dd), "BCR_")]]
keep <- rowSums(m) > 0 & !is.na(dd2$slope) & !is.na(dd2$AHM) & !is.na(dd2$nalc)
keep[is.na(dd2$Structure_Volume_Total_v1) & !is.na(dd2$LandCover_NonVeg_v1.1)] <- FALSE
cnCan <- colnames(dd2001)[colnames(dd2001) %in% colnames(dd2)]
cnUS <- colnames(vi)[colnames(vi) %in% colnames(dd2)]
sum(keep)
for (i in u) {
    cat("------------------------------ BCR", i, "\n")
    if (i %% 100 == 0) {
        tmp <- dd2[dd[[paste0("BCR_", i)]] > 0, cnUS]
        print(table(rowSums(is.na(tmp))))
        keep[rownames(tmp)[rowSums(is.na(tmp)) > 0]] <- FALSE
    } else {
        tmp <- dd2[dd[[paste0("BCR_", i)]] > 0, cnCan]
        print(table(rowSums(is.na(tmp))))
        keep[rownames(tmp)[rowSums(is.na(tmp)) > 0]] <- FALSE
    }
}
sum(keep)
# gap analysis data: 1km^2 res sampling intensity
if (FALSE) {
    xy <- st_coordinates(sf)
    bb <- bbox(xy)
    r <- raster(
        xmn=bb["x", "min"],
        xmx=bb["x", "max"],
        ymn=bb["y", "min"],
        ymx=bb["y", "max"],
        res=1000,
        crs=st_crs(sf)$proj4string)
    sr1k <- rasterize(xy, r, field=1, fun='sum')
    xyz <- cbind(xyFromCell(sr1k, 1:length(sr1k)), z=values(sr1k))
    summary(xyz)
    table(xyz[,3], useNA="a")
    xyz <- xyz[!is.na(xyz[,3]),]
    save(xyz, file=file.path(ROOT, "data", "BAMdb-sampling-intensity-1km-2020-02-05.RData"))
}
#' Subset
ss <- rownames(dd2)[keep]
dd <- dd[ss,]
dd2 <- dd2[ss,]
yy <- yy[ss,]
off <- off[ss,]
sf <- sf[ss,]
#' Check NA patterns
c(sum(is.na(dd)), sum(is.na(dd2)), sum(is.na(off)), sum(is.na(yy)))
plot(table(rowSums(is.na(dd2))))
data.frame(n=colSums(is.na(dd)))
data.frame(n=colSums(is.na(dd2)))
#' Drop species with 0 detections after subsets
any(colSums(yy) == 0)
yy <- yy[,colSums(yy) > 0]
sort(colSums(yy))
#' Evaluate predictor sets based on hist, SD, etc
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
#' `CN` collects column names to be used in BRTs
#' includes factors and ROAD regardless of SD and correlation
CN <- list()
for (i in u) {
    cat("BCR", i, "\n")
    BCR <- paste0("BCR_", i)
    if (i %% 100 == 0) {
        z <- as.matrix(dd2[dd[,BCR] == 1L, cnUS[cnUS != "nalc"]])
        SD <- apply(z, 2, sd)
        CN[[BCR]] <- unique(c("nalc", "ROAD", cnUS))
    } else {
        z <- as.matrix(dd2[dd[,BCR] == 1L, cnCan[cnCan != "nalc"]])
        SD <- apply(z, 2, sd)
        CN[[BCR]] <- unique(c("nalc", "ROAD", get_cn(z[,SD > 0])))
    }
}
sapply(CN, length)
#'
#' Calculating weights
xy <- st_coordinates(sf)
bb <- bbox(xy)
r <- raster(
    xmn=bb["x", "min"],
    xmx=bb["x", "max"],
    ymn=bb["y", "min"],
    ymx=bb["y", "max"],
    res=250,
    crs=st_crs(sf)$proj4string)
sr <- rasterize(xy, r, field=1, fun='sum')
W <- matrix(1,nrow=5,ncol=5)
W[c(1,5,21,25)] <- 0
sr25 <- focal(sr, w=W, na.rm=TRUE)
ni <- extract(sr25, xy)
ni[is.na(ni)] <- 1
wi <- 1/ni
nsub <- ceiling(sum(sapply(sort(unique(ni)), function(z) sum(ni == z)/z)))
dd$ni <- ni
dd$wi <- wi
#' spatial grid IDs
r10 <- aggregate(r, fact=10)
values(r10) <- seq_len(ncell(r10))
dd$cid <- extract(r10, xy)
#' Calculating set of species-subunit combination with >0 sums
fullBCRlist <- names(sort(colSums(dd[,paste0("BCR_", u)])))
b <- as(as.matrix(dd[,paste0("BCR_", u)]), "dgCMatrix")
b <- b[,fullBCRlist]
yy01 <- yy
yy01[yy>0] <- 1
detbcr <- t(sapply(colnames(yy), function(z) colSums(yy01[,z] * b)))
attr(detbcr, "nbcr") <- colSums(b)
SPPBCR <- NULL
for (i in colnames(detbcr)) {
    for (j in rownames(detbcr)) {
        if (detbcr[j,i] > 0)
            SPPBCR <- c(SPPBCR, paste0(j, "-", i))
    }
}
#'
data.frame(n=colSums(is.na(dd)))
data.frame(n=colSums(is.na(dd2)))
#' cluster + year IDs
dd$cyid <- interaction(dd$cid, dd$YEAR, sep="_", drop=TRUE)
#' Save
save(dd, dd2, yy, off, spt, u, CN, nsub, detbcr, SPPBCR, cnCan, cnUS,
    file=file.path(ROOT, "data", "BAMdb-GNMsubset-2020-01-08.RData"))
