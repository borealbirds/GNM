## var importance

library( mefa4)

vi0 <- read.csv("~/repos/api/docs/v4/BAMv4-importance-2020-02-20.csv")
var <- read.csv("~/repos/api/docs/v4/BAMv4-variables-2020-02-20.csv")


non <- c(
    "SpeciesGroups_Unknown_Spp_v1",
    "Structure_Biomass_Branch_v1",
    "Structure_Biomass_Foliage_v1",
    "Structure_Biomass_StemBark_v1",
    "Structure_Biomass_StemWood_v1",
    "Structure_Biomass_TotalDead_v1",
    "Structure_Stand_CrownClosure_v1",
    "Structure_Stand_Height_v1",
    "Structure_Volume_Merch_v1",
    "Structure_Volume_Total_v1",
    "Landsc750_Unknown_Spp_v1",
    "Landsc750_Biomass_Branch_v1",
    "Landsc750_Biomass_Foliage_v1",
    "Landsc750_Biomass_StemBark_v1",
    "Landsc750_Biomass_StemWood_v1",
    "Landsc750_Biomass_TotalDead_v1",
    "Landsc750_Stand_CrownClosure_v1",
    "Landsc750_Stand_Height_v1",
    "Landsc750_Volume_Merch_v1",
    "Landsc750_Volume_Total_v1")
var$non <- var$variable %in% non

u <- c(
    "4 Northwestern Interior Forest",
    "6-0 Boreal Taiga Plains, South",
    "6-1 Boreal Taiga Plains, North",
    "7-0 Taiga Shield & Hudson Plains, West",
    "8-0 Boreal Softwood Shield, West")
U <- c(4, 60, 61, 70, 80)

vi <- droplevels(vi0[vi0$region %in% u & !(vi0$variable %in% non),])



LIMIT <- 0.01
NMAX <- 4
L <- list()
for (spp in levels(vi$id)) {
    for (i in 1:5) {
        reg <- u[i]
        sppreg <- paste0(spp, "-", U[i])
        h <- vi[vi$id == spp & vi$region == reg,]
        h <- h[order(h$importance, decreasing = TRUE),]
        if (nrow(h) > NMAX)
            h <- h[1:NMAX,]
        L[[sppreg]] <- as.character(h$variable)
    }
}
length(L)
str(unique(unlist(L)))
table(sapply(L, length))

d <- data.frame(sppreg=names(L), n=sapply(L, length))
tmp <- strsplit(names(L), "-")
d$spp <- sapply(tmp, "[[", 1)
d$reg <- sapply(tmp, "[[", 2)
table(d$n==0,d$reg)
table(d$n==0,d$spp)

Xtab(ifelse(n==0, 0, 1) ~ spp + reg, d)

cnnew <- sort(unique(unlist(L)))
compare_sets(cnnew, non)

save(L, file="regions/wbi/subsets4.RData")




## BAM v6 data & WT for the WBI study area

library(mefa4)
library(sf)
library(sp)
library(raster)
library(rgdal)

lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
wb <- readOGR("d:/bam/2021/wbi/study-area")
wb <- st_as_sf(wb)
wb <- st_transform(wb, lcc_crs)
levels(wb$HASC_1) <- gsub("CA\\.", "WB", levels(wb$HASC_1))

load("d:/bam/BAM_data_v2019/gnm/data/BAMdb-GNMsubset-2020-01-08.RData")
u <- c(4, 60, 61, 70, 80)
load("d:/bam/2021/rof/BAMv6_ONBBS.RData")
rm(o1, x1, x1off, y1) # ON BBS
load("regions/wbi/subsets.RData")


eee <- new.env()
load("d:/bam/2021/wbi/WT-data-processed-2021-06-23.RData", envir=eee)


#i1 <- dd$bcr %in% c(4, 6, 7, 8)
#dd$BCR_4 > 0 | dd$BCR_60 > 0 | dd$BCR_61 > 0 | dd$BCR_70 > 0 | dd$BCR_80 > 0
pk1 <- rownames(dd)#[i1]

dd1 <- data.frame(
    PKEY=dd$PKEY,
    SS=dd$SS,
    X=dd$X,
    Y=dd$Y,
    YEAR=dd$YEAR,
    ARU=dd$ARU)
rownames(dd1) <- dd1$PKEY
dd1 <- dd1[pk1,]

cnnew <- sort(unique(unlist(L)))

dd1var <- dd2[rownames(dd1), cnnew[!(cnnew %in% c("ARU", "YEAR"))]]

pk2 <- setdiff(rownames(x2), rownames(dd))
pk2 <- intersect(pk2, rownames(y2))
pk2 <- intersect(pk2, rownames(o2))

dd2 <- data.frame(PKEY=rownames(x2),
    SS=x2$SS_SITE2_V4,
    X=x2$longitude,
    Y=x2$latitude,
    YEAR=x2$survey_year,
    ARU=0)
rownames(dd2) <- dd2$PKEY
dd2 <- dd2[pk2,]
dd2 <- dd2[rowSums(is.na(dd2))==0,]


SPP <- sort(unique(substr(names(L), 1, 4)))
Y1 <- yy[pk1, SPP]
colnames(y2)[colnames(y2)=="CAJA"] <- "GRAJ"
Y2 <- y2[pk2, SPP]

O1 <- off[pk1,SPP]


library(mefa4)
library(QPAD)
library(maptools)
library(intrval)
library(raster)

## change this path according your recurring project location
od <- setwd("~/repos/recurring/offset")
load_BAM_QPAD(version = 3)
if (getBAMversion() != "3")
  stop("This script requires BAM version 3")
rlcc <- raster("./data/lcc.tif")
rtree <- raster("./data/tree.tif")
rtz <- raster("./data/utcoffset.tif")
rd1 <- raster("./data/seedgrow.tif")
crs <- proj4string(rtree)
source("functions.R")
setwd(od)

g <- make_off("GRAJ", x2off)$offset
o2 <- cbind(o2, GRAJ=g)

O2 <- o2[pk2,SPP]

table(dd1$ARU)
table(dd2$ARU)
dd1 <- droplevels(dd1[dd1$ARU==0,])
Y1 <- Y1[rownames(dd1),]
O1 <- O1[rownames(dd1),]
dd1var <- dd1var[rownames(dd1),]


#save(dd1, dd2, Y1, Y2, O1, O2, dd1var,
#    file="d:/bam/2021/wbi/WBI-data_BAMv4v6combo_2021-06-10.RData")

## WT data addition

SPPfin <- c("ALFL", "AMCR", "AMGO", "AMRE", "AMRO", "ATSP", "ATTW", "BAOR",
  "BARS", "BAWW", "BBCU", "BBMA", "BBWA", "BBWO", "BCCH", "BEKI",
  "BHCO", "BHVI", "BLBW", "BLJA", "BLPW", "BOBO", "BOCH", "BOWA",
  "BRBL", "BRCR", "BTNW", "CAWA", "CCSP", "CEDW", "CHSP", "CMWA",
  "COGR", "CONW", "CORA", "COYE", "CSWA", "DEJU", "DOWO", "EAKI",
  "EAPH", "EUST", "EVGR", "FOSP", "GCFL", "GCKI", "GCSP", "GCTH",
  "GRAJ", "GRCA", "GRYE", "HAFL", "HAWO", "HETH", "HOLA", "HOSP",
  "HOWR", "KILL", "LALO", "LCSP", "LEFL", "LEYE", "LISP", "MAWA",
  "MAWR", "MODO", "MOWA", "NAWA", "NOFL", "NOWA", "OCWA", "OSFL",
  "OVEN", "PAWA", "PHVI", "PIGR", "PISI", "PIWO", "PUFI", "RBGR",
  "RBNU", "RCKI", "RECR", "REVI", "RUBL", "RUGR", "RWBL", "SAVS",
  "SEWR", "SOSA", "SOSP", "SPSA", "SWSP", "SWTH", "TEWA", "TOSO",
  "TOWA", "TRES", "VATH", "VEER", "VESP", "WAVI", "WBNU", "WCSP",
  "WETA", "WEWP", "WIPT", "WISN", "WIWA", "WIWR", "WTSP", "WWCR",
  "YBFL", "YBSA", "YEWA", "YHBL", "YRWA")

dd3 <- data.frame(PKEY=eee$pk$PKEY,
    SS=eee$pk$SS,
    X=eee$pk$longitude,
    Y=eee$pk$latitude,
    YEAR=eee$pk$YEAR,
    ARU=1)
rownames(dd3) <- dd3$PKEY
dd3 <- dd3[!is.na(dd3$X),]
Y3 <- eee$y[rownames(dd3),SPPfin]
O3 <- eee$OFF[rownames(dd3),SPPfin]


## need to merge poly for the 5 units
## buffer each and see what fall inside

sf <- sf::st_as_sf(dd1, coords = c("X","Y"))
sf <- st_set_crs(sf, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sf <- st_transform(sf, lcc_crs)
sf1 <- sf
rm(sf)

sf <- sf::st_as_sf(dd2, coords = c("X","Y"))
sf <- st_set_crs(sf, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sf <- st_transform(sf, lcc_crs)
sf2 <- sf
rm(sf)

sf <- sf::st_as_sf(dd3, coords = c("X","Y"))
sf <- st_set_crs(sf, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sf <- st_transform(sf, lcc_crs)
sf3 <- sf
rm(sf)

#st_crs(bcrsu) <- st_crs(sf1)
o <- st_join(sf1, wb, join = st_intersects)
o <- o[rownames(dd1),]
dd1$reg <- o$HASC_1
table(dd1$reg, useNA="a")


o <- st_join(sf2, wb, join = st_intersects)
o <- o[rownames(dd2),]
dd2$reg <- o$HASC_1
table(dd2$reg, useNA="a")

o <- st_join(sf3, wb, join = st_intersects)
o <- o[rownames(dd3),]
dd3$reg <- o$HASC_1
table(dd3$reg, useNA="a")

dd1 <- droplevels(dd1[!is.na(dd1$reg),])
dd2 <- droplevels(dd2[!is.na(dd2$reg),])
dd3 <- droplevels(dd3[!is.na(dd3$reg),])


## need to intersect new points with mosaiced WBI layers
ROOT <- "d:/bam/BAM_data_v2019/gnm"

sf <- sf2[rownames(dd2),]
dd2var <- NULL
for (i in 1:222) {
    cat(i, "\n")
    r4 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr4all_1km.grd")), i)
    if (names(r4) %in% colnames(dd1var)) {
        r60 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr60all_1km.grd")), i)
        r61 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr61all_1km.grd")), i)
        r70 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr70all_1km.grd")), i)
        r80 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr80all_1km.grd")), i)
        #rm <- mosaic(r4, r60, r61, r70, r80, fun=mean)
        o4 <- extract(r4, sf)
        o60 <- extract(r60, sf)
        o61 <- extract(r61, sf)
        o70 <- extract(r70, sf)
        o80 <- extract(r80, sf)
        o4[is.na(o4)] <- o60[is.na(o4)]
        o4[is.na(o4)] <- o61[is.na(o4)]
        o4[is.na(o4)] <- o70[is.na(o4)]
        o4[is.na(o4)] <- o80[is.na(o4)]
        o4[is.na(o4)] <- mean(o4, na.rm=TRUE)
        #summary(o4)

        dd2var <- cbind(dd2var, o4)
        colnames(dd2var)[ncol(dd2var)] <- names(r4)
    }
}
dd2var <- data.frame(dd2var)
rownames(dd2var) <- rownames(dd2)
sum(is.na(dd2var))
dd2var <- dd2var[,colnames(dd1var)]
dd2var$nalc <- factor(as.character(dd2var$nalc), levels(dd1var$nalc))


sf <- sf3[rownames(dd3),]
dd3var <- NULL
for (i in 1:222) {
    cat(i, "\n")
    r4 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr4all_1km.grd")), i)
    if (names(r4) %in% colnames(dd1var)) {
        r60 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr60all_1km.grd")), i)
        r61 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr61all_1km.grd")), i)
        r70 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr70all_1km.grd")), i)
        r80 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr80all_1km.grd")), i)
        #rm <- mosaic(r4, r60, r61, r70, r80, fun=mean)
        o4 <- extract(r4, sf)
        o60 <- extract(r60, sf)
        o61 <- extract(r61, sf)
        o70 <- extract(r70, sf)
        o80 <- extract(r80, sf)
        o4[is.na(o4)] <- o60[is.na(o4)]
        o4[is.na(o4)] <- o61[is.na(o4)]
        o4[is.na(o4)] <- o70[is.na(o4)]
        o4[is.na(o4)] <- o80[is.na(o4)]
        o4[is.na(o4)] <- mean(o4, na.rm=TRUE)
        #summary(o4)

        dd3var <- cbind(dd3var, o4)
        colnames(dd3var)[ncol(dd3var)] <- names(r4)
    }
}
dd3var <- data.frame(dd3var)
rownames(dd3var) <- rownames(dd3)
sum(is.na(dd3var))
dd3var <- dd3var[,colnames(dd1var)]
dd3var$nalc <- factor(as.character(dd3var$nalc), levels(dd1var$nalc))

#save(dd1, dd2, Y1, Y2, O1, O2, dd1var, dd2var,
#    file="d:/bam/2021/wbi/WBI-data_BAMv4v6combo_2021-06-10.RData")

dd <- rbind(dd1, dd2, dd3)
ddvar <- rbind(
  dd1var[rownames(dd1),],
  dd2var[rownames(dd2),],
  dd3var[rownames(dd3),])
SPP <- SPPfin
Y <- rbind(
  Y1[rownames(dd1),SPP],
  Y2[rownames(dd2),SPP],
  Y3[rownames(dd3),SPP])
#SPP <- colnames(Y)[colSums(Y>0) >= 200]
#Y <- Y[,SPP]
O <- rbind(
  O1[rownames(dd1),SPP],
  O2[rownames(dd2),SPP],
  O3[rownames(dd3),SPP])

sum(is.na(ddvar))
data.frame(nna=sort(colSums(is.na(ddvar))))
sum(is.na(O))
data.frame(nna=sort(colSums(is.na(O))))
table(rowSums(is.na(ddvar)), rowSums(is.na(O)))

pk <- rownames(dd)[rowSums(is.na(ddvar)) == 0 & rowSums(is.na(O)) == 0]
dd <- dd[pk,]
ddvar <- ddvar[pk,]
O <- O[pk,]
Y <- Y[pk,]

dd <- droplevels(dd)
#SPP <- colnames(Y)[colSums(Y>0) >= 200]
Y <- Y[,SPP]
O <- O[,SPP]


L2 <- L[substr(names(L), 1, 4) %in% SPP]

NN <- data.frame(spp_reg=names(L2),
    spp=substr(names(L2), 1, 4),
    reg=substr(names(L2), 6, nchar(names(L2))))

#system.time(save(dd, Y, O, ddvar, SPP, u, NN, L2,
#    file="d:/bam/2021/wbi/WBI-data_BAMv4v6comboFinal_2021-06-10.RData"))

library(qs)

## write is 10x faster
# read back with qs::qload



#' Calculating weights
library(sf)
library(sp)
library(raster)

lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
sf <- sf::st_as_sf(dd, coords = c("X","Y"))
sf <- st_set_crs(sf, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sf <- st_transform(sf, lcc_crs)

xy <- st_coordinates(sf)
bb <- bbox(xy)
r <- raster(
    xmn=bb["x", "min"],
    xmx=bb["x", "max"],
    ymn=bb["y", "min"],
    ymx=bb["y", "max"],
    res=250,
    crs=st_crs(sf)$proj4string)
r10 <- aggregate(r, fact=10)
values(r10) <- seq_len(ncell(r10))
dd$cid <- extract(r10, xy)
#' cluster + year IDs
dd$cyid <- interaction(dd$cid, dd$YEAR, sep="_", drop=TRUE)



system.time(qsavem(dd, Y, O, ddvar, SPP, u, NN, L2,
    file="d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-25.qRData"))


## some filtering
library(mefa4)
library(qs)

qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-25.qRData")
eee <- new.env()
load("d:/bam/2021/wbi/WT-data-processed-2021-06-23.RData", envir=eee)

pk <- eee$pk
rownames(pk) <- pk$PKEY
compare_sets(rownames(pk), rownames(dd))
z <- pk[intersect(rownames(pk), rownames(dd)),]
OKD <- z$ToD > 2 & z$ToD < 12
OKY <- z$ToY > 125 & z$ToY < 200
OK <- OKD & OKY
NotOK <- rownames(z)[!OK]
str(NotOK)

keep <- !(rownames(dd) %in% NotOK)
table(keep)

dd <- droplevels(dd[keep,])
ddvar <- ddvar[keep,]
Y <- Y[keep,]
O <- O[keep,]

qsavem(dd, Y, O, ddvar, SPP,
    file="d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-30.qRData")

