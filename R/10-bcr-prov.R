## mosaic nalc

library(raster)
ROOT <- "d:/bam/BAM_data_v2019/gnm"

fn <- "BAMdb-GNMsubset-2020-01-08.RData"
PROJ <- "boot"
load(file.path(ROOT, "data", fn))

for (i in u[u < 200]) {
    cat("\nLoading stack for BCR", i, "\n")
    flush.console()
    gc()
    pr <- stack(file.path(ROOT, "data", "subunits", paste0("bcr", i, "all_1km.grd")))
    r <- if (i == u[u < 200][1])
        pr[["nalc"]] else mosaic(r, pr[["nalc"]], fun=min)
}

writeRaster(r, "d:/bam/2021/gnm/regions/nalc.tif")

## processing boundaries to create BCR/PROV classification for pixels

library(sf)
library(rgeos)
library(raster)
library(jsonlite)
library(ggplot2)


p_bcr <- st_read("d:/bam/2021/gnm/regions/bcr/BCR_Terrestrial_master_International.shp")
p_prov <- st_read("d:/bam/gis/jurisdictions/Political Boundaries (Area).shp")

r_reg <- raster("d:/bam/2021/gnm/regions/BCRSubunits.tif")
p_bcr <- st_transform(p_bcr, st_crs(r_reg))
p_prov <- st_transform(p_prov, st_crs(r_reg))
st_crs(r_reg)
st_crs(p_prov)
st_crs(p_bcr)

table(st_is_valid(p_prov))
table(st_is_valid(p_bcr))
p_prov <- st_make_valid(p_prov)
p_bcr <- st_make_valid(p_bcr)

table(p_prov$COUNTRY)
p_prov <- p_prov[p_prov$COUNTRY=="CAN",]
a <- st_intersects(p_bcr, p_prov)
p_bcr <- p_bcr[sapply(a, length) > 0,]

p_bcr <-  p_bcr[,"BCR"]
p_prov <- p_prov[,"NAME"]
p_prov$NAME <- droplevels(p_prov$NAME)

## cell coordinates from raster object
xy <- data.frame(id=1:length(values(r_reg)), subunit=values(r_reg), coordinates(r_reg))
xy1 <- xy[!is.na(xy$subunit),]
xy1 <- st_as_sf(xy1, coords = c("x", "y"), crs = st_crs(r_reg))

# this takes too long
#i1 <- st_intersection(xy1, p_bcr)
#i2 <- st_intersection(xy1, p_prov)

xy1$BCR <- 100
for (j in 1:nrow(p_bcr)) {
    cat(j, "\n")
    flush.console()
    u <- p_bcr$BCR[j]
    i <- st_intersects(xy1, p_bcr[j,])
    k <- sapply(i, length)
    xy1$BCR[k > 0] <- u
}

xy1$PROV <- factor(NA, levels(p_prov$NAME))
for (j in levels(p_prov$NAME)) {
    cat(j, "\n")
    flush.console()
    i <- st_intersects(xy1, p_prov[p_prov$NAME == j,])
    k <- sapply(i, length)
    xy1$PROV[k > 0] <- j
}

summary(xy1)

save(xy, xy1, file="d:/bam/2021/gnm/regions/BCRPROV.RData")

## example species processing
library(raster)
library(mefa4)
library(jsonlite)

api_root <- "https://borealbirds.github.io/api/v4"
tab <- fromJSON(file.path(api_root, "species"))
rownames(tab) <- tab$id

ROOT <- "d:/bam/BAM_data_v2019/gnm/artifacts"
load("d:/bam/2021/gnm/regions/BCRPROV.RData")
xy$bcr <- xy1$BCR[match(xy$id, xy1$id)]
xy$prov <- xy1$PROV[match(xy$id, xy1$id)]
r_reg <- raster("d:/bam/2021/gnm/regions/BCRSubunits.tif")
r_lc <- raster("d:/bam/2021/gnm/regions/nalc.tif")
levels(xy$prov) <- c("Alberta", "British Columbia", "Manitoba",
"New Brunswick", "Newfoundland and Labrador",
"Northwest Territories", "Nova Scotia",
"Nunavut", "Ontario", "Prince Edward Island",
"Quebec", "Saskatchewan", "water", "Yukon Territory")
xy$prov[xy$prov=="water"] <- NA
xy$prov <- droplevels(xy$prov)
xy$bcr[is.na(xy$prov)] <- NA

## subunit x bcr x prov
xy$val <- paste0(xy$subunit, "_", xy$bcr, "_", xy$prov)
xy$val[is.na(xy$prov)] <- NA
xy$val <- as.factor(xy$val)
s <- !is.na(xy$val)

ninfo <- table(xy$val)
save(ninfo, file=paste0("d:/bam/2021/gnm/ninfo.RData"))
v <- xy$val[s]

#spp <- "OVEN"
for (spp in tab$id) {
    cat(spp, "\n")
    flush.console()

    z <- matrix(NA_real_, length(values(r_reg)), 32)

    for (i in 1:32) {
        ri <- raster(file.path(ROOT, spp, paste0("pred-", spp, "-CAN-boot-", i, ".tif")))
        z[,i] <- values(ri)
    }
    z <- z[s,]
    z[is.na(z)] <- 0

    dmean <- groupMeans(z, 1, v)
    save(dmean, file=paste0("d:/bam/2021/gnm/species/", spp, ".RData"))
}


xy$lc <- values(r_lc)
xy$val1 <- paste0(xy$prov, "_", xy$lc)
xy$val2 <- paste0(xy$subunit, "_", xy$lc)
xy$val1[is.na(xy$prov)] <- NA
xy$val1[is.na(xy$lc)] <- NA
xy$val1[is.na(xy$subunit)] <- NA
xy$val2[is.na(xy$val1)] <- NA
xy$val1 <- as.factor(xy$val1)
xy$val2 <- as.factor(xy$val2)
s1 <- !is.na(xy$val1)

ninfo2 <- table(xy$val1)
ninfo3 <- table(xy$val2)
save(ninfo2, ninfo3, file=paste0("d:/bam/2021/gnm/ninfo2.RData"))
v1 <- xy$val1[s1]
v2 <- xy$val2[s1]

#spp <- "OVEN"
for (spp in tab$id) {
    cat(spp, "\n")
    flush.console()

    z <- matrix(NA_real_, length(values(r_reg)), 32)

    for (i in 1:32) {
        ri <- raster(file.path(ROOT, spp, paste0("pred-", spp, "-CAN-boot-", i, ".tif")))
        z[,i] <- values(ri)
    }
    z <- z[s1,]
    z[is.na(z)] <- 0

    dmean1 <- groupMeans(z, 1, v1)
    dmean2 <- groupMeans(z, 1, v2)
    save(dmean1, dmean2, file=paste0("d:/bam/2021/gnm/species-lc/", spp, ".RData"))
}

## Summarize the summaries

library(jsonlite)
library(mefa4)
api_root <- "https://borealbirds.github.io/api/v4"
tab <- fromJSON(file.path(api_root, "species"))
rownames(tab) <- tab$id

load("d:/bam/2021/gnm/ninfo.RData")
load("d:/bam/2021/gnm/ninfo2.RData")

d <- data.frame(ninfo)
colnames(d) <- c("id", "A")
d$id <- as.character(d$id)
tmp <- strsplit(d$id, "_")
d$subunit <- sapply(tmp, "[[", 1)
d$bcr <- sapply(tmp, "[[", 2)
d$prov <- sapply(tmp, "[[", 3)
rownames(d) <- d$id
d$subunit <- as.integer(d$subunit)
d$subunit <- as.character(
    c(4, 5, 9, 10, 11, 12, 13, 14, 60, 61, 70, 71, 80, 81, 82, 83)[d$subunit])

d$bcrprov <- paste(d$bcr, d$prov)
d$canada <- "Canada"

load(paste0("d:/bam/2021/gnm/species/OVEN.RData"))
d <- d[rownames(dmean),]
all(d$id == rownames(dmean))

p <- read.csv("d:/bam/2021/gnm/PopEsts Global 2020.04.29.csv")
p <- p[,c("Common.Name", "Pair.Adjust.Category")]
tab$pifspp <- tab$english
tab$pifspp[tab$pifspp == "Gray Jay"] <- "Canada Jay"
tab$pifspp[tab$pifspp == "Le Conte's Sparrow"] <- "LeConte's Sparrow"
## all other shorebirds are assumed to have pair adj = 2
compare_sets(tab$pifspp, p$Common.Name)
setdiff(tab$pifspp, p$Common.Name)
tab$pair <- p$Pair.Adjust.Category[match(tab$pifspp, p$Common.Name)]
tab$pair[is.na(tab$pair)] <- 2
table(tab$pair, useNA="a")


# N=dmean*A*100*pair_adj
# dmean from matrix, A from table, 100 is km2->ha, pair adj is ~2 but sometimes not.
f <- function(dmean, u, pair) {
    dd <- d[d$bcr != "100",]
    v <- dd[,u]
    A <- dd$A
    Akm <- sum_by(A, v)
    N <- groupSums(dmean[rownames(dd),] * A * 100 * pair, 1, v)[rownames(Akm),,drop=FALSE]
    D <- groupSums(dmean[rownames(dd),] * A/sum(A), 1, v)[rownames(Akm),,drop=FALSE]
    Nq <- t(apply(N, 1, quantile, c(0.5, 0.05, 0.95)))
    Dq <- t(apply(D, 1, quantile, c(0.5, 0.05, 0.95)))
    data.frame(
        region=rownames(Akm),
        Akm2=Akm[,"x"],
        Ninds_mean=rowMeans(N),
        Ninds_median=Nq[,1],
        Ninds_lower=Nq[,2],
        Ninds_upper=Nq[,3],
        Dvoc_mean=rowMeans(D),
        Dvoc_median=Dq[,1],
        Dvoc_lower=Dq[,2],
        Dvoc_upper=Dq[,3]
    )
}

g <- function(out) {
    rbind(
        data.frame(
            SpeciesID=out$species$id,
            CommonName=out$species$english,
            ScientificName=out$species$scientific,
            PairAdjustment=out$species$pair,
            SummaryBy="CANADA",
            BCR=NA,
            JURS=NA,
            BCRJURS=NA,
            out$canada[,-1]),
        data.frame(
            SpeciesID=out$species$id,
            CommonName=out$species$english,
            ScientificName=out$species$scientific,
            PairAdjustment=out$species$pair,
            SummaryBy="BCR",
            BCR=out$bcr$region,
            JURS=NA,
            BCRJURS=NA,
            out$bcr[,-1]),
        data.frame(
            SpeciesID=out$species$id,
            CommonName=out$species$english,
            ScientificName=out$species$scientific,
            PairAdjustment=out$species$pair,
            SummaryBy="JURS",
            BCR=NA,
            JURS=out$prov$region,
            BCRJURS=NA,
            out$prov[,-1]),
        data.frame(
            SpeciesID=out$species$id,
            CommonName=out$species$english,
            ScientificName=out$species$scientific,
            PairAdjustment=out$species$pair,
            SummaryBy="BCRJURS",
            BCR=NA,
            JURS=NA,
            BCRJURS=out$bcrprov$region,
            out$bcrprov[,-1])
    )
}

h <- function(spp) {

    d1 <- data.frame(ninfo2)
    colnames(d1) <- c("id", "A")
    tmp1 <- strsplit(as.character(d1$id), "_")
    d1$prov <- sapply(tmp1, "[[", 1)
    d1$lc <- sapply(tmp1, "[[", 2)
    d1$id <- as.character(d1$id)

    d2 <- data.frame(ninfo3)
    colnames(d2) <- c("id", "A")
    tmp2 <- strsplit(as.character(d2$id), "_")
    d2$subunit <- sapply(tmp2, "[[", 1)
    d2$lc <- sapply(tmp2, "[[", 2)
    d2$id <- as.character(d2$id)

    dmean1 <- dmean1[d1$id,]
    q1 <- t(apply(dmean1, 1, quantile, c(0.5, 0.05, 0.95)))
    dmean2 <- dmean2[d2$id,]
    q2 <- t(apply(dmean2, 1, quantile, c(0.5, 0.05, 0.95)))

    df1 <- data.frame(
        SpeciesID=tab[spp, "id"],
        CommonName=tab[spp, "english"],
        ScientificName=tab[spp, "scientific"],
        PairAdjustment=tab[spp, "pair"],
        SummaryBy="JURSNALC",
        SUBUNIT=NA,
        JURS=d1$prov,
        NALC=d1$lc,
        Akm2=d1$A,
        Dvoc_mean=rowMeans(dmean1),
        Dvoc_median=q1[,1],
        Dvoc_lower=q1[,2],
        Dvoc_upper=q1[,3])
    df2 <- data.frame(
        SpeciesID=tab[spp, "id"],
        CommonName=tab[spp, "english"],
        ScientificName=tab[spp, "scientific"],
        PairAdjustment=tab[spp, "pair"],
        SummaryBy="JURSNALC",
        SUBUNIT=d2$subunit,
        JURS=NA,
        NALC=d2$lc,
        Akm2=d2$A,
        Dvoc_mean=rowMeans(dmean2),
        Dvoc_median=q2[,1],
        Dvoc_lower=q2[,2],
        Dvoc_upper=q2[,3])
    rbind(df1, df2)
}

res <- NULL
res1 <- NULL
for (spp in tab$id) {
    cat(spp, "\n")
    flush.console()

    load(paste0("d:/bam/2021/gnm/species/", spp, ".RData"))
    load(paste0("d:/bam/2021/gnm/species-lc/", spp, ".RData"))
    stopifnot(all(d$id == rownames(dmean)))

    out <- list(
        species=as.list(tab[spp,]),
        canada=f(dmean, "canada", tab[spp, "pair"]),
        subunit=f(dmean, "subunit", tab[spp, "pair"]),
        bcr=f(dmean, "bcr", tab[spp, "pair"]),
        prov=f(dmean, "prov", tab[spp, "pair"]),
        bcrprov=f(dmean, "bcrprov", tab[spp, "pair"]))
    res <- rbind(res, g(out))
    res1 <- rbind(res1, h(spp))
}

write.csv(res, row.names=FALSE, file="d:/bam/2021/gnm/BAM-GNM-Summaries-2021-03-01.csv")
write.csv(res1, row.names=FALSE, file="d:/bam/2021/gnm/BAM-GNM-Densities-2021-03-01.csv")

LCCs <- c(
    "Conifer"=1,
    "Taiga Conifer"=2,
    "Deciduous"=5,
    "Mixedwood"=6,
    "Shrub"=8,
    "Grass"=10,
    "Arctic Shrub"=11,
    "Arctic Grass"=12,
    "Wetland"=14,
    "Cropland"=15,
    ## exclude ---
    "Barren Lands"=16, # 0
    "Urban and Built-up"=17, # ?
    "Water"=18, # 0
    "Snow and Ice"=19) # 0
BCRs0 <- c(
    "Arctic Plains & Mountains"=3,
    "Northwestern Interior Forest"=4,
    "Northern Pacific Rainforest"=5,
    "Boreal Taiga Plains"=6,
    "Taiga Shield & Hudson Plains"=7,
    "Boreal Softwood Shield"=8,
    "Great Basin"=9,
    "Northern Rockies"=10,
    "Prairie Potholes"=11,
    "Boreal Hardwood Transition"=12,
    "Lower Great Lakes/St. Lawrence Plain"=13,
    "Atlantic Northern Forest"=14)
BCRs0x <- c(
    "Arctic Plains & Mountains"=3,
    "Northwestern Interior Forest"=4,
    "Northern Pacific Rainforest"=5,
    "Boreal Taiga Plains, South"=60,
    "Boreal Taiga Plains, North"=61,
    "Taiga Shield & Hudson Plains, West"=70,
    "Taiga Shield & Hudson Plains, East"=71,
    "Boreal Softwood Shield, West"=80,
    "Boreal Softwood Shield, Ontario"=81,
    "Boreal Softwood Shield, East"=82,
    "Boreal Softwood Shield, Newfoundland"=83,
    "Great Basin"=9,
    "Northern Rockies"=10,
    "Prairie Potholes"=11,
    "Boreal Hardwood Transition"=12,
    "Lower Great Lakes/St. Lawrence Plain"=13,
    "Atlantic Northern Forest"=14)
