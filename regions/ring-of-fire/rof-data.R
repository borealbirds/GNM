library(mefa4)
library(jsonlite)
library(sf)
library(rgeos)
library(raster)

tab <- fromJSON("https://borealbirds.github.io/ring-of-fire/species/index.json")
rownames(tab) <- tab$id

# Attached is the file with the locations of ontario bam points north of bcr-12 (SS_V4)
# with designations for whether they are in the study area (study_areaYN),
# which ecozone they are in (ecozone), and whether or not they are in the area of the
# undertaking (aouYN).

# I assume you have the updated habitat data from Lionel,
# where the Ontario land cover goes from columns 101 to 154.
# use SScombo.final.distadded

## rof_bam_pts_ss
load("d:/bam/2021/rof/predictors-data/rof_bam_pts_ss.Rdata")
xx1 <- rof_bam_pts_ss
rm(rof_bam_pts_ss)
e <- new.env()
load("d:/bam/2021/rof/predictors-data/RoF_SpaDESspatialpointcountdata_Feb27.RData", envir=e)
names(e)
xx2 <- e$SScombo.final.distadded
rm(e)
gc()

# eskers
load("d:/bam/2021/rof/predictors-data/RoF_SpaDESspatialpointcountdata_April17.RData")
stopifnot(all(xx2$PKEY_V4 == SScombo.eskerdata.added$PKEY_V4))
# eskerLCC250      prop.eskerLCC250   prob.esker.missed    eskerpoint
summary(SScombo.eskerdata.added[,setdiff(colnames(SScombo.eskerdata.added), colnames(xx2))])
with(SScombo.eskerdata.added, table(lc=eskerLCC250, pt=eskerpoint))

xx2 <- data.frame(xx2, eskerpoint=SScombo.eskerdata.added$eskerpoint)

# V6 data
load("d:/bam/2021/rof/BAMv6_ONBBS.RData")

SPP <- intersect(intersect(tab$id, colnames(y1)), colnames(y2))

# RoF boundary
poly <- st_read("~/GoogleWork/bam/RoF/boundary")
# BAM xy is NAD83, EPSG:4269
xy <- st_as_sf(xx1, coords = c("X", "Y"), crs = 4269)
xy <- st_transform(xy, st_crs(poly))

plot(xy$geometry, pch=".", col=ifelse(xx1$ecozone=="hudson_plain", 2, 4))
plot(poly$geometry, add=TRUE, border="grey", col="#00000022")

di <- st_distance(poly, xy) # units m because of crs
di <- di[1,] / 1000 # unit in km
xx1$dist <- di

rownames(xx1) <- xx1$PKEY_V4
rownames(xx2) <- xx2$PKEY_V4
PKEY <- rownames(xx1)

off <- rbind(o1[,SPP], o2[,SPP])[PKEY,]
y <- rbind(y1[,SPP], y2[,SPP])[PKEY,]
xx2 <- xx2[PKEY,22:ncol(xx2)]

keep <- !is.na(xx2$treecover) & !is.na(xx2$elev) & rowSums(is.na(off)) == 0
xx1 <- xx1[keep,]
xx2 <- xx2[keep,]
y <- y[keep,]
off <- off[keep,]
xx2 <- xx2[,colSums(is.na(xx2))==0]


## high correlations & constant variables
get_cn <- function(z, rmax=0.9) {
    SD <- apply(z, 2, sd, na.rm=TRUE)
    COR <- cor(z[,SD > 0.0001], use="pairwise.complete.obs")
    cr <- mefa:::stack.dist(as.dist(COR), dim.names = TRUE)
    cr <- cr[!is.na(cr$dist),]
    cr <- cr[order(abs(cr$dist), decreasing = TRUE),]
    cr[,1] <- as.character(cr[,1])
    cr[,2] <- as.character(cr[,2])
    while(any(abs(cr$dist) > rmax)) {
        i <- cr[1,2]
        j <- cr[,1] == i | cr[,2] == i
        cr <- cr[!j,]
    }
    union(as.character(cr[,1]), as.character(cr[,2]))
}
cn2 <- get_cn(xx2[,colnames(xx2) != "NEAR_DIST"])

## year matching for surveys ??? not sure what is multi-year, leave it for now

SPP <- colnames(y)[colSums(y) >= 20]
y <- y[,SPP]
off <- off[,SPP]


# resample to have the same n in Boreal Shield than in Hudson Plain
n0 <- sum(xx1$ecozone == "hudson_plain")
i1 <- which(xx1$ecozone == "hudson_plain")
i2 <- which(xx1$ecozone != "hudson_plain")

resample_fun <- function(i) {
    if (i == 1) {
        j1 <- i1
        j2 <- sample(i2, n0, replace=FALSE)
    } else {
        j1 <- sample(i1, n0, replace=TRUE)
        j2 <- sample(i2, n0, replace=TRUE)

    }
    c(j1, j2)
}

B <- 250
set.seed(1)
BB <- sapply(1:B, resample_fun)

save(y, off, xx1, xx2, cn2, SPP, BB, file="d:/bam/2021/rof/BAMv6_RoFpackage_April20.RData")

