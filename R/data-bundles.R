ROOT <- "d:/bam/BAM_data_v2019/gnm"

library(mefa4)

## BAM BBS
e1 <- new.env()
load("d:/bam/Apr2016/out/data_package_2016-12-01.Rdata", envir=e1)
load("d:/bam/Apr2016/out/offsets-v3_2017-04-19.Rdata", envir=e1)
load("d:/bam/Apr2016/out/offsets-v3data_2016-12-01.Rdata", envir=e1)
names(e1)

## BU (AB)
e2 <- new.env()
load("d:/bam/BAM_data_v2019/nwt/BU-nonABMI-offsets-2019-02-04.RData", envir=e2)
names(e2)

## WindTrax (AB)
e3 <- new.env()
load("d:/bam/BAM_data_v2019/nwt/nwt-wildtrax-offsets-2019-01-16.RData", envir=e3)
names(e3)

## Atlas updates
e4 <- new.env()
load("d:/bam/2018/atlas_data/atlas_data_processed-20181018.RData", envir=e4)
names(e4)

## useful attributes

cn <- c("PKEY", "SS", "PCODE","DATE","DATI", "YEAR", "X", "Y", "TSSR", "JDAY", "TREE", "LCC4", "MAXDUR", "MAXDIS")

setdiff(cn, colnames(e2$dd))
e2$dd$PCODE <- interaction("BU", e2$dd$ProjectID, sep="_")
dd2 <- e2$dd[,cn]
dd2$BCR <- NA

setdiff(cn, colnames(e3$dd))
#e3$dd$project.name <- make.names(e3$dd$project.name, unique = TRUE)
e3$dd$PCODE <- interaction("BU", e3$dd$project.name, sep="_")
dd3 <- e3$dd[,cn]
dd3$BCR <- NA

dd4 <- data.frame(e4$PKEY, e4$SS[match(e4$PKEY$SS, e4$SS$SS),])
dd4$YEAR <- dd4$YearCollected
dd4$DATI <- dd4$DATE
dd4$DATE <- as.Date(dd4$DATE)
setdiff(cn, colnames(dd4))
dd4 <- dd4[,cn]
dd4$BCR <- NA

dd1 <- data.frame(e1$PKEY, e1$SS[match(e1$PKEY$SS, e1$SS$SS),], e1$offdat[match(e1$PKEY$PKEY, e1$offdat$PKEY),])
dd1$DATI <- dd1$DATE
dd1$DATE <- as.Date(dd1$DATE)
setdiff(cn, colnames(dd1))
dd1 <- dd1[,c(cn, "BCR")]

dd <- rbind(dd1, dd2, dd3, dd4)
dd <- nonDuplicated(dd, PKEY, TRUE)
dd <- dd[!is.na(dd$X),]
dd <- dd[dd$YEAR <= 2018,]
dd <- dd[dd$YEAR >= 1991,]
dd <- dd[dd$X < 0,]
dd <- dd[dd$Y > 30,]

dd <- dd[,c("PKEY", "SS", "PCODE","DATI","YEAR", "X", "Y","MAXDUR", "MAXDIS")]

## species data and offsets

SPP <- Reduce(intersect, list(colnames(e1$OFF), colnames(e2$off), colnames(e3$off), colnames(e4$OFF)))

off <- rbind(e1$OFF[,SPP], e2$off[,SPP], e3$off[,SPP], e4$OFF[,SPP])
#compare_sets(rownames(off), rownames(dd))
off <- off[rownames(dd),]

YY1 <- Xtab(ABUND ~ PKEY + SPECIES, e1$PCTBL)
YY4 <- e4$YY[,SPP]
YY4 <- YY4[setdiff(rownames(YY4), rownames(YY1)),]
yy <- rbind(YY1[,SPP], e2$y[,SPP], e3$y[,SPP], YY4)
#compare_sets(rownames(yy), rownames(dd))

rn <- intersect(rownames(dd), rownames(yy))
yy <- yy[rn,]
dd <- droplevels(dd[rn,])
off <- off[rn,]
spt <- droplevels(nonDuplicated(e1$TAX, Species_ID, TRUE)[SPP,])

save(dd, yy, off, spt, file=file.path(ROOT, "data", "BAMdb-patched-2019-02-04.RData"))

## making bundle with predictors and buffered BCR indicators

library(mefa4)
library(sf)

load(file.path(ROOT, "data", "BAMdb-patched-2019-02-04.RData"))
lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

library(sp)

sf <- sf::st_as_sf(dd, coords = c("X","Y"))
sf <- st_set_crs(sf, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sf <- st_transform(sf, lcc_crs)

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

dd$YR2 <- ifelse(dd$YEAR < 2006, 1, 2)
dd$SSYR2 <- interaction(dd$SS, dd$YR2, sep="_", drop=TRUE)

dd$o <- seq_len(nrow(dd))
set.seed(1)
dd <- dd[sample(dd$o),]
dd <- dd[!duplicated(dd$SSYR2),]
dd <- dd[order(dd$o),]
dd$o <- NULL

e <- new.env()
load(file.path(ROOT, "data", "predictors", "ss_2001attributes.RData"), envir=e)
dd2001 <- as.data.frame(e$ss)[match(dd$SS, e$ss$SS),]
rownames(dd2001) <- rownames(dd)
rm(e)
dd2001$SS <- dd2001$PCODE <- dd2001$lat <- dd2001$lon <- NULL
dd2001$geometry <- NULL

e <- new.env()
load(file.path(ROOT, "data", "predictors", "ss_2011attributes.RData"), envir=e)
dd2011 <- as.data.frame(e$ss)[match(dd$SS, e$ss$SS),]
rownames(dd2011) <- rownames(dd)
rm(e)
dd2011$SS <- dd2011$PCODE <- dd2011$lat <- dd2011$lon <- NULL
dd2011$geometry <- NULL

stopifnot(all(colnames(dd2001)==colnames(dd2011)))

dd2 <- dd2001 # use 2001 version for -2005 data
dd2[dd$YEAR >= 2006,] <- dd2011[dd$YEAR >= 2006,] # use 2011 version for 2006- data

dd2$lf[dd2$lf < 1] <- 0
dd2$lf <- as.factor(as.integer(dd2$lf))
dd2$nalc <- as.factor(dd2$nalc)

sf <- sf[rownames(dd),]
yy <- yy[rownames(dd),]
off <- off[rownames(dd),]

#i <- 4
for (i in 4:14) {
    cat("\n\nBCR", i, "\n\n")
    bcri <- st_read(file.path(ROOT, "data", "bcr", paste0("bcr",i,"_100km.shp")))
    ddi <- st_intersection(sf, bcri)
    dd[[paste0("BCR_", i)]] <- ifelse(rownames(dd) %in% rownames(ddi), 1L, 0L)
}


ss <- rowSums(is.na(dd2)) == 0
dd <- dd[ss,]
dd2 <- dd2[ss,]
yy <- yy[ss,]
off <- off[ss,]
sf <- sf[ss,]

any(colSums(yy) == 0)

c(sum(is.na(dd)), sum(is.na(dd2)), sum(is.na(off)), sum(is.na(yy)))
colSums(is.na(dd)) # some date missing -- we let it go

## evaluate predictor set based on hist, SD, etc

get_cn <- function(z, rmax=0.9) {
    COR <- cor(z)
    cr <- mefa:::stack.dist(as.dist(COR), dim.names = TRUE)
    cr <- cr[order(abs(cr$dist), decreasing = TRUE),]
    while(any(abs(cr$dist) > rmax)) {
        i <- as.character(cr[1,2])
        j <- cr[,1] == i | cr[,2] == i
        cr <- cr[!j,]
    }
    union(as.character(cr[,1]), as.character(cr[,2]))
}

## factor
cnf <- colnames(dd2)[!sapply(dd2, is.numeric)]

z <- as.matrix(dd2[,sapply(dd2, is.numeric)])
cnn <- get_cn(z)

CN <- list()
for (i in 4:14) {
    cat("BCR", i, "\n")
    BCR <- paste0("BCR_", i)
    z <- as.matrix(dd2[dd[,BCR] == 1L, cnn])
    #MEAN <- colMeans(z)
    SD <- apply(z, 2, sd)
    CN[[BCR]] <- get_cn(z[,SD > 0])
}

sapply(CN, length)

## weights

xy <- st_coordinates(sf)
bb <- bbox(xy)
r <- raster(
    xmn=bb["x", "min"],
    xmx=bb["x", "max"],
    ymn=bb["y", "min"],
    ymx=bb["y", "max"],
    res=250,
    crs=st_crs(sf))
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

save(dd, dd2, yy, off, spt, cnf, cnn, CN, nsub,
    file=file.path(ROOT, "data", "BAMdb-GNMsubset-2019-03-01.RData"))
#load(file.path(ROOT, "data", "BAMdb-GNMsubset-2019-03-01.RData"))
