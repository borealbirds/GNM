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

## making BCR specific bundles

library(mefa4)
library(sf)

load(file.path(ROOT, "data", "BAMdb-patched-2019-02-04.RData"))
lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

library(sp)

sf <- sf::st_as_sf(dd, coords = c("X","Y"))
sf <- st_set_crs(sf, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sf <- st_transform(sf, lcc_crs)

#i <- 4
for (i in 4:14) {
    cat("\n\nBCR", i, "\n\n")
    bcri <- st_read(file.path(ROOT, "data", "bcr", paste0("bcr",i,"_100km.shp")))
    ddi <- st_intersection(sf, bcri)
    ddi$PKEY <- droplevels(ddi$PKEY)
    ddi$SS <- droplevels(ddi$SS)
    ddi$PCODE <- droplevels(ddi$PCODE)
    yyi <- yy[rownames(ddi),]
    yyi <- yyi[,colSums(yyi>0)>0]
    offi <- off[rownames(yyi), colnames(yyi)]
    spti <- droplevels(spt[colnames(yyi),])
    BCR <- i

    save(BCR, ddi, yyi, offi, spti, file=file.path(ROOT, "data", paste0("BAMdb-bcr", i, "-2019-02-04.RData")))
}




