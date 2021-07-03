## load WT data

library(mefa4)

fl <- list.files("d:/bam/2021/wbi/wt")

fl <- fl[!(fl %in% c("english_column_definitions.csv",
  "french_column_definitions.csv"))]

#z <- list()
WT <- list()
for (i in seq_along(fl)) {
  cat(fl[i], "\n")
  flush.console()
  x <- read.csv(paste0("d:/bam/2021/wbi/wt/", fl[i]))
  #z[[i]] <- list(dim=dim(x), cn=colnames(x))
  if (nrow(x)) {
    x$csv_file <- fl[i]
    WT <- rbind(WT, x)
  }
}
#t(sapply(z, function(zz) zz$dim))


for (i in 1:ncol(WT))
  if (is.factor(WT[[i]]))
    WT[[i]] <- as.character(WT[[i]])

## filter out 2nd detections
## this removes all tag time stamps after the 1st one
f <- function(v) {
    v <- strsplit(as.character(v), ",")
    v <- sapply(v, function(z) if (length(z) < 1) NA else z[1])
    v <- as.numeric(v)
    v
}
WT$min0_start_num <- f(WT$min0_start)
WT$min1_start_num <- f(WT$min1_start)
WT$min2_start_num <- f(WT$min2_start)
WT$min3_start_num <- f(WT$min3_start)
WT$min4_start_num <- f(WT$min4_start)
WT$min5_start_num <- f(WT$min5_start)
WT$min6_start_num <- f(WT$min6_start)
WT$min7_start_num <- f(WT$min7_start)
WT$min8_start_num <- f(WT$min8_start)
WT$min9_start_num <- f(WT$min9_start)

WT$method[WT$method == "7mVS+3m1SPM"] <- "3m1SPM+7mVS"
WT$method <- gsub(" ", "", WT$method)
data.frame(table(WT$method))
WT$maxdur <- as.integer(sapply(strsplit(as.character(WT$method), "m"), "[[", 1))
WT$maxdis <- Inf
table(WT$maxdur, WT$maxdis)

WT$Start <- strptime(paste(WT$recording_date, WT$recording_time),  "%Y-%m-%d %H:%M:%S")
WT$ARU <- 1
WT$YEAR <- WT$Start$year + 1900

WT$ToY <- WT$Start$yday
WT$ToD <- WT$Start$hour + WT$Start$min / 60

hist(WT$ToD)
hist(WT$ToY)

data.frame(na=colSums(is.na(WT)))

save(WT, file="d:/bam/2021/wbi/WT-Data-fromErin-2021-06-23.RData")

## --------- further process dump

library(mefa4)

load("d:/bam/2021/wbi/WT-Data-fromErin-2021-06-23.RData")

## keep birds and empty detections
table(WT$species_code[WT$species_class==""])
dt <- WT[WT$species_class %in% c("", "AVES"),]

## add SS and PKEY
dt$SS <- dt$location
dt$SSYR <- paste0(dt$location, "-", dt$YEAR)
dt$PKEY <- paste0(dt$location, "-", dt$YEAR,
  "-", dt$recording_date, "-", dt$recording_time)

table(dt$species_code)
data.frame(table(dt$abundance, useNA="a"))
dt$abundance[dt$abundance=="N/A"] <- 1
dt$abundance[dt$abundance==""] <- 1
dt$abundance[dt$abundance=="CI 1 (1 Frog)"] <- 1
dt$abundance[dt$abundance=="CI 2 (>10 Frogs)"] <- 10
dt$abundance[dt$abundance=="CI 3 (>100 Frogs)"] <- 100
dt$abundance[dt$abundance=="TMTT"] <- 100
dt$abundance <- as.integer(dt$abundance)

y <- Xtab(abundance ~ PKEY + species_code, dt)
pk <- nonDuplicated(dt, PKEY, TRUE)
ss <- nonDuplicated(dt, SS, TRUE)
tax <- nonDuplicated(dt, species_code, TRUE)

xy <- ss[,c("SS", "latitude", "longitude")]
#write.csv(xy, row.names=FALSE, file="d:/bam/2021/wbi/WT-XY-2021-06-23.csv")

## offsets



## check SPP names
## make XT
## save bundle

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

xoff <- make_x(
  dt=pk$recording_date,
  tm=substr(pk$recording_time, 1, 5),
  lon=pk$longitude,
  lat=pmin(pk$latitude, 69),
  dur=pk$maxdur,
  dis=pk$maxdis)
rownames(xoff) <- rownames(pk)

setdiff(colnames(y), getBAMspecieslist())
setdiff(getBAMspecieslist(), colnames(y))

SPP <- intersect(colnames(y), getBAMspecieslist())


OFF <- matrix(0, nrow(xoff), length(SPP))
rownames(OFF) <- rownames(xoff)
colnames(OFF) <- SPP

for (spp in SPP) {
  cat(spp, "\n")
  flush.console()
  o <- make_off(spp, xoff)
  OFF[,spp] <- o$offset
}
str(OFF)

save(pk, y, pk, ss, tax, OFF,
  file="d:/bam/2021/wbi/WT-data-processed-2021-06-23.RData")
