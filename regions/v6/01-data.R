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


## load data for ON BBS

load("d:\\bam\\2021\\rof\\ON-BBS\\BBS.91.19.ON.stopsFeb2021.RData")
#load("d:\\bam\\2021\\rof\\ON-BBS\\BBS.ON.wide.RData")

m5$PKEY <- paste0(m5$SS, ":", m5$Year)
#compare_sets(m5$PKEY, spp_count$PKEY)

y1 <- Xtab(Abund ~ PKEY+Species_ID, m5)
x1 <- nonDuplicated(m5, PKEY, TRUE)[rownames(y1),]

head(x1$date)
head(x1$StartTime)
x1off <- make_x(
    dt=x1$date,
    tm=paste0(x1$StartTime.Hour, ":", x1$StartTime.Min),
    lon=x1$POINT_X,
    lat=x1$POINT_Y,
    dur=3,
    dis=Inf)
rownames(x1off) <- rownames(x1)
summary(x1off)


s1 <- intersect(colnames(y1), getBAMspecieslist())

o1 <- matrix(0, nrow(x1off), length(s1))
rownames(o1) <- rownames(x1)
colnames(o1) <- s1

for (spp in s1) {
  cat(spp, "\n")
  flush.console()
  o1[,spp] <- make_off(spp, x1off)$offset
}
str(o1)
sum(is.na(o1))

## BAM V6

ss <- read.csv("d:\\bam\\2021\\rof\\BAM-V6\\location.txt")
vi <- read.csv("d:\\bam\\2021\\rof\\BAM-V6\\pc_visit.txt")
dt <- read.csv("d:\\bam\\2021\\rof\\BAM-V6\\pc_detection.txt")
l1 <- read.csv("d:\\bam\\2021\\rof\\BAM-V6\\lu_pc_protocol_duration.txt")
l1$maxdur <- as.numeric(sapply(strsplit(as.character(l1[,2]), "-"), function(z) z[length(z)]))
l2 <- read.csv("d:\\bam\\2021\\rof\\BAM-V6\\lu_pc_protocol_distance.txt")
l2$maxdis <- as.numeric( sapply(strsplit(as.character(l2[,2]), "-"), function(z) z[length(z)]))

## AK data
if (FALSE) {
str(ss)
si <- as.character(ss$SS_V4)
si <- sapply(strsplit(si, ":"), "[[", 1)
p <- unique(si)
p[startsWith(p, "AK")]

ss_ak <- droplevels(ss[startsWith(as.character(ss$SS_V4), "AK"),])
vi_ak <- droplevels(vi[startsWith(as.character(vi$PKEY_V4), "AK"),])
dt_ak <- droplevels(dt[startsWith(as.character(dt$PKEY), "AK"),])

i1 <- read.csv("~/repos/bamanalytics/lookup/durint.csv")
i2 <- read.csv("~/repos/bamanalytics/lookup/disint.csv")

vi_ak <- data.frame(vi_ak,
  l1[match(vi_ak$protocol_duration, l1$protocol_duration_id),],
  l2[match(vi_ak$protocol_distance, l2$protocol_distance_numid),])
dt_ak <- data.frame(dt_ak,
  i1[match(dt_ak$duration_interval, i1$DurationID),],
  i2[match(dt_ak$distance_band, i2$DistanceID),])


ss_ak <- ss_ak[,c("SS_V4",
"latitude", "longitude",
"location_comments", "OnRoad_V4", "timezone_V4")]
vi_ak <- vi_ak[,c("location_name_4", "PKEY_V4",
"protocol_duration", "protocol_distance", "ROUND", "survey_year",
"survey_date", "survey_time", "observer_id", "origin_database",
"time_zone", "visit_comments", "METHOD_V4", "DurationMethod_V4",
"DistanceMethod_V4", "MM_V4", "DD_V4", "HR_V4", "MIN_V4", "mm.dd.yyy_V4",
"StartTime_V4", "QC_V4", "Version_V4", "obs_V4", "protocol_duration_id",
"protocol_duration_range", "maxdur", "protocol_distance_numid",
"protocol_distance_range", "maxdis")]
dt_ak <- dt_ak[,c( "PKEY",  "species_code",
"original_species", "duration_interval", "distance_band", "abundance",
"detected", "heard", "seen", "pc_vt", "pc_vt_detail", "age",
"fm", "group", "flyover", "displaytype", "nestevidence", "behaviourother",
"atlas_breeding_code", "detection_comments", "BehCodeBAMV4",
"DurationID", "DUR_Descrip", "TimeEquiv", "Dur_Start", "DUR_end",
"DistanceID", "DIST_START", "DIST_END", "DistanceDescrip")]

save(ss_ak, vi_ak, dt_ak, file="d:/bam/2021/rof/BAM-V6/AK-data.RData")

}

compare_sets(vi$PKEY_V4, dt$PKEY)
compare_sets(vi$PKEY_V6, dt$PKEY_V6)

y2 <- Xtab(abundance ~ PKEY + species_code, dt)
x2 <- data.frame(vi, ss[match(vi$location_name_4, ss$SS_V4),])
rownames(x2) <- x2$PKEY_V4
compare_sets(rownames(y2), rownames(x2))

x2 <- x2[rownames(y2),]
x2$maxdur <- l1$maxdur[match(x2$protocol_duration, l1$protocol_duration_id)]
x2$maxdis <- l2$maxdis[match(x2$protocol_distance, l2$protocol_distance_numid)]

tmp <- strsplit(as.character(x2$survey_time), " ")
table(sapply(tmp, length))
tt <- as.character(x2$survey_time)
tt[sapply(tmp, length) < 2] <- NA
tt <- as.POSIXlt(tt)
hr <- tt$hour
mn <- tt$min
tm <- paste0(
  ifelse(hr < 10, "0", ""),
  hr,
  ":",
  ifelse(mn < 10, "0", ""),
  mn
)
tm[is.na(tt)] <- NA
summary(nchar(tm))

x2$longitude[x2$longitude > -40] <- NA
x2$latitude[x2$latitude < 30] <- NA

x2$longitude[x2$longitude < -163.89] <- -163.89
x2$latitude[x2$latitude > 68.98] <- 68.98
summary(x2$longitude)
summary(x2$latitude)


tmp <- data.frame(
    dt=as.Date(x2$survey_date),
    tm=tm,
    lon=x2$longitude,
    lat=x2$latitude,
    dur=x2$maxdur,
    dis=x2$maxdis)
rownames(tmp) <- rownames(x2)
#tmp <- tmp[!is.na(tmp$lat) & !is.na(tmp$lon),]

x2off <- make_x(
    dt=tmp$dt,
    tm=tmp$tm,
    lon=tmp$lon,
    lat=tmp$lat,
    dur=tmp$dur,
    dis=tmp$dis)
rownames(x2off) <- rownames(tmp)
summary(x2off)


s2 <- intersect(colnames(y2), getBAMspecieslist())

o2 <- matrix(0, nrow(x2off), length(s2))
rownames(o2) <- rownames(x2off)
colnames(o2) <- s2

for (spp in s2) {
  cat(spp, "\n")
  flush.console()
  o2[,spp] <- make_off(spp, x2off)$offset
}
str(o2)
sum(is.na(o2))

save(y1, x1, x1off, o1,
  y2, x2, x2off, o2,
  file="d:\\bam\\2021\\rof\\BAMv6_ONBBS.RData")



