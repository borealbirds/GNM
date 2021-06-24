library(mefa4)
library(dismo)
library(gbm)
library(qs)
library(sf)
library(sp)
library(raster)
library(rgdal)
library(nngeo)

qs::qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-11.qRData")

load("~/repos/GNM/regions/wbi/subsets4.RData")

SU <- list("WBAB"=60,
    "WBBC"=c(4,60),
    "WBYT"=4,
    "WBMB"=c(60,70,80),
    "WBSK"=c(60,80),
    "WBNT"=c(61, 70))

lcc_crs <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
wb <- readOGR("d:/bam/2021/wbi/study-area")
wb <- st_as_sf(wb)
wb <- st_transform(wb, lcc_crs)
levels(wb$HASC_1) <- gsub("CA\\.", "WB", levels(wb$HASC_1))
st_area(wb)
wbd <- st_union(wb)
wbd <- nngeo::st_remove_holes(wbd) # remove holes
wbd <- as_Spatial(wbd)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
for (i in 1:222) {
    r4 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr4all_1km.grd")), i)
    if (names(r4) %in% colnames(ddvar)) {
        cat(names(r4), "\n")
        r60 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr60all_1km.grd")), i)
        r61 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr61all_1km.grd")), i)
        r70 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr70all_1km.grd")), i)
        r80 <- raster(file.path(ROOT, "data", "subunits", paste0("bcr80all_1km.grd")), i)
        rm <- mosaic(r4, r60, r61, r70, r80, fun=mean)
        rm <- mask(rm, wbd)
        writeRaster(rm, paste0("d:/bam/2021/wbi/mosaics/", names(r4), ".tif"), overwrite=TRUE)
    }
}

fl <- list.files("d:/bam/2021/wbi/mosaics/")
Nam <- gsub("\\.tif", "", fl)
wb$RID <- 1:6
r <- raster(paste0("d:/bam/2021/wbi/mosaics/", Nam[1], ".tif"))
ii <- rasterize(as_Spatial(wb), r, "RID")

v <- values(ii)
j <- !is.na(v)
prd <- data.frame(index=which(j), region=v[j])

for (i in 1:length(fl)) {
    cat(Nam[i], "\n")
    flush.console()
    r <- raster(paste0("d:/bam/2021/wbi/mosaics/", Nam[i], ".tif"))
    v <- values(r)
    prd[[Nam[i]]] <- v[j]
}
prd$ARU <- 0
prd$YEAR <- 2011

qsavem(prd, file="d:/bam/2021/wbi/pred.qRData")


