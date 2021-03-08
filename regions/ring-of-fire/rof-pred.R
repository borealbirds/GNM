library(mefa4)
library(jsonlite)
library(sf)
library(rgeos)
library(raster)

#pl <- st_read("d:/bam/2021/rof/predictors-layers/ecoregions-clipped/BCR7and8Ontario.shp")

fl <- list.files("d:/bam/2021/rof/predictors-layers/Clip to BCR 7 and 8 LCC")
names(fl) <- gsub("\\.tif", "", fl)

#load("d:/bam/2021/rof/BAMv6_RoFpackage.RData")

cn2 <- c("agriculture_G750.O", "bedrock_G750.O", "biomass2015.ntems",
    "bog_G750.O", "communities_G750.O", "coniftreed_G750.O", "decidtreed_G750.O",
    "disturbance_G750.O", "elev", "fen_G750.O", "G750LandCover_Veg_v1.grd",
    "G750LandCover_VegNonTreed_v1.grd", "G750LandCover_VegTreed_v1.grd",
    "G750Species_Abie_Bal_v1.grd", "G750Species_Acer_Neg_v1.grd",
    "G750Species_Acer_Pen_v1.grd", "G750Species_Acer_Rub_v1.grd",
    "G750Species_Acer_Sac_v1.grd", "G750Species_Acer_Sah_v1.grd",
    "G750Species_Acer_Spi_v1.grd", "G750Species_Acer_Spp_v1.grd",
    "G750Species_Alnu_Spp_v1.grd", "G750Species_Betu_All_v1.grd",
    "G750Species_Betu_Pap_v1.grd", "G750Species_Betu_Pop_v1.grd",
    "G750Species_Fagu_Gra_v1.grd", "G750Species_Frax_Ame_v1.grd",
    "G750Species_Frax_Nig_v1.grd", "G750Species_Frax_Pen_v1.grd",
    "G750Species_Genc_Spp_v1.grd", "G750Species_Genh_Spp_v1.grd",
    "G750Species_Lari_Lar_v1.grd", "G750Species_Pice_Abi_v1.grd",
    "G750Species_Pice_Gla_v1.grd", "G750Species_Pice_Mar_v1.grd",
    "G750Species_Pice_Rub_v1.grd", "G750Species_Pinu_Ban_v1.grd",
    "G750Species_Pinu_Con_v1.grd", "G750Species_Pinu_Res_v1.grd",
    "G750Species_Pinu_Str_v1.grd", "G750Species_Popu_Bal_v1.grd",
    "G750Species_Popu_Gra_v1.grd", "G750Species_Popu_Spp_v1.grd",
    "G750Species_Popu_Tre_v1.grd", "G750Species_Prun_Pen_v1.grd",
    "G750Species_Quer_Mac_v1.grd", "G750Species_Quer_Rub_v1.grd",
    "G750Species_Sali_Spp_v1.grd", "G750Species_Sorb_Ame_v1.grd",
    "G750Species_Thuj_Occ_v1.grd", "G750Species_Tili_Ame_v1.grd",
    "G750Species_Tsug_Can_v1.grd", "G750Species_Ulmu_Ame_v1.grd",
    "G750SpeciesGroups_Broadleaf_Spp_v1.grd", "G750SpeciesGroups_Needleleaf_Spp_v1.grd",
    "G750SpeciesGroups_Unknown_Spp_v1.grd", "G750Structure_Biomass_TotalDead_v1.grd",
    "G750Structure_Stand_Age_v1.grd", "G750Structure_Volume_Total_v1.grd",
    "heath_G750.O", "height2015.ntems", "LIDARheight", "marsh_G750.O",
    "mixedtreed_G750.O", "mudflat_G750.O", "openwater_G750.O", "road_yesno",
    "slope", "sparsetreed_G750.O", "swamp_G750.O", "TPI", "treecover",
    "turbidwater_G750.O", "volume2015.ntems")
compare_sets(names(fl), cn2)

r <- raster(paste0(
        "d:/bam/2021/rof/predictors-layers/Clip to BCR 7 and 8 LCC/",
       "elev.tif"))
re <- r
values(re)[!is.na(values(re))] <- 1
pl <- rgdal::readOGR("d:/bam/2021/rof/predictors-layers/ecoregions-clipped/BCR7and8Ontario.shp")
pl <- spTransform(pl, proj4string(re))
re <- mask(re, pl[1,])
values(re)[is.na(values(re))] <- 0
values(re)[is.na(values(r))] <- NA

writeRaster(re, paste0(
        "d:/bam/2021/rof/predictors-layers/Clip to BCR 7 and 8 LCC/",
       "ecozone.tif"))

n <- length(values(r))
Chunks <- sort(sample(1:10, n, replace = TRUE))
table(Chunks)
save(Chunks, file="d:/bam/2021/rof/predictors-layers/chunks/chunks.RData")

#j <- 1
cn2x <- c("ecozone", cn2)
for (j in 1:10) {
    cat("Preparing chunk", j, "\n")
    flush.console()
    ss <- Chunks==j
    M <- matrix(0, sum(ss), length(cn2x))
    dimnames(M) <- list(which(ss), cn2x)
    M <- data.frame(M)

    for (k in cn2x) {
        cat(j, k, "\n")
        flush.console()

        r <- raster(paste0(
            "d:/bam/2021/rof/predictors-layers/Clip to BCR 7 and 8 LCC/",
            k, ".tif"))
        M[,k] <- values(r)[ss]

    }

    save(M, file=paste0("d:/bam/2021/rof/predictors-layers/chunks/variables-", j, ".RData"))
}


## making chunked predictions for species


library(dismo)
library(gbm)
library(raster)

r <- raster(paste0(
        "d:/bam/2021/rof/predictors-layers/Clip to BCR 7 and 8 LCC/",
       "elev.tif"))
s <- !is.na(values(r))

load("d:/bam/2021/rof/BAMv6_RoFpackage.RData")

#spp <- "OVEN"
#j <- 1
CHUNKS <- 1:10
for (j in CHUNKS) {
    gc()
    load(paste0("d:/bam/2021/rof/predictors-layers/chunks/variables-", j, ".RData"))
    for (spp in SPP) {
        cat(j, spp, "\n")
        flush.console()
        load(paste0("d:/bam/2021/rof/brt-xv/", spp, ".RData"))
        if (inherits(res, "gbm")) {
            if (!dir.exists(paste0("d:/bam/2021/rof/brt-xv-pred/", spp)))
                dir.create(paste0("d:/bam/2021/rof/brt-xv-pred/", spp))
            p <- predict.gbm(res, M, res$n.trees, type="response")
        save(p,
            file=paste0("d:/bam/2021/rof/brt-xv-pred/", spp, "/", spp, "-chunk-", j, ".RData"))
        }
    }
}


## putting together the pieces
library(raster)
SPP <- list.dirs("d:/bam/2021/rof/brt-xv-pred", full.names=FALSE)[-1]
load("d:/bam/2021/rof/predictors-layers/chunks/chunks.RData")
r <- raster(paste0(
        "d:/bam/2021/rof/predictors-layers/Clip to BCR 7 and 8 LCC/",
       "elev.tif"))
s <- !is.na(values(r))

for (spp in SPP) {
    spppred <- list()
    for (j in 1:10) {
        f <- paste0("d:/bam/2021/rof/brt-xv-pred/",
            spp, "/", spp, "-chunk-", j, ".RData")
        load(f)
        spppred[[j]] <- p
    }
    v <- unlist(spppred)
    ri <- r
    values(ri)[s] <- v[s]
    writeRaster(ri,
        paste0("d:/bam/2021/rof/brt-xv-pred-mosaic/", spp, ".tif"),
        overwrite=TRUE)
}

## make png maps

library(sf)
library(rgeos)
library(raster)

# RoF boundary
r <- raster(paste0("d:/bam/2021/rof/brt-xv-pred-mosaic/OVEN.tif"))
pl <- rgdal::readOGR("~/GoogleWork/bam/RoF/boundary")
pl <- spTransform(pl, proj4string(r))
#u <- mask(r, pl)

for (spp in SPP) {
    cat(spp, "\n")
    flush.console()

    ri <- raster(paste0("d:/bam/2021/rof/brt-xv-pred-mosaic/", spp, ".tif"))
    q <- quantile(ri, 0.999)
    values(ri)[!is.na(values(ri)) & values(ri) > q] <- q
    u <- mask(ri, pl)
    N <- round(2 * sum(values(u), na.rm=TRUE) * 2.5^2 / 10^6, 3)

    png(paste0("d:/bam/2021/rof/brt-xv-pred-mosaic/", spp, ".png"), width=500, height=500)
    op <- par(mfrow=c(1,1), mar=c(1,1,2,6))
    plot(ri, axes=FALSE, box=FALSE, col=hcl.colors(100, "Lajolla"),
        main=paste(spp, "Mean Density (males/ha)\nPopulation size =", N, "M inds."))
    plot(pl, add=TRUE, border=4)
    par(op)
    dev.off()

}

## density by landcover

lc <- cn2[endsWith(cn2, ".O")]
lcc <- sapply(strsplit(lc, "_"), "[[", 1)

LC <- list()

for (k in seq_along(lc)) {
    cat(lc[k], "\n")
    flush.console()

    r <- raster(paste0(
        "d:/bam/2021/rof/predictors-layers/Clip to BCR 7 and 8 LCC/",
        lc[k], ".tif"))
    LC[[lcc[k]]] <- values(r)
}
LC <- do.call(cbind, LC)
rs <- rowSums(LC, na.rm=TRUE)
summary(rs)
range(LC, na.rm=TRUE)
LC[is.na(LC)] <- 0
wm <- find_max(LC)
r <- raster(paste0("d:/bam/2021/rof/brt-xv-pred-mosaic/OVEN.tif"))
wm$index[is.na(values(r))] <- NA

rlcDom <- r
values(rlcDom) <- wm$index
#writeRaster(rlcDom, file="d:/bam/2021/rof/dominant-landcover.tif")
levels(wm$index)

u0 <- mask(rlcDom, pl)
s <- !is.na(values(u0))
DD <- NULL
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()

    ri <- raster(paste0("d:/bam/2021/rof/brt-xv-pred-mosaic/", spp, ".tif"))
    q <- quantile(ri, 0.999)
    ri <- mask(ri, pl)
    values(ri)[!is.na(values(ri)) & values(ri) > q] <- q
    ag <- aggregate(data.frame(D=values(ri)[s]),
        list(LCC=levels(wm$index)[values(u0)[s]]), summary)
    ag <- data.frame(Species=spp, ag)

    DD <- rbind(DD, ag)
}

write.csv(DD, row.names=FALSE, file="d:/bam/2021/rof/SppDensityByOLCC.csv")
