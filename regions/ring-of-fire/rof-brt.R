library(mefa4)
library(gbm)
library(dismo)


load("d:/bam/2021/rof/BAMv6_RoFpackage.RData")
#xx2$NEAR_DIST <- NULL
#cn2 <- cn2[cn2 != "NEAR_DIST"]
cn1 <- c("elev", "treecover", "LIDARheight", "road_yesno", "TPI", "TRI",
    "slope", "agriculture_G750.O", "bedrock_G750.O", "bog_G750.O",
    "communities_G750.O", "coniftreed_G750.O", "decidtreed_G750.O",
    "disturbance_G750.O", "fen_G750.O", "heath_G750.O", "marsh_G750.O",
    "mixedtreed_G750.O", "mudflat_G750.O", "openwater_G750.O", "shoreline_G750.O",
    "sparsetreed_G750.O", "swamp_G750.O", "treedupland_G750.O", "turbidwater_G750.O",
    "G750LandCover_NonVeg_v1.grd", "G750LandCover_Veg_v1.grd", "G750LandCover_VegNonTreed_v1.grd",
    "G750LandCover_VegTreed_v1.grd", "G750Species_Abie_Bal_v1.grd",
    "G750Species_Abie_Spp_v1.grd", "G750Species_Acer_Neg_v1.grd",
    "G750Species_Acer_Pen_v1.grd", "G750Species_Acer_Rub_v1.grd",
    "G750Species_Acer_Sac_v1.grd", "G750Species_Acer_Sah_v1.grd",
    "G750Species_Acer_Spi_v1.grd", "G750Species_Acer_Spp_v1.grd",
    "G750Species_Alnu_Rub_v1.grd", "G750Species_Alnu_Spp_v1.grd",
    "G750Species_Betu_All_v1.grd", "G750Species_Betu_Pap_v1.grd",
    "G750Species_Betu_Pop_v1.grd", "G750Species_Betu_Spp_v1.grd",
    "G750Species_Carp_Car_v1.grd", "G750Species_Cary_Cor_v1.grd",
    "G750Species_Fagu_Gra_v1.grd", "G750Species_Frax_Ame_v1.grd",
    "G750Species_Frax_Nig_v1.grd", "G750Species_Frax_Pen_v1.grd",
    "G750Species_Genc_Spp_v1.grd", "G750Species_Genh_Spp_v1.grd",
    "G750Species_Jugl_Cin_v1.grd", "G750Species_Jugl_Nig_v1.grd",
    "G750Species_Juni_Vir_v1.grd", "G750Species_Lari_Lar_v1.grd",
    "G750Species_Lari_Spp_v1.grd", "G750Species_Malu_Spp_v1.grd",
    "G750Species_Ostr_Vir_v1.grd", "G750Species_Pice_Abi_v1.grd",
    "G750Species_Pice_Eng_v1.grd", "G750Species_Pice_Gla_v1.grd",
    "G750Species_Pice_Mar_v1.grd", "G750Species_Pice_Rub_v1.grd",
    "G750Species_Pice_Spp_v1.grd", "G750Species_Pinu_Alb_v1.grd",
    "G750Species_Pinu_Ban_v1.grd", "G750Species_Pinu_Con_v1.grd",
    "G750Species_Pinu_Mon_v1.grd", "G750Species_Pinu_Pon_v1.grd",
    "G750Species_Pinu_Res_v1.grd", "G750Species_Pinu_Spp_v1.grd",
    "G750Species_Pinu_Str_v1.grd", "G750Species_Pinu_Syl_v1.grd",
    "G750Species_Popu_Bal_v1.grd", "G750Species_Popu_Gra_v1.grd",
    "G750Species_Popu_Spp_v1.grd", "G750Species_Popu_Tre_v1.grd",
    "G750Species_Popu_Tri_v1.grd", "G750Species_Prun_Pen_v1.grd",
    "G750Species_Prun_Ser_v1.grd", "G750Species_Quer_Alb_v1.grd",
    "G750Species_Quer_Mac_v1.grd", "G750Species_Quer_Rub_v1.grd",
    "G750Species_Quer_Spp_v1.grd", "G750Species_Sali_Spp_v1.grd",
    "G750Species_Sorb_Ame_v1.grd", "G750Species_Thuj_Occ_v1.grd",
    "G750Species_Tili_Ame_v1.grd", "G750Species_Tsug_Can_v1.grd",
    "G750Species_Tsug_Spp_v1.grd", "G750Species_Ulmu_Ame_v1.grd",
    "G750SpeciesGroups_Broadleaf_Spp_v1.grd", "G750SpeciesGroups_Needleleaf_Spp_v1.grd",
    "G750SpeciesGroups_Unknown_Spp_v1.grd", "G750Structure_Biomass_Branch_v1.grd",
    "G750Structure_Biomass_Foliage_v1.grd", "G750Structure_Biomass_StemBark_v1.grd",
    "G750Structure_Biomass_StemWood_v1.grd", "G750Structure_Biomass_TotalDead_v1.grd",
    "G750Structure_Biomass_TotalLiveAboveGround_v1.grd", "G750Structure_Stand_Age_v1.grd",
    "G750Structure_Stand_CrownClosure_v1.grd", "G750Structure_Stand_Height_v1.grd",
    "G750Structure_Volume_Merch_v1.grd", "G750Structure_Volume_Total_v1.grd",
    "biomass2015.ntems", "volume2015.ntems", "height2015.ntems")

#cn2 <- get_cn(xx2[,cn1])
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

## run BRT with xv

run_brt_xv <- function(spp, RATE=0.001) {
	i <- 1
    si <- BB[,i]
    if (sum(y[si, spp]) < 1)
        return(structure(sprintf("0 detections for %s", spp), class="try-error"))
    xi <- data.frame(
        count=as.numeric(y[si, spp]),
        offset=off[si, spp],
        ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0)[si],
        xx2[si, cn2])
    out <- try(gbm.step(xi,
        gbm.y = 1,
        gbm.x = 3:ncol(xi),
        offset = xi$offset,
        family = "poisson",
        tree.complexity = 3,
        learning.rate = RATE,
        bag.fraction = 0.5))
    if (!inherits(out, "try-error"))
        out$rof_settings <- list(RATE=RATE, spp=spp, i=i)
    out
}


for (spp in SPP) {
    cat("\n\n------------------------------", spp, "------------------------------\n\n")
    res <- run_brt_xv(spp)
    save(res, file=paste0("d:/bam/2021/rof/brt-xv/", spp, ".RData"))
}

## check variable importance

library(mefa4)
library(gbm)
library(dismo)

rel_inf <- function(res) {
    rel.inf <- relative.influence(res, res$n.trees)
    rel.inf[rel.inf < 0] <- 0
    i <- order(-rel.inf)
    rel.inf <- 100 * rel.inf/sum(rel.inf)
    out <- data.frame(var = res$var.names[i], rel.inf = rel.inf[i])
    attr(out, "n.trees") <- res$n.trees
    out
}

RIall <-NULL
for (spp in SPP) {
    cat(spp, "\n")
    load(paste0("d:/bam/2021/rof/brt-xv/", spp, ".RData"))
    if (inherits(res, "gbm")) {
        u <- rel_inf(res)
        u$spp <- spp
        RIall <- rbind(RIall, u)
    }
}
write.csv(RIall, row.names=FALSE, file="d:/bam/2021/rof/SppBRTVarImp.csv")

library(ggplot2)

vi <- read.csv("d:/bam/2021/rof/SppBRTVarImp.csv")
dd <- read.csv("d:/bam/2021/rof/SppDensityByOLCC.csv")

spp <- "OVEN"
for (spp in levels(vi$spp)) {
    cat(spp, "\n")

    vii <- vi[vi$spp == spp,]

    p1 <- ggplot(data=vii[1:10,], aes(x=var, y=rel.inf)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title=spp) +
        xlab("Variables") + ylab("% importance") +
        theme_minimal()
    #p1
    ggsave(sprintf("d:/bam/2021/rof/brt-xv-pred-mosaic/%s-varimp.png", spp), p1)

    ddi <- dd[dd$Species == spp,]
    p2 <- ggplot(data=ddi, aes(x=LCC, y=D.Mean)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title=spp) +
        xlab("Land cover") + ylab("Density (males/ha)") +
        theme_minimal()
    #p2
    ggsave(sprintf("d:/bam/2021/rof/brt-xv-pred-mosaic/%s-dens-by-lcc.png", spp), p2)

}

##--

g <- paste0(xx1$SS_V4, "_", xx1$survey_year)
## looks like there are within year revisits in this data set
length(g)
length(unique(g))



spp <- "OVEN"
i <- 1
RATE=0.001
ntree=1000

si <- BB[,i]

if (sum(y[si, spp]) < 1)
    return(structure(sprintf("0 detections for %s", spp), class="try-error"))
xi <- data.frame(
    count=as.numeric(y[si, spp]),
    offset=off[si, spp],
    ecozone=xx1$ecozone[si],
    xx2[si, cn2])

t1 <- system.time({
    out1 <- try(gbm::gbm(xi$count ~ . + offset(xi$offset),
        data=xi[,-(1:2)],
        n.trees = ntree,
        interaction.depth = 3,
        shrinkage = RATE,
        bag.fraction = 0.5,
        distribution = "poisson",
        var.monotone = NULL,#rep(0, length(4:ncol(DAT))),
        keep.data = FALSE,
        verbose = TRUE,
        n.cores = 1))
})
t2 <- system.time({
    out2 <- try(gbm.step(xi,
        gbm.y = 1,
        gbm.x = 3:ncol(xi),
        offset = xi$offset,
        family = "poisson",
        tree.complexity = 3,
        learning.rate = RATE,
        bag.fraction = 0.5))
})


spp <- "OVEN"
i <- 1
RATE=0.001
ntree=1000


t1 <- system.time({
    out1 <- try(gbm::gbm(xi$count ~ . + offset(xi$offset),
        data=xi[,-(1:2)],
        n.trees = ntree,
        interaction.depth = 3,
        shrinkage = RATE,
        bag.fraction = 0.5,
        distribution = "poisson",
        var.monotone = NULL,#rep(0, length(4:ncol(DAT))),
        keep.data = FALSE,
        verbose = TRUE,
        n.cores = 1))
})
t2 <- system.time({
    out2 <- try(gbm.step(xi,
        gbm.y = 1,
        gbm.x = 3:ncol(xi),
        offset = xi$offset,
        family = "poisson",
        tree.complexity = 3,
        learning.rate = RATE,
        bag.fraction = 0.5))
})



run_brt1 <- function(RUN, SUB=NULL, RATE=0.001) {
    ## parse input
    tmp <- strsplit(RUN, "-")[[1L]]
    spp <- tmp[1L]
    BCR <- tmp[2L]
    ## create data subset
    ss <- dd[,BCR] == 1L
    if (sum(yy[ss, spp]) < 1)
        return(structure(
            sprintf("0 detections for %s in %s", spp, BCR),
            class="try-error"))
    DAT <- data.frame(
        count=as.numeric(yy[ss, spp]),
        offset=off[ss, spp],
        weights=dd$wi[ss],
        ARU=dd$ARU[ss], # ARU added here, but not as layer
        dd2[ss, CN[[BCR]]])
    if (!is.null(SUB)) {
        DAT$weights <- 1
        SUB <- min(SUB, nrow(DAT))
        ## this is not bootstrap resampling
        ## just subset to reduce memory footprint
        DAT <- DAT[sample(nrow(DAT), SUB, prob=DAT$weights),]
    }
    cat("\nFitting gbm.step for", spp, "in", BCR, "\n")
    cat("    Sample size =", nrow(DAT), "\n\n")
    out <- try(gbm.step(DAT,
        gbm.y = 1,
        gbm.x = 4:ncol(DAT),
        offset = DAT$offset,
        site.weights = DAT$weights,
        family = "poisson",
        tree.complexity = 3,
        learning.rate = RATE,
        bag.fraction = 0.5))
    out
}
run_brt2 <- function(RUN, SUB=NULL, RATE=0.001, ntree=1000) {
    ## parse input
    tmp <- strsplit(RUN, "-")[[1L]]
    spp <- tmp[1L]
    BCR <- tmp[2L]
    ## create data subset
    ss <- dd[,BCR] == 1L
    if (sum(yy[ss, spp]) < 1)
        return(structure(
            sprintf("0 detections for %s in %s", spp, BCR),
            class="try-error"))
    DAT <- data.frame(
        count=as.numeric(yy[ss, spp]),
        offset=off[ss, spp],
        weights=dd$wi[ss],
        ARU=dd$ARU[ss], # ARU added here, but not as layer
        dd2[ss, CN[[BCR]]])
    if (!is.null(SUB)) {
        DAT$weights <- 1
        SUB <- min(SUB, nrow(DAT))
        ## this is not bootstrap resampling
        ## just subset to reduce memory footprint
        DAT <- DAT[sample(nrow(DAT), SUB, prob=DAT$weights),]
    }
    cat("\nFitting gbm for", spp, "in", BCR, "\n")
    cat("    Sample size =", nrow(DAT), "\n")
    cat("    N trees     =", ntree, "\n\n")
    out <- try(gbm::gbm(DAT$count ~ . + offset(DAT$offset),
        data=DAT[,-(1:3)],
        n.trees = ntree,
        interaction.depth = 3,
        shrinkage = RATE,
        bag.fraction = 0.5,
        weights = DAT$weights,
        distribution = "poisson",
        var.monotone = NULL,#rep(0, length(4:ncol(DAT))),
        keep.data = FALSE,
        verbose = TRUE,
        n.cores = 1))
    out
}
