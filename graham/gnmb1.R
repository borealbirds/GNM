## --- settings ---
## file name for data bundle, need to be in /data/ dir
fn <- "BAMdb-GNMsubset-2020-01-08.RData"
## project name for storing the output
PROJ <- "all2"
## output directory
OUTDIR <- if (interactive())
    paste0("d:/bam/BAM_data_v2019/gnm/out/", PROJ) else paste0("/scratch/psolymos/out/", PROJ)

cat("* Loading packages and sourcing functions:")
library(parallel)
library(mefa4)
library(gbm)

run_brt_boot <- function(b, spp) {
    #sppbcr <- SPPBCR[grep(spp, SPPBCR)]
    sppbcr <- paste0(spp, "-BCR_", u)
    if (!dir.exists(paste0(OUTDIR, "/", spp)))
        dir.create(paste0(OUTDIR, "/", spp))
    for (RUN in sppbcr) {
        bcr <- strsplit(RUN, "-")[[1]][2]
        if (!dir.exists(paste0(OUTDIR, "/", spp, "/", bcr)))
            dir.create(paste0(OUTDIR, "/", spp, "/", bcr))
        fout <- paste0(OUTDIR, "/", spp, "/", bcr, "/gnmboot-",
            spp, "-", bcr, "-", b, ".RData")
        if (!file.exists(fout)) {
            out <- .run_brt_boot(b, RUN)
            save(out, file=fout)
        }
    }
    invisible(TRUE)
}

## b: bootstrap id (>= 1)
## RUN: SPP-BCR_NO tag
.run_brt_boot <- function(b, RUN, verbose=interactive()) {
    t0 <- proc.time()["elapsed"]
    ## parse input
    tmp <- strsplit(RUN, "-")[[1L]]
    spp <- tmp[1L]
    BCR <- tmp[2L]
    bcr <- as.integer(strsplit(BCR, "_")[[1L]][2L])
    if (bcr %% 100 == 0) {
        ## US all clim+topo
        cn <- CN[[BCR]]
        nt <- ntmax[spp]
    } else {
        ## Canada: based on XV
        cn <- names(xvinfo[[RUN]]$varimp)[xvinfo[[RUN]]$varimp > 0]
        nt <- xvinfo[[RUN]]$ntree
    }

    ## create data subset for BCR unit
    ss <- dd[,BCR] == 1L
    DAT <- data.frame(
        count=as.numeric(yy[ss, spp]),
        offset=off[ss, spp],
        #weights=dd$wi[ss],
        cyid=dd$cyid[ss],
        YEAR=dd$YEAR[ss],
        ARU=dd$ARU[ss], # ARU added here, but not as layer
        dd2[ss, cn])
    ## subsample based on 2.5x2.5km^2 cell x year units
    DAT <- DAT[sample.int(nrow(DAT)),]
    DAT <- DAT[!duplicated(DAT$cyid),]
    if (b > 1)
        DAT <- DAT[sample.int(nrow(DAT), replace=TRUE),]
    ## 0 detection output
    if (sum(DAT$count) < 1) {
        out <- structure(
            sprintf("0 detections for %s in %s", spp, BCR),
            class="try-error")
    } else {
        if (verbose)
            cat("\nFitting gbm::gbm for", spp, "in", BCR, "/ b =", b, "/ n =",nrow(DAT), "\n")
        out <- try(gbm::gbm(DAT$count ~ . + offset(DAT$offset),
            data=DAT[,-(1:3)],
            n.trees = nt,
            interaction.depth = 3,
            shrinkage = 1/nt,
            bag.fraction = 0.5,
            #weights = DAT$weights,
            distribution = "poisson",
            var.monotone = NULL,#rep(0, length(4:ncol(DAT))),
            keep.data = FALSE,
            verbose = verbose,
            n.cores = 1))
    }
    attr(out, "__settings__") <- list(
        species=spp, region=BCR, iteration=b, elapsed=proc.time()["elapsed"]-t0)
    out
}
run_brt_boot_all <- function(spp, b=1, nmax=50000, nt=10000, clim_only=FALSE) {
    if (!dir.exists(paste0(OUTDIR, "/", spp)))
        dir.create(paste0(OUTDIR, "/", spp))
    bcr <- "ALL"
    if (!dir.exists(paste0(OUTDIR, "/", spp, "/", bcr)))
        dir.create(paste0(OUTDIR, "/", spp, "/", bcr))
    fout <- paste0(OUTDIR, "/", spp, "/", bcr, "/gnmboot-",
        spp, "-", bcr, "-", b,
        if (clim_only) "-climonly" else "",
        ".RData")
    if (!file.exists(fout)) {
        out <- .run_brt_boot_all(b, spp, nmax=nmax, nt=nt, clim_only = clim_only)
        save(out, file=fout)
    }
    invisible(TRUE)
}
.run_brt_boot_all <- function(b, spp, nmax=50000, nt=10000, clim_only=FALSE, verbose=interactive()) {
    t0 <- proc.time()["elapsed"]
    ## create data subset for BCR unit
    cn1 <- c("nalc", "ROAD", "Landsc750_Tsug_Spp_v1", "Landsc750_Alnu_Rub_v1",
        "LandCover_VegTreed_v1.1", "Species_Fagu_Gra_v1", "Landsc750_Betu_All_v1",
        "Structure_Volume_Total_v1", "Landsc750_Acer_Sah_v1", "Landsc750_Tsug_Het_v1",
        "Landsc750_Popu_Tri_v1", "Landsc750_Abie_Las_v1", "Landsc750_Jugl_Nig_v1",
        "Landsc750_Popu_Spp_v1", "Landsc750_Betu_Pap_v1", "Species_Jugl_Nig_v1",
        "Landsc750_Acer_Spp_v1", "Landsc750_Ulmu_Ame_v1", "Landsc750_Pinu_Str_v1",
        "Landsc750_Arbu_Men_v1", "Landsc750_Pinu_Ban_v1", "Landsc750_Prun_Pen_v1",
        "Landsc750_Tili_Ame_v1", "Landsc750_Acer_Mac_v1", "Landsc750_Pice_Spp_v1",
        "Structure_Stand_Age_v1", "Landsc750_Genc_Spp_v1", "Landsc750_Cham_Noo_v1",
        "Landsc750_Abie_Spp_v1", "Landsc750_Pice_Gla_v1", "Landsc750_Ostr_Vir_v1",
        "Landsc750_Frax_Ame_v1", "Landsc750_Popu_Bal_v1", "Landsc750_Pinu_Pon_v1",
        "Species_Betu_All_v1", "Landsc750_Popu_Gra_v1", "TD", "Landsc750_Thuj_Occ_v1",
        "Landsc750_Pice_Eng_v1", "Structure_Stand_CrownClosure_v1", "Landsc750_Pice_Sit_v1",
        "Landsc750_Abie_Ama_v1", "Landsc750_Lari_Occ_v1", "Species_Quer_Alb_v1",
        "Landsc750_Tsug_Mer_v1", "SpeciesGroups_Broadleaf_Spp_v1", "Landsc750_Pinu_Mon_v1",
        "Species_Genh_Spp_v1", "PPT_wt", "Landsc750_Carp_Car_v1", "Landsc750_Sorb_Ame_v1",
        "Landsc750_Pinu_Res_v1", "Landsc750_Acer_Sac_v1", "Landsc750_Betu_Spp_v1",
        "Landsc750_Lari_Spp_v1", "Species_Quer_Rub_v1", "Landsc750_Lari_Lar_v1",
        "Landsc750_Frax_Nig_v1", "Landsc750_Frax_Pen_v1", "Species_Tsug_Can_v1",
        "Species_Alnu_Rub_v1", "Structure_Biomass_TotalDead_v1", "Landsc750_Pinu_Spp_v1",
        "Landsc750_Prun_Ser_v1", "Species_Tsug_Het_v1", "Species_Prun_Pen_v1",
        "Landsc750_Quer_Mac_v1", "NFFD", "Species_Tili_Ame_v1", "Landsc750_Betu_Pop_v1",
        "Species_Tsug_Spp_v1", "Tave_sm", "Landsc750_Alnu_Spp_v1", "Landsc750_Pinu_Syl_v1",
        "SpeciesGroups_Unknown_Spp_v1", "Landsc750_Acer_Spi_v1", "Species_Popu_Tri_v1",
        "Species_Popu_Tre_v1", "Species_Acer_Spp_v1", "Landsc750_Cary_Cor_v1",
        "Landsc750_Pinu_Alb_v1", "SHM", "Species_Pice_Spp_v1", "Species_Betu_Pap_v1",
        "Landsc750_Pice_Abi_v1", "Landsc750_Malu_Spp_v1", "Species_Frax_Ame_v1",
        "Species_Pice_Rub_v1", "Landsc750_Lari_Lya_v1", "Species_Pinu_Ban_v1",
        "Landsc750_Juni_Vir_v1", "roughness", "Species_Pinu_Str_v1",
        "Landsc750_Sali_Spp_v1", "LandCover_Veg_v1.1", "Landsc750_Acer_Neg_v1",
        "Species_Acer_Rub_v1", "Species_Ostr_Vir_v1", "Landsc750_Jugl_Cin_v1",
        "Species_Popu_Spp_v1", "Species_Thuj_Pli_v1", "Landsc750_Acer_Pen_v1",
        "Species_Carp_Car_v1", "Species_Ulmu_Ame_v1", "Species_Popu_Gra_v1",
        "Species_Thuj_Occ_v1", "Species_Pice_Mar_v1", "Species_Sorb_Ame_v1",
        "Species_Pseu_Men_v1", "Species_Pice_Eng_v1", "Species_Pinu_Con_v1",
        "Species_Prun_Ser_v1", "Species_Frax_Nig_v1", "Species_Acer_Sah_v1",
        "Species_Betu_Pop_v1", "Species_Quer_Spp_v1", "Species_Pinu_Mon_v1",
        "Species_Abie_Bal_v1", "Species_Genc_Spp_v1", "Species_Pice_Gla_v1",
        "led750", "Species_Pinu_Res_v1", "Species_Cham_Noo_v1", "Species_Lari_Spp_v1",
        "Species_Pice_Sit_v1", "Species_Popu_Bal_v1", "Species_Lari_Lar_v1",
        "Species_Abie_Spp_v1", "Species_Pinu_Pon_v1", "Species_Alnu_Spp_v1",
        "Species_Tsug_Mer_v1", "TPI", "Species_Quer_Mac_v1", "Species_Acer_Spi_v1",
        "Species_Betu_Spp_v1", "Species_Pinu_Alb_v1", "Species_Acer_Pen_v1",
        "Species_Cary_Cor_v1", "Species_Malu_Spp_v1", "Species_Pice_Abi_v1",
        "Species_Pinu_Syl_v1", "Species_Jugl_Cin_v1", "Species_Acer_Sac_v1",
        "Species_Abie_Las_v1", "Species_Lari_Occ_v1", "Species_Pinu_Spp_v1",
        "Species_Sali_Spp_v1", "Species_Juni_Vir_v1", "Species_Abie_Ama_v1",
        "Species_Lari_Lya_v1", "Species_Acer_Mac_v1", "Species_Frax_Pen_v1",
        "Species_Acer_Neg_v1", "dev750")
    cn2 <- c(#"TPI", "TRI", "slope", "roughness",
        "AHM", "bFFP", "CMD", "DD_0",
        "DD_18", "DD18", "DD5", "eFFP", "EMT", "EXT", "FFP", "MAP", "MAT",
        "MCMT", "MSP", "MWMT", "NFFD", "PPT_sm", "PPT_wt", "SHM", "Tave_sm",
        "Tave_wt", "TD", "ROAD")
    cn <- if (clim_only)
        cn2 else cn1
    DAT <- data.frame(
        count=as.numeric(yy[, spp]),
        offset=off[, spp],
        cyid=dd$cyid,
        YEAR=dd$YEAR,
        ARU=dd$ARU, # ARU added here, but not as layer
        dd2[,cn]) # no cn subset here
    if (clim_only)
        DAT$count <- ifelse(DAT$count > 0, 1, 0)
    ## subsample based on 2.5x2.5km^2 cell x year units
    DAT <- DAT[sample.int(nrow(DAT)),]
    DAT <- DAT[!duplicated(DAT$cyid),]
    DAT <- DAT[sample.int(nrow(DAT), size=nmax, replace=b>1),]
    ## 0 detection output
    if (sum(DAT$count) < 1) {
        out <- structure(
            sprintf("0 detections for %s in %s", spp),
            class="try-error")
    } else {
        if (verbose)
            cat("\nFitting gbm::gbm for", spp, "/ b =", b, "/ n =",nrow(DAT), "\n")
        out <- try(gbm::gbm(DAT$count ~ . + offset(DAT$offset),
            data=DAT[,-(1:3)],
            n.trees = nt,
            interaction.depth = 3,
            shrinkage = 1/nt,
            bag.fraction = 0.5,
            distribution = if (clim_only) "bernoulli" else "poisson",
            var.monotone = NULL,
            keep.data = FALSE,
            verbose = verbose,
            n.cores = 1))
    }
    attr(out, "__settings__") <- list(
        species=spp, region="ALL", iteration=b, elapsed=proc.time()["elapsed"]-t0)
    out
}

## Create an array from the NODESLIST environnement variable
if (interactive()) {
    nodeslist <- 2
    setwd("d:/bam/BAM_data_v2019/gnm")
} else {
    cat("OK\n* Getting nodes list ... ")
    nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))
    cat("OK\n  Nodes list:\n")
    print(nodeslist)
}

cat("OK\n* Loading data on master ... ")
load(file.path("data", fn))
load(file.path("data", "xvinfo.RData"))
## get max tree size for each species to inform US models
tmp <- strsplit(names(xvinfo), "-")
x <- data.frame(sppbcr=names(xvinfo),
    spp=sapply(tmp, "[[", 1),
    bcr=sapply(tmp, "[[", 2),
    nt=sapply(xvinfo, "[[", "ntree"),
    vi=sapply(xvinfo, function(z)
        length(z$varimp[z$varimp > 0])))
x <- x[x$vi>0 & x$nt>0,]
x$pocc <- colMeans(yy>0)[as.character(x$spp)]
x$pbcr <- 0
for (i in levels(x$bcr)) {
    cm <- colMeans(yy[dd[[i]] > 0,]>0)
    x$pbcr[x$bcr == i] <- cm[as.character(x$spp[x$bcr == i])]
}
ntmax <- aggregate(x$nt, list(spp=x$spp), max)
ntmax <- structure(ntmax$x, names=as.character(ntmax$spp))

## Create the cluster with the nodes name.
## One process per count of node name.
## nodeslist = node1 node1 node2 node2, means we are starting 2 processes
## on node1, likewise on node2.
cat("* Spawning workers...")
cl <- makePSOCKcluster(nodeslist, type = "PSOCK")

cat("OK\n* Loading packages on workers ... ")
tmpcl <- clusterEvalQ(cl, library(mefa4))
tmpcl <- clusterEvalQ(cl, library(gbm))

cat("OK\n* Exporting and data loading on workers ... ")
if (interactive())
    tmpcl <- clusterEvalQ(cl, setwd("d:/bam/BAM_data_v2019/gnm"))
clusterExport(cl, c("dd", "dd2", "off", "yy", "CN",
    "PROJ", "xvinfo", "ntmax", "OUTDIR", "SPPBCR", "u",
    ".run_brt_boot", ".run_brt_boot_all"))

cat("OK\n* Establishing checkpoint ... ")
SPP <- colnames(yy)
DONE <- list.files(OUTDIR)
TOGO <- setdiff(SPP, DONE)

cat("OK\n* Start running models:")
set.seed(as.integer(Sys.time()))
ncl <- if (interactive())
    nodeslist else length(nodeslist)
#while (length(TOGO) > 0) {
for (counter in 1:10) { # run only 5 species (~10hrs)
    spp <- sample(TOGO, min(32, length(TOGO)))
    cat("\n  -", length(DONE), "done,", length(TOGO), "more to go, doing", spp, "on", date(), "... ")
    #res <- lapply(X=seq_len(ncl), fun=run_brt_boot, spp=spp)
#    parLapply(cl=cl, X=seq_len(ncl), fun=run_brt_boot, spp=spp)
    parLapply(cl=cl, X=spp, fun=run_brt_boot_all, clim_only=TRUE)
    DONE <- list.files(OUTDIR)
    TOGO <- setdiff(SPP, DONE)
}

## Releaseing resources.
cat("\n* Shutting down ... ")
stopCluster(cl)
cat("OK\nDONE!\n")
q("no")
