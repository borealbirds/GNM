library(mefa4)
library(gbm)
library(dismo)
library(ggplot2)
library(segmented)

load("d:/bam/2021/rof/BAMv6_RoFpackage_April20.RData")
cn2 <- c("eskerpoint",
    "agriculture_G750.O", "bedrock_G750.O", "biomass2015.ntems",
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

rel_inf <- function(res) {
    rel.inf <- relative.influence(res, res$n.trees)
    rel.inf[rel.inf < 0] <- 0
    i <- order(-rel.inf)
    rel.inf <- 100 * rel.inf/sum(rel.inf)
    out <- data.frame(
        var = factor(res$var.names[i], res$var.names[i]),
        rel.inf = rel.inf[i])
    attr(out, "n.trees") <- res$n.trees
    out$cinfl <- cumsum(out$rel.inf)
    out$q <- seq_len(nrow(out))/nrow(out)
    out
}


get_fff <- function(res, thr=5) {
    u <- rel_inf(res)
    vars <- c("count", "offset", as.character(u$var[u$rel.inf > 0]))

    itop <- max(2, max(which(u$rel.inf >= thr)))
    utop <- as.character(u$var[seq_len(itop)])
    xtop <- as.character(u$var[-seq_len(itop)])
    xtop <- xtop[xtop %in% vars]

    paste0(
        "count ~ ",
        paste0(vars[-(1:2)], collapse=" + "), " + ",
        paste0(paste0("I(", utop, "^2)"), collapse=" + "))
}


run_glm1 <- function(i, spp, res) {
    fff <- get_fff(res)
    si <- BB[,i]
    if (sum(y[si, spp]) < 1)
        return(structure(sprintf("0 detections for %s", spp), class="try-error"))
    xi <- data.frame(
        count=as.numeric(y[si, spp]),
        offset=off[si, spp],
        ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0)[si],
        xx2[si, cn2])
    coef(glm(as.formula(fff), data=xi,
        offset = xi$offset,
        family = "poisson", model=FALSE, y=FALSE))
}


run_gam1 <- function(i, spp, res) {
    thr <- 5
    si <- BB[,i]
    if (sum(y[si, spp]) < 1)
        return(structure(sprintf("0 detections for %s", spp), class="try-error"))
    xi <- data.frame(
        count=as.numeric(y[si, spp]),
        offset=off[si, spp],
        ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0)[si],
        xx2[si, cn2])
    u <- rel_inf(res)
    vars <- c("count", "offset", as.character(u$var[u$rel.inf > 0]))

    itop <- max(2, max(which(u$rel.inf >= thr)))
    utop <- as.character(u$var[seq_len(itop)])
    xtop <- as.character(u$var[-seq_len(itop)])
    xtop <- xtop[xtop %in% vars]

    ## GAM fit
    ff <- as.formula(paste0(
        "count ~ ",
        paste0("s(", utop, ")", collapse=" + "), " + ",
        paste0(xtop, collapse=" + ")))
    try(mgcv::gam(ff,
        data=xi,
        offset = xi$offset,
        family = "poisson"))
}

pred_glm1 <- function(cc, fff, x) {
    ff0 <- as.formula(gsub("count ", "", fff))
    X <- model.matrix(formula(ff0), x)[,rownames(cc)]
    drop(X %*% cc)
}

run_brt1 <- function(i, spp, res) {
     RATE=0.001
    u <- rel_inf(res)
    vars <- c("count", "offset", as.character(u$var[u$rel.inf > 0]))
    ntree <- res$n.trees
    si <- BB[,i]
    if (sum(y[si, spp]) < 1)
        return(structure(sprintf("0 detections for %s", spp), class="try-error"))
    xi <- data.frame(
        count=as.numeric(y[si, spp]),
        offset=off[si, spp],
        ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0)[si],
        xx2[si, cn2])[,vars]
    out <- try(gbm::gbm(xi$count ~ . + offset(xi$offset),
        data=xi[,-(1:2)],
        n.trees = ntree,
        interaction.depth = 3,
        shrinkage = RATE,
        bag.fraction = 0.5,
        distribution = "poisson",
        var.monotone = NULL,
        keep.data = FALSE,
        verbose = FALSE,
        n.cores = 1))
    if (!inherits(out, "try-error"))
        out$rof_settings <- list(RATE=RATE, spp=spp, i=i)
    out
}


do_all <- function(i, spp, res) {
    # fresh models
    m1 <- run_brt1(i, spp, res)
    m2 <- run_gam1(i, spp, res)
    m3 <- run_glm1(i, spp, res)
    # predictions
    ff0 <- as.formula(gsub("count ", "", get_fff(res)))
    x <- data.frame(
        count=as.numeric(y[, spp]),
        offset=off[, spp],
        ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0),
        xx2[, cn2])
    X <- model.matrix(formula(ff0), x)[,names(m3)]
    pr1 <- suppressWarnings(predict.gbm(m1, x, res$n.trees, type="link"))
    pr2 <- predict(m2, newdata=x, type="link")
    pr3 <- drop(X %*% m3)
    # smoothers
    thr <- 5
    u <- rel_inf(res)
    vars <- c("count", "offset", as.character(u$var[u$rel.inf > 0]))
    itop <- max(2, max(which(u$rel.inf >= thr)))
    utop <- as.character(u$var[seq_len(itop)])
    xtop <- as.character(u$var[-seq_len(itop)])
    xtop <- xtop[xtop %in% vars]
    si <- BB[,i]
    xi <- data.frame(
        count=as.numeric(y[si, spp]),
        offset=off[si, spp],
        ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0)[si],
        xx2[si, cn2])[,vars]
    ff4 <- as.formula(paste0(
        "pr1[si] ~ ",
        paste0(vars[-(1:2)], collapse=" + "), " + ",
        paste0(paste0("I(", utop, "^2)"), collapse=" + ")))
    m4 <- lm(ff4, xi, model=FALSE)
    pr4 <- predict(m4, newdata=x)
    ## GAM smoother
    ff5 <- as.formula(paste0(
        "pr1[si] ~ ",
        paste0("s(", utop, ")", collapse=" + "), " + ",
        paste0(xtop, collapse=" + ")))
    m5 <- mgcv::gam(ff5, data=xi)
    pr5 <- predict(m5, newdata=x)

    out <- cbind(brt=pr1, gam=pr2, glm=pr3, gams=pr5, glms=pr4)
    attr(out, "settings") <- list(i=i, spp=spp)
    out
}


## BRT run #1 --> variable screening
spp <- "OVEN"
load(paste0("d:/bam/2021/rof/brt2-xv/", spp, ".RData"))
xiall <- data.frame(
    ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0),
    xx2[, cn2])
pr0 <- suppressWarnings(predict.gbm(res, xiall, res$n.trees, type="link"))

B <- 100

Res <- array(NA, c(nrow(xx1), 5, B))
#i <- 1
for (i in 1:B) {
    cat(spp, i, "\n")
    flush.console()
    z <- try(do_all(i, spp, res))
    if (!inherits(z, "try-error"))
        Res[,,i] <- z
}

save(Res1, file=paste0("d:/bam/2021/rof/", spp, "-boot1.RData"))
save(Res2, file=paste0("d:/bam/2021/rof/", spp, "-boot2.RData"))


## calculate boot avg
## calculate deviance and weight
## use the xx1$dist == 0 subset

ps_brt <- z[,1]
ps_gam <- z[,2]
ps_glm <- z[,3]
ps_gams <- z[,4]
ps_glms <- z[,5]

# gof

load(paste0("d:/bam/2021/rof/", spp, "-boot1.RData"))
load(paste0("d:/bam/2021/rof/", spp, "-boot2.RData"))
ps_brt <- cbind(Res1[,1,], Res2[,1,])
ps_gam <- cbind(Res1[,2,], Res2[,2,])
ps_glm <- cbind(Res1[,3,], Res2[,3,])
ps_gams <- cbind(Res1[,4,], Res2[,4,])
ps_glms <- cbind(Res1[,5,], Res2[,5,])

#ps_brt <- rowMeans(Res[,1,], na.rm=TRUE)
#ps_gam <- rowMeans(Res[,2,], na.rm=TRUE)
#ps_glm <- rowMeans(Res[,3,], na.rm=TRUE)
#ps_gams <- rowMeans(Res[,4,], na.rm=TRUE)
#ps_glms <- rowMeans(Res[,5,], na.rm=TRUE)

    thr <- 5
    u <- rel_inf(res)
    vars <- c("count", "offset", as.character(u$var[u$rel.inf > 0]))
    itop <- max(2, max(which(u$rel.inf >= thr)))
    utop <- as.character(u$var[seq_len(itop)])
    xtop <- as.character(u$var[-seq_len(itop)])
    xtop <- xtop[xtop %in% vars]
    ff4 <- as.formula(paste0(
        "rowMeans(ps_brt) ~ ",
        paste0(vars[-(1:2)], collapse=" + "), " + ",
        paste0(paste0("I(", utop, "^2)"), collapse=" + ")))
    m4 <- lm(ff4, xiall, model=FALSE)
ps_glmss <- predict(m4, newdata=xiall)
    ## GAM smoother
    ff5 <- as.formula(paste0(
        "rowMeans(ps_brt) ~ ",
        paste0("s(", utop, ")", collapse=" + "), " + ",
        paste0(xtop, collapse=" + ")))
    m5 <- mgcv::gam(ff5, data=xiall)
ps_gamss <- predict(m5, newdata=xiall)



si <- which(xx1$dist == 0)
#si <- BB[,2]
#si <- 1:nrow(xx1)

mu <- cbind(
#    saturated=0,
#    null=off[si,spp],
    brtxv=pr0[si] + off[si,spp],
    brt=rowMeans(ps_brt[si,] + off[si,spp]),
    gam=rowMeans(ps_gam[si,]+ off[si,spp]),
    glm=rowMeans(ps_glm[si,]+ off[si,spp]),
    gams=rowMeans(ps_gams[si,]+ off[si,spp]),
    glms=rowMeans(ps_glms[si,]+ off[si,spp]))
mu <- cbind(mu, ens=rowMeans(mu[,c("brt", "gam", "glm")]))
#    gamss=ps_gamss[si],
#    glmss=ps_glmss[si])
lam <- exp(mu)
#lam[,"null"] <- mean(y[si,spp])
#lam[,"saturated"] <- y[si,spp]
yy <- lam
yy[] <- as.numeric(y[si,spp])
ll <- dpois(as.matrix(yy), as.matrix(lam), log=TRUE)

# delta deviance
sort(colSums(-2*ll)-min(colSums(-2*ll)))
dev <- colSums(-2*ll)-min(colSums(-2*ll))

# IC weights
round(exp(-dev/2)/sum(exp(-dev/2)), 4)

# Pseudo R2
ll0 <- sum(dpois(y[si,spp], mean(y[si,spp]), log=TRUE))
lls <- sum(dpois(y[si,spp], y[si,spp], log=TRUE))
llf <- colSums(ll)
R2 <- 1 - (lls - llf) / (lls - ll0)
R2
sort(R2)


plot(0:6, rep(1.5,7), type="b", col=1, xlab="observed count", ylab="predicted", ylim=c(0,12))
points(0:6, 0:6, col=2, type="b")
points(jitter(rep(0:6,100)), jitter(rpois(700, rep(0:6,100))), col=4)

print(plot(res,"G750Species_Popu_Spp_v1.grd"))
zz=.plot_fun(1,res,u)
plot(zz$x, zz$y, xlab="G750Species_Popu_Spp_v1.grd", ylab="OVEN density", type="l")
lines(zz$x, fitted(lm(y ~ x+I(x^2), zz)), col=2)
lines(zz$x, fitted(mgcv::gam(y ~ s(x), data=zz)), col=4)

dd=data.frame(
    Density=c(lam[,"brt"], lam[,"gam"], lam[,"glm"]),
    Counts=as.factor(y[si,spp]),
    Model=rep(c("BRT (36%)", "GAM (27%)", "GLM (26%)"), each=length(si)))
tab <- table(y[si,spp])
levels(dd$Counts) <- paste0(names(tab), "\n(n=", tab, ")")
ggplot(data=dd, aes(x=Counts, y=Density, fill=Model)) +
    geom_boxplot() +
    theme_minimal()


ia <- gbm.interactions(res)

gbm.perspec(res,
    which(colnames(res$gbm.call$dataframe)[-(1:2)]==rownames(u)[1]),
    which(colnames(res$gbm.call$dataframe)[-(1:2)]==rownames(u)[2]),
    theta=300)


library(GGally)


limitRange <- function(data, mapping, ...) {
    Mx <- max(data)
    ggplot(data = data, mapping = mapping, ...) +
        geom_point(alpha=0.1, ...) +
#        geom_smooth(method = "lm", se = FALSE) +
        scale_y_continuous(limits = c(0, Mx)) +
        scale_x_continuous(limits = c(0, Mx)) +
        geom_abline(slope=1, intercept=0, col="grey")
}

ggpairs(data.frame(exp(mu[,c("brt", "gams", "glms", "gam", "glm")])),
    title=paste(spp, "density"),
    progress = FALSE,
    lower = list(continuous = limitRange)) +
    theme_minimal()


## R2

xiall <- data.frame(
    ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0),
    xx2[, cn2])
B <- 100

si <- which(xx1$dist == 0)
#si <- BB[,2]
#si <- 1:nrow(xx1)


RR <- list()
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    load(paste0("d:/bam/2021/rof/brt2-xv/", spp, ".RData"))
    if (inherits(res, "gbm")) {
        load(paste0("d:/bam/2021/rof/brt2-xv/", spp, ".RData"))
        pr0 <- suppressWarnings(predict.gbm(res, xiall, res$n.trees, type="link"))
        ff0 <- as.formula(gsub("count ", "", get_fff(res)))
        x <- data.frame(
            count=as.numeric(y[, spp]),
            offset=off[, spp],
            xiall)
        #m3 <- run_glm1(i, spp, res)
        m3 <- try(pbapply::pbsapply(1:B, run_glm1, spp=spp, res=res))
        if (!inherits(m3, "try-error")) {
            X <- model.matrix(formula(ff0), x)[,rownames(m3)]
            pr3 <- X %*% m3

            mu <- cbind(
                brtxv=pr0[si] + off[si,spp],
                glm=rowMeans(pr3[si,]+ off[si,spp]))
            lam <- exp(mu)
            yy <- lam
            yy[] <- as.numeric(y[si,spp])
            ll <- dpois(as.matrix(yy), as.matrix(lam), log=TRUE)
            ll0 <- sum(dpois(y[si,spp], mean(y[si,spp]), log=TRUE))
            lls <- sum(dpois(y[si,spp], y[si,spp], log=TRUE))
            llf <- colSums(ll)
            R2 <- 1 - (lls - llf) / (lls - ll0)

            #RR[[spp]] <- list(dev=colSums(-2*ll), R2=R2, lam=cbind(y=y[si,spp], lam))
            RR[[spp]] <- list(dev=colSums(-2*ll), R2=R2)
        }

    }
}

r <- t(sapply(RR, "[[", "R2"))
r[r < 0] <- NA
summary(r[,"glm"]/r[,"brtxv"])

plot(r, xlim=c(0,1), ylim=c(0,1))
abline(0,1, lty=2)
abline(lm(glm ~ brtxv-1, data.frame(r)))
