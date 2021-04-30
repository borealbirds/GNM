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

.plot_fun <- function(i, res, u) {
    j <- as.character(u$var[i])
    x <- plot.gbm(res, j,
        n.trees = res$n.trees,
        return.grid=TRUE,
        type="response",
        ylab=paste(res$rof_settings$spp, "density (males/ha)"),
        xlab=paste0(j, " (", round(u$rel.inf[i], 2), "%)"))
    colnames(x) <- c("x", "y")
    x$var <- paste0(j, " (", round(u$rel.inf[i], 2), "%)")
    attr(x, "out.attrs") <- NULL
    mo <- lm(y ~ x, x)
    s <- try(segmented(mo, npsi=2), silent = TRUE)
    if (inherits(s, "try-error"))
        s <- segmented(mo, npsi=1)
    x$segm <- predict(s)
    x$segm[x$segm <= 0] <- min(x$segm[x$segm > 0])
    x
}

plot_fun <- function(res, m=12, log=FALSE) {
    u <- rel_inf(res)
    xx <- do.call(rbind, lapply(seq_len(m), .plot_fun, res, u))
    xx$var <- factor(xx$var, unique(xx$var))
    if (log) {
        xx$y <- log(xx$y)
        xx$segm <- log(xx$segm)
    }
    p <- ggplot(xx, aes(x=x, y=y)) +
        geom_line() +
        geom_line(aes(x=x, y=segm), col=4) +
        facet_wrap(vars(var), scales="free_x") +
        ylab(paste(res$rof_settings$spp, if (log) "log" else "", "density (males/ha)")) +
        xlab("Predictor values") +
        theme_minimal()
}

glm_smooth <- function(res, log=TRUE) {
    pr <- suppressWarnings(predict.gbm(res, xi, res$n.trees,
        type=if (log) "link" else "response"))
    u <- rel_inf(res)
    vars <- as.character(u$var[u$rel.inf > 0])
    xxi <- xi[,vars,drop=FALSE]
    xxi2 <- xi2[,colnames(xi2) %in% vars,drop=FALSE]
    colnames(xxi2) <- paste0(colnames(xxi2), "2")
    xxi <- cbind(xxi, xxi2)
    m <- lm(pr ~ ., xxi, model=FALSE)
    m$brt <- list(pred = pr, log=log, vars=vars, influence=u)
    m
}
run_glm <- function(res) {
    spp <- res$rof_settings$spp
    i <- res$rof_settings$i
    si <- BB[,i]
    if (sum(y[si, spp]) < 1)
        return(structure(sprintf("0 detections for %s", spp), class="try-error"))
    xi <- data.frame(
        count=as.numeric(y[si, spp]),
        offset=off[si, spp],
        ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0)[si],
        xx2[si, cn2])
    xu <- apply(xi, 2, function(z) length(unique(z))/length(z))
    xi2 <- xi^2
    xi2 <- xi2[,xu>0.5]

    u <- rel_inf(res)
    vars <- c("count", "offset", as.character(u$var[u$rel.inf > 0]))
    xxi <- xi[,vars,drop=FALSE]
    xxi2 <- xi2[,colnames(xi2) %in% vars,drop=FALSE]
    colnames(xxi2) <- paste0(colnames(xxi2), "2")
    xxi <- cbind(xxi, xxi2)

    out <- try(glm(count ~ .-offset, data=xxi,
        offset = xi$offset,
        family = "poisson", model=FALSE, y=FALSE))
    if (!inherits(out, "try-error")) {
        out$rof_settings <- list(spp=spp, i=i)
        out$pr_brt <- suppressWarnings(predict.gbm(res, xi, res$n.trees,
            type="link"))
        X <- model.matrix(formula(out), xxi)[,names(coef(out))]
        out$pr_glm <- drop(X %*% coef(out))
    }
    out
}

run_all <- function(res, thr=5) {
    spp <- res$rof_settings$spp
    i <- res$rof_settings$i
    si <- BB[,i]
    if (sum(y[si, spp]) < 1)
        return(structure(sprintf("0 detections for %s", spp), class="try-error"))
    xiall <- data.frame(
        count=as.numeric(y[, spp]),
        offset=off[, spp],
        ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0),
        xx2[, cn2])
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

    ## BRT fit
    out <- list()
    out$pr_brt <- suppressWarnings(predict.gbm(res, xiall, res$n.trees, type="link"))

    ## GLM fit
    ff0 <- as.formula(paste0(
        "count ~ ",
        paste0(vars[-(1:2)], collapse=" + "), " + ",
        paste0(paste0("I(", utop, "^2)"), collapse=" + ")))
    m0 <- try(glm(ff0, data=xi,
        offset = xi$offset,
        family = "poisson", model=FALSE, y=FALSE))
    X <- model.matrix(formula(m0), xiall)[,names(coef(m0))]
    out$pr_glm <- drop(X %*% coef(m0))

    ## GAM fit
    ff <- as.formula(paste0(
        "count ~ ",
        paste0("s(", utop, ")", collapse=" + "), " + ",
        paste0(xtop, collapse=" + ")))
    m1 <- try(mgcv::gam(ff,
        data=xi,
        offset = xi$offset,
        family = "poisson"))
    out$pr_gam <- as.numeric(predict(m1, newdata=xiall, type="link"))

    ## GLM smoother
    ff2 <- as.formula(paste0(
        "out$pr_brt[si] ~ ",
        paste0(vars[-(1:2)], collapse=" + "), " + ",
        paste0(paste0("I(", utop, "^2)"), collapse=" + ")))
    m2 <- lm(ff2, xi, model=FALSE)
    out$smooth_glm <- predict(m2, newdata=xiall)
    ## GAM smoother
    ff3 <- as.formula(paste0(
        "out$pr_brt[si] ~ ",
        paste0("s(", utop, ")", collapse=" + "), " + ",
        paste0(xtop, collapse=" + ")))
    m3 <- mgcv::gam(ff3, data=xi)
    out$smooth_gam <- predict(m3, newdata=xiall)

    as.data.frame(out)
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

glm0 <- function(i, spp, fff) {
    s1 <- BB[,1]
    si <- sample(s1, replace=TRUE)
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

glm1 <- function(i, spp, fff) {
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

pred_glm1 <- function(cc, fff, x) {
    ff0 <- as.formula(gsub("count ", "", fff))
    X <- model.matrix(formula(ff0), x)[,rownames(cc)]
    drop(X %*% cc)
}


## BRT run #1 --> variable screening
spp <- "OVEN"
load(paste0("d:/bam/2021/rof/brt2-xv/", spp, ".RData"))
xiall <- data.frame(
    ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0),
    xx2[, cn2])
fff <- get_fff(res)

z=run_all(res)

#cc <- pbapply::pbsapply(1:250, glm1, spp=spp, fff=fff)
cc <- pbapply::pbsapply(1:100, glm0, spp=spp, fff=fff)
z$pr_glm2 <- rowMeans(pred_glm1(cc, fff, xiall))



# gof
i <- 1
si <- BB[,i]
nsi <- which(!(seq_len(nrow(xx1)) %in% si))
lamIn <- exp(z[si,] + off[si,spp])
lamOut <- exp(z[nsi,] + off[nsi,spp])
yIn <- lamIn
yIn[] <- as.numeric(y[si,spp])

yOut <- lamOut
yOut[] <- as.numeric(y[nsi,spp])

llIn <- dpois(as.matrix(yIn), as.matrix(lamIn), log=TRUE)
llOut <- dpois(as.matrix(yOut), as.matrix(lamOut), log=TRUE)

# delta deviance
dIn <- colSums(-2*llIn)-min(colSums(-2*llIn))
dOut <- colSums(-2*llOut)-min(colSums(-2*llOut))

# IC weights
round(exp(-dIn/2)/sum(exp(-dIn/2)), 4)
round(exp(-dOut/2)/sum(exp(-dOut/2)), 4)


library(GGally)


limitRange <- function(data, mapping, ...) {
    Mx <- max(data)
    ggplot(data = data, mapping = mapping, ...) +
        geom_point(alpha=0.1, ...) +
        geom_smooth(method = "lm", se = FALSE) +
        scale_y_continuous(limits = c(0, Mx)) +
        scale_x_continuous(limits = c(0, Mx)) +
        geom_abline(slope=1, intercept=0, col="grey")
}

for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    load(paste0("d:/bam/2021/rof/brt2-xv/", spp, ".RData"))
    if (inherits(res, "gbm")) {
        z=run_all(res)
        ez=exp(z)
        p <- ggpairs(ez,
            title=paste(spp, "density"),
            progress = FALSE,
            lower = list(continuous = limitRange)) +
            theme_minimal()
        ggsave(sprintf("d:/bam/2021/rof/brt2-xv-pred-mosaic/%s-glm2.png", spp), p)
    }
}


## GLM based on RelImp is not so terrible: can we use bootstrap?

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

glm1 <- function(i, spp, fff) {
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

pred_glm1 <- function(cc, fff, xi) {
    ff0 <- as.formula(gsub("count ", "", fff))
    X <- model.matrix(formula(ff0), xi)[,rownames(cc)]
    drop(X %*% cc)
}

spp <- "OVEN"
load(paste0("d:/bam/2021/rof/brt2-xv/", spp, ".RData"))
i <- 1
si <- BB[,i]
xi <- data.frame(
    ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0)[si],
    xx2[si, cn2])
fff <- get_fff(res)
cc <- pbapply::pbsapply(1:100, glm1, spp=spp, fff=fff)
pr_glms <- pred_glm1(cc, fff, xi)
pr_brt <- suppressWarnings(predict.gbm(res, xi, res$n.trees, type="link"))

plot(pr_brt, rowMeans(pr_glms))

plot(pr_brt, pr_glms[,1], col="grey")
points(pr_brt, rowMeans(pr_glms))

i <- 1
si <- BB[,i]
xi <- data.frame(
#        count=as.numeric(y[si, spp]),
#        offset=off[si, spp],
        ecozone=ifelse(xx1$ecozone=="hudson_plain", 1, 0)[si],
        xx2[si, cn2])
xu <- apply(xi, 2, function(z) length(unique(z))/length(z))
plot(apply(xi,2,function(z) mean(z/max(z))) ~ xu)
xi2 <- xi^2
xi2 <- xi2[,xu>0.5]

spp <- "OVEN"
load(paste0("d:/bam/2021/rof/brt2-xv/", spp, ".RData"))
p <- plot_fun(res, log=TRUE)
p

u <- rel_inf(res)
xx <- do.call(rbind, lapply(1:12, .plot_fun, res, u))

p <- ggplot(xx, aes(x=x, y=y)) +
        geom_line() +
        geom_line(aes(x=x, y=segm), col=4) +
        facet_wrap(vars(var), scales="free_x") +
        ylab(paste(res$rof_settings$spp, "density (males/ha)")) +
        xlab("Predictor values") +
        theme_minimal()

u <- rel_inf(res)
plot(c(0,cinfl) ~ c(0,q),u,type="l")
abline(h=c(80,90,95), lty=2)

o1 <- run_glm(res)
o2 <- run_gam(res)
plot(o1$pr_brt, o1$pr_glm)
abline(0,1,col=2)
plot(o2$pr_brt, o2$pr_gam)
abline(0,1,col=2)
plot(o1$pr_glm, o2$pr_gam)
abline(0,1,col=2)
sd(o1$pr_brt)
sd(o1$pr_glm)
sd(o2$pr_gam)


RIall <-NULL
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    load(paste0("d:/bam/2021/rof/brt2-xv/", spp, ".RData"))
    if (inherits(res, "gbm")) {

        u <- rel_inf(res)
        u$spp <- spp
        RIall <- rbind(RIall, u)
        p <- plot_fun(res, log=TRUE)
        ggsave(sprintf("d:/bam/2021/rof/brt2-xv-pred-mosaic/%s-effects12.png", spp), p)

        m <- glm_smooth(res, log=TRUE)
        png(sprintf("d:/bam/2021/rof/brt2-xv-pred-mosaic/%s-glm2.png", spp),
            width=800, height=400)
        op <- par(mfrow=c(1,2))
        ra <- range(fitted(m), m$brt$pred)
        plot(y=fitted(m), x=m$brt$pred, col="#00000022",pch=19,xlab="logD BRT", ylab="logD GLM",
            main=spp, xlim=ra, ylim=ra)
        abline(0,1,col=2)
        ra <- range(exp(fitted(m)), exp(m$brt$pred))
        plot(y=exp(fitted(m)), x=exp(m$brt$pred), col="#00000022",pch=19,xlab="D BRT",
            ylab="D GLM", xlim=ra, ylim=ra)
        abline(0,1,col=2)
        par(op)
        dev.off()
    }
}
write.csv(RIall, row.names=FALSE, file="d:/bam/2021/rof/SppBRTVarImp_v2.csv")



for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    load(paste0("d:/bam/2021/rof/brt2-xv/", spp, ".RData"))
    if (inherits(res, "gbm")) {

        png(sprintf("d:/bam/2021/rof/brt2-xv-pred-mosaic/%s-glm2.png", spp),
            width=800, height=800)
        op <- par(mfrow=c(2,2))

        m <- run_glm(res)
        ra <- range(m$pr_brt, m$pr_glm)
        plot(y=m$pr_glm, x=m$pr_brt, col="#00000022",pch=19,xlab="logD BRT", ylab="logD GLM",
            main=paste(spp, "(real GLM)"), xlim=ra, ylim=ra)
        abline(0,1,col=2)
        ra <- range(exp(m$pr_brt), exp(m$pr_glm))
        plot(y=exp(m$pr_glm), x=exp(m$pr_brt), col="#00000022",pch=19,xlab="D BRT",
            ylab="D GLM", xlim=ra, ylim=ra)
        abline(0,1,col=2)

        m <- glm_smooth(res, log=TRUE)
        ra <- range(fitted(m), m$brt$pred)
        plot(y=fitted(m), x=m$brt$pred, col="#00000022",pch=19,xlab="logD BRT", ylab="logD GLM",
            main=paste(spp, "(smooth)"), xlim=ra, ylim=ra)
        abline(0,1,col=2)
        ra <- range(exp(fitted(m)), exp(m$brt$pred))
        plot(y=exp(fitted(m)), x=exp(m$brt$pred), col="#00000022",pch=19,xlab="D BRT",
            ylab="D GLM", xlim=ra, ylim=ra)
        abline(0,1,col=2)

        par(op)
        dev.off()
    }
}


# make prediction for the whole data set
#pr <- predict.gbm(res, xi, res$n.trees, type="response")
pr <- predict.gbm(res, xi, res$n.trees, type="link")

fun <- function(v) {
    vars <- as.character(u$var[u$cinfl <= 100*v])
    m <- lm(pr ~ ., xi[,vars])
    c(v=v,
        cor=cor(fitted(m), pr),
        sig2=sigma(m)^2,
        R2=summary(m)$r.squared)
}

fun(0.5)

a=t(sapply(seq(0.5, 1, 0.05), fun))
plot(sig2 ~ v,a, type="l")
plot(R2 ~ v,a, type="l")


v <- 0.9
vars <- as.character(u$var[u$cinfl <= 100*v])

xxi <- xi[,vars]
xxi2 <- xi2[,colnames(xi2) %in% vars]
colnames(xxi2) <- paste0(colnames(xxi2), "2")
xxi <- cbind(xxi, xxi2)

m <- lm(pr ~ ., xxi)
plot(fitted(m), pr, col="#00000022",pch=19)
abline(0,1,col=2)
cor(fitted(m), pr)
sigma(m)^2
summary(m)$r.squared

par(mfrow=c(1,2))
plot(fitted(m), pr, col="#00000022",pch=19)
abline(0,1,col=2)
plot(exp(fitted(m)), exp(pr), col="#00000022",pch=19)
abline(0,1,col=2)
par(mfrow=c(1,1))

print(object.size(res), units="Mb")
print(object.size(m), units="Mb")
print(object.size(coef(m)), units="Mb")


RIall <-NULL
for (spp in SPP) {
    cat(spp, "\n")
    load(paste0("d:/bam/2021/rof/brt2-xv/", spp, ".RData"))
    if (inherits(res, "gbm")) {
        u <- rel_inf(res)
        u$spp <- spp
        RIall <- rbind(RIall, u)
        p <- plot_fun(res, log=TRUE)
        ggsave(sprintf("d:/bam/2021/rof/brt2-xv-pred-mosaic/%s-effects12.png", spp), p)
    }
}
write.csv(RIall, row.names=FALSE, file="d:/bam/2021/rof/SppBRTVarImp_v2.csv")

RIall$n0 <- ifelse(RIall$rel.inf > 0, 1, 0)
z <- xtabs(~var+n0,RIall)
z <- z[,"1"]/rowSums(z)
sort(z)
z["eskerpoint"]
which(names(sort(z)) == "eskerpoint")
# list species with importance matrics


z <- xtabs(rel.inf~var+n0,RIall)[,2]
z <- rev(sort(z))
o <- seq_along(z)/length(z)
l <- cumsum(z)/sum(z)
plot(o, l)
str(z[l<.8])

sort(names(z)[l<.95])
