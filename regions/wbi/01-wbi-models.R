library(mefa4)
library(dismo)
library(gbm)
library(qs)
library(fastglm)
library(eflm)

qs::qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-25.qRData")

load("~/repos/GNM/regions/wbi/subsets4.RData")

SU <- list("WBAB"=60,
    "WBBC"=c(4,60),
    "WBYT"=4,
    "WBMB"=c(60,70,80),
    "WBSK"=c(60,80),
    "WBNT"=c(61, 70))
## polygon area in M kmsq
AA <- c(
    WBAB = 440142,
    WBBC = 283930,
    WBYT = 425409,
    WBMB = 542250,
    WBSK = 394355,
    WBNT = 993122
)
PP <- AA / sum(AA)

tmp <- table(dd$reg)/sum(table(dd$reg))
dx <- data.frame(Area=PP, Pts=as.numeric(tmp[names(PP)]))
dx$Ratio <- dx$Pts / dx$Area

get_data_all <- function(spp, replace=FALSE, nmax=10^4) {
    nn <- round(PP * nmax)
    cn0 <- NULL
    for (i in names(SU)) {
        regs <- paste0(spp, "-", SU[[i]])
        cn <- L[regs]
        cn <- sort(unique(unlist(cn)))
        cn0 <- c(cn0, cn)
    }
    cn0 <- sort(unique(cn0))
    DAT <- NULL
    for (i in names(SU)) {
        d <- get_data_by_reg(spp, i, cn=cn0, replace=replace,
            nmax=nmax)#nmax=nn[i])
        d$reg <- factor(i, levels(dd$reg))
        DAT <- rbind(DAT, d)
    }
    DAT
}

get_data_by_reg <- function(spp, reg, cn=NULL, replace=FALSE, nmax=10^4) {
    ss <- dd$reg == reg
    if (is.null(cn)) {
        regs <- paste0(spp, "-", SU[[reg]])
        cn <- L[regs]
        cn <- sort(unique(unlist(cn)))
    }
    cn <- unique(c("nalc", cn))
    cn <- cn[cn != "YEAR"]
    cn <- cn[cn != "ARU"]
    DAT <- data.frame(
        count=as.numeric(Y[ss, spp]),
        offset=O[ss, spp],
        cyid=dd$cyid[ss],
        YEAR=dd$YEAR[ss],
        ARU=dd$ARU[ss], # ARU added here, but not as layer
        ddvar[ss, cn, drop=FALSE])
    ## subsample based on 2.5x2.5km^2 cell x year units
    DAT <- DAT[sample.int(nrow(DAT)),]
    DAT <- DAT[!duplicated(DAT$cyid),]
    if (replace)
        DAT <- DAT[sample.int(nrow(DAT),  nrow(DAT), replace=TRUE),]
    if (nrow(DAT) > nmax)
        DAT <- DAT[seq_len(nmax),]
    DAT$cyid <- NULL
    DAT
}
simple_roc <- function(labels, scores){
    Labels <- labels[order(scores, decreasing=TRUE)]
    data.frame(
        TPR=cumsum(Labels)/sum(Labels),
        FPR=cumsum(!Labels)/sum(!Labels),
        Labels=Labels)
}
simple_auc <- function(ROC) {
    ROC$inv_spec <- 1-ROC$FPR
    dx <- diff(ROC$inv_spec)
    sum(dx * ROC$TPR[-1]) / sum(dx)
}

fit_fun <- function(i, spp, reg=NULL) {

    d <- if (is.null(reg)) {
        get_data_all(spp, replace=i>1)
    } else {
        get_data_by_reg(spp, reg, replace=i>1)
    }
    d <- droplevels(d)
    d$reg <- NULL
    if (nrow(d) < 1)
        return(structure("no data", class="try-error"))
    if (sum(d$count) < 1)
        return(structure("0 detections", class="try-error"))

    b <- try(gbm::gbm(d$count ~ . + offset(d$offset),
                data=d[,-(1:2)],
                n.trees = 10^4,
                interaction.depth = 3,
                shrinkage = 1/10^4,
                bag.fraction = 0.5,
                distribution = "poisson",
                var.monotone = NULL,#rep(0, length(4:ncol(DAT))),
                keep.data = FALSE,
                n.cores = 1))
    if (inherits(b, "try-error"))
        return(structure("gbm failed: err", class="try-error"))
    if (is.null(b))
        return(structure("gbm failed: null", class="try-error"))
    nd <- d
    logDgbm <- suppressWarnings(predict(b, newdata=nd, n.trees=b$n.trees))
    m <- glm(count ~ .-offset, data=d, family="poisson", offset=d$offset)
    logDglm <- model.matrix(m) %*% coef(m)

    s1 <- lm(logDgbm ~ nalc, data=d)
    s2 <- lm(logDgbm ~ ., data=d[,-(1:2)])
    p1 <- fitted(s1)
    p2 <- fitted(s2)

    v <- data.frame(y=d$count, off=d$offset,
        glm=logDglm, gbm=logDgbm, ssn=p1, ssa=p2)
    rownames(v) <- rownames(d)

    AUC <- apply(v, 2, function(z)
        simple_auc(simple_roc(ifelse(v$y>0, 1, 0), exp(z+v$off))))

    ## add here glm based smooth

    out <- list(i=i, spp=spp,
        reg=if (is.null(reg)) "ALL" else reg,
        vars=colnames(d)[-(1:2)],
        pk=rownames(v),
        #obs=v,
        AUC=AUC,
        gbm=b,
        glm=coef(m),
        ssn=coef(s1), ssa=coef(s2))
    class(out) <- "wb_fit"
    out
}

print.wb_fit <- function(x, ...) {
    cat("WB fit:", x$spp, "/", x$reg, "/", x$i, "\n")
    print(x$AUC)
    cat("\n")
    invisible(x)
}


dim(get_data_by_reg("OVEN", "WBAB"))
dim(get_data_by_reg("OVEN", "WBBC"))
dim(get_data_by_reg("OVEN", "WBMB"))
dim(get_data_by_reg("OVEN", "WBNT"))
dim(get_data_by_reg("OVEN", "WBSK"))
dim(get_data_by_reg("OVEN", "WBYT"))
dim(get_data_all("OVEN"))

## one run for each spp
## use all data
i <- 1
SPPx <- SPP[1:30]
for (spp in SPPx) {
    gc()
    cat(i, spp, "\n")
    flush.console()
    tmp <- try(fit_fun(i, spp, reg=NULL))
    if (inherits(tmp, "try-error"))
        tmp <- structure(as.character(tmp), class="try-error")
    RES <- tmp
    dir.create(paste0("d:/bam/2021/wbi/out/", spp))
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qsavem(RES, file=fn)
}

# 100 run for 1 spp
spp <- "OVEN"
for (i in 2:100) {
    gc()
    cat(i, spp, "\n")
    flush.console()
    tmp <- try(fit_fun(i, spp, reg=NULL))
    if (inherits(tmp, "try-error"))
        tmp <- structure(as.character(tmp), class="try-error")
    RES <- tmp
    dir.create(paste0("d:/bam/2021/wbi/out/", spp))
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qsavem(RES, file=fn)
}

## some boot
SPPx <- SPP[1:30]
SPPx <- SPP[31:60]
SPPx <- SPP[61:90]
SPPx <- SPP[91:117]

B <- 2:10
for (i in Bv) {
    for (spp in SPPx) {
        gc()
        if (!dir.exists(paste0("d:/bam/2021/wbi/out/", spp)))
            dir.create(paste0("d:/bam/2021/wbi/out/", spp))
        fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
        #if (!file.exists(fn)) {
        if (TRUE) {
            RES <- list()
            for (reg in names(SU)) {
                cat(i, spp, reg, "\n")
                flush.console()
                tmp <- try(fit_fun(i, spp, reg))
                if (inherits(tmp, "try-error"))
                    tmp <- structure(as.character(tmp), class="try-error")
                RES[[reg]] <- tmp
            }
            qsavem(RES,file=fn)
        } else {
            cat(i, spp, "- skipping\n")
        }
    }
}

