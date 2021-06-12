library(mefa4)
library(dismo)
library(gbm)
library(qs)
library(fastglm)
library(eflm)

qs::qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-11.qRData")

load("~/repos/GNM/regions/wbi/subsets4.RData")

SU <- list("WBAB"=60,
    "WBBC"=c(4,60),
    "WBYT"=4,
    "WBMT"=c(60,70,80),
    "WBSK"=c(60,80),
    "WBNT"=c(61, 70))


get_data_by_reg <- function(spp, reg, cn=NULL, replace=FALSE, nmax=10^4) {
    ss <- dd$reg == reg
    if (is.null(cn)) {
        regs <- paste0(spp, "-", SU[[reg]])
        cn <- L[regs]
        cn <- sort(unique(unlist(cn)))
    }
    cn <- cn[cn != "YEAR"]
    cn <- unique(c("nalc", cn))
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

fit_fun <- function(i, spp, reg) {

    d <- droplevels(get_data_by_reg(spp, reg, replace=i>1))
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

    out <- list(i=i, spp=spp, reg=reg, vars=colnames(d)[-(1:2)],
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

## one run for each spp
i <- 1
for (spp in SPP) {
    RES <- list()
    for (reg in names(SU)) {
        cat(i, spp, reg, "\n")
        flush.console()
        tmp <- try(fit_fun(i, spp, reg))
        if (inherits(tmp, "try-error"))
            tmp <- structure(as.character(tmp), class="try-error")
        RES[[reg]] <- tmp
    }
    dir.create(paste0("d:/bam/2021/wbi/out/", spp))
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qsave(RES,file=fn)
}

# 100 run for 1 spp
spp <- "OVEN"
dir.create(paste0("d:/bam/2021/wbi/out/", spp))
for (i in 2:100) {
    RES <- list()
    for (reg in names(SU)) {
        cat(i, spp, reg, "\n")
        flush.console()
        tmp <- try(fit_fun(i, spp, reg))
        if (inherits(tmp, "try-error"))
            tmp <- structure(as.character(tmp), class="try-error")
        RES[[reg]] <- tmp
    }
    fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
    qsave(RES,file=fn)
}

