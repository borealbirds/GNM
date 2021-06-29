get_data_all <- function(spp, replace=FALSE, nmax=10^4, cn0=NULL) {
    #nn <- round(PP * nmax)
    if (is.null(cn0)) {
        cn0 <- NULL
        for (i in names(SU)) {
            regs <- paste0(spp, "-", SU[[i]])
            cn <- L[regs]
            cn <- sort(unique(unlist(cn)))
            cn0 <- c(cn0, cn)
        }
        cn0 <- sort(unique(cn0))
    }
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
    ## need to deal with TMTT
    yyy <- as.numeric(Y[, spp])
    q <- ceiling(quantile(yyy[yyy>0 & yyy<100], 0.99))
    DAT$count[DAT$count>q] <- q
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

fit_fun <- function(i, spp, reg=NULL, cn=NULL) {
    if (is.null(cn))
        cn <- unique(c("ROAD", names(DAIC[[spp]][1:10])))

    d <- if (is.null(reg)) {
        get_data_all(spp, replace=i>1, cn0=cn)
    } else {
        get_data_by_reg(spp, reg, replace=i>1, cn=cn)
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

recast_pred <- function(pr) {
    rs <- r
    u <- pr
    q <- quantile(u, 0.9999, na.rm=TRUE)
    u[!is.na(u) & u > q] <- q
    u <- u[match(1:length(values(rs)),prd$index)]
    values(rs) <- u
    rs
}

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

plot_fun <- function(res, log=FALSE) {
    u <- rel_inf(res)
    u <- u[rownames(u)!="nalc",]
    m <- nrow(u)
    xx <- do.call(rbind, lapply(seq_len(m), .plot_fun, res, u))
    xx$var <- factor(xx$var, unique(xx$var))
    if (log) {
        xx$y <- log(xx$y)
        xx$segm <- log(xx$segm)
    }
    p <- ggplot(xx, aes(x=x, y=y)) +
        geom_line() +
        geom_line(aes(x=x, y=segm), col=4) +
#        facet_wrap(vars(var), scales="free_x") +
        facet_wrap(vars(var), scales="free") +
        ylab(paste(res$rof_settings$spp, if (log) "log" else "", "density (males/ha)")) +
        xlab("Predictor values") +
        theme_minimal()
    p
}
