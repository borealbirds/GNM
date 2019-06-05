#' ---
#' title: "Survey data for BAM Generalized National Models"
#' author: "Peter Solymos <solymos@ualberta.ca>"
#' date: "`r as.Date(Sys.time())`"
#' output: pdf_document
#' ---
#+ echo=FALSE
knitr::opts_chunk$set(eval=FALSE)
ROOT <- "d:/bam/BAM_data_v2019/gnm"
#par(las=1, scipen=999)
#'
#' This script pulls together the BAM pached db and the analysis bundles.
#'
#' # Preamble
library(mefa4)
library(intrval)
#'
#' # Pieces
#'
#' BAM BBS
e1 <- new.env()
load("d:/bam/Apr2016/out/data_package_2016-12-01.Rdata", envir=e1)
load("d:/bam/Apr2016/out/offsets-v3_2017-04-19.Rdata", envir=e1)
load("d:/bam/Apr2016/out/offsets-v3data_2016-12-01.Rdata", envir=e1)
names(e1)
#' BU (AB)
e2 <- new.env()
load("d:/bam/BAM_data_v2019/nwt/BU-nonABMI-offsets-2019-02-04.RData", envir=e2)
names(e2)
#' WildTrax (AB)
e3 <- new.env()
load("d:/bam/BAM_data_v2019/nwt/nwt-wildtrax-offsets-2019-01-16.RData", envir=e3)
names(e3)
#' Atlas updates
e4 <- new.env()
load("d:/bam/2018/atlas_data/atlas_data_processed-20181018.RData", envir=e4)
names(e4)
#'
#' # PKEY level attributes
#'
#' Useful attributes to keep
cn <- c("PKEY", "SS", "PCODE","DATE","DATI", "YEAR", "X", "Y",
    "TSSR", "JDAY", "TREE", "LCC4", "MAXDUR", "MAXDIS", "ROAD")
#' BAM+BBS: joining tables
dd1 <- data.frame(e1$PKEY, e1$SS[match(e1$PKEY$SS, e1$SS$SS),], e1$offdat[match(e1$PKEY$PKEY, e1$offdat$PKEY),])
dd1$DATI <- dd1$DATE
dd1$DATE <- as.Date(dd1$DATE)
setdiff(cn, colnames(dd1))
dd1 <- dd1[,cn]
#dd1 <- dd1[,c(cn, "BCR")]
#' BU data
setdiff(cn, colnames(e2$dd))
e2$dd$PCODE <- interaction("BU", e2$dd$ProjectID, sep="_")
e2$dd$ROAD <- 0
dd2 <- e2$dd[,cn]
#dd2$BCR <- NA
#' WildTrax
setdiff(cn, colnames(e3$dd))
#e3$dd$project.name <- make.names(e3$dd$project.name, unique = TRUE)
e3$dd$PCODE <- interaction("BU", e3$dd$project.name, sep="_")
e3$dd$ROAD <- 0
dd3 <- e3$dd[,cn]
#dd3$BCR <- NA
#' Atlas updates
dd4 <- data.frame(e4$PKEY, e4$SS[match(e4$PKEY$SS, e4$SS$SS),])
dd4$YEAR <- dd4$YearCollected
dd4$DATI <- dd4$DATE
dd4$DATE <- as.Date(dd4$DATE)
dd4$ROAD <- 0
setdiff(cn, colnames(dd4))
dd4 <- dd4[,cn]
#dd4$BCR <- NA
#' Binding tables together
dd <- rbind(dd1, dd2, dd3, dd4)
#' No duplicated PKEYS allowed
sum(duplicated(dd$PKEY))
dd <- nonDuplicated(dd, PKEY, TRUE)
#' Coordinates, dates, times are needed
dd <- dd[!is.na(dd$X),]
dd <- dd[!is.na(dd$DATI),]
#' Very few surveys outside of 1991-2018
dd <- dd[dd$YEAR <= 2018,]
dd <- dd[dd$YEAR >= 1991,]
#' Obviously wrong coordinates dropped
dd <- dd[dd$X < 0,]
dd <- dd[dd$Y > 30,]
#' Inspect day of year
quantile(dd$DATI$yday, c(0, 0.01, 0.025, 0.05, 0.95, 0.975, 0.99, 1))
dd <- dd[dd$DATI$yday %[]% c(150, 200),]
hist(dd$DATI$yday)
#' Inspect time of day
quantile(dd$DATI$hour, c(0, 0.01, 0.025, 0.05, 0.95, 0.975, 0.99, 1))
dd <- dd[dd$DATI$hour %[]% c(2, 10),]
hist(dd$DATI$hour)
#' Keep most important columns
dd <- dd[,c("PKEY", "SS", "PCODE","DATI","YEAR", "X", "Y","MAXDUR", "MAXDIS", "ROAD")]
#'
#' # Species data and offsets
#'
#' Union of species codes
SPP <- sort(Reduce(union, list(colnames(e1$OFF), colnames(e2$off), colnames(e3$off), colnames(e4$OFF))))
#' Species cross tables, appended with 0s for missing species
#' BAM+BBS
YY1 <- Xtab(ABUND ~ PKEY + SPECIES, e1$PCTBL)
compare_sets(colnames(YY1), SPP)
YY1 <- YY1[,SPP]
#' BU (function `f` appends the missing columns)
YY2 <- e2$y
compare_sets(colnames(YY2), SPP)
f <- function(y) {
    yx <- Matrix(0, nrow(y), length(setdiff(SPP, colnames(y))))
    colnames(yx) <- setdiff(SPP, colnames(y))
    cbind(y, yx)[,SPP]
}
YY2 <- f(YY2)
#' Windtrax
YY3 <- e3$y
compare_sets(colnames(YY3), SPP)
YY3 <- f(YY3)
#' New Atlas data: drop possible duplicates
YY4 <- e4$YY
compare_sets(rownames(YY4), rownames(YY1))
YY4 <- YY4[setdiff(rownames(YY4), rownames(YY1)),]
compare_sets(colnames(YY4), SPP)
YY4 <- f(YY4)
yy <- rbind(YY1, YY2, YY3, YY4)
#' Offsets
off <- rbind(e1$OFF[,SPP], f(e2$off), f(e3$off), f(e4$OFF))
#' Standardizing rows across tables
rn <- intersect(rownames(dd), rownames(yy))
yy <- yy[rn,]
dd <- droplevels(dd[rn,])
off <- off[rn,]
spt <- droplevels(nonDuplicated(e1$TAX, Species_ID, TRUE)[SPP,])
#' Saving patched objects
if (interactive())
    save(dd, yy, off, spt,
        file=file.path(ROOT, "data", "BAMdb-patched-2019-06-04.RData"))
