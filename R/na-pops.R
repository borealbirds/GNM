library(mefa4)
library(intrval)
e1 <- new.env()
load("d:/bam/Apr2016/out/data_package_2016-12-01.Rdata", envir=e1)
load("d:/bam/Apr2016/out/offsets-v3_2017-04-19.Rdata", envir=e1)
load("d:/bam/Apr2016/out/offsets-v3data_2016-12-01.Rdata", envir=e1)
names(e1)

PC <- e1$PCTBL
PK <- e1$PKEY[match(PC$PKEY, e1$PKEY$PKEY),]
SS <- e1$SS[match(PC$SS, e1$SS$SS),]
TX <- e1$TAX[match(PC$SPECIES, e1$TAX$Species_ID),]

DT <- data.frame(
        PC[,c(c("PCODE", "SS", "PKEY", "SPECIES", "ABUND", "dur", "dis"))],
        PK[,c("YEAR", "MONTH", "DAY", "HOUR", "MIN", "MAXDUR", "MAXDIS", "DURMETH", "DISMETH", "ROAD")],
        SS[,c("JURS", "COUNTRY", "TZONE", "BOREALLOC", "BCR")],
        TX[,c("English_Name", "French_Name", "Spanish_Name", "Scientific_Name")])
MN <- droplevels(DT[startsWith(as.character(DT$PCODE), "MN"), ])
table(DT$PCODE, DT$COUNTRY, useNA="a")
DT <- DT[is.na(DT$COUNTRY) | DT$COUNTRY=="CAN",]
DT <- DT[DT$PCODE != "BBS",]
DT <- DT[DT$PCODE != "MNNFB",]
DT <- DT[!is.na(DT$YEAR),]
DT <- droplevels(DT)
str(DT)
summary(DT)

## little bit of magic for dur/dis factor levels
load("s:/Peter/bam/NA-POPS/BAM-ltdur-ltdis.RData")
MN$dur <- factor(as.character(MN$dur), colnames(ltdur$x))
MN$dis <- factor(as.character(MN$dis), colnames(ltdis$x))
table(MN$dur, useNA="a")
table(MN$dis, useNA="a")

DT$dur <- as.character(DT$dur)
DT$dur[DT$dur=="10-10"] <- "6.66-10"
DT$dur <- factor(DT$dur, colnames(ltdur$x))
DT$dis <- factor(as.character(DT$dis), colnames(ltdis$x))
table(DT$dur, useNA="a")
table(DT$dis, useNA="a")

compare_sets(DT$dur, colnames(ltdur$x))
compare_sets(DT$dis, colnames(ltdis$x))
compare_sets(MN$dur, colnames(ltdur$x))
compare_sets(MN$dis, colnames(ltdis$x))

save(DT,file="s:/Peter/bam/NA-POPS/BAM-CAN-longform-2020-08-13.RData")
save(MN,file="s:/Peter/bam/NA-POPS/BAM-MN-longform-2020-08-13.RData")

## -- processing data --

library(mefa4)
library(detect)

ROOT <- "s:/Peter/bam/NA-POPS"
load(file.path(ROOT, "BAM-CAN-longform-2020-08-13.RData"))
load(file.path(ROOT, "BAM-MN-longform-2020-08-13.RData"))
load(file.path(ROOT, "BAM-ltdur-ltdis.RData"))


## change data here: MN or DT (all BAM Canada)
DAT <- DT

## a list of pkey x duration matrices, list elements for each species
YDU <- Xtab(ABUND ~ PKEY + dur + SPECIES, DAT)

## a list of pkey x distance matrices, list elements for each species
YDI <- Xtab(ABUND ~ PKEY + dis + SPECIES, DAT)

## pkey level attributes (space/time, no XY)
X <- nonDuplicated(
    DAT[,c("PCODE", "SS", "PKEY", "YEAR",
        "MONTH", "DAY", "HOUR", "MIN", "MAXDUR", "MAXDIS", "DURMETH",
        "DISMETH", "ROAD", "JURS", "COUNTRY", "TZONE", "BOREALLOC", "BCR")],
    PKEY,
    TRUE)
## species level attributes
Z <- nonDuplicated(
    DAT[,c("SPECIES", "English_Name", "French_Name", "Spanish_Name", "Scientific_Name")],
    SPECIES,
    TRUE)


spp <- "OVEN"

## removal

f <- Y | D ~ 1

Y0 <- as.matrix(YDU[[spp]])

## make sure that columns (intervals) match up
stopifnot(all(colnames(Y0) == colnames(ltdur$x)))
## make sure we are not missing any methodologies
stopifnot(length(setdiff(levels(X$DURMETH), rownames(ltdur$end))) == 0)

## interval end matrix
D <- ltdur$end[match(X$DURMETH, rownames(ltdur$end)),]
## exclude 0 sum and <1 interval rows
iob <- rowSums(Y0) > 0 & rowSums(!is.na(D)) > 1
X0 <- droplevels(X[iob,])
Y0 <- Y0[iob,]
D <- D[iob,]
n <- nrow(D)
## arranging counts into position
Yid <- ltdur$id[match(X0$DURMETH, rownames(ltdur$id)),]
Y <- matrix(NA, n, ncol(ltdur$id))
for (i in seq_len(n)) {
    w <- Yid[i,]
    w <- w[!is.na(w)]
    Y[i,seq_len(length(w))] <- Y0[i,w]
}
mod <- cmulti(f, X, type="rem")



## distance

f <- Y | D ~ 1

Y0 <- as.matrix(YDI[[spp]])

## make sure that columns (intervals) match up
stopifnot(all(colnames(Y0) == colnames(ltdis$x)))
## make sure we are not missing any methodologies
stopifnot(length(setdiff(levels(X$DisMETH), rownames(ltdis$end))) == 0)

## interval end matrix
D <- ltdis$end[match(X$DISMETH, rownames(ltdis$end)),]
## exclude 0 sum and <1 interval rows
iob <- rowSums(Y0) > 0 & rowSums(!is.na(D)) > 1
X0 <- droplevels(X[iob,])
Y0 <- Y0[iob,]
D <- D[iob,]
n <- nrow(D)
## arranging counts into position
Yid <- ltdis$id[match(X0$DISMETH, rownames(ltdis$id)),]
Y <- matrix(NA, n, ncol(ltdis$id))
for (i in seq_len(n)) {
    w <- Yid[i,]
    w <- w[!is.na(w)]
    Y[i,seq_len(length(w))] <- Y0[i,w]
}
mod <- cmulti(f, X, type="dis")



