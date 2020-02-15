library(mefa4)
library(jsonlite)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
PROJ <- "run3"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2019-10-29.RData"))

str(spt)

SPP <- rownames(spt)

#{"id":"ALFL","common":"Alder Flycatcher",
#    "scientific":"Empidonax alnorum",
#    "family":"Tyrannidae",
#    "combo":"Alder Flycatcher (Empidonax alnorum)"}
fr <- as.character(spt$French_Name)
tools::showNonASCII(fr)
fr <- gsub("é", "&eacute;", fr)
fr <- gsub("É", "&Eacute;", fr)
fr <- gsub("è", "&egrave;", fr)
fr <- gsub("à", "&agrave;", fr)
fr <- gsub("ê", "&ecirc;", fr)
fr <- gsub("â", "&acirc;", fr)
fr <- gsub("ô", "&ocirc;", fr)
tools::showNonASCII(fr)

x <- data.frame(id=spt$Species_ID,
    idnext=c(SPP[-1], SPP[1]),
    idprevious=c(SPP[length(SPP)], SPP[-length(SPP)]),
    scientific=spt$Scientific_Name,
    english=spt$English_Name,
    #french=spt$French_Name,
    french=fr,
    family=spt$Family_Sci,
    show=rep(TRUE, nrow(spt)))
head(x)
writeLines(toJSON(x), "~/repos/borealbirds.github.io/src/data/species.json")
writeLines(toJSON(x), "~/repos/api/docs/v4/species/index.json")

xxx <- x
rownames(xxx) <- xxx$id
as.list(droplevels(xxx["ALFL",]))

vv <- data.frame(title=xxx$english, path=paste0("/species/", rownames(xxx)),
    summary=xxx$scientific)
writeLines(toJSON(vv), "~/repos/borealbirds.github.io/static/search.json")

## species stuff

#bcrs <- 4:14
#lccs <- c(1L, 2L, 5L, 6L, 8L, 10L, 11L, 12L, 14L, 15L)

LCCs <- c(
    "Conifer"=1,
    "Taiga Conifer"=2,
    "Deciduous"=5,
    "Mixedwood"=6,
    "Shrub"=8,
    "Grass"=10,
    "Arctic Shrub"=11,
    "Arctic Grass"=12,
    "Wetland"=14,
    "Cropland"=15)
BCRs0 <- c(
    "Arctic Plains & Mountains"=3,
    "Northwestern Interior Forest"=4,
    "Northern Pacific Rainforest"=5,
    "Boreal Taiga Plains"=6,
    "Taiga Shield & Hudson Plains"=7,
    "Boreal Softwood Shield"=8,
    "Great Basin"=9,
    "Northern Rockies"=10,
    "Prairie Potholes"=11,
    "Boreal Hardwood Transition"=12,
    "Lower Great Lakes/St. Lawrence Plain"=13,
    "Atlantic Northern Forest"=14)

fround <- function(x) {
    x <- round(x, 4)
    if (all(x > 0.1))
        x <- round(x, 3)
    if (all(x > 1))
        x <- round(x, 2)
    if (all(x > 10))
        x <- round(x, 1)
    if (all(x > 100))
        x <- round(x, 0)
    x
}
fabu <- function(tmp, atmp, reg, id) {
    list(
        #id = id,
        region = reg,
        abundance = list(
            estimate=fround(median(tmp)/10^6),
            lower=fround(unname(quantile(tmp, 0.05)/10^6)),
            upper=fround(unname(quantile(tmp, 0.95)/10^6))),
        density = list(
            estimate=fround(median(0.01*tmp/atmp)),
            lower=fround(unname(quantile(0.01*tmp/atmp, 0.05))),
            upper=fround(unname(quantile(0.01*tmp/atmp, 0.95)))),
        areakmsq=fround(atmp/10^6)
    )
}
fden <- function(tmp, atmp, reg, id) {
    v <- 0.01*tmp/atmp
    lcc <- names(LCCs)[match(rownames(tmp), as.character(LCCs))]
    q <- t(apply(v, 1, quantile, c(0.5, 0.05, 0.95)))
    list(
        #id = id,
        region = reg,
        data = list(
            landcover=lcc,
            estimate=fround(q[,1]),
            lower=fround(q[,2]),
            upper=fround(q[,3]))
    )
}

#spp <- "ALFL"
for (spp in SPP) {
cat(spp, "\n")
flush.console()
x <- read.csv(
    sprintf(
        "c:/Users/Peter Solymos/GoogleWork/bam/website/%s_densities.csv", spp))

#str(x)
#hist(x$mean)
#plot(x$mean, x$Q50)
#abline(0,1)

tmp <- substr(as.character(x$spec), 9, nchar(as.character(x$spec)))
tmp[tmp=="NA"] <- NA
tmp <- as.integer(as.character(tmp))
table(tmp,useNA="a")
x$iter <- tmp

#spp <- substr(as.character(x$spec)[1], 1, 4)
z <- data.frame(
    nalc=x$nalc,
    landcover=factor(names(LCCs)[match(x$nalc, LCCs)], names(LCCs)),
    bcr=x$BCR,
    region=factor(names(BCRs)[match(x$BCR, BCRs)], names(BCRs)),
    density=x$mean,
    areakmsq=x$area,
    run=x$iter)
z <- z[!is.na(z$run),]
z <- droplevels(z[z$bcr != 3,])
z$density[is.na(z$density)] <- 0
summary(z)
rownames(z) <- paste0("b", z$bcr, "_l", z$nalc, "_r", z$run)
head(z)
bl <- as.matrix(Xtab(~ bcr + nalc, z))
bl[bl > 0] <- 1

BCRs <- BCRs0[BCRs0 > 3]

## bcr x nalc x 32 density
D <- Xtab(density ~ bcr + nalc + run, z)
## bcr x nalc area
A <- Xtab(areakmsq ~ bcr + nalc + run, z)[[1]]
## bcr x nalc x 32 abundance
N <- list(D[[1]] * A * 100)
for (i in 2:32)
    N[[i]] <- D[[i]] * A * 100


## Total Canada
tmp <- sapply(N, sum)
atmp <- sum(A)
Total <- fabu(tmp, atmp, "Canada", spp)

#toJSON(Total,pretty=TRUE,auto_unbox=TRUE)

## BCR's
tmp <- sapply(N, rowSums)
atmp <- rowSums(A)
Bcr <- lapply(seq_along(atmp), function(i) {
    reg <- paste(names(atmp)[i], names(BCRs)[as.character(BCRs) == names(atmp)[i]])
    fabu(tmp[,i], atmp[i], reg, spp)
})

#toJSON(c(list(Total), Bcr), pretty=TRUE,auto_unbox=TRUE)

## Densities

## Total Canada
tmp <- sapply(N, colSums)
atmp <- colSums(A)
Total2 <- fden(tmp, atmp, "Canada", spp)
#toJSON(Total2,pretty=TRUE,auto_unbox=TRUE)

## BCR's
Bcr2 <- lapply(seq_along(BCRs), function(i) {
    reg <- as.character(BCRs[i])
    tmp <- sapply(1:32, function(i) as.matrix(N[[i]])[reg,])
    atmp <- A[reg,]
    s <- bl[reg,] > 0
    fden(tmp[s,], atmp[s], reg, spp)
})
#toJSON(c(list(Total2), Bcr2), pretty=TRUE,auto_unbox=TRUE)


out <- list(
    species=as.list(droplevels(xxx[spp,])),
    popsize=c(list(Total), Bcr),
    densplot=c(list(Total2), Bcr2))
if (!dir.exists(sprintf("~/repos/api/docs/v4/species/%s", spp)))
    dir.create(sprintf("~/repos/api/docs/v4/species/%s", spp))
writeLines(toJSON(out, auto_unbox=TRUE),
    sprintf("~/repos/api/docs/v4/species/%s/index.json", spp))

}


for (spp in SPP) {
cat(spp, "\n")
flush.console()
if (!dir.exists(sprintf("~/repos/api/docs/v4/species/%s/images", spp)))
    dir.create(sprintf("~/repos/api/docs/v4/species/%s/images", spp))
fin <- sprintf("~/GoogleWork/bam/website/%s_pred1km2.png", spp)
fout <- sprintf("~/repos/api/docs/v4/species/%s/images/mean-pred.png", spp)
file.copy(fin, fout)
}

