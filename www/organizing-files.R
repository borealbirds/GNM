library(jsonlite)
library(mefa4)

dir <- "d:/bam/BAM_data_v2019/gnm/www/"
#dir <- "~/GoogleWork/tmp/"

out <- "~/repos/api/docs/v4/"

st <- read.csv("d:/bam/2020/speclist.csv")
rownames(st) <- st$species
spp_keep <- as.character(st$species)[st$drop==0]
spp_drop <- as.character(st$species)[st$drop==1]

s <- fromJSON(readLines("~/repos/borealbirds.github.io/static/search.json"))
writeLines(toJSON(s, pretty=TRUE), "~/repos/borealbirds.github.io/static/search.json")

## species info
if (FALSE) {
e <- new.env()
load("d:/bam/BAM_data_v2019/gnm/data/BAMdb-GNMsubset-2020-01-08.RData", envir=e)

sp <- data.frame(id=rownames(e$spt),
    e$spt[,c("English_Name", "Scientific_Name", "Family_Sci")])
rownames(sp) <- sp$id
for (i in 1:ncol(sp))
    sp[,i] <- as.character(sp[,i])
colnames(sp) <- c("id", "common", "scientific", "family")
sp$combo <- paste0(sp$common, " (", sp$scientific, ")")
sp <- sp[spp_keep,]

writeLines(toJSON(list(data=sp),
    pretty=FALSE, rownames=FALSE),
    paste0(out, "species/index.json"))
}

sp <- fromJSON(readLines(paste0(out, "species/index.json")))$data
rownames(sp) <- sp$id

#GRYE, LEYE, WISN, SOSA, SPSA, KILL, RUGR, MODO, ROPI, SPGR

lcl <- c("Temperate or sub-polar needleleaf forest"=1,
    "Sub-polar taiga needleleaf fores"=2,
    "Temperate or sub-polar broadleaf deciduous"=5,
    "Mixed Forest"=6,
    "Temperate or sub-polar shrubland"=8,
    "Temperate or sub-polar grassland"=10,
    "Sub-polar or polar shrubland-lichen-moss"=11,
    "Sub-polar or polar grassland-lichen-moss"=12,
    "Wetland"=14,
    "Cropland"=15)
lc <- c("Conifer"=1,
    "Taiga Conifer"=2,
    "Deciduous"=5,
    "Mixedwood"=6,
    "Shrub"=8,
    "Grass"=10,
    "Arctic Shrub"=11,
    "Arctic Grass"=12,
    "Wetland"=14,
    "Cropland"=15)

bc <- c(
    "Northwestern Interior Forest"=4,
    "Northern Pacific Rainforest"=5,
    "Boreal Taiga Plains"=6,
    "Taiga Shield And Hudson Plains"=7,
    "Boreal Softwood Shield"=8,
    "Great Basin"=9,
    "Northern Rockies"=10,
    "Prarie Potholes"=11,
    "Boreal Hardwood Transition"=12,
    "Lower Great Lakes/ St. Lawrence Plain"=13,
    "Atlantic Northern Forest"=14)


bcr <- read.csv("~/GoogleWork/tmp/bcrlandcov.csv")
bcr <- droplevels(bcr[bcr$BCR %in% bc,])
bcr$BCR_NALC <- as.factor(paste0(bcr$BCR, "_", bcr$nalc))

CN <- c("Area", "Nmean", "NQ5", "NQ50", "NQ95")
CN2 <- c("Area", "Nmean")
sp2 <- sp
sp2$previous <- c(NA, rownames(sp)[-nrow(sp2)])
sp2[["next"]] <- c(rownames(sp)[-1], NA)


spp <- "ALFL"

## figures
for (spp in rownames(sp2)) {

fmap1 <- paste0(dir, "/ds3/", spp, "_pred1km1.png")
fmap2 <- paste0(dir, "/ds3/", spp, "_pred1km2.png")
fplot <- paste0(dir, "/ds2/", spp, "_boxplot.png")


omap1 <- paste0(out, "species/", spp, "/map.png")
omap2 <- paste0(out, "species/", spp, "/mapwdet.png")
oplot <- paste0(out, "species/", spp, "/dbynalc.png")

if (!dir.exists(paste0(out, "species/", spp)))
    dir.create(paste0(out, "species/", spp))
file.copy(fmap1, omap1, overwrite=TRUE)
file.copy(fmap2, omap2, overwrite=TRUE)
file.copy(fplot, oplot, overwrite=TRUE)

}

## json data
for (spp in rownames(sp2)) {

fdat <- paste0(dir, "/ds2/", spp, "_densities.csv")


odat <- paste0(out, "species/", spp, "/index.json")

x <- read.csv(fdat)
x <- droplevels(x[x$BCR %in% 4:14,])
x$BCR_NALC <- factor(paste0(x$BCR, "_", x$nalc), levels(bcr$BCR_NALC))
x$NAME <- as.factor(names(bc))[match(x$BCR, bc)]
x$NALC <- names(lc)[match(x$nalc, lc)]
x$Area <- round(bcr$area[match(x$BCR_NALC, bcr$BCR_NALC)])
x$Nmean <- round(2 * 100 * x$mean * x$Area)
x$NQ5 <- round(2 * 100 * x$Q5 * x$Area)
x$NQ50 <- round(2 * 100 * x$Q50 * x$Area)
x$NQ95 <- round(2 * 100 * x$Q95 * x$Area)


CAN <- aggregate(x[,CN2], list(bcr=rep("Canada", nrow(x))), sum)
CAN$Dmean <- CAN$Nmean / (100 * CAN$Area)
CAN
BCR <- aggregate(x[,CN2], list(bcr=x$BCR), sum)
BCR$Dmean <- BCR$Nmean / (100 * BCR$Area)
BCR
LCC <- aggregate(x[,CN2], list(nalc=x$nalc), sum)
LCC$Dmean <- LCC$Nmean / (100 * LCC$Area)
LCC

nbybcr <- as.matrix(Xtab(Nmean ~ BCR + nalc, x))
abybcr <- as.matrix(Xtab(Area ~ BCR + nalc, x))
dbybcr <- nbybcr / (100 * abybcr)
dbybcr[is.na(dbybcr)] <- 0

table(x$nalc,x$BCR)

div <- 10^6
CAN$Area <- CAN$Area / div
CAN$Nmean <- CAN$Nmean / div
BCR$Area <- BCR$Area / div
BCR$Nmean <- BCR$Nmean / div
LCC$Area <- LCC$Area / div
LCC$Nmean <- LCC$Nmean / div


L <- list(data=list(
    species=as.list(sp2[spp,]),
    bcr=paste(bc, names(bc)),
    nalc=paste(lc, names(lc)),
    popsize=list(
        columns=c("BCR", "Area", "Nmean", "Dmean"),
        ntotal=c(0, unlist(CAN[,-1])),
        nbybcr=as.matrix(BCR)
    ),
    habitat=list(
        columns=colnames(dbybcr),
        rows=rownames(dbybcr),
        dtotal=as.list(LCC),
        dbynalc=dbybcr
    )
))


toJSON(L, pretty=TRUE, rownames=FALSE)

if (!dir.exists(paste0(out, "species/", spp)))
    dir.create(paste0(out, "species/", spp))
writeLines(toJSON(L, pretty=FALSE, rownames=FALSE, dataframe="rows"),
    odat)

}


## updating maps for species

st["SOSA", "mapnum2"] <- 6
st["SOSA", "mapnum1"] <- 5

tmp <- read.csv("~/repos/api/docs/v4/BAMv4-abundances-2020-02-20.csv")
tmp <- tmp[tmp$id %in% spp_keep,]
write.csv(tmp, row.names = FALSE, file="~/repos/api/docs/v4/BAMv4-abundances-2020-02-20.csv")

tmp <- read.csv("~/repos/api/docs/v4/BAMv4-densities-2020-02-20.csv")
tmp <- tmp[tmp$id %in% spp_keep,]
write.csv(tmp, row.names = FALSE, file="~/repos/api/docs/v4/BAMv4-densities-2020-02-20.csv")

for (spp in spp_keep) {

    unlink(paste0("~/repos/api/docs/v4/species/", spp, "/mean-pred.png"))
    unlink(paste0("~/repos/api/docs/v4/species/", spp, "/mean-det.png"))

    frpred <- paste0("d:/bam/2020/map-images/", spp, "_pred1km", st[spp, "mapnum1"], ".png")
    topred <- paste0("~/repos/api/docs/v4/species/", spp, "/images/mean-pred.png")
    file.copy(frpred, topred, overwrite=TRUE)

    frdet <- paste0("d:/bam/2020/map-images/", spp, "_pred1km", st[spp, "mapnum2"], ".png")
    todet <- paste0("~/repos/api/docs/v4/species/", spp, "/images/mean-det.png")
    file.copy(frdet, todet, overwrite=TRUE)
}



