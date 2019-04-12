library(jsonlite)
library(mefa4)

#dir <- "d:/bam/BAM_data_v2019/gnm/www/"
dir <- "~/GoogleWork/tmp/"

out <- "~/repos/api/v4/"

e <- new.env()
load("d:/bam/BAM_data_v2019/gnm/data/BAMdb-GNMsubset-2019-03-01.RData", envir=e)

sp <- data.frame(id=rownames(e$spt),
    e$spt[,c("English_Name", "Scientific_Name", "Family_Sci")])
rownames(sp) <- sp$id
for (i in 1:ncol(sp))
    sp[,i] <- as.character(sp[,i])
colnames(sp) <- c("id", "common", "scientific", "family")
sp$combo <- paste0(sp$common, " (", sp$scientific, ")")

writeLines(toJSON(list(data=sp),
    pretty=FALSE, rownames=FALSE),
    paste0(out, "species/index.json"))

sp <- fromJSON(readLines(paste0(out, "species/index.json")))$data
rownames(sp) <- sp$id

lc <- c("Temperate or sub-polar needleleaf forest"=1,
    "Sub-polar taiga needleleaf fores"=2,
    "Temperate or sub-polar broadleaf deciduous"=5,
    "Mixed Forest"=6,
    "Temperate or sub-polar shrubland"=8,
    "Temperate or sub-polar grassland"=10,
    "Sub-polar or polar shrubland-lichen-moss"=11,
    "Sub-polar or polar grassland-lichen-moss"=12,
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

#writeLines(toJSON(list(data=data.frame(id=as.character(bc), name=names(bc))),
#    pretty=FALSE, rownames=FALSE),
#    paste0(out, "bcr/index.json"))

#writeLines(toJSON(list(data=data.frame(id=as.character(lc), name=names(lc))),
#    pretty=TRUE, rownames=FALSE),
#    paste0(out, "nalc/index.json"))


#bcr <- read.csv("d:/bam/BAM_data_v2019/gnm/www/bcrlandcov.csv")
bcr <- read.csv("~/GoogleWork/tmp/bcrlandcov.csv")
bcr <- droplevels(bcr[bcr$BCR %in% bc,])
bcr$BCR_NALC <- as.factor(paste0(bcr$BCR, "_", bcr$nalc))
#A <- as.matrix(Xtab(area ~ BCR + nalc, bcr))

CN <- c("Area", "Nmean", "NQ5", "NQ50", "NQ95")
df <- function(x) {
    x$Dmean <- x$Nmean / x$Area
    x$DQ5 <- x$NQ5 / (100 * x$Area)
    x$DQ50 <- x$NQ50 / (100 * x$Area)
    x$DQ95 <- x$NQ95 / (100 * x$Area)
    x
}
sp2 <- sp
sp2$previous <- c(NA, rownames(sp)[-nrow(sp2)])
sp2[["next"]] <- c(rownames(sp)[-1], NA)


fl <- list.files(dir)

spp <- "AMCR"

fmap <- paste0(dir, "/ds/", spp, "_pred1km1.png")
fdat <- paste0(dir, "/ds/", spp, "_densities.csv")

omap <- paste0(out, "species/", spp, "/map.png")
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

CAN <- df(aggregate(x[,CN], list(id=x$spec), sum))
LCC <- df(aggregate(x[,CN], list(nalc=as.character(x$nalc)), sum))
BCR <- df(aggregate(x[,CN], list(bcr=as.character(x$BCR)), sum))

L <- list(data=list(
    species=as.list(sp2[spp,]),
    national=list(total=as.list(CAN),
    nalc=as.list(LCC),
    bcr=as.list(BCR))))
#toJSON(L, pretty=TRUE, rownames=FALSE)

if (!dir.exists(paste0(out, "species/", spp)))
    dir.create(paste0(out, "species/", spp))
file.copy(fmap, omap)
writeLines(toJSON(L, pretty=FALSE, rownames=FALSE),
    odat)


z <- x[!duplicated(x$BCR_NAME),]
z <- z[order(z$BCR),]

cat(paste(z$BCR, z$BCR_NAME), sep="\n")


"_boxplot.png"
"_densities.csv"
"_densityplot.png"
"_pred1km1.png"



