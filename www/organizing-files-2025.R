library(jsonlite)
library(mefa4)

# TODO
# - DONE update density and detection maps in /api
# - DONE update bcr/lcc level svg files in /api -- ISSUE with Canadian results
# - DONE update bcr level pop size json data blob in /api
# - DONE add new xlsx file to /api and link to it in html
# - DONE update gridsome site
# - update maps (maps missing)
# - update figures (include all LC classes even 0s)
# - update zenodo and github site with bugfix blurb

dir <- "/Users/Peter/Dropbox/collab/bam-gnm-v4/qpad-bugfix/"

sp <- openxlsx::read.xlsx(paste0(dir, "BAMv4-results-2025-05-21.xlsx"), 2) # species
pop <- openxlsx::read.xlsx(paste0(dir, "BAMv4-results-2025-05-21.xlsx"), 6) # density & population by region
drl <- openxlsx::read.xlsx(paste0(dir, "BAMv4-results-2025-05-21.xlsx"), 7) # density by region/lcc

drl$landcover[drl$landcover == "Arctic grass"] <- "Arctic Grass"
drl$landcover[drl$landcover == "Arctic shrub"] <- "Arctic Shrub"
drl$landcover[drl$landcover == "Taiga conifer"] <- "Taiga Conifer"
drl$reg2 <- sapply(strsplit(drl$region, " "), "[[", 1)

z <- read.csv(paste0(dir, "speclist_update.csv"))
z[z$mapnum2==5,]
z$mapnum1[z$species %in% c("CLSW", "SOSA")] <- 5
z$mapnum2[z$species == "SOSA"] <- 6
z$mapnum2[z$mapnum2 == 5] <- 6 # CLSW, LEYE, WISN
z[z$mapnum2==5,]
table(n1=z$mapnum1[z$drop==0], n2=z$mapnum2[z$drop==0])

# overwrite the density maps
# The naming system is: SPEC_pred1kmX, where SPEC = 4-letter code and X = map number.
# Maps 1, 2, and 5 are the versions without occurrence points and range map boundaries.
# Maps 3, 4, and 6 are the versions with occurrence points and range boundaries.
#They are coupled, with (1,3), (2,4), and (5,6) having the same legend breaks.
# The best map in terms of legend breaks varies by species, so please use the following 
# look-up table to identify which map number to use for which species (also attached, 
# and in the root of the “website” directory)

out <- "/Users/Peter/git/github.com/borealbirds/api/docs/v4/"
SPP <- sp$id

zzz <- NULL
for (spp in SPP) {

    mn1 <- z$mapnum1[z$species == spp]
    mn2 <- z$mapnum2[z$species == spp]
    if (spp %in% c("LEYE", "WISN"))
        mn2 <- "6b"
    f_in_1 <- paste0(dir, "NewMosaicApproach/", spp, "_pred1km", mn1, ".png")
    f_in_2 <- paste0(dir, "NewMosaicApproach/", spp, "_pred1km", mn2, ".png")
    if (!file.exists(f_in_1))
        cat(spp, "1\n")
    if (!file.exists(f_in_2))
        cat(spp, "2\n")
    if (!file.exists(f_in_1) || !file.exists(f_in_2))
        zzz <- c(zzz, spp)
}
zzz
z[z$species %in% zzz,]

for (spp in SPP) {

    message(spp)

    mn1 <- z$mapnum1[z$species == spp]
    mn2 <- z$mapnum2[z$species == spp]
    if (spp %in% c("LEYE", "WISN"))
        mn2 <- "6b"

    f_in_1 <- paste0(dir, "NewMosaicApproach/", spp, "_pred1km", mn1, ".png")
    f_in_2 <- paste0(dir, "NewMosaicApproach/", spp, "_pred1km", mn2, ".png")

    f_out_1 <- paste0(out, "species/", spp, "/images/mean-pred.png")
    f_out_2 <- paste0(out, "species/", spp, "/images/mean-det.png")

    file.copy(f_in_1, f_out_1, overwrite=TRUE)
    file.copy(f_in_2, f_out_2, overwrite=TRUE)

}

# update json blob

for (spp in SPP) {
    j <- fromJSON(readLines(paste0(out, "species/", spp, "/index.json")))
    s <- sp[sp$id == spp,]
    p <- pop[pop$id == spp,-(1:3)]
    l <- drl[drl$id == spp,-(1:3)]

    # popsize
    jp <- j$popsize
    jp <- jp[!(jp$region %in% c("6 Boreal Taiga Plains", "7 Taiga Shield & Hudson Plains", "8 Boreal Softwood Shield")),]
    p <- p[match(jp$region, p$region),]

    jp$abundance$estimate <- p$abundance_estimate
    jp$abundance$lower <- p$abundance_lower
    jp$abundance$upper <- p$abundance_upper

    jp$density$estimate <- p$density_estimate
    jp$density$lower <- p$density_lower
    jp$density$upper <- p$density_upper

    # land cover
    jl <- j$densplot
    jl <- jl[jl$region %in% unique(l$reg2),]
    for (k in 1:nrow(jl)) {
        r <- jl$region[k]
        for (i in seq_along(jl$data$landcover[[k]])) {
            lcc <- jl$data$landcover[[k]][i]
            ll <- l[l$reg2 == r & l$landcover == lcc,]
            if (nrow(ll) < 1) {
                jl$data$estimate[[k]][i] <- 0
                jl$data$lower[[k]][i] <- 0
                jl$data$upper[[k]][i] <- 0
            } else {
                jl$data$estimate[[k]][i] <- ll$estimate[1]
                jl$data$lower[[k]][i] <- ll$lower[1]
                jl$data$upper[[k]][i] <- ll$upper[1]
            }
            # if (is.na(jl$data$estimate[[k]][i])) {
            #     stop(paste(r, lcc))
            # }
        }
    }

    jj <- j
    jj$popsize <- jp
    jj$densplot <- jl
    jjj <- toJSON(jj)
    writeLines(jjj, paste0(out, "species/", spp, "/index.json"))

}

# update landcover to include CANADA as well
# rg <- unique(drl$reg2)
rg <- c("Can",
    "4", "5", "6-0", "6-1", "7-0", "7-1", "8-0", "8-1", "8-2", "8-3", 
    "9", "10", "11", "12", "13", "14")

for (spp in SPP) {
    j <- fromJSON(readLines(paste0(out, "species/", spp, "/index.json")))
    l <- drl[drl$id == spp,-(1:3)]

    # land cover
    jl <- j$densplot

    jl2 <- data.frame(region=rg, data=NA)
    q <- lapply(rg, \(r) {
        ll <- l[l$reg2 == r, c("landcover", "estimate", "lower", "upper")]
        ll <- ll[match(LC, ll$landcover),]
        ll$landcover <- LC
        ll[is.na(ll)] <- 0
        as.list(ll)
    })

    d2 <- data.frame(landcover = rep(NA, length(rg)), 
        estimate = rep(NA, length(rg)), 
        lower = rep(NA, length(rg)), 
        upper = rep(NA, length(rg)))
    d2$landcover <- lapply(q, \(x) x[["landcover"]])
    d2$estimate <- lapply(q, \(x) x[["estimate"]])
    d2$lower <- lapply(q, \(x) x[["lower"]])
    d2$upper <- lapply(q, \(x) x[["upper"]])
    jl2$data <- d2

    jj <- j
    jj$densplot <- jl2
    jjj <- toJSON(jj)
    writeLines(jjj, paste0(out, "species/", spp, "/index.json"))

}


# create svg's
LC <- c(
    "Deciduous", 
    "Mixedwood", 
    "Conifer", 
    "Taiga Conifer", 
    "Shrub", 
    "Arctic Shrub",
    "Grass", 
    "Arctic Grass", 
    "Cropland", 
    "Wetland")

for (spp in SPP) {

    l <- drl[drl$id == spp,]
    for (rv in unique(l$reg2)) {

        cat(spp, rv, "\n")

        v <- l[l$reg2 == rv,]
        v <- v[match(rev(LC), v$landcover),]
        # v0 <- Tab2[Tab2$id==spp,]

        est <- v$estimate
        names(est) <- rev(LC)

        g <- if (rv == "Canada") "can" else rv
        fout <- paste0(out, "species/", spp, "/images/dbylc-", g, ".svg")

        svg(fout, width=8, height=5)
        op <- par(mar=c(4,8,4,2))
        tck <- barplot(est, horiz=TRUE, las=1, main=v$region[v$reg2 == rv][1], xlab="Density (males/ha)",
            xlim=c(0, max(l[,c("estimate", "lower", "upper")])), col="#007a7c", border=NA)
        segments(v$lower, tck, est, lwd=3, lend=1, col="white")
        segments(est, tck, v$upper, lwd=3, lend=1, col="#007a7c")
        par(op)
        dev.off()
    }

    if (file.exists(paste0(out, "species/", spp, "/images/dbylc-6.svg")))
        unlink(paste0(out, "species/", spp, "/images/dbylc-6.svg"))
    if (file.exists(paste0(out, "species/", spp, "/images/dbylc-7.svg")))
        unlink(paste0(out, "species/", spp, "/images/dbylc-7.svg"))
    if (file.exists(paste0(out, "species/", spp, "/images/dbylc-8.svg")))
        unlink(paste0(out, "species/", spp, "/images/dbylc-8.svg"))
}



## -----------


# dir <- "d:/bam/BAM_data_v2019/gnm/www/"
#dir <- "~/GoogleWork/tmp/"

out <- "/Users/Peter/git/github.com/borealbirds/api/docs/v4/"

# st <- read.csv("d:/bam/2020/speclist.csv")
# rownames(st) <- st$species
# spp_keep <- as.character(st$species)[st$drop==0]
# spp_drop <- as.character(st$species)[st$drop==1]

# s <- fromJSON(readLines("~/repos/borealbirds.github.io/static/search.json"))
# writeLines(toJSON(s, pretty=TRUE), "~/repos/borealbirds.github.io/static/search.json")

# ## species info
# if (FALSE) {
# e <- new.env()
# load("d:/bam/BAM_data_v2019/gnm/data/BAMdb-GNMsubset-2020-01-08.RData", envir=e)

# sp <- data.frame(id=rownames(e$spt),
#     e$spt[,c("English_Name", "Scientific_Name", "Family_Sci")])
# rownames(sp) <- sp$id
# for (i in 1:ncol(sp))
#     sp[,i] <- as.character(sp[,i])
# colnames(sp) <- c("id", "common", "scientific", "family")
# sp$combo <- paste0(sp$common, " (", sp$scientific, ")")
# sp <- sp[spp_keep,]

# writeLines(toJSON(list(data=sp),
#     pretty=FALSE, rownames=FALSE),
#     paste0(out, "species/index.json"))
# }

sp <- fromJSON(readLines(paste0(out, "species/index.json")))
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

tmpl <- fromJSON(readLines(paste0(out, "species/ALFL/index.json")))


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



