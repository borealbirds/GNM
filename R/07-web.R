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
#writeLines(toJSON(x), "~/repos/borealbirds.github.io/src/data/species.json")
#writeLines(toJSON(x), "~/repos/api/docs/v4/species/index.json")

xxx <- x
rownames(xxx) <- xxx$id
as.list(droplevels(xxx["ALFL",]))

vv <- data.frame(title=xxx$english, path=paste0("/species/", rownames(xxx)),
    summary=xxx$scientific)
#writeLines(toJSON(vv), "~/repos/borealbirds.github.io/static/search.json")

## species stuff

xxx <- jsonlite::fromJSON("https://borealbirds.github.io/api/v4/species/")
SPP <- as.character(xxx$id)
rownames(xxx) <- SPP

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
BCRs0x <- c(
    "Arctic Plains & Mountains"=3,
    "Northwestern Interior Forest"=4,
    "Northern Pacific Rainforest"=5,
    "Boreal Taiga Plains, South"=60,
    "Boreal Taiga Plains, North"=61,
    "Taiga Shield & Hudson Plains, West"=70,
    "Taiga Shield & Hudson Plains, East"=71,
    "Boreal Softwood Shield, West"=80,
    "Boreal Softwood Shield, Ontario"=81,
    "Boreal Softwood Shield, East"=82,
    "Boreal Softwood Shield, Newfoundland"=83,
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

RES <- list()
#spp <- "ALFL"
for (spp in SPP) {
    cat(spp, "\n")
    flush.console()
    x1 <- read.csv(
        sprintf(
            "d:/bam/BAM_data_v2019/gnm/data/summaries/%s_densities.csv", spp))
    x2 <- read.csv(
        sprintf(
            "d:/bam/BAM_data_v2019/gnm/data/summaries/%s_densities_subunit.csv", spp))

    OUT <- list()
    for (j in 1:2) {
        x <- if (j == 1)
            x1 else x2
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
            region=factor(names(BCRs0)[match(x$BCR, BCRs0)], names(BCRs0)),
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

        BCRs <- if (j == 1)
            BCRs0[BCRs0 > 3] else BCRs0x[!(BCRs0x %in% BCRs0)]

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
        tmp <- sapply(N, rowSums)[as.character(BCRs),]
        atmp <- rowSums(A)[as.character(BCRs)]
        Bcr <- lapply(seq_along(atmp), function(i) {
            reg <- paste(names(atmp)[i], names(BCRs)[as.character(BCRs) == names(atmp)[i]])
            if (j==2)
                reg <- paste0(substr(reg, 1, 1), "-", substr(reg, 2, nchar(reg)), collapse="")
            fabu(tmp[i,], atmp[i], reg, spp)
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
            if (j==2)
                reg <- paste0(substr(reg, 1, 1), "-", substr(reg, 2, nchar(reg)), collapse="")
            fden(tmp[s,], atmp[s], reg, spp)
        })
        #toJSON(c(list(Total2), Bcr2), pretty=TRUE,auto_unbox=TRUE)


        OUT[[j]] <- list(
            species=as.list(droplevels(xxx[spp,])),
            popsize=c(list(Total), Bcr),
            densplot=c(list(Total2), Bcr2))
    }

    out <- OUT[[1]]
    out$popsize <- c(out$popsize, OUT[[2]]$popsize[-1])
    out$densplot <- c(out$densplot, OUT[[2]]$densplot[-1])

    RES[[spp]] <- out
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
fin <- sprintf("~/GoogleWork/bam/website/map-images/%s_pred1km2.png", spp)
fout <- sprintf("~/repos/api/docs/v4/species/%s/images/mean-pred.png", spp)
file.copy(fin, fout, overwrite=TRUE)

fin <- sprintf("~/GoogleWork/bam/website/map-images/%s_pred1km3.png", spp)
fout <- sprintf("~/repos/api/docs/v4/species/%s/images/mean-det.png", spp)
file.copy(fin, fout, overwrite=TRUE)
}

## making density figure plots ~ grrrrr plotly/mapbox issues in gridsome !@#$*%

Tab1 <- NULL
Tab2 <- NULL
for (x in RES) {
#x <- RES[[1]]
    xx <- cbind(
        as.data.frame(x$species),
        do.call(rbind, lapply(x$popsize, as.data.frame)))
    xx$idnext <- xx$idprevious <- xx$french <- xx$show <- xx$family <- NULL
    Tab1 <- rbind(Tab1, xx)
    xx <- cbind(
        as.data.frame(x$species),
        do.call(rbind, lapply(x$densplot, as.data.frame)))
    xx$idnext <- xx$idprevious <- xx$french <- xx$show <- xx$family <- NULL
    Tab2 <- rbind(Tab2, xx)
}
colnames(Tab1) <- gsub("\\.", "_", colnames(Tab1))
colnames(Tab2) <- gsub("data\\.", "", colnames(Tab2))

levs <- c(
    "Canada"="Canada",
    "4 Northwestern Interior Forest"=4,
    "5 Northern Pacific Rainforest"=5,
    "6 Boreal Taiga Plains"=6,
    "7 Taiga Shield & Hudson Plains"=7,
    "8 Boreal Softwood Shield"=8,
    "9 Great Basin"=9,
    "10 Northern Rockies"=10,
    "11 Prairie Potholes"=11,
    "12 Boreal Hardwood Transition"=12,
    "13 Lower Great Lakes/St. Lawrence Plain"=13,
    "14 Atlantic Northern Forest"=14,
    "6-0 Boreal Taiga Plains, South"="6-0",
    "6-1 Boreal Taiga Plains, North"="6-1",
    "7-0 Taiga Shield & Hudson Plains, West"="7-0",
    "7-1 Taiga Shield & Hudson Plains, East"="7-1",
    "8-0 Boreal Softwood Shield, West"="8-0",
    "8-1 Boreal Softwood Shield, Ontario"="8-1",
    "8-2 Boreal Softwood Shield, East"="8-2",
    "8-3 Boreal Softwood Shield, Newfoundland"="8-3")


Tab2$region <- as.character(Tab2$region)
Tab2$region <- names(levs)[match(Tab2$region, levs)]
#Tab2$region[is.na(Tab2$region)] <- "Canada"
Tab2$region <- as.factor(Tab2$region)

write.csv(Tab1, row.names = FALSE, file="~/repos/api/docs/v4/BAMv4-abundances-2020-02-20.csv")
write.csv(Tab2, row.names = FALSE, file="~/repos/api/docs/v4/BAMv4-densities-2020-02-20.csv")

## saving csv files with population size estimates

#spp <- "ALFL"
r <- levels(Tab2$region)
#rv <- r[12]
for (spp in SPP) {
    for (rv in r) {
        cat(spp, rv, "\n")
        v <- Tab2[Tab2$id==spp & Tab2$region == rv,]
        v0 <- Tab2[Tab2$id==spp,]

        est <- v$estimate
        names(est) <- as.character(v$landcover)
        g <- strsplit(rv, " ")[[1]]
        g <- if (length(g) == 1)
            "can" else g[1]
        fout <- sprintf("~/repos/api/docs/v4/species/%s/images/dbylc-%s.svg", spp, g)
#        fout <- sprintf("dbylc-%s.svg", g)
        svg(fout, width=8, height=5)
        op <- par(mar=c(4,8,4,2))
        tck <- barplot(est, horiz=TRUE, las=1, main=rv, xlab="Density (males/ha)",
            xlim=c(0, max(v0[,6:8])), col="#007a7c", border=NA)
        segments(v$lower, tck, est, lwd=3, lend=1, col="white")
        segments(est, tck, v$upper, lwd=3, lend=1, col="#007a7c")
        par(op)
        dev.off()
    }
}

lf <- list()
for (spp in SPP) {
    lf[[spp]] <- list.files(
        sprintf("~/repos/api/docs/v4/species/%s/images/", spp))
}

table(sapply(lf, length))
names(lf)[sapply(lf, length) < 14]

for (spp in names(lf)[sapply(lf, length) < 14]) {
fin <- sprintf("~/GoogleWork/bam/website/map-images/%s_pred1km2.png", spp)
fout <- sprintf("~/repos/api/docs/v4/species/%s/images/mean-det.png", spp)
file.copy(fin, fout, overwrite=TRUE)
}



xxx <- jsonlite::fromJSON("https://borealbirds.github.io/api/v4/species/")
SPP <- as.character(xxx$id)
rownames(xxx) <- SPP

tmp <- list.files("~/repos/api/docs/v4/species/")
