# Ring of Fire summaries
library(sf)
library(rgeos)
library(raster)
library(jsonlite)
library(ggplot2)


ROOT <- "d:/bam/BAM_data_v2019/gnm/artifacts"
OUT <- "~/repos/ring-of-fire/docs"

api_root <- "https://borealbirds.github.io/api/v4"
tab <- fromJSON(file.path(api_root, "species"))
rownames(tab) <- tab$id


# RoF boundary
poly <- st_read("~/GoogleWork/bam/RoF/boundary")

# Canadian LCC
lc <- raster("~/repos/recurring/offset/data/lcc.tif")

## species raster template
rt <- raster(file.path(ROOT, "ALFL/pred-ALFL-CAN-Mean.tif"))
poly <- st_transform(poly, st_crs(rt))

lc <- trim(mask(lc, poly))
rt <- trim(mask(rt, poly))


s <- !is.na(values(lc)) & !is.na(values(rt))
LC <- data.frame(
    lc=values(lc),      # land cover classes, integer
    area=1)             # km^2
## land cover classes
labs <- read.csv("~/repos/bamanalytics/lookup/lcc05.csv")
LC$label <- labs$BAMLCC05V2_label2[match(LC$lc, labs$lcc05v1_2)]
LC <- LC[s,]
table(LC$label)
LC$label <- droplevels(LC$label)


RES <- list()
SPP <- tab$id

#spp <- "ALFL"
for (spp in SPP) {
    NN <- numeric(32)
    DD <- matrix(0, 32, nlevels(LC$label))
    colnames(DD) <- levels(LC$label)

    rm <- raster(file.path(ROOT, spp, paste0("pred-", spp, "-CAN-Mean.tif")))
    rm <- trim(mask(rm, poly))
    rs <- raster(file.path(ROOT, spp, paste0("pred-", spp, "-CAN-SD.tif")))
    rs <- trim(mask(rs, poly))

    #i <- 1
    for (i in 1:32) {
        cat(spp, i, "\n")
        flush.console()

        ri <- raster(file.path(ROOT, spp, paste0("pred-", spp, "-CAN-boot-", i, ".tif")))
        ri <- trim(mask(ri, poly))

        Di <- values(ri)[s]
        Di[is.na(Di)] <- 0
        # million individuals with pair adjustment
        Ni <- sum(Di, na.rm=TRUE) * 100 * 2 / 10^6
        PSi <- aggregate(list(density=Di), list(landcover=LC$label), mean)

        NN[i] <- Ni
        DD[i,as.character(PSi$landcover)] <- PSi$density
    }

    if (!dir.exists(paste0(OUT, "/species/", spp)))
        dir.create(paste0(OUT, "/species/", spp))
    png(paste0(OUT, "/species/", spp, "/map.png"), width=800, height=400)
    op <- par(mfrow=c(1,2), mar=c(1,1,2,6))
    plot(rm, axes=FALSE, box=FALSE, col=hcl.colors(100, "Lajolla"),
        main="Mean Density (males/ha)")
    plot(rs, axes=FALSE, box=FALSE, col=hcl.colors(100, "Blue-Red"),
        main="Std. Dev.")
    par(op)
    dev.off()

    CI <- t(apply(DD, 2, quantile, c(0.05, 0.95)))
    DF <- data.frame(landcover=factor(colnames(DD),
            c("Conif", "Mixed", "Decid", "Grass", "Wet")),
        estimate=colMeans(DD),
        lower=CI[,1], upper=CI[,2])

    p <- ggplot(DF,
        aes(x=landcover, y=estimate)) +
        geom_bar(stat="identity", fill="#95B6C1") +
        coord_flip() +
        geom_errorbar(aes(ymin=lower, ymax=upper),
                      width=0.2, color="#105A73") +
        ylab("Mean Density (males/ha)") +
        xlab("Landcover") +
        theme_minimal()
    ggsave(paste0(OUT, "/species/", spp, "/density.png"), p)

    o <- list(
        abundance=list(
            estimate=mean(NN),
            lower=unname(quantile(NN, 0.05)),
            upper=unname(quantile(NN, 0.95)),
            boot=NN),
        density=list(
            landcover=as.character(DF$landcover),
            estimate=DF$estimate,
            lower=DF$lower,
            upper=DF$upper,
            boot=lapply(structure(
                c("Conif", "Mixed", "Decid", "Grass", "Wet"),
                names=c("Conif", "Mixed", "Decid", "Grass", "Wet")),
                function(z) DD[,z])))

    RES[[spp]] <- o

}

for (spp in SPP) {

    o <- RES[[spp]]
    o$species <- list(id=spp,
        english=tab[spp, "english"],
        french=tab[spp, "french"],
        scientific=tab[spp, "scientific"])
    o <- o[c(3,1,2)]
    RES[[spp]] <- o

    writeLines(
        toJSON(o, pretty=TRUE, auto_unbox = TRUE),
        paste0(OUT, "/species/", spp, "/index.json")
    )

}


fN <- function(o) {
    data.frame(Species_ID=o$species$id,
        English_name=o$species$english,
        French_name=o$species$french,
        Scientific_name=o$species$scientific,
        N_estimate=o$abundance$estimate,
        N_lower=o$abundance$lower,
        N_upper=o$abundance$upper,
        Unit="M individuals")
}
fD <- function(o) {
    data.frame(Species_ID=o$species$id,
        English_name=o$species$english,
        French_name=o$species$french,
        Scientific_name=o$species$scientific,
        Landcover=o$density$landcover,
        D_estimate=o$density$estimate,
        D_lower=o$density$lower,
        D_upper=o$density$upper,
        Unit="males per ha")
}

NNN <- do.call("rbind", lapply(RES, fN))
DDD <- do.call("rbind", lapply(RES, fD))

write.csv(NNN, row.names = FALSE,
    file=paste0(OUT, "/species/rof-abundance.csv"))
write.csv(DDD, row.names = FALSE,
    file=paste0(OUT, "/species/rof-density.csv"))

writeLines(toJSON(tab), paste0(OUT, "/species/index.json"))

save(RES, tab, NNN, DDD, file="~/GoogleWork/bam/RoF/rof-results.RData")

## writing pages

load("~/GoogleWork/bam/RoF/rof-results.RData")

TOC <-
'---
title: List of Species
---

Click on a species to see the summaries.
'

for (spp in SPP) {
    o <- RES[[spp]]
    NAM <- paste0(o$species$english, " &ndash; ", o$species$french,
        " (_", o$species$scientific, "_)")
    LINK <- paste0("{{ site.baseurl }}/species/", spp, "/")
    TOC <- c(TOC, paste0("- [", NAM, "](", LINK, ")"))
}
cat(TOC, sep="\n")
writeLines(TOC, paste0(OUT, "/species/index.md"))

library(whisker)

TMP <-
'---
title: {{ENG}} - {{FRE}}
subtitle: (_{{SCI}}_)
---

The population size of {{ENG}} in the Ring of Fire region was {{EST}} ({{LO}}, {{HI}}) million individuals based on the [BAM National Models](https://dx.doi.org/10.5281/zenodo.4018335).

# Maps

Mean density (males per ha) is the average of 32 bootstrap based prediction maps.

![Distribution map]({{B}}/species/{{ID}}/map.png)

# Density

Mean densities and bootstrap based confidence intervals by land cover type.
Density summaries were calculated using post-hoc binning ([BAM 2020](https://dx.doi.org/10.5281/zenodo.4018335)), we used the 2005 Canadian Land Cover layer.

![Density by land cover type]({{B}}/species/{{ID}}/density.png)

# Download

The following summaries are awailable for download (all species combined):

- Species abundances in the Ring of Fire region ([CSV file]({{B}}/species/rof-abundance.csv))
with estimates and 90% confidence interval (in million individuals)
- Species densities in land cover classes within the Ring of Fire region ([CSV file]({{B}}/species/rof-density.csv))
with estimates and 90% confidence interval (in males per ha)

'

SPP <- rownames(tab)
for (spp in SPP) {
    o <- RES[[spp]]

    LIST <- list(
        B="{{ site.baseurl }}",
        ID=o$species$id,
        ENG=o$species$english,
        FRE=o$species$french,
        SCI=o$species$scientific,
        EST=round(o$abundance$estimate, 3),
        LO=round(o$abundance$lower, 3),
        HI=round(o$abundance$upper, 3))

    writeLines(whisker.render(TMP, LIST),
        paste0(OUT, "/species/", spp, "/index.md"))

}
