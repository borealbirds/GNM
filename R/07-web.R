library(mefa4)
library(jsonlite)

ROOT <- "d:/bam/BAM_data_v2019/gnm"
PROJ <- "run3"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2019-10-29.RData"))

str(spt)

#{"id":"ALFL","common":"Alder Flycatcher",
#    "scientific":"Empidonax alnorum",
#    "family":"Tyrannidae",
#    "combo":"Alder Flycatcher (Empidonax alnorum)"}
x <- data.frame(id=spt$Species_ID,
    scientific=spt$Scientific_Name,
    english=spt$English_Name,
    french=spt$French_Name,
    #french=utils::URLencode(as.character(spt$French_Name)),
    family=spt$Family_Sci,
    show=rep(TRUE, nrow(spt)))
writeLines(toJSON(x), "~/repos/borealbirds.github.io/static/species.json")

x <- read.csv("d:/bam/BAM_data_v2019/gnm/AMRE_densities.csv")

str(x)
hist(x$mean)
plot(x$mean, x$Q50)
abline(0,1)

tmp <- strsplit(as.character(x$spec), "-")
spp <- substr(tmp[[1]][1], 1, 4)
x$iter <- sapply(tmp, function(z) if (length(x) < 2) NA else as.integer(z[2]))
z <- data.frame(lancover=x$nalc, region=x$BCR, density=x$mean, areakmsq=x$area, run=x$iter)
z <- z[!is.na(z$run),]
z$density[is.na(z$density)] <- 0
summary(z)
