library(mefa4)
Tab1 <- read.csv("~/repos/api/docs/v4/BAMv4-abundances-2020-02-20.csv")
Tab2 <- read.csv("~/repos/api/docs/v4/BAMv4-densities-2020-02-20.csv")
pif <- read.csv("~/GoogleWork/bam/PIF-AB/popBCR-CAN_v3_27-Feb-2019.csv")

compare_sets(Tab1$english, pif$English.Name)
pif$id1 <- as.character(Tab1$id[match(pif$English.Name, Tab1$english)])
pif$id2 <- as.character(Tab1$id[match(pif$Scientific.Name, Tab1$scientific)])
pif$id <- pif$id2
pif$id[is.na(pif$id2)] <- pif$id1[is.na(pif$id2)]
compare_sets(Tab1$id, pif$id)
pif <- pif[!is.na(pif$id) & pif$Country == "CAN",]

SPP <- sort(unique(pif$id))

pif$BCR <- as.factor(pif$BCR)
pif$N <- as.numeric(gsub(",", "", as.character(pif$Population.Estimate..unrounded.)))
summary(pif$N)
pif$N <- pif$N / pif$Pair.Adjust.Category
pif$N <- pif$N / 10^6

a1 <- aggregate(pif$N, list(id=pif$id, region=pif$Country), sum, na.rm=TRUE)
a2 <- aggregate(pif$N, list(id=pif$id, region=pif$BCR), sum, na.rm=TRUE)
head(a1)
head(a2)
a <- rbind(a1,a2)
summary(a)

Tab1$reg <- Tab1$region
levels(Tab1$reg) <- sapply(strsplit(levels(Tab1$reg), " "), function(z) z[1])

rownames(Tab1) <- paste0(Tab1$id, "_", Tab1$reg)
rownames(a) <- paste0(a$id, "_", a$region)
compare_sets(rownames(Tab1), rownames(a))

i <- intersect(rownames(Tab1), rownames(a))

x <- data.frame(id=Tab1[i,"id"], region=Tab1[i,"reg"],
    Npix=Tab1[i,"abundance_estimate"],
    Npif=a[i,"x"])
summary(x)
x$logNpix <- log(x$Npix)
x$logNpif <- log(x$Npif)
x$region <- droplevels(x$region)

plot(log(Npif) ~ log(Npix), x)
abline(0,1)
plot(Npif ~ Npix, x)
abline(0,1)

lm(Npix ~ Npif-1, x)
b <- list(Canada=coef(lm(Npix ~ Npif-1, x)))
for (i in levels(x$region))
    b[[paste("BCR", i)]] <- coef(lm(Npix ~ Npif-1, x[x$region == i,]))
bb <- unlist(b)
names(bb) <- names(b)

data.frame(b=bb)

boxplot(log(Npix/Npif) ~ region, x)

library(ggplot2)

x$logRatio <- log10(x$Npix/x$Npif)
x <- droplevels(x[!is.na(x$logRatio) & is.finite(x$logRatio),])
x2 <- x
levels(x2$region)[] <- "Canada"
x <- rbind(x, x2)

p <- ggplot(x, aes(x=region, y=logRatio)) +
    geom_violin(draw_quantiles=0.5, scale="width", fill="gold") +
    geom_hline(yintercept=0, col=2, lty=2) +
    geom_hline(yintercept=median(x$logRatio), col=2) +
    xlab("Bird Conservation Regions") +
    ylab(expression(log[10](N[PIX]/N[PIF]))) +
    theme_minimal() +
    theme(legend.position = "none")
ggsave("~/GoogleWork/bam/pif-pix-gnm.png", p, width=8, height=6)

## calculate BCR/prov level abundances

library(sp)
library(raster)
library(rgdal)

g <- readOGR("d:/bam/BAM_data_v2019/gnm/data/regions/BCR_BAMSubunits.shp")
g <- spTransform(g, "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")
g@data$BCRCode[g@data$BCRCode == 2] <- 200
g@data$BCRPROV <- paste(g@data$BCRCode, g@data$PROVINCE_S)

gc <- g[g@data$BCRCode < 200,]
table(gc@data$BCRPROV)

SPP <- list.dirs("d:/bam/BAM_data_v2019/gnm/artifacts",
    full.names = FALSE, recursive=FALSE)

N <- matrix(0, length(SPP), length(unique(gc@data$BCRPROV)))
dimnames(N) <- list(SPP, unique(gc@data$BCRPROV))

for (spp in SPP) {
    fn <- paste0("d:/bam/BAM_data_v2019/gnm/artifacts/",
        spp, "/pred-", spp, "-CAN-Mean.tif")
    r <- raster(fn)
    for (i in unique(gc@data$BCRPROV)) {
        cat(spp, i)
        flush.console()
        gci <- gc[gc@data$BCRPROV == i,]
        ri <- trim(mask(r, gci))
        N[spp, i] <- round(100 * sum(values(ri), na.rm=TRUE))
        cat("\t", round(N[spp, i]/10^6, 2), "\n")
    }
}

save(N, file="d:/bam/BAM_data_v2019/gnm/data/N-by-bcr-prov.RData")

write.csv(N, file="~/GoogleWork/bam/PIF-AB/N-by-bcr-prov-2020-03-16.csv")

## compare BCR/prov results: BAM v4 vs PIF v3

library(mefa4)
library(lhreg)
N <- read.csv("~/GoogleWork/bam/PIF-AB/N-by-bcr-prov-2020-03-16.csv")
Tab1 <- read.csv("~/repos/api/docs/v4/BAMv4-abundances-2020-02-20.csv")
Tab2 <- read.csv("~/repos/api/docs/v4/BAMv4-densities-2020-02-20.csv")
pif <- read.csv("~/GoogleWork/bam/PIF-AB/popBCR-CAN_v3_27-Feb-2019.csv")
lhr <- lhreg::lhreg_data

compare_sets(Tab1$english, pif$English.Name)
pif$id1 <- as.character(Tab1$id[match(pif$English.Name, Tab1$english)])
pif$id2 <- as.character(Tab1$id[match(pif$Scientific.Name, Tab1$scientific)])
pif$id <- pif$id2
pif$id[is.na(pif$id2)] <- pif$id1[is.na(pif$id2)]
compare_sets(Tab1$id, pif$id)
pif <- pif[!is.na(pif$id) & pif$Country == "CAN",]

SPP <- sort(unique(pif$id))

pif$BCR <- as.factor(pif$BCR)
pif$N <- as.numeric(gsub(",", "", as.character(pif$Population.Estimate..unrounded.)))
summary(pif$N)
pif$N <- pif$N / pif$Pair.Adjust.Category
pif$N <- pif$N / 10^6
#Lower.95..bound
#Upper.95..bound

x <- data.frame(
    species=pif$id,
    MDD=pif$Detection.Distance.Category..m.,
    Tadj=pif$Time.Adjust.Mean,
    PairAdj=pif$Pair.Adjust.Category,
    bcr=pif$BCR,
    prov=pif$Province...State...Territory,
    Npif=pif$N)
x <- droplevels(x[!is.na(x$Npif),])
x$bcr <- as.integer(as.character(x$bcr))
rownames(x) <- paste(x$species, x$bcr, x$prov, sep="-")

rownames(N) <- N[,1]
N[[1]] <- NULL
N <- as.matrix(N)
pix <- Melt(N)
pix$cols <- as.character(pix$cols)

cn <- data.frame(cn=colnames(N), stringsAsFactors = FALSE)
tmp <- strsplit(cn$cn, "\\.")
bcr <- sapply(tmp, "[[", 1)
bcr <- substr(bcr, 2, nchar(bcr))
cn$bcr <- as.integer(bcr)
prov <- sapply(tmp, function(z) paste(z[-1], collapse=" "))
pr <- c(
    ALBERTA = "AB",
    `BRITISH COLUMBIA` = "BC",
    MANITOBA = "MB",
    `NEW BRUNSWICK` = "NB",
    NEWFOUNDLAND = "NL",
    `NORTHWEST TERRITORIES` = "NT",
    `NOVA SCOTIA` = "NS",
    NUNAVUT = "NU",
    ONTARIO = "ON",
    `PRINCE EDWARD ISLAND` = "PE",
    QUEBEC = "QC",
    SASKATCHEWAN = "SK",
    YUKON = "YT")
cn$prov <- prov
cn$pr <- pr[match(prov, names(pr))]

pix <- data.frame(pix, cn[match(pix$cols, cn$cn),])
rownames(pix) <- with(pix, paste(rows, bcr, pr, sep="-"))

compare_sets(rownames(x), rownames(pix))
x$Npix <- pix$value[match(rownames(x), rownames(pix))]/10^6
x <- droplevels(x[!is.na(x$Npix),])

x$logNpix <- log(x$Npix)
x$logNpif <- log(x$Npif)

op <- par(mfrow=c(1,3))
plot(log(Npif) ~ log(Npix), x)
abline(0,1)
plot(Npif ~ Npix, x)
abline(0,1)
hist(x$logNpix-x$logNpif)
par(op)

summary(x$Npix/x$Npif)
exp(summary(x$logNpix-x$logNpif))

## species info

compare_sets(x$species, lhr$spp)
x <- droplevels(x[x$species %in% lhr$spp,])
x <- data.frame(x, lhr[match(x$species, lhr$spp),])
x$p3min <- 1-exp(-3*exp(x$logphi))
x <- data.frame(uid=rownames(x), x)

write.csv(x, row.names=FALSE,
    file="~/GoogleWork/bam/PIF-CAN/pifpix-can-estimates-2021-05-05.csv")


