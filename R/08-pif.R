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
levels(Tab1$reg) <- c("10", "11", "12", "13", "14", "4", "5", "6", "7", "8", "9", "CAN")

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

plot(log(Npif) ~ log(Npix), x)
abline(0,1)
plot(Npif ~ Npix, x)
abline(0,1)
