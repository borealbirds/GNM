library(jsonlite)
library(mefa4)
library(lhreg)

## spp info
tab <- fromJSON("https://borealbirds.github.io/api/v4/species")
rownames(tab) <- tab$id

p <- read.csv("d:/bam/2021/gnm/PopEsts Global 2020.04.29.csv")
tab$pifspp <- tab$english
tab$pifspp[tab$pifspp == "Gray Jay"] <- "Canada Jay"
tab$pifspp[tab$pifspp == "Le Conte's Sparrow"] <- "LeConte's Sparrow"
## all other shorebirds are assumed to have pair adj = 2
compare_sets(tab$pifspp, p$Common.Name)
setdiff(tab$pifspp, p$Common.Name)
tab$pair <- p$Pair.Adjust.Category[match(tab$pifspp, p$Common.Name)]
tab$pair[is.na(tab$pair)] <- 2
table(tab$pair, useNA="a")
tab$tadj <- p$Time.Adjust.Mean[match(tab$pifspp, p$Common.Name)]
tab$mdd <- p$Detection.Distance.Category[match(tab$pifspp, p$Common.Name)]

lhr <- lhreg::lhreg_data
compare_sets(tab$id, lhr$spp)
tab <- data.frame(tab, lhr[match(tab$id, lhr$spp),])
tab$p3min <- 1-exp(-3*exp(tab$logphi))
tabs <- tab
rm(tab, p, lhr)
tabs <- droplevels(tabs[!is.na(tabs$mdd),])

## BAM estimates

# N=dmean*A*100*pair_adj
# D is vocalizing inds / ha
tabn <- read.csv("d:/bam/2021/gnm/BAM-GNM-Summaries-2021-03-01-CI95.csv")
#tabd <- read.csv("d:/bam/2021/gnm/BAM-GNM-Densities-2021-03-01-CI95.csv")

## PIF estimates

tabp <- read.csv("d:/bam/2021/gnm/PopEsts BCR x ProvState 2020.04.24.csv")
tabp <- droplevels(tabp[tabp$Country=="CAN", c("English.Name", "Scientific.Name", "Introduced",
    "BCR", "Province...State...Territory", "Country",
    "Population.Estimate..unrounded.",
    "Lower.80..bound..unrounded.", "Upper.80..bound..unrounded.",
    "Lower.95..bound..unrounded.", "Upper.95..bound..unrounded.")])
colnames(tabp) <- c("English.Name", "Scientific.Name", "Introduced",
    "BCR", "Jurs", "Country",
    "Npif",
    "NpifLo80", "NpifUp80",
    "NpifLo95", "NpifUp95")
for (i in c("Npif", "NpifLo80", "NpifUp80", "NpifLo95", "NpifUp95"))
    tabp[[i]] <- as.numeric(gsub(",", "", as.character(tabp[[i]])))

compare_sets(tabs$pifspp, tabp$English.Name)
tabp$id <- as.factor(tabs$id[match(tabp$English.Name, tabs$pifspp)])

compare_sets(tabn$SpeciesID, tabp$id)

table(tabn$BCR)
table(tabp$BCR)

table(tabn$JURS)
table(tabp$Jurs)

tabp$Jurs2 <- tabp$Jurs
levels(tabp$Jurs2) <- c("Alberta", "British Columbia", "Manitoba", "New Brunswick",
    "Newfoundland and Labrador",  "Nova Scotia", "Northwest Territories","Nunavut",
    "Ontario", "Prince Edward Island",
    "Quebec", "Saskatchewan", "Yukon Territory")
tabp$sbj <- paste0(tabp$id, " ", tabp$BCR, " ", tabp$Jurs2)

tabn <- tabn[!is.na(tabn$BCRJURS),]
tabn$sbj <- paste0(tabn$SpeciesID, " ", tabn$BCRJURS)

tabp <- tabp[!is.na(tabp$id),]

compare_sets(tabn$sbj, tabp$sbj)


## regional lc
LCCs <- c(
    "Conifer"=1,
    "Taiga Conifer"=2,
    "Deciduous"=5,
    "Mixedwood"=6,
    "Shrub"=8,
    "Grass"=10,
    "Arctic Shrub"=11,
    "Arctic Grass"=12,
    "Arctic Barren"=13,
    "Wetland"=14,
    "Cropland"=15,
    ## exclude ---
    "Barren Lands"=16, # 0
    "Urban and Built-up"=17, # ?
    "Water"=18, # 0
    "Snow and Ice"=19) # 0
lc <- read.csv("d:/bam/2021/gnm/BAM-GNM-NALCbyBCRJURSwithCI95-2021-05-05.csv")
lc$LCC <- as.factor(names(LCCs)[match(lc$nalc, LCCs)])
lc$bcrjurs <- paste(lc$bcr, lc$jurs)
lc$bcrjurs[is.na(lc$bcr) | is.na(lc$jurs)] <- NA
lc$bcrjurs <- as.factor(lc$bcrjurs)
lc$sbj <- paste(lc$spp, lc$bcr, lc$jurs)
lc$sbj[is.na(lc$bcrjurs)] <- NA
lc$sbj <- as.factor(lc$sbj)

if (FALSE) {
library(ggplot2)
z <- lc[lc$spp=="OVEN" &
        !is.na(lc$jurs) & lc$jurs=="Alberta" &
        !is.na(lc$bcr) & lc$bcr == 6,]
z <- lc[lc$spp=="OVEN" &
        is.na(lc$jurs) &
        !is.na(lc$bcr) & lc$bcr == 6,]
ggplot(data=z, aes(x=LCC, y=Mean)) + geom_col()
}

z1 <- droplevels(lc[lc$spp=="OVEN" &
        !is.na(lc$bcrjurs) &
        lc$bcr != 100,])
str(z1)
ncell <- as.matrix(mefa4::Xtab(ncell ~ bcrjurs + LCC, z1))
dcell <- as.matrix(mefa4::Xtab(Mean ~ sbj + LCC, lc))


compare_sets(tabn$sbj, rownames(dcell))

