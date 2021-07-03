library(mefa4)
library(dismo)
library(gbm)
library(qs)

## define the bootstrap run (1-10)
i <- 1

## define the species code (see the wbi-stats.csv for a list)
spp <- "OVEN"

## load the GBM object
## adjust the root path as needed
## the object to be loaded is called RES
fn <- paste0("d:/bam/2021/wbi/out/", spp, "/", "WB-", spp, "-ALL-", i, ".qRData")
qload(fn)

## the pkeys that were used to train the model
str(RES$pk)

## the variables that were used to train the model
## this is the same across bootstrap runs for the same species
RES$vars

## organize your predictors into a data frame
## see the columns in the training data, I am using that here
qload("d:/bam/2021/wbi/WBI-data_BAMv4v6-WBAreas_2021-06-30.qRData")

## call the predictor data frame pr (this is just an example)
pr <- cbind(dd, ddvar)[1:100,RES$vars]

## some data manipulation is needed
## make sure 'nalc' is factor with levels
## "1"  "2"  "5"  "6"  "8"  "10" "11" "12" "13" "14" "15" "16" "17"
pr$nalc <- factor(as.character(pr$nalc), levels(ddvar$nalc))
## unused levels will end up being NA, so might get dropped later
## densities corresponding to these NAs can be set to 0 later

## set yeat to 2011, year of the latest veg/biomass data used in training
pr$YEAR <- 2011

## set ARU to 0, this is how we account for methods differences
pr$ARU <- 0

## finally the prediction
estD <- predict(RES$gbm, pr, type="response")
## you can safely ignore the warning:
##  predict.gbm does not add the offset to the predicted values
## this is exactly what we want!
