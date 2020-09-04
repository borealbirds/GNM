---
title: "Applications"
output:
  md_document:
    preserve_yaml: true
---

The BAM National Model results are meant to be used in various
applications. We provide some R scripts here to facilitate the use of
the results.

Prerequisites
-------------

We will assume that you use R version &gt;= 3.6 with the following
packages installed: [raster](https://CRAN.R-project.org/package=raster),
[sf](https://CRAN.R-project.org/package=sf),
[jsonlite](https://CRAN.R-project.org/package=jsonlite),
[readxl](https://CRAN.R-project.org/package=readxl).

Working with the JSON API
-------------------------

The JSON API uses JSON as data exchange format. Let’s define the url for
the BAM v4 API:

    library(jsonlite)
    api_root <- "https://borealbirds.github.io/api/v4"

### Get the list of species

First, we need to get the list of species from the JSON API, so that we
know what species codes to use:

    library(jsonlite)
    api_root <- "https://borealbirds.github.io/api/v4"
    tab <- fromJSON(file.path(api_root, "species"))
    str(tab)

    ## 'data.frame':    143 obs. of  8 variables:
    ##  $ id        : chr  "ALFL" "AMCR" "AMGO" "AMPI" ...
    ##  $ idnext    : chr  "AMCR" "AMGO" "AMPI" "AMRE" ...
    ##  $ idprevious: chr  "YTVI" "ALFL" "AMCR" "AMGO" ...
    ##  $ scientific: chr  "Empidonax alnorum" "Corvus brachyrhynchos" "Spinus tristis" "Anthus rubescens" ...
    ##  $ english   : chr  "Alder Flycatcher" "American Crow" "American Goldfinch" "American Pipit" ...
    ##  $ french    : chr  "Moucherolle des aulnes" "Corneille d'Am&eacute;rique" "Chardonneret jaune" "Pipit d'Am&eacute;rique" ...
    ##  $ family    : chr  "Tyrannidae" "Corvidae" "Fringillidae" "Motacillidae" ...
    ##  $ show      : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...

    head(tab[,c("id", "english")])

    ##     id            english
    ## 1 ALFL   Alder Flycatcher
    ## 2 AMCR      American Crow
    ## 3 AMGO American Goldfinch
    ## 4 AMPI     American Pipit
    ## 5 AMRE  American Redstart
    ## 6 AMRO     American Robin

### Get estimates for a species

We can use the `id` column if we need to loop over multiple species. Now
we’ll only use one species, Canada Warbler (`"CAWA"`):

    spp <- "CAWA"
    results <- fromJSON(file.path(api_root, "species", spp))
    str(results, max.level=2)

    ## List of 3
    ##  $ species :List of 8
    ##   ..$ id        : chr "CAWA"
    ##   ..$ idnext    : chr "CCSP"
    ##   ..$ idprevious: chr "BWWA"
    ##   ..$ scientific: chr "Cardellina canadensis"
    ##   ..$ english   : chr "Canada Warbler"
    ##   ..$ french    : chr "Paruline du Canada"
    ##   ..$ family    : chr "Parulidae"
    ##   ..$ show      : logi TRUE
    ##  $ popsize :'data.frame':    12 obs. of  4 variables:
    ##   ..$ region   : chr [1:12] "Canada" "4 Northwestern Interior Forest" "5 Northern Pacific Rainforest" "6 Boreal Taiga Plains" ...
    ##   ..$ abundance:'data.frame':    12 obs. of  3 variables:
    ##   ..$ density  :'data.frame':    12 obs. of  3 variables:
    ##   ..$ areakmsq : num [1:12] 6.21 0.546 0.136 1.18 1.49 1.32 0.0552 0.372 0.444 0.373 ...
    ##  $ densplot:'data.frame':    12 obs. of  2 variables:
    ##   ..$ region: chr [1:12] "Canada" "4" "5" "6" ...
    ##   ..$ data  :'data.frame':   12 obs. of  4 variables:

The `species` element in the list contain the species info saw already
in the species `tab`le.

The `popsize` element contains population sizes, densities, and areas
for the various regions:

    do.call(cbind, results$popsize)

    ##                                     region abundance.estimate abundance.lower
    ## 1                                   Canada             4.8100          4.5900
    ## 2           4 Northwestern Interior Forest             0.3090          0.2070
    ## 3            5 Northern Pacific Rainforest             0.0014          0.0009
    ## 4                    6 Boreal Taiga Plains             1.0600          0.9430
    ## 5           7 Taiga Shield & Hudson Plains             0.2060          0.1560
    ## 6                 8 Boreal Softwood Shield             1.4600          1.3400
    ## 7                            9 Great Basin             0.0002          0.0000
    ## 8                      10 Northern Rockies             0.0353          0.0208
    ## 9                      11 Prairie Potholes             0.1750          0.1440
    ## 10           12 Boreal Hardwood Transition             1.0400          0.9920
    ## 11 13 Lower Great Lakes/St. Lawrence Plain             0.0986          0.0919
    ## 12             14 Atlantic Northern Forest             0.3880          0.3510
    ##    abundance.upper density.estimate density.lower density.upper areakmsq
    ## 1           5.2100           0.0077        0.0074        0.0084   6.2100
    ## 2           0.4010           0.0057        0.0038        0.0073   0.5460
    ## 3           0.0023           0.0001        0.0001        0.0002   0.1360
    ## 4           1.2700           0.0090        0.0080        0.0107   1.1800
    ## 5           0.2870           0.0014        0.0010        0.0019   1.4900
    ## 6           1.7000           0.0111        0.0102        0.0129   1.3200
    ## 7           0.0011           0.0000        0.0000        0.0002   0.0552
    ## 8           0.0534           0.0009        0.0006        0.0014   0.3720
    ## 9           0.2000           0.0040        0.0032        0.0045   0.4440
    ## 10          1.1100           0.0280        0.0266        0.0297   0.3730
    ## 11          0.1030           0.0097        0.0090        0.0101   0.1020
    ## 12          0.4320           0.0199        0.0180        0.0222   0.1950

The `densplot` element contains the regional landcover specific
densities:

    ## pick a region
    results$densplot$region

    ##  [1] "Canada" "4"      "5"      "6"      "7"      "8"      "9"      "10"    
    ##  [9] "11"     "12"     "13"     "14"

    region <- 1
    results$densplot$region[region]

    ## [1] "Canada"

    ## densities for this region
    as.data.frame(lapply(results$densplot$data, "[[", 1))

    ##        landcover estimate  lower  upper
    ## 1        Conifer   0.0066 0.0062 0.0073
    ## 2  Taiga Conifer   0.0032 0.0021 0.0040
    ## 3      Deciduous   0.0158 0.0150 0.0167
    ## 4      Mixedwood   0.0152 0.0146 0.0159
    ## 5          Shrub   0.0049 0.0043 0.0058
    ## 6          Grass   0.0050 0.0045 0.0053
    ## 7   Arctic Shrub   0.0014 0.0009 0.0020
    ## 8   Arctic Grass   0.0020 0.0014 0.0026
    ## 9        Wetland   0.0059 0.0051 0.0068
    ## 10      Cropland   0.0052 0.0044 0.0056

These results reflect population and density as males only (million
males, males/ha, respectively). To apply a pair adjustment, the numbers
have to be multiplied by 2.

### Density map images

The density map images can be accessed as
`https://borealbirds.github.io/api/v4/species/CAWA/images/mean-pred.png`
and
`https://borealbirds.github.io/api/v4/species/CAWA/images/mean-det.png`.

![](https://borealbirds.github.io/api/v4/species/CAWA/images/mean-pred.png)

![](https://borealbirds.github.io/api/v4/species/CAWA/images/mean-det.png)

Assessing results
-----------------

Accessing the `BAMv4-results-2020-02-20.xlsx` file gives us the
following tables (sheet names in parenthesis):

-   metadata explaining sheets and columns (metadata)
-   list species of species (species)
-   list of variables (variables)
-   variable importance (importance)
-   model validation results (validation)
-   population estimates (abundances)
-   lend cover based density estimates (densities)

Download the Excel file in a temporary file:

    library(readxl)
    tmp <- tempfile()
    download.file(file.path(api_root, "BAMv4-results-2020-02-20.xlsx"), tmp)

Now we can read in different sheets:

    vars <- read_xlsx(tmp, sheet="variables")
    str(vars)

    ## tibble [219 × 4] (S3: tbl_df/tbl/data.frame)
    ##  $ variable  : chr [1:219] "YEAR" "ARU" "AHM" "bFFP" ...
    ##  $ definition: chr [1:219] "Year of survey" "ARU (1) or human point count (0)" "Annual heat:moisture" "Beginning of the frost free period" ...
    ##  $ resolution: chr [1:219] NA NA "1 km" "1 km" ...
    ##  $ source    : chr [1:219] NA NA "Wang T., Hamann A., Spittlehouse D., & Carroll C. (2016) Locally Downscaled and Spatially Customizable Climate "| __truncated__ "Wang T., Hamann A., Spittlehouse D., & Carroll C. (2016) Locally Downscaled and Spatially Customizable Climate "| __truncated__ ...

Let’s read in all the tables and delete the temp file:

    sheets <- c(
        "metadata",
        "species",
        "variables",
        "importance",
        "validation",
        "abundances",
        "densities")
    tabs <- list()
    for (sheet in sheets)
        tabs[[sheet]] <- read_xlsx(tmp, sheet)
    unlink(tmp)

    hist(rnorm(1000))

![](https://borealbirds.github.io/GNM/applications/index_files/figure-markdown_strict/unnamed-chunk-9-1.png)
