# Explaining what, where, how

## Step 1: organizing training data

*Input*: various databases, flat files, previously bundled results
Input data sets include: BAM db, BBS db, WindTrax data dump, patches (individual projects not yet merged into BAM).

- project info (methodology)
- location info (long/lat): important for attribution when not year specific, often joined by survey table
- survey level info (location + date/time): these include revisits, yearly visits et to same location, this is the unit of analysis, a single survey event
- count table: this lists species abundances by survey ID, optionally subdivided into distance and time intervals in a removal sense (i.e. adding up the counts should be the total number of individuals counted by species)

How pivot tables are done:

- <https://peter.solymos.org/qpad-workshop/day1-2-data-processing.html>

Derive predictors.

*Output*: versioned (by date) data bundle ready for modeling. The goals is to organize the data in a common format, no matter what the inputs were.

- counts: PKEY x SPP matrix, cell values are the total counts
- offsets: PKEY x SPP matrix
- design variables: PKEY x P
- predictors: PKEY x CN data frame
- taxonomy (optional): SPP x P data frame (english, scientific, etc.)

PKEY is the unique ID for sampling events. SPP is the vector os species codes. CN is the column names for predictors.

I differentiate between predictor variables (these can change) and design variables. Design variables (long/lat, **year**/date/time, max duration/distance, **roadside**, **ARU**) are important to account for methods differences, i.e. used in models as fixed effects and to calculate offsets.

When you have only locations influencing attribution: you have lat/long.

When you have yearly attrib then you need to loop through the years.

### Species pivot tables

Get from long to wide format

Note about TMTT

### Offsets

<https://github.com/borealbirds/qpad-offsets>

## Step 2: BRT modeling

replication, local and HPC cluster

## Step 3: BRT model diagnostics

AUC, boxplots, OCCC, etc

## Step 4: organizing data for prediction

## Step 5: prediction

density & aggregate summaries

## Step 6: organizing results

this gets displayed on website

<https://github.com/borealbirds/mini-website-template>
