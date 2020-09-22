---
title: Generalized National Models
subtitle: Documentation
---

Reliable information on species' population sizes, trends, habitat associations, and distributions is important for conservation and land-use planning, as well as status assessment and recovery planning for species at risk. However, the development of such estimates at a national scale is challenged by a variety of factors, including sparse data coverage in remote regions ([Stralberg et al. 2015](https://dx.doi.org/10.1890/13-2289.1)), differential habitat selection across large geographies ([Crosby et al. 2019](https://doi.org/10.1111/ddi.12991)), and variation in survey protocols ([Sólymos et al. 2013](https://doi.org/10.1111/2041-210X.12106)).

With these factors in mind, we developed a generalized analytical approach to model species density in relation to environmental covariates, using the Boreal Avian Modelling Project database of point-count surveys (through 2018) and widely available spatial predictors ([Cumming et al. 2010](https://doi.org/10.7939/R3Z31NW3X), [Barker et al. 2015](http://dx.doi.org/10.1002/wsb.567)). We developed separate models for each geographic region (bird conservation regions intersected by jurisdiction boundaries) based on covariates such as tree species biomass (local and landscape scale), forest age, topography, land use, and climate. We used machine learning to allow for variable interactions and non-linear responses while avoiding time-consuming species-by-species parameterization. We applied cross-validation to avoid overfitting and bootstrap resampling to estimate uncertainty associated with our density estimates.

## Contact

Please contact us if you have questions or suggestions via these channels:

* [borealbirds.ualberta.ca](https://borealbirds.ualberta.ca/)
* [Twitter](https://twitter.com/borealbirds)
* [GitHub](https://github.com/borealbirds)

## Citing the models

Please cite this document when using the BAM National Model results as:

Boreal Avian Modelling Project, 2020.
*BAM Generalized National Models Documentation, Version 4.0*. Available at 
[https://borealbirds.github.io/](https://borealbirds.github.io/).
DOI: [10.5281/zenodo.4018335](https://dx.doi.org/10.5281/zenodo.4018335).
