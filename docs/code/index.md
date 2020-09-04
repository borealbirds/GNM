---
title: Code
---

## Generating results

> The code for the BAM Generalized National Models (GNM) is hosted in the [borealbirds/GNM](https://github.com/borealbirds/GNM) GitHub (git) repository.

The code in the `/R` folder of the repository includes R scripts for
processing the data (observations, offsets, and predictors)
and mosaicing together the predictive maps pieces.

The code in the `/graham` folder contains
code to run on [Compute Canada](https://www.computecanada.ca/)'s
[Graham cluster](https://docs.computecanada.ca/wiki/Graham).
The code includes regional boosted regression models and predictions.

The `/www` folder of the repository includes scripts to summarize
the outputs and organize the results into a presentable format
hosted as part of the API repository.

## Storing the results

> The results from the BAM Generalized National Models (GNM) are hosted in the [borealbirds/api](https://github.com/borealbirds/api) GitHub (git) repository.

Assets are served from the `/docs` folder of the git master branch via [GitHub pages](https://pages.github.com/).

The species mean density raster layers are available in GeoTIFF format from [here](https://drive.google.com/drive/folders/1exWa6vfhGo1DNUL4ei2baDz77as7jYzY?usp=sharing).

## Presenting the results

> The results from the BAM Generalized National Models (GNM) are presented via GitHub pages based on the [borealbirds/borealbirds.github.io](https://github.com/borealbirds/borealbirds.github.io/) GitHub (git) repository.

The website's features include:

- Uses [Vue](https://vuejs.org/) and [Gridsome](https://gridsome.org/)
- [Tailwind CSS v1](https://tailwindcss.com) (with PurgeCSS)
- Based on [this](https://github.com/drehimself/gridsome-portfolio-starter) Gridsome template with light/dark theme.
- Search among species with [Fuse.js](https://fusejs.io) and [vue-fuse](https://github.com/shayneo/vue-fuse)
- 404 Page
- RSS Feed
- Sitemap in XML
- Google search crawling is allowed
- Comments via [Disqus](https://disqus.com/)
- Bird images from [Unsplash](https://unsplash.com/collections/9507373/birds) for the 404 page

To install the tools to build the website, follow these steps:

```
# Install Gridsome CLI tool
npm install --global @gridsome/cli

# Clone the repo
git clone -b dev https://github.com/borealbirds/borealbirds.github.io
cd borealbirds.github.io

# Install dependencies
npm install

# Run development server with hot reloading
gridsome develop

## now look at http://localhost:8080
```

To locally build and deploy, use the `_build.sh` script (you will need write access to the GitHub repo).

Automatic build & deployment enabled via [GitHub Actions](https://github.com/borealbirds/borealbirds.github.io/actions),
see setup in `.github/workflows/build.yml` file.
