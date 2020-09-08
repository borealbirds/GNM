library(rmarkdown)
library(knitr)
library(magick)

## read in files
f <- list(
    "Introduction"="docs/index.md",
    "Methods"="docs/methods/index.md",
    "Results"="docs/results/index.md",
    "Applications"="docs/applications/index.Rmd",
    "Source Code"="docs/code/index.md"
)
tx <- lapply(f, readLines)
hd <- readLines("docs/_src/head.txt")

## remove yaml header
tx <- lapply(tx, function(x) {
    i <- which(startsWith(x, "---"))[2] + 1
    x[i:length(x)]
})

## remove 1st chunk from Rmd
x <- tx[["Applications"]]
i <- which(startsWith(x, "```"))[2] + 1
x <- x[i:length(x)]
tx[["Applications"]] <- x

## sticthing together
md <- c(
    hd,
    ""
)
for (i in 1:length(tx)) {
    md <- c(
        md,
        paste("#", names(tx)[i]),
        "",
        tx[[i]],
        ""
    )
}
md <- gsub("{{ site.baseurl }}", "https://borealbirds.github.io/GNM", md, fixed=TRUE)

## making images local for latex
j <- which(startsWith(md, "!["))
jj <- gsub(")", "", sapply(strsplit(md[j], "](", fixed=TRUE), "[[", 2), fixed=TRUE)
jjj <- paste0("docs/_src/", basename(jj))
jjjj <- paste0("![](", gsub(".svg", ".png", basename(jj), fixed=TRUE), ")")
for (i in 1:length(j)) {
    download.file(jj[i], jjj[i])
}
## convert svg to png
im <- image_read("docs/_src/dbylc-can.svg")
im <- image_convert(im, "png")
image_write(im, "docs/_src/dbylc-can.png")
md[j] <- jjjj

writeLines(md, "docs/_src/gnm-docs.Rmd")

rmarkdown::render("docs/_src/gnm-docs.Rmd")
file.copy("docs/_src/gnm-docs.pdf",
          "docs/bam-national-models-v40.pdf", overwrite=TRUE)

