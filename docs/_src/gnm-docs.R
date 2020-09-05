library(rmarkdown)
library(knitr)

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
jjj <- basename(jj) # paste0("docs/_src/", basename(jj))
jjjj <- paste0("![](", jjj, ")")
for (i in 1:length(j)) {
    download.file(jj[i], jjj[i])
}
md[j] <- jjjj


writeLines(md, "docs/_src/gnm-docs.Rmd")

