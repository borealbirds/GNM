library(rmarkdown)
library(whisker)


tmp1 <- '### {{species}}

Boxplot showing expected (density with QPAD correction) vs. observed counts:

![](d:/bam/2021/wbi/out/{{spp}}/WB-{{spp}}-box-1.png)

Detections:

![](d:/bam/2021/wbi/out/{{spp}}/WB-{{spp}}-det.png)

Map:

![](d:/bam/2021/wbi/out/{{spp}}/WB-{{spp}}-pred-1.png)

Marginal effects for NALC:

![](d:/bam/2021/wbi/out/{{spp}}/WB-{{spp}}-avg-nalc-1.png)

Marginal effects for continuous predictors:

![](d:/bam/2021/wbi/out/{{spp}}/WB-{{spp}}-resp-1.png)

Bootstrap prediction concordance:

![](d:/bam/2021/wbi/out/{{spp}}/WB-{{spp}}-occc-10.png)


'

tmp0 <- readLines("regions/wbi/wbi-model-development.Rmd")

txt <- unname(c(tmp0, sapply(SPP, function(spp) {
    whisker.render(tmp1,
        list(spp=spp, species=tab[spp,"english"]))
})))
writeLines(txt, "regions/wbi/WBI-report.Rmd")

render("regions/wbi/WBI-report.Rmd")
