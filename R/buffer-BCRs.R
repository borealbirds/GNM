library(sf)
w <-"G:/Boreal/NationalModelsV2/"
setwd(w)
bcr <- st_read("F:/GIS/basemaps/BCRs/bcrfinallcc.shp")

for (i in 4:14) {
  b <- bcr[bcr$BCR==i,]
  #br <- rasterize(b,r1k)
  bb <- st_buffer(b, dist=100000)
  bbb <- st_union(bb)
  st_write(bbb,paste("bcr",i,"_100km.shp",sep=""), delete_layer=TRUE)
}
