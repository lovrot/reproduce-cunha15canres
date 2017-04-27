message("Downloading preprocessed data for GSE81954 from GEO")
destdir <- file.path("cache", "geo")
gse81954 <- getGEO("GSE81954", destdir = destdir, GSEMatrix = TRUE)[[1]]

stopifnot(validObject(gse81954))
