message("Downloading preprocessed data for GSE48091 from GEO")
destdir <- file.path("cache", "geo")
gse48091 <- getGEO("GSE48091", destdir = destdir, GSEMatrix = TRUE)[[1]]

stopifnot(validObject(gse48091))
