## Extract characteristics and append
for (id in c("1", "1.1", "1.2", "1.3")) {
  tmp <- pData(gse81954)[, paste0("characteristics_ch", id)]
  varnam <- unique(unlist(lapply(strsplit(tmp, ": "), function(x) x[1])))
  varval <- unlist(lapply(strsplit(tmp, ": "), function(x) x[2]))
  pData(gse81954)[, varnam] <- varval
}
varLabels(gse81954)[varLabels(gse81954) == "case-control status"] <- "casecontstat"

## Munge phenotypic data
pData(gse81954) <- pData(gse81954) %>%
  as_tibble() %>%
  mutate(
    tumorid = ID,
    piece = as.integer(piece),
    setnr = as.integer(setnr),
    casecontstat = as.integer(casecontstat),
    casecontcd = factor(casecontstat,
      levels = c(0, 1),
      labels = c("control", "case"))
  ) %>%
  (function(x) data.frame(x, row.names = x$geo_accession))

## Append annotation data
annotation(gse81954) <- "HuRSTA2a520709"
fData(gse81954) <- fData(gse81954) %>%
  mutate(
    probeid = featureNames(gse81954),
    entrezid = unlist(mget(probeid, HuRSTA2a520709ENTREZID, ifnotfound = NA)),
    symbol = unlist(mget(probeid, HuRSTA2a520709SYMBOL, ifnotfound = NA))
  ) %>%
  data.frame(row.names = featureNames(gse81954))

stopifnot(validObject(gse81954))
