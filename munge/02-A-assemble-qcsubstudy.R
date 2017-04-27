tmpdf <- tibble(geo_accession = colnames(qcsubstudy_exprs)) %>%
  left_join(
    rbind(
      pData(gse48091)[, c("geo_accession", "tumorid", "piece")],
      pData(gse81954)[, c("geo_accession", "tumorid", "piece")]),
    by = "geo_accession") %>%
  left_join(
    pData(gse81954)[, c("tumorid", "setnr", "casecontstat", "casecontcd")],
    by = "tumorid") %>%
  data.frame()
rownames(tmpdf) <- tmpdf$geo_accession

qcsubstudy <- ExpressionSet(qcsubstudy_exprs)
pData(qcsubstudy) <- tmpdf

## Munge phenotypic data
pData(qcsubstudy) <- pData(qcsubstudy) %>%
  as_tibble() %>%
  mutate(
    rna_extract = factor(piece, levels = c(1, 2),
      labels = c("original", "reextract"))
  ) %>%
  (function(x) data.frame(x, row.names = x$geo_accession))

## Append annotation data
annotation(qcsubstudy) <- "HuRSTA2a520709"
fData(qcsubstudy) <- fData(qcsubstudy) %>%
  mutate(
    probeid = featureNames(qcsubstudy),
    entrezid = unlist(mget(probeid, HuRSTA2a520709ENTREZID, ifnotfound = NA)),
    symbol = unlist(mget(probeid, HuRSTA2a520709SYMBOL, ifnotfound = NA))
  ) %>%
  data.frame(row.names = featureNames(qcsubstudy))

stopifnot(validObject(qcsubstudy))
