sel_symbols <- c("ERBB2", "GRB7")

eset <- casecontstudy[fData(casecontstudy)$symbol %in% sel_symbols, ] %>%
  genefilter::featureFilter()
featureNames(eset) <- as.character(fData(eset)$symbol)
pData(eset) <- cbind(pData(eset), t(exprs(eset)))
casecontstudy$HER2_amplicon_metagene <- colMeans(exprs(eset))

eset <- qcsubstudy[fData(qcsubstudy)$symbol %in% sel_symbols, ] %>%
  averageExprsByFDataVar("symbol")
featureNames(eset) <- as.character(fData(eset)$symbol)
qcsubstudy$HER2_amplicon_metagene <- colMeans(exprs(eset))
