setnr <- strsplit(gse48091$setnr, split = ",")
names(setnr) <- gse48091$tumorid
setnr <- lapply(setnr, as.integer)

casecontstat <- strsplit(gse48091$casecontstat, split = ",")
names(casecontstat) <- gse48091$tumorid
casecontstat <- lapply(casecontstat, as.integer)

casecontstudy_design <-
  tibble(
    tumorid = substr(names(unlist(setnr)), 1, 8),
    setnr = unlist(setnr),
    casecontstat = unlist(casecontstat)) %>%
  mutate(
    casecontcd = factor(casecontstat,
      levels = c(0, 1),
      labels = c("control", "case")),
    subjid = paste(tumorid, setnr, casecontstat, sep = "_"),
    trtgrp = factor(round(setnr / 1000) * 1000,
      levels = c(1000, 3000, 2000), labels = c("ET", "CT+ET", "CT"))
  ) %>%
  left_join(pData(gse48091)[, c("tumorid", "geo_accession")], by = "tumorid") %>%
  data.frame()

rownames(casecontstudy_design) <- casecontstudy_design$subjid
