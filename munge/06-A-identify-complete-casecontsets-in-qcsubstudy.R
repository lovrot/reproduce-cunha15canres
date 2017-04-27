pData(qcsubstudy) <- pData(qcsubstudy) %>%
  as_tibble() %>%
  ## Identify sets with one case and at least one control
  group_by(rna_extract, setnr) %>%
  mutate(
    n_cases_in_set = sum(casecontcd == "case"),
    n_controls_in_set = sum(casecontcd == "control"),
    complete_casecontset =
      n_cases_in_set == 1 & n_controls_in_set >= 1) %>%
  (function(x) data.frame(x, row.names = x$geo_accession))

stopifnot(validObject(qcsubstudy))
