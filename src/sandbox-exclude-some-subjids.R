## Objective: Construct a data split into training and test sets such that
## * the QC study (complete) case-control sets are in the test partition
## * study subjects that correspond to indiviudal tumours included in several
##   case-control sets are in either one of the two paritions

ProjectTemplate::reload.project()

select <- dplyr::select  # mask AnnotationDbi::select

n_keep <- 5000

## Non-specific filtering of features

eset <- genefilter::featureFilter(casecontstudy)
eset <- genefilter::varFilter(eset, var.cutoff = 1 - n_keep / nrow(eset))
featureNames(eset) <- make.names(fData(eset)$probeid)

all_data <- cbind(
  pData(eset)[, "casecontcd", drop = FALSE],
  as.data.frame(t(exprs(eset))))

## (Current) Data splitting

all_casecontsets <- unique(casecontstudy$setnr)
casecontsets_not_in_qcsubstudy <- setdiff(all_casecontsets,
  subset(pData(qcsubstudy), complete_casecontset)$setnr)

set.seed(20101216)
train_sets <- sample(casecontsets_not_in_qcsubstudy,
  size = length(all_casecontsets)*2/3)
train_subjects <- subset(pData(casecontstudy), setnr %in% train_sets)$subjid
train_idx <- which(rownames(all_data) %in% train_subjects)

training  <- all_data[train_idx, ]
testing  <- all_data[-train_idx, ]

## Visualise parition (casecontset-centric view)
casecontsets <- unique(casecontstudy$setnr)
tibble(
  casecontset = casecontsets,
  partition = factor(
    casecontset %in% pData(casecontstudy)[rownames(training), "setnr"],
    levels = c(TRUE, FALSE),
    labels = c("Training set", "Test set")),
  in_qcsubstudy = factor(
    casecontset %in% subset(pData(qcsubstudy), complete_casecontset)$setnr,
    levels = c(TRUE, FALSE),
    labels = c("QC substudy", "(rest)"))) %>%
  dplyr::select(in_qcsubstudy, partition) %>%
  table() %>%
  addmargins() %>%
  knitr::kable(caption = "Table. Number of case-control sets in each partition of the full study.")

## Visualise parition (individual tumour-centric view)
tumorids <- gse48091$tumorid
tibble(
  tumorid = tumorids,
  partition = factor(
    tumorid %in% pData(casecontstudy)[rownames(training), "tumorid"],
    levels = c(TRUE, FALSE),
    labels = c("Training set", "Test set")),
  in_qcsubstudy = factor(
    tumorid %in% subset(pData(qcsubstudy), complete_casecontset)$tumorid,
    levels = c(TRUE, FALSE),
    labels = c("QC substudy", "(rest)"))) %>%
  dplyr::select(in_qcsubstudy, partition) %>%
  table() %>%
  addmargins() %>%
  knitr::kable(caption = "Table. Number of individual tumours in each partition of the full study.")

tmp_tbl <- pData(casecontstudy) %>%
  mutate(
    partition = factor(
      subjid %in% rownames(training),
      levels = c(TRUE, FALSE),
      labels = c("Training set", "Test set"))
  ) %>%
  select(tumorid, partition) %>%
  table() %>%
  addmargins(margin = 2) %>%
  as.data.frame.matrix() %>%
  rownames_to_column(var = "tumorid") %>%
  as_tibble()
names(tmp_tbl) <- gsub(" ", "_", names(tmp_tbl), fixed = TRUE)

tmp_tbl  %>%
  arrange(desc(Sum)) %>%
  head()

tmp_tbl %>%
  select(Sum) %>%
  table() %>%
  addmargins()

tmp_tbl2 <- tmp_tbl  %>%
  filter(Sum > 1) %>%
  arrange(desc(Sum))

tmp_tbl2$in_qcsubstudy <-
  tmp_tbl2$tumorid %in%
  subset(pData(qcsubstudy), complete_casecontset)$tumorid
tmp_tbl2$is_case_in_training_set <-
  tmp_tbl2$tumorid %in%
  subset(pData(casecontstudy), casecontstat == 1 & subjid %in% rownames(training))$tumorid
tmp_tbl2$is_case_in_test_set <-
  tmp_tbl2$tumorid %in%
  subset(pData(casecontstudy), casecontstat == 1 & subjid %in% rownames(testing))$tumorid

tmp_tbl3 <- tmp_tbl2 %>%
  mutate(
    keep_in_test_set = !is_case_in_test_set &
      (in_qcsubstudy | is_case_in_test_set | (Test_set >= Training_set))
  )

if (interactive()) try(View(tmp_tbl3))

##

training <- training[
  !(pData(casecontstudy)[rownames(training), "tumorid"] %in%
      filter(tmp_tbl3, keep_in_test_set)$tumorid), ]
testing <- testing[
  !(pData(casecontstudy)[rownames(testing), "tumorid"] %in%
      filter(tmp_tbl3, !keep_in_test_set)$tumorid), ]

## Visualise parition (individual tumour-centric view)
tumorids <- gse48091$tumorid
tibble(
  tumorid = tumorids,
  partition =
    ifelse(tumorid %in% pData(casecontstudy)[rownames(training), "tumorid"], "Training set",
      ifelse(tumorid %in% pData(casecontstudy)[rownames(testing), "tumorid"], "Test set", NA)),
  in_qcsubstudy = factor(
    tumorid %in% subset(pData(qcsubstudy), complete_casecontset)$tumorid,
    levels = c(TRUE, FALSE),
    labels = c("QC substudy", "(rest)"))) %>%
  dplyr::select(in_qcsubstudy, partition) %>%
  table(useNA = "ifany") %>%
  addmargins() %>%
  knitr::kable(caption = "Table. Number of individual tumours in each partition of the full study.")

tmp_tbl <- pData(casecontstudy) %>%
  mutate(
    partition =
      ifelse(tumorid %in% pData(casecontstudy)[rownames(training), "tumorid"], "Training set",
        ifelse(tumorid %in% pData(casecontstudy)[rownames(testing), "tumorid"], "Test set", NA))
  ) %>%
  select(tumorid, partition) %>%
  table() %>%
  addmargins(margin = 2) %>%
  as.data.frame.matrix() %>%
  rownames_to_column(var = "tumorid") %>%
  as_tibble()
names(tmp_tbl) <- gsub(" ", "_", names(tmp_tbl), fixed = TRUE)

tmp_tbl  %>%
  arrange(desc(Sum)) %>%
  head()
