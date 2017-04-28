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

tmp_tbl  %>%
  arrange(desc(Sum)) %>%
  head()

tmp_tbl %>%
  select(Sum) %>%
  table() %>%
  addmargins()

## Identify clusters of casecontsets that share (at least one) individual
## tumour(s), but not with other clusters

casecontsets <- split(casecontstudy_design$tumorid, casecontstudy_design$setnr)

## http://stackoverflow.com/questions/25130462/get-disjoint-sets-from-a-list-in-r
library(igraph)
igraph_method <- function(lst) {
  edg <- do.call("rbind", lapply(lst, function(x) {
    if (length(x) > 1) t(combn(x, 2)) else NULL
  }))
  g <- graph.data.frame(edg)
  split(V(g)$name, clusters(g)$membership)
}
disjoint_sets <- igraph_method(casecontsets)

disjoint_sets_of_casecontsets <- list()  # initialise
for (i in seq(along = disjoint_sets)) {
    y <- lapply(casecontsets, function(x) any(is.element(x, disjoint_sets[[i]])))
    set_of_casecontsets <- as.integer(names(y)[unlist(y)])
    disjoint_sets_of_casecontsets[[i]] <- set_of_casecontsets
}
names(disjoint_sets_of_casecontsets) <-
  lapply(disjoint_sets_of_casecontsets, paste, collapse = "_")

## Revised data partition
all_sets_of_casecontsets <- names(disjoint_sets_of_casecontsets)

y <- lapply(disjoint_sets_of_casecontsets,
    function(x) !any(is.element(x,
      subset(pData(qcsubstudy), complete_casecontset)$setnr)))
sets_of_casecontsets_not_in_qcsubstudy <- names(y)[unlist(y)]

n_all <- unlist(lapply(disjoint_sets_of_casecontsets, length))
n <- n_all[sets_of_casecontsets_not_in_qcsubstudy]

set.seed(20101216)
n_randomly_reordered <- n[sample.int(n = length(n))]
x <- cumsum(n_randomly_reordered)

all_casecontsets <- unique(casecontstudy$setnr)
x[x <= length(all_casecontsets)*2/3]
