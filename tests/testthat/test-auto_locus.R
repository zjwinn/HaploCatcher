# Tests exercise the full pipeline against the bundled example data sets
# (geno_mat, gene_comp, marker_info). Model sizes are kept small so the suite
# runs quickly while still touching every exported function and both the
# default and user-defined case-label code paths.

# shared fixtures: load the example data once and define a train/test split
data("geno_mat")
data("gene_comp")
data("marker_info")
gene <- "sst1_solid_stem"
set.seed(022294)
train_ids <- sample(rownames(geno_mat), size = round(nrow(geno_mat) * 0.8, 0))
test_ids  <- rownames(geno_mat)[!rownames(geno_mat) %in% train_ids]
set.seed(NULL)

test_that("locus_cv returns the expected structure", {
  set.seed(1)
  fit <- locus_cv(geno_mat = geno_mat, gene_file = gene_comp, gene_name = gene,
                  marker_info = marker_info, chromosome = "3B",
                  ncor_markers = 10, n_neighbors = 5,
                  include_hets = FALSE, verbose = FALSE)
  expect_named(fit, c("data_frames", "test_predictions", "confusion_matrices"))
  expect_s3_class(fit$confusion_matrices$rf, "confusionMatrix")
})

test_that("locus_train + locus_pred forward-predict held-out lines", {
  set.seed(2)
  train_gc <- gene_comp[gene_comp$FullSampleName %in% train_ids, ]
  fit <- locus_train(geno_mat = geno_mat, gene_file = train_gc, gene_name = gene,
                     marker_info = marker_info, chromosome = "3B",
                     ncor_markers = 10, n_neighbors = 5,
                     models_request = "rf", set_seed = 2, verbose = FALSE)
  pred <- locus_pred(fit, geno_mat, test_ids)
  expect_true(all(c("FullSampleName", "Prediction_RF") %in% names(pred)))
  expect_gt(nrow(pred), 0)
})

test_that("locus_perm_cv returns the five-element summary object", {
  set.seed(3)
  res <- locus_perm_cv(n_perms = 2, geno_mat = geno_mat, gene_file = gene_comp,
                       gene_name = gene, marker_info = marker_info, chromosome = "3B",
                       ncor_markers = 10, n_neighbors = 5, verbose = FALSE)
  expect_length(res, 5)
  expect_named(res, c("Overall_Parameters", "By_Class_Parameters",
                      "Overall_Summary", "By_Class_Summary", "Raw_Permutation_Info"))
})

test_that("auto_locus single-model pipeline runs end to end", {
  expect_no_error(
    auto_locus(geno_mat = geno_mat, gene_file = gene_comp, gene_name = gene,
               marker_info = marker_info, chromosome = "3B",
               training_genotypes = train_ids, testing_genotypes = test_ids,
               ncor_markers = 10, n_neighbors = 5, n_perms = 2,
               set_seed = 022294, include_hets = TRUE,
               verbose = FALSE, plot_cv_results = FALSE)
  )
})

test_that("auto_locus voting pipeline runs end to end", {
  res <- auto_locus(geno_mat = geno_mat, gene_file = gene_comp, gene_name = gene,
                    marker_info = marker_info, chromosome = "3B",
                    training_genotypes = train_ids, testing_genotypes = test_ids,
                    ncor_markers = 10, n_neighbors = 5, n_perms = 2,
                    predict_by_vote = TRUE, n_votes = 2,
                    verbose = FALSE, plot_cv_results = FALSE)
  expect_true("consensus_predictions" %in% names(res))
  expect_true("Consensus_Call" %in% names(res$consensus_predictions))
})

test_that("user-defined het_label matches the default convention", {
  # the bundled gene uses the het_<gene> convention, so passing the explicit
  # label should select exactly the same individuals as the default
  set.seed(4); a <- locus_cv(geno_mat, gene_comp, gene, marker_info, "3B",
                             ncor_markers = 10, n_neighbors = 5, verbose = FALSE)
  set.seed(4); b <- locus_cv(geno_mat, gene_comp, gene, marker_info, "3B",
                             ncor_markers = 10, n_neighbors = 5, verbose = FALSE,
                             het_label = "het_sst1_solid_stem")
  expect_identical(nrow(a$data_frames$training), nrow(b$data_frames$training))
})

test_that("input checking rejects malformed arguments", {
  # missing required gene_file column
  bad <- gene_comp[, setdiff(names(gene_comp), "Call")]
  expect_error(locus_cv(geno_mat, bad, gene, marker_info, "3B"),
               "Gene", ignore.case = TRUE)
  # out-of-range training proportion
  expect_error(locus_cv(geno_mat, gene_comp, gene, marker_info, "3B",
                        percent_training = 1.5))
  # unknown chromosome
  expect_error(locus_cv(geno_mat, gene_comp, gene, marker_info, "ZZ"))
})
