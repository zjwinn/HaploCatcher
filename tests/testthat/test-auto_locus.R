test_that("auto_locus function works", {

  data("geno_mat")
  data("gene_comp")
  data("marker_info")
  set.seed(022294)
  training_genotypes=sample(rownames(geno_mat), size = round(nrow(geno_mat)*0.8, 0))
  testing_genotypes=rownames(geno_mat)[!rownames(geno_mat) %in% training_genotypes]
  set.seed(NULL)

  expect_no_error(auto_locus(geno_mat = geno_mat,
                             gene_file = gene_comp,
                             gene_name = "sst1_solid_stem",
                             marker_info = marker_info,
                             chromosome = "3B",
                             training_genotypes = training_genotypes,
                             testing_genotypes = testing_genotypes,
                             set_seed = 022294,
                             n_perms = 5,
                             include_hets = TRUE))


})
