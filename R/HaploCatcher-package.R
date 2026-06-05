#' HaploCatcher: A Predictive Haplotyping Package
#'
#' HaploCatcher predicts the allelic state (haplotype) of a genotype at a
#' specific locus / QTL / gene by training k-nearest neighbors (KNN) and random
#' forest (RF) models on genome-wide markers. A training panel that has already
#' been screened for the locus teaches the models, which then predict the
#' haplotype of un-screened, genome-wide genotyped lines. The method follows
#' Winn et al. (2022) <doi:10.1007/s00122-022-04178-w> and the package itself is
#' described in Winn et al. (2023) <doi:10.1002/tpg2.20412>.
#'
#' The public pipeline mirrors Figure 1b of the package paper:
#' permutation cross-validation (\code{\link{locus_perm_cv}} over
#' \code{\link{locus_cv}}), best-model selection by kappa or accuracy, then
#' forward prediction by a single seeded model (\code{\link{locus_train}} +
#' \code{\link{locus_pred}}) or by majority-rule voting. The wrapper
#' \code{\link{auto_locus}} runs the whole pipeline; \code{\link{plot_locus_perm_cv}}
#' visualizes cross-validation. All shared work lives in hidden \code{.hc_*}
#' helpers so the exported functions stay short and consistent.
#'
#' @keywords internal
#' @importFrom foreach %dopar% foreach
#' @importFrom randomForest randomForest
#' @importFrom ggplot2 ggplot
#' @importFrom patchwork plot_layout
#' @importFrom caret train knn3 createFolds trainControl confusionMatrix
"_PACKAGE"
