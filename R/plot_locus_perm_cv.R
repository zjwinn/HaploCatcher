#' Visualize Permutation CV Results
#'
#' Takes the result of [locus_perm_cv()] and draws a composite of accuracy,
#' kappa, sensitivity, and specificity across permutations. When more than one
#' call class is present (heterozygotes retained), sensitivity and specificity
#' are faceted by class.
#'
#' @param results A list produced by [locus_perm_cv()].
#' @param individual_images Logical; also print each panel on its own. Default FALSE.
#' @param het_label Optional character vector of class labels to treat as heterozygous when relabeling facets. When NULL (default), the "het_" prefix is used.
#' @param neg_label Optional character vector of class labels to treat as the negative/wild-type case when relabeling facets. When NULL (default), the "non_" prefix is used.
#'
#' @return Invisibly returns NULL; called for its plotting side effect.
#'
#' @export
#'
#' @examples
#'
#' #refer to vignette for an in depth look at the plot_locus_perm_cv function
#' vignette("An_Intro_to_HaploCatcher", package = "HaploCatcher")
#'
#' @importFrom patchwork plot_layout
#' @importFrom ggplot2 ggplot

plot_locus_perm_cv <- function(results, individual_images = FALSE,
                               het_label = NULL, neg_label = NULL) {

  # the result must be the five-element object from locus_perm_cv()
  if (!is.list(results) || length(results) != 5L ||
      is.null(results$Overall_Parameters) || is.null(results$By_Class_Parameters))
    stop("'results' is not a list produced by 'locus_perm_cv()'!", call. = FALSE)

  # bind the column names referenced in aes() to satisfy R CMD check
  Accuracy <- Kappa <- Sensitivity <- Specificity <- Model <- Class <- NULL

  # a small constructor for the recurring boxplot panels
  panel <- function(df, mapping, title, legend = FALSE, facet = FALSE) {
    # base boxplot of the requested metric by model
    p <- ggplot2::ggplot(df, mapping) +
      ggplot2::geom_boxplot() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::labs(title = title)
    # optionally facet sensitivity/specificity by class
    if (facet) p <- p + ggplot2::facet_grid(~Class)
    # drop the legend on secondary panels for a cleaner composite
    if (!legend) p <- p + ggplot2::theme(legend.position = "none")
    p
  }

  # overall accuracy and kappa come from the overall-parameter table
  ov <- results$Overall_Parameters
  b <- panel(ov, ggplot2::aes(y = Accuracy, x = Model, fill = Model), "Overall Accuracy", legend = TRUE)
  c <- panel(ov, ggplot2::aes(y = Kappa, x = Model, fill = Model), "Overall Kappa")

  # sensitivity and specificity come from the by-class table
  bc <- results$By_Class_Parameters
  # decide whether the locus is multiclass (heterozygotes retained)
  multiclass <- length(unique(stats::na.omit(bc$Class))) > 1
  if (multiclass) {
    # relabel raw class names to standard haplotype notation for the facets,
    # honoring any user-supplied negative/heterozygous label sets
    bc$Class <- ifelse(.hc_is_neg(bc$Class, neg_label), "-/-",
                       ifelse(.hc_is_het(bc$Class, het_label), "+/-", "+/+"))
    bc$Class <- factor(bc$Class, levels = c("+/+", "+/-", "-/-"))
    d <- panel(bc, ggplot2::aes(y = Sensitivity, x = Model, fill = Model), "By-Class Sensitivity", facet = TRUE)
    e <- panel(bc, ggplot2::aes(y = Specificity, x = Model, fill = Model), "By-Class Specificity", facet = TRUE)
  } else {
    d <- panel(bc, ggplot2::aes(y = Sensitivity, x = Model, fill = Model), "Overall Sensitivity")
    e <- panel(bc, ggplot2::aes(y = Specificity, x = Model, fill = Model), "Overall Specificity")
  }

  # optionally print each panel on its own
  if (individual_images) { print(b); print(c); print(d); print(e) }

  # print the 2x2 composite with a shared legend
  print(((b + c) / (d + e)) + patchwork::plot_layout(guides = "collect"))
  # nothing meaningful to return
  invisible(NULL)
}
