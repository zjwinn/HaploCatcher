# HaploCatcher 2.0.0

## Major restructure

* All shared logic was extracted into hidden internal helpers (the `.hc_*`
  family) so the exported functions are short and consistent. The data
  preparation, het filtering, marker selection, model fitting, evaluation,
  permutation summarization, and vote tabulation are now each defined once
  instead of being copy across `locus_cv()`, `locus_train()`,
  `locus_perm_cv()`, and `auto_locus()`.
* Marker selection is dramatically faster. The previous code built a full
  marker-by-marker correlation matrix (O(m^2)) and discarded all but one
  column; it now correlates the call against every marker in a single pass
  (O(m)). On the bundled example data this selects an identical set of top
  markers with no change in the correlation values.
* Cross-validation control was simplified. `trainControl()` previously
  requested `method = "repeatedcv"` with `repeats = 1000` while also supplying
  an explicit `index`, which caret silently ignores. It now uses a plain
  5-fold CV over the same folds — same behavior, far less overhead.
* A bug in the random-forest leave-one-out fallback (which retrained a KNN
  model instead of RF) was fixed.

## New features

* Heterozygous and negative cases can now be defined explicitly via the new
  `het_label` (and, for plotting, `neg_label`) arguments on `locus_cv()`,
  `locus_train()`, `locus_perm_cv()`, `plot_locus_perm_cv()`, and
  `auto_locus()`. When these are left `NULL` the package keeps its original
  `gene` / `het_gene` / `non_gene` naming convention, so existing workflows are
  unaffected.
* Input checking is centralized and more robust, with clearer error messages
  for malformed genotype matrices, gene files, marker info, proportions, and
  case labels.

## Other

* A single, directly source-able copy of the package is provided as
  `HaploCatcher-source.R` (regenerate with `tools/build_single_source.R`).
* `lattice` was dropped from Imports; `stats` was added.
