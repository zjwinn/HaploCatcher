## Resubmission

This is a resubmission. In response to the previous submission, I have fixed
the (possibly) invalid file URI reported by CRAN: the README linked to NEWS.md
with a relative file path. The README now links to an absolute, resolvable URL
(the project's GitHub release notes) instead.

## Submission

This is an update of HaploCatcher (version 2.0.1).

It is a major internal restructure for speed and maintainability. User-facing
highlights:

* Shared logic was consolidated into hidden internal helpers; marker selection
  is now O(m) instead of O(m^2) (identical results on the example data).
* New optional `het_label` / `neg_label` arguments let users define the
  positive / heterozygous / negative call labels instead of relying on the
  gene / het_gene / non_gene naming convention (default behaviour unchanged).
* The forward-prediction step now falls back to the alternate model (with a
  warning) when the selected model fails (e.g. k-nearest neighbors "too many
  ties"), and only errors if both models fail.

## Test environments

* local: Windows 11, R 4.6.0 -- 0 errors | 0 warnings | 0 notes
* win-builder: R-devel and R-release -- 0 errors | 0 warnings | 0 notes

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN (confirmed with revdepcheck).
