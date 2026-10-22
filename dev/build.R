# Helper: generate man/ from roxygen comments, build vignette, run checks.
# Run from the package root. Requires: devtools, roxygen2, knitr, rmarkdown.
#   install.packages(c("devtools", "roxygen2", "knitr", "rmarkdown"))
#   (and optionally "Synth", "testthat")

devtools::document()   # (re)generates man/*.Rd and NAMESPACE from roxygen
devtools::test()       # runs tests/testthat/
devtools::build_vignettes()
devtools::check()      # full R CMD check
# devtools::install()  # install into your library
