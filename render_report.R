required_packages <- c("rmarkdown", "knitr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}

rmarkdown::render(
  input = "report/analysis_report.Rmd",
  output_file = "analysis_report.html",
  output_dir = "report",
  knit_root_dir = getwd(),
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)
