doc_stem <- "gsea_rotation_vs_permutation"

output_dir <- Sys.getenv("QUARTO_OUTPUT_DIR", unset = "")
if (!nzchar(output_dir)) {
  output_dir <- "report"
}

cleanup_if_duplicated <- function(root_path, out_path) {
  if (!file.exists(root_path)) {
    return(invisible(FALSE))
  }
  if (!file.exists(out_path)) {
    return(invisible(FALSE))
  }

  root_norm <- normalizePath(root_path, winslash = "/", mustWork = FALSE)
  out_norm <- normalizePath(out_path, winslash = "/", mustWork = FALSE)
  if (identical(root_norm, out_norm)) {
    return(invisible(FALSE))
  }

  unlink(root_path, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

cleanup_if_duplicated(
  root_path = paste0(doc_stem, "_files"),
  out_path = file.path(output_dir, paste0(doc_stem, "_files"))
)

cleanup_if_duplicated(
  root_path = paste0(doc_stem, "_cache"),
  out_path = file.path(output_dir, paste0(doc_stem, "_cache"))
)
