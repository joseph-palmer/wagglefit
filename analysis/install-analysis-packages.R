library_calls <- purrr::map(
  list.files("analysis/", pattern = "*.Rmd|.R"),
  ~ {
    stringr::str_match(readLines(paste0("analysis/", .x)), "library\\((.*?)\\)")[,2]
  }
)

library_calls <- unlist(library_calls)
packages <- library_calls[!is.na(library_calls)]
new.packages <- packages[
  !(packages %in% installed.packages()[,"Package"])
]
if( length(new.packages)) {
  install.packages(new.packages)
}
