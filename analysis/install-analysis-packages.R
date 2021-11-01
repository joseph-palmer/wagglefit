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
new.packages <- new.packages[new.packages != "wagglefit"]
if( length(new.packages)) {
  message("Installing the following packages:")
  print(new.packages)
  install.packages(new.packages)
} else {
  message("Nothing to install, all dependencies satisfied")
}
