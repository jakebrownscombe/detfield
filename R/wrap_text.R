# ---- Messages and printed lines wrapped to the console width ----

# Wraps `text` to 8 characters under getOption("width"), so a message still
# fits after the "Error : " or "Warning: " that R prints before it.
# Continuation lines are indented by `exdent` spaces. A trailing newline is
# dropped before wrapping.
wrap_text <- function(text, exdent = 2) {
  text <- sub("\n$", "", text)
  width <- getOption("width", 80) - 8
  paste(strwrap(text, width = width, exdent = exdent), collapse = "\n")
}

# cat() of one wrapped line, ending in a newline. Leading spaces of `text`
# are kept as the indent of the first line, and continuation lines are
# indented two spaces further. Lines stop 4 characters short of the width,
# so a line still fits after the "#> " that knitr prints before output.
cat_wrapped <- function(text) {
  lead <- nchar(text) - nchar(sub("^ +", "", text))
  body <- sub("^ +", "", sub("\n$", "", text))
  out <- strwrap(body, width = getOption("width", 80) - 4 - lead, exdent = 2)
  cat(paste0(strrep(" ", lead), out), sep = "\n")
}
