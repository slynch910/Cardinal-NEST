# Cardinal NEST - shared color palette (bslib theme, CSS vars, plotting)

cardinal_palette <- list(
  primary    = "#B22222",  # Cardinal Red
  secondary  = "#4A5568",  # Slate Gray
  accent     = "#D4A017",  # Gold
  background = "#F8F9FA",  # Off White
  surface    = "#FFFFFF",  # White
  text       = "#2D3748",  # Charcoal
  success    = "#38A169",  # Sage Green
  warning    = "#DD6B20"   # Amber
)

# bslib Bootstrap theme built from the palette.
cardinal_bs_theme <- function(palette = cardinal_palette) {
  bslib::bs_theme(
    version    = 5,
    bg         = palette$surface,
    fg         = palette$text,
    primary    = palette$primary,
    secondary  = palette$secondary,
    success    = palette$success,
    warning    = palette$warning,
    "body-bg"  = palette$background
  )
}

# Palette as CSS custom properties (var(--cn-primary), ...) for use in
# inline styles and custom classes.
cardinal_palette_css <- function(palette = cardinal_palette) {
  vars <- paste0(
    sprintf("  --cn-%s: %s;", names(palette), unlist(palette)),
    collapse = "\n"
  )
  shiny::tags$head(
    shiny::tags$style(shiny::HTML(sprintf(":root {\n%s\n}", vars)))
  )
}
