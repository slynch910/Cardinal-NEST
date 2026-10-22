# Cardinal NEST - application entry point

# Runs the app from a source checkout without installing the package:
#   shiny::runApp(".")
# Installed users should call CardinalNEST::runCardinalNEST() instead.

library(shiny)
library(Cardinal)
library(ggplot2)
library(DT)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

cardinal_nest_app(options = list(port = 8080, host = "127.0.0.1"))
