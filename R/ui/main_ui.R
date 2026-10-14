# Cardinal NEST - main UI

library(shiny)

main_ui <- function() {
  navbarPage(
    title = "Cardinal NEST",
    id = "main_tabs",
    theme = cardinal_bs_theme(),
    header = tagList(
      cardinal_palette_css(),
      # DataTables rendered while their tab is hidden get zero column widths;
      # re-adjust visible tables each time a tab is shown.
      tags$script(HTML(
        "$(document).on('shown.bs.tab', function() {",
        "  if (window.jQuery && $.fn && $.fn.dataTable) {",
        "    $.fn.dataTable.tables({ visible: true, api: true }).columns.adjust();",
        "  }",
        "});"
      ))
    ),
    tabPanel("Import", icon = icon("upload"),
             dataImportUI("import")),
    tabPanel("Explore", icon = icon("search"),
             metadataUI("explore")),
    tabPanel("Preprocess", icon = icon("cogs"),
             preprocessingUI("preprocess")),
    tabPanel("Spatial", icon = icon("draw-polygon"),
             spatialUI("spatial")),
    tabPanel("Visualize", icon = icon("image"),
             visualizationUI("visualize_content")),
    tabPanel("Analyze", icon = icon("chart-bar"),
             analysisUI("analyze")),
    tabPanel("Export", icon = icon("download"),
             exportUI("export")),
    tabPanel("Workflow", icon = icon("diagram-project"),
             workflowUI("workflow")),
    tabPanel("Performance", icon = icon("gauge-high"),
             performanceUI("performance")),
    tabPanel("Help", icon = icon("question-circle"),
             helpUI("help"))
  )
}
