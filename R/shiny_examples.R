
#' Seasonality Shiny
#'
#' \code{run_seasonality_shiny} runs a shiny app that displays the outcomes 
#'   of the seasonality investigation. The shiny shows how the bias in the 
#'   observed prevalence of pfhrp2 deletions changes for each level 1 admin 
#'   region changes throughout the transmission season
#'
#' @export
run_seasonality_shiny <- function() {
  appDir <- system.file("shiny-seasonality", "app", package = "hrp2malaRia")
  if (appDir == "") {
    stop("Could not find shiny directory. Try re-installing `hrp2malaRia`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}