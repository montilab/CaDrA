library(devtools)
load_all(".")

library(shiny)

port <- Sys.getenv("PORT")
app <- CaDrA::CaDrA_App()

shiny::runApp(app, host='0.0.0.0', port=ifelse(is.na(as.numeric(port)), 3838, port))
