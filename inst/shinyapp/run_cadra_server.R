
library(devtools)
load_all("/srv/shiny-server/")

library(shiny)

app <- CaDrA::CaDrA_App()

shiny::runApp(app, host='0.0.0.0', port=3838)
