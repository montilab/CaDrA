
library(devtools)
load_all("/srv/shiny-server/")

library(plumber)

shiny::runApp("cadra_api.R", host='0.0.0.0', port=3838)
