
library(devtools)
load_all(".")

library(plumber)

cadra_api <- plumber::plumb("./inst/plumber/cadra_api.R")
cadra_api$run(host='0.0.0.0', port=3838)
