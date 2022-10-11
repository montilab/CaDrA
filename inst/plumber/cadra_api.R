
##Devtools packages####
library(devtools)
load_all(".")

##R packages####
library(tidyverse)
library(Biobase)
library(jsonlite)

# Increase the memory limit for reading and downloading data from API
library(unix)
unix::rlimit_as(1e12, 1e12)

# Global expression sets
dataset_choices <- list(
  "CCLE_MUT_SCNA" =  system.file("data/CCLE_MUT_SCNA.rda", package = "CaDrA"),
  "sim.ES" = system.file("data/sim.ES.RData", package = "CaDrA"),
  "BRCA_GISTIC_MUT_SIG" = system.file("data/BRCA_GISTIC_MUT_SIG.rda", package = "CaDrA")
)

# Global input scores
score_choices <- list(
  "CTNBB1_reporter" = system.file("data/CTNBB1_reporter.rda", package = "CaDrA"),
  "sim.Scores" =  system.file("data/sim.Scores.rda", package = "CaDrA"),
  "TAZYAP_BRCA_ACTIVITY" = system.file("data/TAZYAP_BRCA_ACTIVITY.rda", package = "CaDrA")
)

# Obtain the external data
get_extdata <- function(dataset_choices, score_choices){
  
  # Check if external data exists in package
  if(file.exists(system.file("extdata", "datalist.csv", package = "CaDrA"))){
    
    # Read in a list of files in datalist.csv 
    datalist <- read.csv(system.file("extdata/datalist.csv", package = "CaDrA"), header=TRUE)
    datalist <- datalist[which(datalist$eset_paths != "" & !is.na(datalist$eset_paths) & datalist$eset_names != "" & !is.na(datalist$eset_names)),]
    
    # Obtain external expression sets
    eset_paths <- datalist$eset_paths
    eset_names <- datalist$eset_names
    
    # Create a labels for each file
    names(eset_paths) <- eset_names
    
    if(length(eset_paths) > 0){
      dataset_choices <- c(dataset_choices, eset_paths)
    }
    
    # Obtain external input scores
    score_paths <- datalist$score_paths
    score_names <- datalist$score_names
    
    # Create a labels for each file
    names(score_paths) <- score_names
    
    if(length(score_paths) > 0){
      score_choices <- c(score_choices, score_paths)
    }
    
  }
  
  return(list(eset_choices=dataset_choices, input_score_choices=score_choices))
  
}

## Get a list of feature set available on database
#* @param orderby
#' @get /feature_set
#' @post /feature_set
feature_set <- function(orderby="asc"){
  
  # Combine extdata with global expression set and scores dataset if it was provided
  eset_choices <- get_extdata(dataset_choices, score_choices)[["eset_choices"]] %>% unlist() %>% names() 
  
  # Sort projects by ascending or descending order
  if(orderby == "asc"){
    eset_choices <- sort(eset_choices, decreasing = FALSE)
  }else if(orderby == "desc"){
    eset_choices <- sort(eset_choices, decreasing = TRUE)
  }
  
  return(toJSON(eset_choices, pretty=TRUE))
  
}

## Get expression set of the selected feature set
#* @param feature_set
#' @serializer rds
#' @get /expression_set
#' @post /expression_set
expression_set <- function(res, req, feature_set, include_scores=FALSE){
  
  # Remove multiple spaces
  feature_set <- gsub("\\s+", " ", feature_set)
  
  # Combine extdata with global expression set and scores dataset if it was provided
  extdata <- get_extdata(dataset_choices, score_choices)
  eset_choices <- extdata[["eset_choices"]] %>% unlist()
  input_score_choices <- extdata[["input_score_choices"]] %>% unlist()
  eset_names <- names(eset_choices)
  
  # Check if the project is available on GeneHive, if not, return warning message
  if(toupper(feature_set) %in% toupper(eset_names)){
    
    dl_data <- list()
    
    dataset <- eset_choices[which(toupper(eset_names) == toupper(feature_set))]
    
    if(tools::file_ext(dataset) == "rda" | tools::file_ext(dataset) == "RData"){
      envir_name <- load(dataset)
      ES <- get(envir_name)
    }else{
      ES <- readRDS(system.file("extdata", "eset", dataset, package = "CaDrA"))
    }
    
    dl_data <- c(dl_data, ES=ES)
    
    if(include_scores){
      
      scores <- input_score_choices[which(eset_choices == dataset)] %>% unlist()
      
      if(!is.na(names(scores))){
        
        if(tools::file_ext(scores) == "rda" | tools::file_ext(scores) == "RData"){
          envir_name <- load(scores)
          input_score <- get(envir_name)
        }else{
          input_score <- readRDS(system.file("extdata", "input_score", scores, package = "CaDrA"))
        }
        
        dl_data <- c(dl_data, input_score=list(input_score))

      }else{
        
        dl_data <- c(dl_data, input_score=list("Currently, there are no input scores associated with this feature set on CaDrA portal."))
        
      }
      
    }
    
    return(as_attachment(dl_data, paste0(feature_set, ".rds")))
    
  }else{
    
    res$status <- 404  
    return(list(error=paste0(feature_set, " feature set is not available on the CaDrA portal")))
    
  }
  
}




