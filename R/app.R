
#' Shiny UI modules 
#' 
#' @param id A unique namespace identifier
#' @return Shiny ui elements
#'
#' @import shiny
#' 
#' @export 
CaDrA_UI <- function(id){
  
  ns <- NS(id)
  
  fluidRow(
    style = "padding: 5px 10px 50px 10px;",
    
    column(
      width = 4, 
      style = "border: 1px solid gray; border-radius: 3px; background: lightgrey; padding: 5px 10px 10px 10px;",
      
      h2("CaDrA Options", style="text-align: center;"),
      br(),
      
      tagList(
        fileInput(inputId = ns("ES_file"), label = strong(span(style = "color: red;", "*"), "Choose a binary feature file:"), width = "100%"),
        radioButtons(inputId = ns("ES_file_type"), label = "File type:", choices=c(".csv", ".rds"), selected = ".csv", inline = TRUE),
        helpText("Note: The binary feature file must contain unique features (rownames) which are used to search for best features"),
        
        br(),
        
        fileInput(inputId = ns("input_score_file"), label = strong(span(style = "color: red;", "*"), "Choose an input score file:"), width = "100%"),
        radioButtons(inputId = ns("input_score_file_type"), label = "File type:", choices=c(".csv", ".rds"), selected = ".csv", inline = TRUE),
        helpText("Note: The input score file must have names or labels that match the colnames of the binary feature file\n"),
        
        br(),
        
        radioButtons(inputId = ns("method"), label = strong(span(style="color:red;", "*"), "Scoring method:"), choices = c("ks", "wilcox", "revealer", "custom"), selected = "ks", inline = TRUE),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'ks'", ns("method")),
          checkboxInput(inputId = ns("weighted_ks"), label = "Compute the weighted ks?", value = FALSE), 
          conditionalPanel(
            condition = sprintf("input['%s'] == true", ns("weighted_ks")),
            fileInput(inputId = ns("weights_file"), label = strong(span(style = "color: red;", "*"), "Choose a weight file:"), width = "100%"),
            radioButtons(inputId = ns("weights_file_type"), label = "File type:", choices=c(".csv", ".rds"), selected = ".csv", inline = TRUE),        
            helpText("Note: Weights must have names or labels that match the names of the input score file"),
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'ks' | input['%s'] == 'wilcox'", ns("method"), ns("method")),
          selectInput(inputId = ns("alternative"), label = strong(span(style="color:red;", "*"), "Alternative:"), choices = c("less", "two.sided", "greater"), selected = "less", width = "100%"),
        ),
        fluidRow(
          column(
            width = 6,
            radioButtons(inputId = ns("metric"), label = strong(span(style="color:red;", "*"), "Type of metric:"), choices=c("pval", "stat"), selected = "pval", inline = FALSE)
          ),
          column(
            width = 6,
            radioButtons(inputId = ns("search_method"), label = strong(span(style="color:red;", "*"), "Search method:"), choices=c("forward"="forward", "forward and backward"="both"), selected = "forward", inline = FALSE)
          )
        ),
        numericInput(inputId = ns("max_size"), label = strong(span(style="color:red;", "*"), "Select a maximum size that a meta-feature can extend to do for a given search"), min = 1, max = 100, step = 1, value = 7, width = "100%"),
        radioButtons(inputId = ns("initial_seed"), label = strong(span(style="color:red;", "*"), "How to start the search?"), choices = c("Top N seeds"="top_N_seeds", "Known seeds"="search_start_seeds"), selected = "top_N_seeds", inline = TRUE),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'top_N_seeds'", ns("initial_seed")),
          numericInput(inputId = ns("top_N"), label = strong(span(style="color:red;", "*"), "Select a number to evaluate over the top N features"), min = 1, max = 100, step = 1, value = 10, width = "100%"),
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'search_start_seeds'", ns("initial_seed")),
          textAreaInput(inputId = ns("search_start"), label = strong(span(style = "color:red;", "*"), "Enter a list of integers/character strings (separated by commas) which specifies an index or feature name within the expression set object to start the search with"), value="", width="100%")
        ),
        checkboxInput(inputId = ns("permutation_test"), label = strong("Perform permutation testing?"), value = FALSE), 
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("permutation_test")),
          numericInput(inputId = ns("n_perm"), label = strong(span(style="color:red;", "*"), "Number of permutations to perform"), min = 1, max = 1000, step = 1, value = 10),
          textInput(inputId = ns("seed"), label = strong(span(style="color:red;", "*"), "A seed set for random permutation"), value = 123),
          numericInput(inputId = ns("ncores"), label = strong(span(style="color:red;", "*"), "Number of cores to perform parallelization for permutation testing"), min = 1, max = 10, step = 1, value = 1)
        ),
        
        br(),
        
        uiOutput(outputId = ns("error_message")),
        
        actionButton(inputId = ns("run_cadra"), label = strong("RUN"))
      )
    ),
    
    column(
      width = 8,
      style = "padding: 5px 10px 10px 10px;",

      tabsetPanel(
        id = "tabs",
        
        tabPanel(
          title = "Home", 
          style = "padding: 5px 10px 10px 10px;",
          
          div(
            h2("CaDrA"),
            p('CaDrA is an R package that supports a heuristic search framework aimed at identifying candidate drivers of a molecular phenotype of interest. The main function takes two inputs: i) a binary multi-omics dataset (where the rows are 1/0 vectors indicating the presence/absence of "omics" features such as somatic mutations, copy number alterations, epigenetic marks, etc.); and ii) and a molecular phenotype represented as a vector of continuous scores (sample-specific scores representing a phenotypic readout of interest, such as protein expression, pathway activity, etc.). Based on this input, CaDrA implements a forward/backward search algorithm to find the set of features that together is maximally associated with the observed input scores, based on one of several scoring functions (Kolmogorov-Smirnov, Conditional Mutual Information, Wilcoxon, custom-defined scoring function), making it useful to find complementary omics features likely driving the input molecular phenotype.'),
            p('For more information, please see the associated manuscript', a(target="_blank", href="https://www.frontiersin.org/articles/10.3389/fgene.2019.00121/full", 'Kartha et al. (2019)'))
          ) 
        ),
        
        tabPanel(
          title = "Data Formats", 
          style = "padding: 5px 10px 10px 10px;",
          
          h2("Datasets Required in CaDrA:"),
          tags$ol(
            tags$li("A binary features dataset (such as somatic mutations, copy number alterations, chromosomal translocations, etc.) The 1/0 vectors indicating the presence/absence of 'omics' features in the samples. The binary features matrix must be an object of class ExpressionSet from Biobase package)"),
            tags$li("A vector of continuous scores (or input_score) represents a functional response of interest (such as protein expression, pathway activity, etc.)")
          )
        ),
        
        tabPanel(
          title = "Run CaDrA", 
          style = "padding: 5px 10px 10px 10px;",
          
          div(
            uiOutput(outputId = ns("instructions"))
          ),
          
          div(
            uiOutput(outputId = ns("featureData_title")),
            DT::dataTableOutput(outputId = ns("featureData"))
          ),
          
          div(
            uiOutput(outputId = ns("inputScoreData_title")),
            DT::dataTableOutput(outputId = ns("inputScoreData"))
          ),
          
          div(
            uiOutput(outputId = ns("meta_plot_title")),
            plotOutput(outputId = ns("meta_plot"))
          ),
          
          div(
            uiOutput(outputId = ns("topn_plot_title")),
            plotOutput(outputId = ns("topn_plot"))
          ),
          
          div(
            uiOutput(outputId = ns("permutation_plot_title")),
            plotOutput(outputId = ns("permutation_plot"))
          )
          
        ),
        
        tabPanel(
          title = "Documentation", 
          style = "padding: 5px 10px 10px 10px;",
          
          h2("Installation"),
          p("You can install the development version of CaDrA from GitHub (Recommended)"), 
          
          tags$pre(
            tags$code(
            style="text-align: left; padding: 10px 10px 10px 10px;",
              "
              library(devtools)
              devtools::install_github('montilab/CaDrA', ref='dev')
              "
            )
          ),
          
          h2("Quickstart"),
          
          tags$pre(
            tags$code(
            style="text-align: left; padding: 10px 10px 10px 10px;",
              "
              library(CaDrA)
              library(Biobase)
              "
            )
          )
          
        ),
        
        tabPanel(
          title = "Publication", 
          style = "padding: 5px 10px 10px 10px;",
          
          h2("Citation"),
          p("Kartha VK, Kern JG, Sebastiani P, Zhang L, Varelas X, Monti S (2017) CaDrA: A computational framework for performing candidate driver analyses using binary genomic features.", a(href="https://www.biorxiv.org/content/early/2017/11/23/221846", "{bioRxiv}")),
          
          h2("Github"),
          p("Kartha V, Monti S, Chau R, Bulekova K (2022). CaDrA: Candidate Driver Analysis. R package version 2.0.0, ", a(target="_blank", href="https://github.com/montilab/CaDrA/", "https://github.com/montilab/CaDrA/"), "."),
          
          tags$pre(
            tags$code(
              "
              @Manual{,
                title = {CaDrA: Candidate Driver Analysis},
                author = {Vinay Kartha and Stefano Monti and Reina Chau and Katia Bulekova},
                year = {2022},
                note = {R package version 2.0.0},
                url = {https://github.com/montilab/CaDrA/},
              }
              "
            )
          )
          
        ),
        
        tabPanel(
          title = "Contract Us", 
          style = "padding: 5px 10px 10px 10px;",
          
          h2("Developers"),
          
          p(a(href="mailto:vkartha@bu.edu", strong("Vinay Kartha")), em(". Author.")),
          
          p(a(href="mailto:smonti@bu.edu", strong("Stefano Monti")), em(". Author.")),
          
          p(a(href="mailto:rchau88@bu.edu", strong("Reina Chau")), em(". Author, maintainer.")),
          
          p(a(href="mailto:ktrn@bu.edu", strong("Katia Bulekova")), em(". Author"))
          
        )
      )
    )
  )
  
}

#' Shiny Server modules 
#' 
#' @param id A unique namespace identifier
#' @return Shiny ui elements
#'
#' @import shiny dplyr readr tibble DT methods
#' 
#' @export 
CaDrA_Server <- function(id){
  
  moduleServer(
    id,
    function(input, output, session) {

      # Create reative values
      candidate_search_result <- reactiveVal()
      permutation_result <- reactiveVal()  
      instructions_message <- reactiveVal(TRUE)
      error_message <- reactiveVal()
      
      # Output instructions message
      output$instructions <- renderUI({
        
        req(instructions_message())
        
        div(
          h2("Instructions"),
          
          tags$pre(
            tags$code(
              "
                Select the 'CaDrA Options' on the right and Click 'RUN' button on the bottom
                "
            )
          )
        )
        
      })

      observeEvent(input$run_cadra, {
        
        candidate_search_result(NULL)
        permutation_result(NULL) 
        instructions_message(FALSE)
        error_message(NULL)
        
        inputfile <- input$ES_file;
        inputtype <- input$ES_file_type;
        
        if(is.null(inputfile)){
          error_message("Please choose a file to import.")
          return(NULL)
        }
        
        csv_ext <-  grep(toupper(".csv"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
        rds_ext <-  grep(toupper(".rds"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
        
        if(inputtype %in% ".csv" & length(csv_ext) > 0){
          
          # read in the Eset file
          Eset <- readr::read_csv(inputfile$datapath) %>% column_to_rownames(var="Features") %>% mutate_all(as.integer)
          
          # convert Eset to matrix
          Eset <- as.matrix(Eset, nrow=nrow(Eset), ncol=ncol(Eset), byrow=TRUE, dimnames=list(rownames(Eset), colnames(Eset)))
          
          #create phenotypic data
          pData <- data.frame(Samples = colnames(Eset), stringsAsFactors=TRUE)
          rownames(pData) <- pData$Samples
          phenoData <- new("AnnotatedDataFrame", data=pData)
          
          #create feature data
          fData <- data.frame(Features = rownames(Eset), stringsAsFactors = TRUE)
          rownames(fData) <- fData$Features
          featureData <- new("AnnotatedDataFrame", data=fData)
          
          #create expression set
          ES <- ExpressionSet(assayData=Eset, phenoData=phenoData, featureData=featureData)
          
        }else if(inputtype %in% ".rds" & length(rds_ext) > 0){
          
          ES <- readRDS(inputfile$datapath)
          
        }else{
          
          error_message("Incorrect file format. Please check your file again.")
          return(NULL)
          
        }
        
        inputfile <- input$input_score_file;
        inputtype <- input$input_score_file_type;
        
        if(is.null(inputfile)){
          error_message("Please choose a file to import.")
          return(NULL)
        }
        
        csv_ext <-  grep(toupper(".csv"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
        rds_ext <-  grep(toupper(".rds"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
        
        if(inputtype %in% ".csv" & length(csv_ext) > 0){
          
          dat <- readr::read_csv(inputfile$datapath) %>% column_to_rownames(var="Samples") %>% mutate_all(as.numeric)
          input_score <- as.numeric(dat$Scores)
          names(input_score) <- rownames(dat)
          
        }else if(inputtype %in% ".rds" & length(rds_ext) > 0){
          
          dat <- readRDS(inputfile$datapath) %>% column_to_rownames(var="Samples") %>% mutate_all(as.numeric)
          input_score <- as.numeric(dat$Scores)
          names(input_score) <- rownames(dat)
          
        }else{
          
          error_message("Incorrect file format. Please check your file again.")
          return(NULL)
          
        }
        
        req(ES, input_score)
        
        # Check if the ES is provided and is a BioBase ExpressionSet object
        if(length(ES) == 0 || class(ES)[1] != "ExpressionSet") 
          error_message("'ES' must be an ExpressionSet class argument (required).")
        
        # Check if the dataset has only binary 0 or 1 values 
        if(!all(exprs(ES) %in% c(0,1))){
          error_message("The expression matrix (ES) must contain only binary values (0/1) with no NAs.\n")
        }
        
        # Make sure the input ES has rownames for features tracking
        if(is.null(rownames(ES)))
          error_message("The ES object does not have rownames or featureData to track the features by. Please provide unique features or rownames for the expression matrix.\n")
        
        # Check input_score is provided and are continuous values with no NAs
        if(length(input_score) == 0 || any(!is.numeric(input_score)) || any(is.na(input_score)))
          error_message("input_score must be a vector of continous values (with no NAs) where the vector names match the colnames of the expression matrix.\n")
        
        # Make sure the input_score has names or labels that are the same as the colnames of ES
        if(is.null(names(input_score)))
          error_message("The input_score object must have names or labels to track the samples by. Please provide unique sample names or labels that matches the colnames of the expression matrix.\n")
        
        # Make sure the input_score has the same length as number of samples in ES
        if(length(input_score) != ncol(ES)){
          error_message("The input_score must have the same length as the number of columns in ES.\n")
        }else{
          if(any(!names(input_score) %in% colnames(ES))){
            error_message("The input_score object must have names or labels that match the colnames of the expression matrix.\n")
          }
          # match colnames of expression matrix with names of provided input_score
          ES <- ES[,names(input_score)]
        }
        
        # Get method
        method = input$method; 
        
        if(method == "ks"){
          
          if(input$weighted_ks == TRUE){
            
            inputfile <- input$weights_file;
            inputtype <- input$weights_file_type;
            
            if(is.null(inputfile)){
              error_message("Please choose a file to import.")
              return(NULL)
            }
            
            csv_ext <-  grep(toupper(".csv"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
            rds_ext <-  grep(toupper(".rds"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
            
            if(inputtype %in% ".csv" & length(csv_ext) > 0){
              
              dat <- read_csv(inputfile$datapath) %>% column_to_rownames(var="Samples") %>% mutate_all(as.numeric)
              weights <- as.numeric(dat$Weights)
              names(weights) <- rownames(dat)
              
            }else if(inputtype %in% ".rds" & length(rds_ext) > 0){
              
              dat <- readRDS(inputfile$datapath) %>% column_to_rownames(var="Samples") %>% mutate_all(as.numeric)
              weights <- as.numeric(dat$Weights)
              names(weights) <- rownames(dat)
              
            }else{
              
              error_message("Incorrect file format. Please check your file again.")
              return(NULL)
              
            }
            
            # Check weights is provided and are continuous values with no NAs
            if(length(weights) == 0 || any(!is.numeric(weights)) || any(is.na(weights)))
              error_message("weights must be a vector of continous values (with no NAs) where the vector names match the colnames of the expression matrix.\n")
            
            # Make sure the weights has names or labels that are the same as the colnames of ES
            if(is.null(names(weights)))
              error_message("The weights object must have names or labels to track the samples by. Please provide unique sample names or labels that matches the colnames of the expression matrix.\n")
            
            # Make sure the weights has the same length as number of samples in ES
            if(length(weights) != ncol(ES)){
              error_message("The weights must have the same length as the number of columns in ES.\n")
            }else{
              if(any(!names(weights) %in% colnames(ES))){
                error_message("The weights object must have names or labels that match the colnames of the expression matrix.\n")
              }
              # match colnames of expression matrix with names of provided weights
              weights <- weights[colnames(ES)]
            }
            
          }else{
            
            weights <- NULL
            
          }
          
        }
        
        if(method %in% c("ks", "wilcox")){
          alternative <- input$alternative
        }
        
        metric = input$metric;
        
        search_method = input$search_method;
        
        max_size = as.numeric(input$max_size)
        
        if(is.na(max_size) || length(max_size)==0){
          error_message("Please specify an integer value specifies a maximum size that a meta-feature can extend to do for a given search.\n")
          return(NULL)
        }
        
        initial_seed = input$initial_seed
        
        if(initial_seed == "top_N_seeds"){
          
          search_start = NULL
          top_N = as.numeric(input$top_N)
          
          if(is.na(top_N) || length(top_N)==0){
            error_message("Please specify a NUMERIC top_N value to evaluate over top N features.\n")
            return(NULL)
          }
          
          
        }else{
          
          search_start = input$search_start
          top_N = NULL
          
          if(is.numeric(search_start)){ 
            # User-specified feature index (has to be an integer from 1:nrow(ES))
            verbose("Starting with specified sorted feature indexes...\n")
            
            if(search_start > nrow(ES)){ # Index out of range
              error_message("Invalid starting index specified... Please specify a valid starting index within the range of the existing ES...\n")
              return(NULL)
            }
          }
          
          if(is.character(search_start)){
            # User-specified feature name (has to be a character from rownames(1:nrow(ES)))
            verbose("Starting with specified feature names...\n")
            
            if(!(search_start %in% rownames(ES))){ #provided feature name not in rownames
              error_message("Provided starting feature does not exist among ES's rownames.\n\n")
              return(NULL)
            }
          }
          
        }
        
        candidate_search_result(
          CaDrA::candidate_search(
            ES = ES,
            input_score = input_score,
            method = method,  
            custom_function = NULL,
            custom_parameters = NULL,
            alternative = alternative,        
            metric = metric,             
            weights = weights,   
            search_start = search_start,
            top_N = top_N,
            max_size = max_size,               
            best_score_only = FALSE,
            do_plot = FALSE,
            verbose = FALSE
          )
        ) 
        
        permutation = input$permutation_test
        
        if(permutation == TRUE){
          
          n_perm = as.numeric(input$n_perm)
          
          if(is.na(top_N) || length(top_N)==0){
            error_message("Please specify a NUMERIC top_N value to evaluate over top N features.\n")
            return(NULL)
          }
          
          seed = as.numeric(seed)
          
          if(is.na(top_N) || length(top_N)==0){
            error_message("Please specify a NUMERIC top_N value to evaluate over top N features.\n")
            return(NULL)
          }
          
          ncores = as.numeric(ncores)
          
          if(is.na(top_N) || length(top_N)==0){
            error_message("Please specify a NUMERIC top_N value to evaluate over top N features.\n")
            return(NULL)
          }
          
          permutation_result(
            CaDrA::CaDrA(
              ES = ES, 
              input_score = input_score, 
              method = method,                
              custom_function = NULL,
              custom_parameters = NULL,
              alternative = alternative,        
              metric = metric,             
              weights = weights,  
              top_N = top_N,
              search_start = search_start,
              search_method = search_method,      
              max_size = max_size,               
              n_perm = n_perm,               
              seed = seed,                  
              smooth = TRUE,
              obs_best_score = NULL,
              plot = FALSE, 
              ncores = ncores,
              cache_path = NULL,
              verbose = FALSE
            )
          )
          
        }
        
        error_message("Analysis is Completed! See 'Run CaDrA'")
        
      })

      output$error_message <- renderUI({
        
        req(error_message())
        
        p(style="color: red; font-weight: bold;", error_message())
        
      })
      
      output$featureData_title <- renderUI({
        
        req(candidate_search_result())
        
        h2("Binary Meta-Feature")
        
      })
      
      output$featureData <- DT::renderDataTable({
        
        req(candidate_search_result())
        
        Eset <- candidate_search_result()[[1]][["ESet"]]
        
        table <- exprs(Eset) %>% as.data.frame() %>% 
          datatable(
            rownames = TRUE,
            extensions = 'Buttons',
            selection = "single",
            options = list(
              deferRender = FALSE,
              paging = TRUE,
              searching = TRUE,
              ordering = TRUE,
              pageLength = 20,
              scrollX = TRUE,
              scrollY = 400,
              scrollCollapse = TRUE,
              dom = 'T<"clear">Blfrtip',
              buttons=c('copy','csv','print')
            )
          )
        
        return(table)
        
      })
      
      output$inputScoreData_title <- renderUI({
        
        req(candidate_search_result())
        
        h2("Observed Input Scores")
        
      })
      
      output$inputScoreData <- DT::renderDataTable({
        
        req(candidate_search_result())
        
        Score <- candidate_search_result()[[1]][["input_score"]]
        
        table <- matrix(Score, nrow=1, ncol=length(Score), byrow=T, dimnames=list("input_score", names(Score))) %>% 
          datatable(
            rownames = TRUE,
            extensions = 'Buttons',
            selection = "single",
            options = list(
              deferRender = FALSE,
              paging = TRUE,
              searching = TRUE,
              ordering = TRUE,
              pageLength = 20,
              scrollX = TRUE,
              scrollY = 400,
              scrollCollapse = TRUE,
              dom = 'T<"clear">Blfrtip',
              buttons=c('copy','csv','print')
            )
          )
        
        return(table)
        
      })
      
      output$meta_plot_title <- renderUI({
        
        req(candidate_search_result())
        
        h2("Meta-Features Plot")
        
      })
      
      output$meta_plot <- renderPlot({
        
        req(candidate_search_result())
        
        topn_res <- candidate_search_result()
        topn_best_meta <- CaDrA::topn_best(topn_res)
        
        CaDrA::meta_plot(topn_best_list = topn_best_meta)
        
      })
      
      output$topn_plot_title <- renderUI({
        
        req(candidate_search_result())
        
        h2("Top N Overlapping Heatmap")
        
      })
      
      output$topn_plot <- renderPlot({
        
        req(candidate_search_result())

        topn_res <- candidate_search_result()
        CaDrA::topn_plot(topn_res)
        
      })      
      
      output$permutation_plot_title <- renderUI({
        
        req(permutation_result())
        
        h2("Permutation-Based Testing")
        
      })
      
      output$permutation_plot <- renderPlot({
        
        req(permutation_result())

        perm_res <- permutation_result()
        permutation_plot(perm_res)
        
      })
      
      candidate_search_result
      
    }
  )
  
}

#' Run both Shiny UI and Server Modules 
#' 
#' @return Shiny application 
#'
#' @import shiny
#' 
#' @export
CaDrA_App <- function() {
  
  ui <- fluidPage(
    
    shinyjs::useShinyjs(),
    
    titlePanel("CaDrA: Candidate Drivers Analysis"),
    
    helpText("Multi-Omic Search for Candidate Drivers of Functional Signatures"),

    CaDrA_UI(id = "CaDrA")
    
  )
  
  server <- function(input, output, session) {
    
    CaDrA_Server(id = "CaDrA")
    
  }
  
  shinyApp(ui=ui, server=server)
  
}







