
# function to hover columns in data table
create_hover_txt <- function(table){
  
  column_names <- colnames(table)
  
  th_tr <- lapply(1:length(column_names), function(l){ 
    title = column_names[l]
    name = ifelse(nchar(title) > 4, paste0(substr(title, 1, 4), "..."), title)
    th <- sprintf('<th title = "%s">%s</th>\n', title, name) 
  }) %>% purrr::flatten_chr() %>% paste0(., collapse = "")
  
  th_tr <- paste0('<th title=""></th>\n', th_tr) %>% HTML()
  
  sketch = htmltools::withTags(
    tags$table(
      class = 'display',
      tags$thead(
        th_tr
      )
    )
  )
  
  return(sketch)
  
}


#' Shiny UI modules 
#' 
#' @param id A unique namespace identifier
#' @return Shiny ui elements
#'
#' @import shiny
#' @importFrom markdown markdownToHTML
#' 
#' @export 
CaDrA_UI <- function(id){
  
  ns <- NS(id)
  
  fluidRow(
    style = "padding: 5px 10px 10px 10px;",
    
    column(
      width = 4, 
      style = "border: 1px solid gray; border-radius: 3px; background: lightgrey; padding: 5px 10px 10px 10px; min-height: 850px;",
      
      h2("CaDrA Options", style="text-align: center;"),
      br(),
      
      tagList(
        selectInput(inputId = ns("dataset"), label="Choose a dataset:", choices=c("Oncogenic YAP/TAZ Activity in Human Breast Cancer"="BRCA_GISTIC_MUT_SIG", "Transcriptional Activation of B-catenin in Cancer"="CCLE_MUT_SCNA", "Simulated Data"="sim.Data", "Import Data"), selected="BRCA_GISTIC_MUT_SIG", width = "100%"),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Import Data'", ns("dataset")),
          
          fileInput(inputId = ns("ES_file"), label = strong(span(style = "color: red;", "*"), "Choose a binary feature file:"), width = "100%"),
          radioButtons(inputId = ns("ES_file_type"), label = "File type:", choices=c(".csv", ".rds"), selected = ".csv", inline = TRUE),
          helpText("Note: The binary feature file must contain unique features (rownames) which are used to search for best features"),
          
          br(),
          
          fileInput(inputId = ns("input_score_file"), label = strong(span(style = "color: red;", "*"), "Choose an input score file:"), width = "100%"),
          radioButtons(inputId = ns("input_score_file_type"), label = "File type:", choices=c(".csv", ".rds"), selected = ".csv", inline = TRUE),
          helpText("Note: The input score file must have names or labels that match the colnames of the binary feature file\n")
        ),
        
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
            radioButtons(inputId = ns("search_method"), label = strong(span(style="color:red;", "*"), "Search method:"), choices=c("forward and backward"="both", "forward"="forward"), selected = "both", inline = FALSE)
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
          numericInput(inputId = ns("n_perm"), label = strong(span(style="color:red;", "*"), "Number of permutations to perform"), min = 1, max = Inf, step = 1, value = 100),
          numericInput(inputId = ns("seed"), label = strong(span(style="color:red;", "*"), "A seed set for random permutation"), min = 1, max = Inf, step = 1, value = 123),
          numericInput(inputId = ns("ncores"), label = strong(span(style="color:red;", "*"), "Number of cores to perform parallelization for permutation testing"), min = 1, max = Inf, step = 1, value = 1)
        ),
        
        br(),
        
        uiOutput(outputId = ns("error_message")),
        
        actionButton(inputId = ns("run_cadra"), label = strong("RUN"), style="background: blue; color: white;"),
        
        br(), br(), br(), br(),
        
        HTML("<p style='text-align: center;'><span class='footer-info'>&copy; Monti Lab &diams; <script>document.write(new Date().getFullYear());</script> &diams; All Rights Reserved.</span></p>")
      )
    ),
    
    column(
      width = 8,
      style = "border: 1px solid gray; padding: 5px 10px 10px 10px; min-height: 850px;",
      
      tabsetPanel(
        id = "tabs",
        
        tabPanel(
          title = "Run CaDrA", 
          style = "padding: 5px 10px 10px 10px;",
          icon = icon(name = "running", lib = "font-awesome"),
          
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
          title = "Help", 
          style = "padding: 5px 10px 10px 10px;",
          icon = icon(name = "question", lib = "font-awesome"),
          
          HTML(markdown::markdownToHTML('README.md', fragment.only=T))
          
        ),
        
        tabPanel(
          title = "Publication", 
          style = "padding: 5px 10px 10px 10px;",
          icon = icon(name = "book", lib = "font-awesome"),
          
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
          title = "Contact Us", 
          style = "padding: 5px 10px 10px 10px;",
          icon = icon(name = "envelope", lib = "font-awesome"),
          
          h2("Developers"),
          
          p(a(href="mailto:rchau88@bu.edu", strong("Reina Chau")), em(". Author, maintainer.")),
          
          p(a(href="mailto:ktrn@bu.edu", strong("Katia Bulekova")), em(". Author")),
          
          p(a(href="mailto:vkartha@bu.edu", strong("Vinay Kartha")), em(". Author.")),
          
          p(a(href="mailto:smonti@bu.edu", strong("Stefano Monti")), em(". Author."))
          
        )
      )
    )
  )
  
}

#' Shiny Server modules 
#' 
#' @param id A unique namespace identifier
#' @return Shiny server elements
#'
#' @import shiny methods htmltools
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom dplyr mutate_all
#' @importFrom stats rnorm 
#' @importFrom utils read.csv write.csv
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
        
        if(input$dataset == "BRCA_GISTIC_MUT_SIG"){
          
          ## Read in BRCA GISTIC+Mutation ESet object
          eset_mut_scna <-  CaDrA::BRCA_GISTIC_MUT_SIG
          
          ## Read in input score
          input_scores <-  CaDrA::TAZYAP_BRCA_ACTIVITY
          
          ## Samples to keep based on the overlap between the two inputs
          overlap <- intersect(names(input_scores), Biobase::sampleNames(eset_mut_scna))
          eset_mut_scna <- eset_mut_scna[,overlap]
          input_scores <- input_scores[overlap]
          
          ## Binarize ES to only have 0's and 1's
          exprs(eset_mut_scna)[exprs(eset_mut_scna) > 1] <- 1.0
          
          ## Pre-filter ESet based on occurrence frequency
          eset_mut_scna_flt <- CaDrA::prefilter_data(
            ES = eset_mut_scna,
            max.cutoff = 0.6, # max frequency (60%)
            min.cutoff = 0.03 # min frequency (3%)
          ) 
          
          ES <- eset_mut_scna_flt
          input_score <- input_scores
          
        }else if(input$dataset == "CCLE_MUT_SCNA"){
          
          ES <- CaDrA::CCLE_MUT_SCNA
          input_score <- CaDrA::CTNBB1_reporter        
          
        }else if(input$dataset == "sim.Data"){
          
          ES <- CaDrA::sim.ES
          
          # set seed
          set.seed(123)
          
          # Provide a vector of continuous scores for a target profile with names to each score value
          input_score = rnorm(n = ncol(ES))
          names(input_score) <- colnames(ES)
          
        }else if(input$dataset == "Import Data"){
          
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
            Eset <- read.csv(inputfile$datapath, header=T) %>% tibble::column_to_rownames(var="Features") %>% dplyr::mutate_all(as.integer)
            
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
            
            dat <- read.csv(inputfile$datapath, header = T) %>% tibble::column_to_rownames(var="Samples") %>% dplyr::mutate_all(as.numeric)
            input_score <- as.numeric(dat$Scores)
            names(input_score) <- rownames(dat)
            
          }else if(inputtype %in% ".rds" & length(rds_ext) > 0){
            
            dat <- readRDS(inputfile$datapath) %>% tibble::column_to_rownames(var="Samples") %>% dplyr::mutate_all(as.numeric)
            input_score <- as.numeric(dat$Scores)
            names(input_score) <- rownames(dat)
            
          }else{
            
            error_message("Incorrect file format. Please check your file again.")
            return(NULL)
            
          }
          
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
              
              dat <- read.csv(inputfile$datapath, header=T) %>% tibble::column_to_rownames(var="Samples") %>% mutate_all(as.numeric)
              weights <- as.numeric(dat$Weights)
              names(weights) <- rownames(dat)
              
            }else if(inputtype %in% ".rds" & length(rds_ext) > 0){
              
              dat <- readRDS(inputfile$datapath) %>% tibble::column_to_rownames(var="Samples") %>% mutate_all(as.numeric)
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
          
          seed = as.numeric(input$seed)
          
          if(is.na(top_N) || length(top_N)==0){
            error_message("Please specify a NUMERIC top_N value to evaluate over top N features.\n")
            return(NULL)
          }
          
          ncores = as.numeric(input$ncores)
          
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
        
        error_message("")
        
      })
      
      output$error_message <- renderUI({
        
        req(error_message())
        
        p(style="color: red; font-weight: bold;", error_message())
        
      })
      
      output$featureData_title <- renderUI({
        
        req(candidate_search_result())
        
        dataset <- isolate({ input$dataset })
        
        if(dataset == "BRCA_GISTIC_MUT_SIG"){
          
          title <- "Dataset: Oncogenic YAP/TAZ Activity in Human Breast Cancer"
          description <- "An ExpressionSet object consists of 16,873 genomic features across 963 samples."
          
        }else if(dataset == "CCLE_MUT_SCNA"){
          
          title <- "Dataset: Transcriptional Activation of B-catenin in Cancer"
          description <- "An ExpressionSet object consists of 17,723 genomic features across 82 samples."
          
        }else if(dataset == "sim.Data"){
          
          title <- "Dataset: Simulated Data"
          description <- "A simulated ExpressionSet that comprises of 1000 genomic features (rows) and 200 sample profiles (columns). Each row feature is represented by a vector of binary number (1/0) indicating the presence/absence of the feature in the samples. This simulated data includes 10 left-skewed (i.e. True Positive or TP) and 990 uniformly-distributed (i.e. True Null or TN) features."
          
        }else if(dataset == "Import Data"){
          
          title <- "Dataset: Imported Data"
          description <- ""
          
        }
        
        div(
          h2(title),
          br(),
          p(description),
          h2("Best Meta-Features Set")
        )
        
      })
      
      output$featureData <- DT::renderDataTable({
        
        req(candidate_search_result())
        
        ns <- session$ns
        
        topn_best_meta <- topn_best(topn_list=candidate_search_result())
        
        ES_table <- topn_best_meta[["ESet"]] %>% exprs(.) %>% as.data.frame(.)
        
        hover_columns <- create_hover_txt(table = ES_table)
        
        table <- ES_table %>% 
          DT::datatable(
            container = hover_columns,
            rownames = TRUE,
            extensions = 'Buttons',
            selection = "single",
            options = list(
              deferRender = FALSE,
              paging = FALSE,
              searching = TRUE,
              ordering = TRUE,
              pageLength = 20,
              scrollX = TRUE,
              scrollY = 400,
              scrollCollapse = TRUE,
              dom = 'T<"clear">Blfrtip',
              buttons = list(
                list(
                  extend = "collection",
                  text = 'Download All Results',
                  action = DT::JS(
                    sprintf(
                      paste0(
                        "function ( e, dt, node, config ) {",
                        "Shiny.setInputValue('%s', true, {priority: 'event'});",
                        "}"
                      ), ns('Download_Eset')
                    )
                  )
                )
              )
            )
          )
        
        return(table)
        
      })
      
      observeEvent(input$Download_Eset, {
        
        ns <- session$ns
        
        shiny::showModal(
          shiny::modalDialog(
            title = "Download Best Meta-Features Set",
            downloadButton(outputId = ns("downloadEsetCSV"), "Download Table as CSV file"),
            br(), br(),
            downloadButton(outputId = ns("downloadEsetRDS"), "Download Table as RDS file"),
          )
        )
        
      })
      
      ## Download CSV File ####
      output$downloadEsetCSV <- downloadHandler(
        
        filename = function() {
          paste0("CaDrA-Best-Meta-Features-Eset.csv")
        },
        
        content = function(file) {
          
          topn_best_meta <- topn_best(topn_list=candidate_search_result())
          
          Eset_table <- topn_best_meta[["ESet"]] %>% exprs(.) %>% as.data.frame(.) %>% tibble::rownames_to_column(., var="Features")
          
          write.csv(Eset_table, file, row.names=F)
          
        }
        
      )
      
      ## Download RDS File ####
      output$downloadEsetRDS <- downloadHandler(
        
        filename = function() {
          paste0("CaDrA-Best-Meta-Features-Eset.rds")
        },
        
        content = function(file) {
          
          topn_best_meta <- topn_best(topn_list=candidate_search_result())
          
          Eset_table <- topn_best_meta[["ESet"]]
          
          saveRDS(Eset_table, file)
          
        }
        
      )
      
      output$inputScoreData_title <- renderUI({
        
        req(candidate_search_result())
        
        h2("Observed Input Scores")
        
      })
      
      output$inputScoreData <- DT::renderDataTable({
        
        req(candidate_search_result())
        
        ns <- session$ns
        
        topn_best_meta <- topn_best(topn_list=candidate_search_result())
        
        Scores <- topn_best_meta[["input_score"]] %>% signif(., digits = 4)
        
        score_table <- matrix(Scores, nrow=1, ncol=length(Scores), byrow=T, dimnames=list("input_score", names(Scores)))
        
        hover_columns <- create_hover_txt(table = score_table)
        
        table <- score_table %>% 
          DT::datatable(
            container = hover_columns,
            rownames = TRUE,
            extensions = 'Buttons',
            selection = "single",
            options = list(
              deferRender = FALSE,
              paging = FALSE,
              searching = TRUE,
              ordering = TRUE,
              pageLength = 20,
              scrollX = TRUE,
              scrollY = 400,
              scrollCollapse = TRUE,
              dom = 'T<"clear">Blfrtip',
              buttons = list(
                list(
                  extend = "collection",
                  text = 'Download All Results',
                  action = DT::JS(
                    sprintf(
                      paste0(
                        "function ( e, dt, node, config ) {",
                        "Shiny.setInputValue('%s', true, {priority: 'event'});",
                        "}"
                      ), ns('Download_InputScore')
                    )
                  )
                )
              )
            )
          )
        
        return(table)
        
      })  
      
      
      observeEvent(input$Download_InputScore, {
        
        ns <- session$ns
        
        shiny::showModal(
          shiny::modalDialog(
            title = "Download Observed Input Scores",
            downloadButton(outputId = ns("downloadScoreCSV"), "Download Table as CSV file"),
            br(), br(),
            downloadButton(outputId = ns("downloadScoreRDS"), "Download Table as RDS file"),
          )
        )
        
      })
      
      ## Download CSV File ####
      output$downloadScoreCSV <- downloadHandler(
        
        filename = function() {
          paste0("CaDrA-Observed-Input-Scores.csv")
        },
        
        content = function(file) {
          
          topn_best_meta <- topn_best(topn_list=candidate_search_result())
          
          Scores <- topn_best_meta[["input_score"]]
          
          score_table <- data.frame(
            Samples = names(Scores),
            Scores = Scores,
            stringsAsFactors = F
          )
          
          write.csv(score_table, file, row.names=F)
          
        }
        
      )
      
      ## Download RDS File ####
      output$downloadScoreRDS <- downloadHandler(
        
        filename = function() {
          paste0("CaDrA-Observed-Input-Scores.rds")
        },
        
        content = function(file) {
          
          topn_best_meta <- topn_best(topn_list=candidate_search_result())
          Scores <- topn_best_meta[["input_score"]]
          
          score_table <- data.frame(
            Samples = names(Scores),
            Scores = Scores,
            stringsAsFactors = F
          )
          
          saveRDS(score_table, file)
          
        }
        
      )
      
      output$meta_plot_title <- renderUI({
        
        req(candidate_search_result())
        
        h2("Best Meta-Features Plot")
        
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
    
    titlePanel("CaDrA: Candidate Drivers Analysis"),
    
    helpText("Multi-Omic Search for Candidate Drivers of Functional Signatures"),
    
    CaDrA_UI(id = "CaDrA")
    
  )
  
  server <- function(input, output, session) {
    
    CaDrA_Server(id = "CaDrA")
    
  }
  
  shinyApp(ui=ui, server=server)
  
}







