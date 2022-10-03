
# function to hover columns in data table
create_hover_txt <- function(table){
  column_names <- colnames(table)
  th_tr <- lapply(seq_along(column_names), function(l){ 
    title <- column_names[l]
    name <- ifelse(nchar(title) > 4, paste0(substr(title, 1, 4), "..."), title)
    th <- sprintf('<th title = "%s">%s</th>\n', title, name) 
  }) %>% purrr::flatten_chr() %>% paste0(., collapse = "")
  th_tr <- paste0('<th title=""></th>\n', th_tr) %>% HTML()
  sketch <- htmltools::withTags(
    tags$table(
      class = 'display',
      tags$thead(
        th_tr
      )
    )
  )
  return(sketch)
}

# Global expression sets
dataset_choices <- list(
  "CCLE MUT + SCNAs in Cancers (CCLE_MUT_SCNA)" =  system.file("data/CCLE_MUT_SCNA.rda", package = "CaDrA"),
  "Simulated Expression Set (sim.ES)" = system.file("data/sim.ES.RData", package = "CaDrA"),
  "TCGA BRCA GISTIC + Mutation Signatures (BRCA_GISTIC_MUT_SIG)" = system.file("data/BRCA_GISTIC_MUT_SIG.rda", package = "CaDrA")
)

# Global input scores
score_choices <- list(
  "Activation of B-catenin in Cancers (CTNBB1_reporter)" = system.file("data/CTNBB1_reporter.rda", package = "CaDrA"),
  "Simulated Input Scores from rnorm(n=length(sim.ES), mean=0, sd=1) (sim.Scores)" =  system.file("data/sim.Scores.rda", package = "CaDrA"),
  "YAP/TAZ Activity in Human Breast Cancers (TAZYAP_BRCA_ACTIVITY)" = system.file("data/TAZYAP_BRCA_ACTIVITY.rda", package = "CaDrA")
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
  
#' Shiny UI modules 
#' 
#' @param id A unique namespace identifier
#' @return Shiny ui elements
#'
#' @examples
#' 
#' # Load R library
#' library(shiny)
#'
#' # Create ui and server for Shiny app
#' id <- "myapp"
#' 
#' ui <- fluidPage(
#'    CaDrA::CaDrA_UI(id = id)
#' )
#' server <- function(input, output, session) {
#'    CaDrA::CaDrA_Server(id = id)
#' }
#' app <-  shinyApp(ui=ui, server=server)
#' 
#' # Launch Shiny app (NOT RUN)
#' # shiny::runApp(app, host='0.0.0.0', port=3838)
#' 
#' @import shiny markdown
#' @importFrom htmltools includeMarkdown
#' 
#' @export 
CaDrA_UI <- function(id) 
{
  
  # Combine extdata with global expression set and scores dataset if it was provided
  eset_choices <- get_extdata(dataset_choices, score_choices)[["eset_choices"]] %>% unlist()
  
  ns <- NS(id)
  
  fluidRow(
    style = "padding: 5px 10px 10px 10px;",
    
    tags$style(
      HTML(
        "
        .side-bar-options {
          border: 1px solid gray; 
          border-radius: 3px; 
          background: lightgrey; 
          padding: 5px 10px 10px 10px; 
          min-height: 850px;
        }
        .tooltip-txt {
          color: red;
          font-weight: bold;
          width: 70px;
        }
        .loading_div {
          display: flex;
          text-align: center;
          align-items: center;
          justify-content: center;
          width: 100%;
          height: 400px;
        }
        .loading_text {
          text-align: center;
          margin-left: -30px;
        }
        .loader {
          border: 16px solid #f3f3f3;
          border-top: 16px solid #3498db;
          border-radius: 50%;
          width: 120px;
          height: 120px;
          animation: spin 2s linear infinite;
        }
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }
        "
      )
    ),
    tags$script(
      HTML(
        paste0(
          "Shiny.addCustomMessageHandler('ToggleOperation', function(message){",
          "var x = document.getElementById(message.id);",
            "if (message.display === 'yes') {",
              "x.style.display = 'flex';",
            "} else {",
              "x.style.display = 'none';",
            "}",
          "})"
        )
      )
    ),
    column(
      width = 4, 
      class = "side-bar-options",
      
      h2("CaDrA Options", style="text-align: center;"),
      
      br(),
      
      tagList(
        selectInput(
          inputId = ns("dataset"), 
          label = "Feature Set", 
          choices = c(eset_choices, "Import Data"),
          selected = eset_choices[1], 
          width = "100%"
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Import Data'", ns("dataset")),
          fileInput(
            inputId = ns("ES_file"), 
            label = strong(span(style = "color: red;", "*"), 
                           "Feature Set file:"), 
            width = "100%"
          ),
          radioButtons(
            inputId = ns("ES_file_type"), 
            label = HTML(paste0(
              'File type ', 
              '<a class="tooltip-txt" data-html="true" ',
              'data-tooltip-toggle="tooltip" data-placement=',
              '"top" title=\"NOTE: If file is in csv format, ',
              'the \'Feature Set\' must be a data ',
              'frame including a \'Features\' column name ',
              'that contains unique names or labels to ',
              'search for best features. Otherwise, \'Feature ',
              'Set\' must be an object of class ExpressionSet ',
              'from BioBase package.\">?</a>')), 
            choices = c(".csv", ".rds"), 
            selected = ".csv", 
            inline = TRUE
          )
        ),
        selectInput(
          inputId = ns("scores"), 
          label = "Input Score", 
          choices = "Import Data",
          width = "100%"
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Import Data'", ns("scores")),
          fileInput(
            inputId = ns("input_score_file"), 
            label = strong(span(style = "color: red;", "*"), 
                           "Input Score file:"), 
            width = "100%"
          ),
          radioButtons(
            inputId = ns("input_score_file_type"), 
            label = HTML(paste0(
              'File type ', 
              '<a class="tooltip-txt" data-html="true" ',
              'data-tooltip-toggle="tooltip" data-placement=',
              '"top" title=\"NOTE: If file is in csv format, ', 
              'then the \'Input Score\' file ',
              'must be a data frame with two columns ',
              '(Samples and Scores) and the \'Samples\' column ',
              'must match the colnames of \'Feature Set\'. ',
              'Otherwise, \'Input Score\' must be a list of  ',
              'vectors and have names or labels that match the ',
              'colnames of the \'Feature Set\'.\">?</a>')), 
            choices = c(".csv", ".rds"), 
            selected = ".csv", 
            inline = TRUE
          )
        ),
        numericInput(
          inputId = ns("min_cutoff"), 
          label = HTML(paste0(
            '<strong>Min Event Frequency (n)</strong> ', 
            '<a class="tooltip-txt" data-html="true" ',
            'data-tooltip-toggle="tooltip" data-placement="top" ',
            'title=\"Minimum number of \'occurrences\' a feature ',
            '(e.g., a mutation) must have to be included in the ',
            '\'Feature Set\`. Features with fewer events than the ',
            'specified number will be removed.\n\nNOTE: \'Min event ',
            'frequency\' must be >= 5.\">?</a>')),
          value = 30,
          min = 5, 
          max = Inf, 
          step = 1,
          width = "100%"
        ),
        numericInput(
          inputId = ns("max_cutoff"), 
          label = HTML(paste0(
            '<strong>Max Event Frequency (%)</strong> ',
            '<a class="tooltip-txt" data-html="true" ',
            'data-tooltip-toggle="tooltip" data-placement="top" ',
            'title=\"Maximum number (expressed as % of total) of ',
            '\'occurrences\' a feature (e.g., a mutation) can have ',
            'to be included in the \'Feature Set\`. Features with a ',
            'higher percentage of events than the specified number ',
            'will be removed.\n\nNOTE: \'Max event frequency\' must ',
            'be <= 90.\">?</a>')),
          value = 60, 
          min = 0, 
          max = 100, 
          step = 1, 
          width = "100%"
        ),
        radioButtons(
          inputId = ns("method"), 
          label = strong(span(style="color:red;", "*"), 
                         "Scoring method:"), 
          choices = c("ks", "wilcox", "revealer"), 
          selected = "ks", inline = TRUE
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'ks'", ns("method")),
          checkboxInput(
            inputId = ns("weighted_ks"), 
            label = "Compute weighted KS?", 
            value = FALSE
          ), 
          conditionalPanel(
            condition = sprintf("input['%s'] == true", ns("weighted_ks")),
            fileInput(
              inputId = ns("weights_file"), 
              label = strong(span(style = "color: red;", "*"), 
                             "Choose a weight file:"), 
              width = "100%"
            ),
            radioButtons(
              inputId = ns("weights_file_type"), 
              label = HTML(paste0(
                'File type ', 
                '<a class="tooltip-txt" data-html="true" ',
                'data-tooltip-toggle="tooltip" ',
                'data-placement="top" title=\"NOTE: ',
                'If file is in csv format, then ',
                'the \'Weights\' file must be a data frame ',
                'with two columns (Samples and Weights) and ',
                'the \'Samples\' column must match the colnames of ',
                '\'Feature Set\'. Otherwise, \'Weights\' file',
                'must contain a list of vectors and have names or ',
                'labels that match the colnames of \'Feature Set\'.\">?</a>')), 
              choices=c(".csv", ".rds"), 
              selected = ".csv", 
              inline = TRUE
            )     
          )
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'ks' | 
                              input['%s'] == 'wilcox'", 
                              ns("method"), 
                              ns("method")),
          selectInput(
            inputId = ns("alternative"), 
            label = strong(span(style="color:red;", "*"), 
                           "Alternative:"), 
            choices = c("less", "two.sided", "greater"), 
            selected = "less", width = "100%"
          ),
        ),
        fluidRow(
          column(
            width = 6,
            radioButtons(
              inputId = ns("metric"), 
              label = strong(span(style="color:red;", "*"), 
                             "Type of metric:"), 
              choices=c("pval", "stat"), 
              selected = "pval", inline = FALSE
            )
          ),
          column(
            width = 6,
            radioButtons(
              inputId = ns("search_method"), 
              label = strong(span(style="color:red;", "*"), 
                             "Search method:"), 
              choices=c("forward and backward"="both", 
                        "forward"="forward"), 
              selected = "both", inline = FALSE
            )
          )
        ),
        numericInput(
          inputId = ns("max_size"), 
          label = HTML(paste0(
            '<span style=\"color:red;\">*</span> Max meta-feature size ', 
            '<a class="tooltip-txt" data-html="true" ',
            'data-tooltip-toggle="tooltip" data-placement="top" ',
            'title=\"Max possible number of features to be ',
            'included in the meta-feature (search will stop ',
            'after max is reached)\">?</a>')), 
          min = 1, 
          max = 100, 
          step = 1, 
          value = 7, 
          width = "100%"
        ),
        radioButtons(
          inputId = ns("initial_seed"), 
          label = HTML(paste0(
            '<span style=\"color:red;\">*</span> Search modality ', 
            '<a class="tooltip-txt" data-html="true" ',
            'data-tooltip-toggle="tooltip" data-placement="top" ',
            'title=\"\'Top N\' repeats the search starting from ',
            'each of the top N scoring features. \'Custom seeds\' ',
            'repeats the search starting from each of the custom ',
            'seeds. WARNING: If number of seeds specified is greater ',
            'than 10, this may result in a longer search time.\">?</a>')), 
          choices = c("Top N seeds"="top_N_seeds", "Custom seeds"="search_start_seeds"), 
          selected = "top_N_seeds", 
          inline = TRUE
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'top_N_seeds'", 
                              ns("initial_seed")),
          numericInput(
            inputId = ns("top_N"), 
            label = strong(span(style = "color:red;", "*"), 
                           paste0("Top N value")), 
            min = 1, 
            max = 100, 
            step = 1, 
            value = 10, 
            width = "100%"
          ),
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'search_start_seeds'", 
                              ns("initial_seed")),
          textAreaInput(
            inputId = ns("search_start"), 
            label = strong(span(style = "color:red;", "*"), 
                           paste0('Enter a list of character strings ',
                                  '(separated by commas) corresponding ', 
                                  'to feature names within the ',
                                  '\'Feature Set\' object')), 
            value="", 
            width="100%"
          )
        ),
        checkboxInput(
          inputId = ns("permutation_test"), 
          label = strong("Perform permutation testing?"), 
          value = FALSE
        ), 
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("permutation_test")),
          numericInput(
            inputId = ns("n_perm"), 
            label = strong(span(style="color:red;", "*"), 
                           paste0("Number of permutations to perform")), 
            min = 1, 
            max = Inf,
            step = 1, 
            value = 100,
            width = "100%"
          ),
          numericInput(
            inputId = ns("ncores"), 
            label = strong(span(style="color:red;", "*"), 
                           paste0("Number of cores to perform parallelization for permutation testing")), 
            min = 1, 
            max = Inf, 
            step = 1, 
            value = 1,
            width = "100%"
          )
        ),
        br(),
        uiOutput(outputId = ns("error_message")),
        actionButton(
          inputId = ns("run_cadra"), 
          label = strong("RUN"), 
          style="background: blue; color: white;"
        ),
        actionButton(
          inputId = ns("stop_cadra"), 
          label = strong("STOP"), 
          style="background: blue; color: white;"
        ),
        br(), br(), br(), br(),
        HTML(
          paste0(
            "<p style='text-align: center;'>",
            "<span class='footer-info'>&copy; Monti Lab &diams; ",
            "<script>document.write(new Date().getFullYear());",
            "</script> &diams; All Rights Reserved.</span></p>"
          )
        )
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
            id = ns("loading_icon"), class = "loading_div", style="display: none;",
            span(
              div(class = "loader"),
              br(),
              p(class = "loading_text", "Running Candidate Search...")
            )
          ),
          div(
            uiOutput(outputId = ns("featureData_title"))
          ),
          div(
            uiOutput(outputId = ns("inputScoreData_title")),
            DT::dataTableOutput(outputId = ns("inputScoreData"))
          ),
          div(
            uiOutput(outputId = ns("bestFeatureData_title")),
            DT::dataTableOutput(outputId = ns("bestFeatureData"))
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
            id = ns("permutation_loading_icon"), class = "loading_div", style="display: none;",
            span(
              div(class = "loader"),
              br(),
              p(class = "loading_text", "Running Permutation Testing...")
            )
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
          
          htmltools::includeMarkdown(file.path(system.file('shinyapp', package="CaDrA"), "README.md"))
          
        ),
        tabPanel(
          title = "Publication", 
          style = "padding: 5px 10px 10px 10px;",
          icon = icon(name = "book", lib = "font-awesome"),
          
          h2("Citation"),
          p("Kartha VK, Kern JG, Sebastiani P, Zhang L, Varelas X, Monti S (2019) CaDrA: A computational framework for performing candidate driver analyses using genomic features.", a(href="https://www.frontiersin.org/articles/10.3389/fgene.2019.00121/full", "{Frontiers in Genetics}")),
          
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
#' @examples
#' 
#' # Load R library
#' library(shiny)
#' 
#' # Create ui and server for Shiny app
#' id <- "myapp"
#' 
#' ui <- fluidPage(
#'    CaDrA::CaDrA_UI(id = id)
#' )
#' 
#' server <- function(input, output, session) {
#'    CaDrA::CaDrA_Server(id = id)
#' }
#' 
#' app <-  shinyApp(ui=ui, server=server)
#' 
#' # Launch and deploy Shiny app (NOT RUN)
#' # shiny::runApp(app, host='0.0.0.0', port=3838)
#'  
#' @import shiny htmltools Biobase parallel
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom dplyr mutate_all
#' @importFrom stats rnorm 
#' @importFrom utils read.csv write.csv data
#' @importFrom methods new
#' @importFrom tools pskill
#' 
#' @export 
CaDrA_Server <- function(id){
  
  moduleServer(
    id,
    function(input, output, session) {

      # Combine extdata with global expression set and scores dataset if it was provided
      extdata <- get_extdata(dataset_choices, score_choices)
      eset_choices <- extdata[["eset_choices"]] %>% unlist()
      input_score_choices <- extdata[["input_score_choices"]] %>% unlist()
      
      # Detect number of cores on machine
      num_of_cores <- detectCores()

      # Create reactive values
      rVal <- reactiveValues()
      rVal$candidate_search_process <- NULL
      rVal$candidate_search_obs <- NULL
      rVal$cadra_permutation_process <- NULL
      rVal$cadra_permutation_obs <- NULL      
      feature_set_description <- reactiveVal()
      feature_set_data <- reactiveVal()
      input_score_data <- reactiveVal()
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
              Select the 'CaDrA Options' on the left and Click 'RUN' at the bottom
              "
            )
          )
        )
      })
      
      observeEvent(input$dataset, {
        
        selected_dataset <- isolate({ input$dataset }) 
        
        if(selected_dataset != "Import Data"){
          
          selection <- input_score_choices[which(eset_choices == selected_dataset)] %>% unlist()
          
          if(!is.na(names(selection))){
            updateSelectInput(session, inputId = "scores", choices = c(selection, "Import Data"), selected = selection[1])
          }else{
            updateSelectInput(session, inputId = "scores", choices = "Import Data")
          }
          
          if(tools::file_ext(selected_dataset) == "rda" | tools::file_ext(selected_dataset) == "RData"){
            selected_dataset <- load(selected_dataset)
          }
          
          if(selected_dataset == "BRCA_GISTIC_MUT_SIG"){
            updateNumericInput(session, inputId = "min_cutoff", value = 30)
          }else {
            updateNumericInput(session, inputId = "min_cutoff", value = 5)
          }
          
        }else{
          
          updateSelectInput(session, inputId = "scores", choices = "Import Data")
          
        }
        
      })
      
      #
      # Start the process
      #    
      observeEvent(input$run_cadra, {
        
        ns <- session$ns
        
        if(!is.null(rVal$candidate_search_process))
          return()
        
        rVal$candidate_search_result <- NULL
        rVal$cadra_permutation_result <- NULL
        feature_set_description(NULL)
        feature_set_data(NULL)
        input_score_data(NULL)
        instructions_message(FALSE)
        error_message(NULL)
        
        ## Show loading icon ####
        session$sendCustomMessage(type = "ToggleOperation", message = list(id=ns("loading_icon"), display="yes"))
        
        dataset <- isolate({ input$dataset })
        scores <- isolate({ input$scores })
        
        if(dataset == "Import Data"){
          
          inputfile <- input$ES_file;
          inputtype <- input$ES_file_type;
          
          if(is.null(inputfile)){
            error_message("Please choose a 'Feature Set' file to import.")
            return(NULL)
          }
          
          csv_ext <-  grep(toupper(".csv"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
          rds_ext <-  grep(toupper(".rds"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)

          if(inputtype %in% ".csv" & length(csv_ext) > 0){
            
            # read in the Eset file
            Eset <- read.csv(inputfile$datapath, header=TRUE, check.names=FALSE) 
            
            if("Features" %in% colnames(Eset)){
              
              Eset <- Eset %>% tibble::column_to_rownames(var="Features") %>% 
                dplyr::mutate_all(as.numeric)
              
              # convert Eset to matrix
              Eset <- as.matrix(Eset, nrow=nrow(Eset), ncol=ncol(Eset), 
                                byrow=TRUE, 
                                dimnames=list(rownames(Eset), colnames(Eset)))
              
              #create phenotypic data
              pData <- data.frame(Samples = colnames(Eset), stringsAsFactors=TRUE)
              rownames(pData) <- pData$Samples
              phenoData <- methods::new("AnnotatedDataFrame", data=pData)
              
              #create feature data
              fData <- data.frame(Features = rownames(Eset), stringsAsFactors=TRUE)
              rownames(fData) <- fData$Features
              featureData <-  methods::new("AnnotatedDataFrame", data=fData)
              
              #create expression set
              ES <- Biobase::ExpressionSet(
                assayData=Eset, 
                phenoData=phenoData, 
                featureData=featureData
              )
              
            }else{
              
              error_message("The 'Feature Set' file must contain a 'Features' column name that contains unique names or labels to search for best features.")
              return(NULL)
              
            }
            
          }else if (inputtype %in% ".rds" & length(rds_ext) > 0){
            
            ES <- readRDS(inputfile$datapath)
            
          }else{
            
            error_message("Incorrect file format. Please check your 'Feature Set' file again.")
            return(NULL)
            
          }
          
        }else{
          
            if(tools::file_ext(dataset) == "rda" | tools::file_ext(dataset) == "RData"){
              envir_name <- load(dataset)
              ES <- get(envir_name)
            }else{
              ES <- readRDS(dataset)
            }
          
        }
        
        if(scores == "Import Data"){
          
          inputfile <- input$input_score_file;
          inputtype <- input$input_score_file_type;
          
          if(is.null(inputfile)){
            error_message("Please choose a 'Input Score' file to import.")
            return(NULL)
          }
          
          csv_ext <-  grep(toupper(".csv"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
          rds_ext <-  grep(toupper(".rds"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
          
          if(inputtype %in% ".csv" & length(csv_ext) > 0) {
            
            dat <- read.csv(inputfile$datapath, header = TRUE, check.names = FALSE)
            
            if(all(c("Samples", "Scores") %in% colnames(dat))){
              input_score <- as.numeric(dat$Scores)
              names(input_score) <- as.character(dat$Samples)
            }else{
              error_message("The 'Input Score' file must be a data frame with two columns: Samples and Scores.")
              return(NULL)
            }
            
          } else if (inputtype %in% ".rds" & length(rds_ext) > 0){
            
            input_score <- readRDS(inputfile$datapath) 
            
          } else {
            
            error_message("Incorrect file format. Please check your 'Input Score' file again.")
            return(NULL)
            
          }
          
        }else{

          if(tools::file_ext(scores) == "rda" | tools::file_ext(scores) == "RData"){
            envir_name <- load(scores)
            input_score <- get(envir_name)
          }else{
            input_score <- readRDS(scores)
          }
          
        }
        
        # Obtain the pre-filter parameters (min_cutoff)
        min_cutoff <- as.integer(input$min_cutoff)
        
        #print(sprintf("minimum cutoff: %s", min_cutoff))
        
        if(is.na(min_cutoff) || length(min_cutoff)==0 || min_cutoff < 5){
          
          error_message("Please specify an integer value for Min Event Frequency >= 5 \n")
          return(NULL)
          
        } else {
          
          if(min_cutoff > ncol(ES)){
            error_message(sprintf("There are not enough samples  in \'Feature Set\' to meet Min Event Frequency = %s \n", min_cutoff))
            return(NULL)
          }

          percent_min_cutoff <- round(min_cutoff/ncol(ES), 2)
          
        }
        
        # Obtain the pre-filter parameters (max_cutoff)
        max_cutoff <- as.numeric(input$max_cutoff)
        
        #print(sprintf("maximum cutoff: %s", max_cutoff))
        
        if(is.na(max_cutoff) || length(max_cutoff)==0 || 
           max_cutoff <= 0 || max_cutoff > 100){
          error_message("Please specify a value for Max Event Frequency between 1 and 100\n")
          return(NULL)
        } else {
          max_cutoff <- max_cutoff/100
        }
        
        ## keep a record of the number of features in original ES
        n_orig_features <- nrow(ES)
        
        ## Pre-filter ESet based on occurrence frequency
        ES <- CaDrA::prefilter_data(
          ES = ES,
          max.cutoff = max_cutoff,
          min.cutoff = percent_min_cutoff 
        ) 
        
        # Make sure matrix is not empty after removing uninformative features
        if(nrow(exprs(ES)) == 0){
          error_message("Feature filtering based on specified 'Min Event Frequency' and 'Max Event Frequency' yielded an empty \'Feature Set\'.\n")
          return(NULL)
        }
        # Check if the ES is provided and is a BioBase ExpressionSet object
        if(length(ES) == 0 || !is(ES, "ExpressionSet")){
          error_message("'ES' must be a Biobase::ExpressionSet.")
          return(NULL)
        }
        # Check if the dataset has only binary 0 or 1 values 
        if(!all(exprs(ES) %in% c(0,1))){
          error_message("The \'Feature Set\' (ES) must contain only binary values (0/1) with no NAs.\n")
          return(NULL)
        }
        # Make sure the input ES has rownames for features tracking
        if(is.null(rownames(ES))){
          error_message("The ES object does not have rownames or featureData to track the features by. Please provide unique features or rownames for the \'Feature Set\'.\n")
          return(NULL)
        }
        # Check input_score is provided and are continuous values with no NAs
        if(length(input_score) == 0 || any(!is.numeric(input_score)) || any(is.na(input_score))){
          error_message("input_score must be a vector of continous values (with no NAs) where the vector names match the colnames of the \'Feature Set\'.\n")
          return(NULL)
        }
        # Make sure the input_score has names or labels that are the same as the colnames of ES
        if(is.null(names(input_score))){
          error_message("The input_score object must have names or labels to track the samples by. Please provide unique sample names or labels that matches the colnames of the \'Feature Set\'.\n")
          return(NULL)
        }
        # Make sure the input_score has the same length as number of samples in ES
        if (length(input_score) != ncol(ES)) {
          
          error_message("The input_score must have the same length as the number of columns in ES.\n")
          return(NULL)
          
        } else {
          
          # pos <- which(!names(input_score) %in% colnames(ES))
          # print(names(input_score)[pos])
          # print(colnames(ES)[pos])
          
          if(any(!names(input_score) %in% colnames(ES))){
            error_message("The input_score object must have names or labels that match the colnames of the \'Feature Set\'.\n")
            return(NULL)
          }
          # match colnames of Feature Set with names of provided input_score
          ES <- ES[,names(input_score)]
        }
        # Get method
        method <- input$method; 
        
        if(method == "ks"){
          
          if(input$weighted_ks == TRUE){
            
            inputfile <- input$weights_file;
            inputtype <- input$weights_file_type;
            
            if(is.null(inputfile)){
              error_message("Please choose a 'Weighted KS' file to import.")
              return(NULL)
            }
            
            csv_ext <-  grep(toupper(".csv"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
            rds_ext <-  grep(toupper(".rds"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
            
            if(inputtype %in% ".csv" & length(csv_ext) > 0){
              
              dat <- read.csv(inputfile$datapath, header=TRUE, check.names=FALSE)
              
              if(all(c("Samples", "Weights") %in% colnames(dat))){
                weights <- as.numeric(dat$Weights)
                names(weights) <- as.character(dat$Samples)
              }else{
                error_message("The 'Weighted KS' file must be a data.frame with two columns: Samples and Weights.")
                return(NULL)
              }
              
            }else if(inputtype %in% ".rds" & length(rds_ext) > 0){
              
              weights <- readRDS(inputfile$datapath)
              
            }else{
              
              error_message("Incorrect file format. Please check your 'Weighted KS' file again.")
              return(NULL)
              
            }
            
            # Check weights is provided and are continuous values with no NAs
            if(length(weights) == 0 || any(!is.numeric(weights)) || any(is.na(weights))){
              error_message("weights must be a vector of continous values (with no NAs) with the vector names matching the colnames of the \'Feature Set\'.\n")
              return(NULL)
            }
            
            # Make sure the weights has names or labels that are the same as the colnames of ES
            if(is.null(names(weights))){
              error_message("The weights object must have names or labels to track the samples by. Please provide unique sample names or labels that match the colnames of the \'Feature Set\'.\n")
              return(NULL)
            }

            # Make sure the weights has the same length as number of samples in ES
            if (length(weights) != ncol(ES)) {

              error_message("The weights must have the same length as the number of columns in the \'Feature Set\'.\n")
              return(NULL)
              
            } else {
              
              if (any(!names(weights) %in% colnames(ES))) {
                error_message("The weights object must have names or labels that match the colnames of the \'Feature Set\'.\n")
                return(NULL)
              }
              # match colnames of Feature Set with names of provided weights
              weights <- weights[colnames(ES)]
            }
          } else {
            
            weights <- NULL
            
          }
        }
        if(method %in% c("ks", "wilcox")){
          alternative <- input$alternative
        }
        metric <- input$metric;
        
        search_method <- input$search_method;
        
        max_size <- as.integer(input$max_size)
        
        if(is.na(max_size) || length(max_size) == 0 || max_size <= 0){
          error_message("Please specify an integer value specifies a maximum size that a meta-feature can extend to do for a given search (max_size must be >= 1).\n")
          return(NULL)
        }
        
        if(max_size > nrow(ES)){
          error_message("Please specify a \'Max meta-feature size\' lesser than the number of features in the Feature Set.\n")
          return(NULL)
        }
        
        initial_seed <- input$initial_seed
        
        if(initial_seed == "top_N_seeds"){
          
          search_start <- NULL
          top_N <- as.integer(input$top_N)
          
          if(is.na(top_N) || length(top_N) == 0 || top_N <= 0){
            error_message("Please specify an INTEGER Top N value to evaluate over Top N features (top_N must be >= 1).\n")
            return(NULL)
          }
          
          if(top_N > nrow(ES)){
            error_message("Please specify a Top N value lesser than the number of features in the Feature Set.\n")
            return(NULL)
          }
          
        } else {
          
          search_start <- strsplit(as.character(input$search_start), ",", fixed=TRUE) %>% unlist() %>% trimws();
          search_start <- search_start[search_start != ""]
          top_N <- NULL
          
          if(length(search_start) == 0 || any(!search_start %in% rownames(ES))){ 
            error_message("The provided \'Custom Seeds\' do not exist among Feature Set's rownames.\n\n")
            return(NULL)
          }
        }
        
        permutation <- input$permutation_test
        
        if(permutation == TRUE){
          
          ## show permutation loading icon ####
          session$sendCustomMessage(type = "ToggleOperation", message = list(id=ns("permutation_loading_icon"), display="yes"))

          n_perm <- as.integer(input$n_perm)
          
          if(is.na(n_perm) || length(n_perm)==0 || n_perm <= 0){
            error_message("Please specify an INTEGER number of permutations (n_perm must be >= 1).\n")
            return(NULL)
          }
          
          ncores <- as.integer(input$ncores)
          
          if(is.na(ncores) || length(ncores)==0 || ncores <= 0){
            error_message("Please specify the number of parallelization cores for permutation testing (ncores must be >= 1).\n")
            return(NULL)
          }

          if(ncores > num_of_cores){
            error_message(paste0("There are ONLY ", num_of_cores, " cores available on the system. Please specify the number of parallelization cores for permutation testing (ncores <= ", num_of_cores, ")."))
            return(NULL)
          }
          
          rVal$cadra_permutation_process <- parallel::mcparallel({
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
              smooth = TRUE,
              obs_best_score = NULL,
              plot = FALSE, 
              ncores = ncores,
              cache_path = NULL,
              verbose = FALSE
            )
          })
          
          #print(paste0("cadra permutation process: ", rVal$cadra_permutation_process$pid, " started"))
          
        }
        
        # Export the ES and input_score to reactive datase
        feature_set_description(
          sprintf(
            paste0(
              'After filtering features with Min Event Frequency = %s ',
              '(or having < %s%% prevalence across all samples) and Max ',
              'Event Frequency = %s (or having > %s%% prevalance across ',
              'all samples), the \'Feature Set\' retained %s genomic ',
              'features out of %s supplied features across %s samples.'), 
            min_cutoff, 
            percent_min_cutoff*100, 
            max_cutoff*100, 
            max_cutoff*100, 
            format(nrow(ES), big.mark = ","), 
            format(n_orig_features, big.mark = ","), 
            format(ncol(ES), big.mark =",")
          )
        )
        
        feature_set_data(exprs(ES))
        input_score_data(input_score)
        
        # Compute the candidate search algorithm
        rVal$candidate_search_process <- parallel::mcparallel({
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
        })
        
        #print(paste0("candidate search process: ", rVal$candidate_search_process$pid, " started"))
        
        error_message("NONE")
        
      })
      #
      # Stop the process
      #
      observeEvent(input$stop_cadra, {
        
        if(!is.null(rVal$candidate_search_process)) {
          #print(paste0("candidate search process: ", rVal$candidate_search_process$pid, " killed"))
          tools::pskill(rVal$candidate_search_process$pid)
          rVal$candidate_search_process <- NULL

          if (!is.null(rVal$candidate_search_process)) {
            rVal$candidate_search_obs$destroy()
          }
        }
        
        if(!is.null(rVal$cadra_permutation_process)) {
          #print(paste0("cadra permutation process: ", rVal$cadra_permutation_process$pid, " killed"))
          tools::pskill(rVal$cadra_permutation_process$pid)
          rVal$cadra_permutation_process <- NULL
          
          if (!is.null(rVal$cadra_permutation_process)) {
            rVal$cadra_permutation_obs$destroy()
          }
        }
        
        ns <- session$ns
        
        rVal$candidate_search_result <- NULL; rVal$cadra_permutation_result <- NULL;
        
        error_message("Your process has been interrupted")
        
        ## Hide loading icon ####
        session$sendCustomMessage(type = "ToggleOperation", message = list(id=ns("loading_icon"), display="no"))

        ## Hide permutation loading icon ####
        session$sendCustomMessage(type = "ToggleOperation", message = list(id=ns("permutation_loading_icon"), display="no"))
        
        # Show instruction message
        instructions_message(TRUE)
        
      }, ignoreInit = TRUE) 
      #
      # Handle candidate search process event
      #
      observeEvent(rVal$candidate_search_process, {
        rVal$candidate_search_obs <- observe({
          shiny::invalidateLater(500, session)
          isolate({
            ns <- session$ns
            result <- parallel::mccollect(rVal$candidate_search_process, wait = FALSE)
            if(!is.null(result)) {
              rVal$candidate_search_result <- result[[1]]
              rVal$candidate_search_obs$destroy()
              rVal$candidate_search_process <- NULL
              ## Hide loading icon ####
              session$sendCustomMessage(type = "ToggleOperation", message = list(id=ns("loading_icon"), display="no"))
            }
          })
        })
      }, ignoreInit = TRUE)
      #
      # Handle cadra permutation process event
      #      
      observeEvent(rVal$cadra_permutation_process, {
        rVal$cadra_permutation_obs <- observe({
          shiny::invalidateLater(500, session)
          isolate({
            ns <- session$ns
            result <- parallel::mccollect(rVal$cadra_permutation_process, wait = FALSE)
            if(!is.null(result)) {
              rVal$cadra_permutation_result <- result[[1]]
              rVal$cadra_permutation_obs$destroy()
              rVal$cadra_permutation_process <- NULL
              ## Hide permutation loading icon ####
              session$sendCustomMessage(type = "ToggleOperation", message = list(id=ns("permutation_loading_icon"), display="no"))
            }
          })
        })
      })
      #
      # Render messages
      #
      output$error_message <- renderUI({

        req(error_message())

        ns <- session$ns

        if(error_message() != "NONE"){
          
          ## Update loading icon ####
          session$sendCustomMessage(type = "ToggleOperation", message = list(id=ns("loading_icon"), display="no"))

          ## Update loading icon for permutation test
          session$sendCustomMessage(type = "ToggleOperation", message = list(id=ns("permutation_loading_icon"), display="no"))
          
          # Show instructions
          instructions_message(TRUE)
          
          # Display error message
          p(style="color: red; font-weight: bold; margin-bottom: 10px;", error_message())
          
        }
        
      })
      output$featureData_title <- renderUI({

        req(rVal$candidate_search_result, feature_set_description())

        ns <- session$ns

        description <- feature_set_description()

        dataset <- isolate({ input$dataset })

        if (dataset == "Import Data"){
          title <- "Dataset: Imported Data"
        }else{
          title <- paste0("Dataset: ", names(eset_choices[which(eset_choices == dataset)]))
        }
        div(
          h2(title),
          br(),
          p(description),
          downloadButton(outputId = ns("download_featureset"), label="Download Filtered Feature Set")
        )
      })
      output$download_featureset <- downloadHandler(
        filename = function() {
          paste0("CaDrA-Filtered-Features-Eset.csv")
        },
        content = function(file) {
          Eset_table <- feature_set_data() %>% as.data.frame(.) %>% rownames_to_column(var="Features")
          write.csv(Eset_table, file, row.names=FALSE)
        }
      )
      output$bestFeatureData_title <- renderUI({
        req(rVal$candidate_search_result)
        h2("Best Meta-Feature Eset")
      })
      output$bestFeatureData <- DT::renderDataTable({

        req(rVal$candidate_search_result)

        ns <- session$ns

        topn_best_meta <- topn_best(topn_list=rVal$candidate_search_result)

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
                  text = 'Download Results',
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
            title = "Download Best Meta-Feature Eset",
            downloadButton(outputId = ns("downloadEsetCSV"), 
                           "Download Table as CSV file"),
            br(), br(),
            downloadButton(outputId = ns("downloadEsetRDS"), 
                           "Download Table as RDS file"),
          )
        )
      })
      ## Download CSV File ####
      output$downloadEsetCSV <- downloadHandler(

        filename = function() {
          paste0("CaDrA-Best-Meta-Feature-Eset.csv")
        },

        content = function(file) {
          topn_best_meta <- topn_best(topn_list=rVal$candidate_search_result)
          Eset_table <- topn_best_meta[["ESet"]] %>% exprs(.) %>% as.data.frame(.) %>% tibble::rownames_to_column(., var="Features")
          write.csv(Eset_table, file, row.names=FALSE)

        }
      )
      ## Download RDS File ####
      output$downloadEsetRDS <- downloadHandler(

        filename = function() {
          paste0("CaDrA-Best-Meta-Feature-Eset.rds")
        },
        content = function(file) {

          topn_best_meta <- topn_best(topn_list=rVal$candidate_search_result)

          Eset_table <- topn_best_meta[["ESet"]]

          saveRDS(Eset_table, file)
        }
      )
      output$inputScoreData_title <- renderUI({

        req(rVal$candidate_search_result, input_score_data())

        selection <- isolate({ input$scores })

        if(selection == "Import Data"){
          title <- "Imported Data"
        }else{
          title <- names(input_score_choices[which(input_score_choices == selection)])
        }
        
        h2("Observed Input Scores:", title)
        
      })
      output$inputScoreData <- DT::renderDataTable({

        req(rVal$candidate_search_result, input_score_data())

        ns <- session$ns

        Scores <- input_score_data() %>% signif(., digits = 4)
        score_table <- matrix(Scores, nrow=1, ncol=length(Scores), byrow=TRUE, dimnames=list("input_score", names(Scores)))

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
                  text = 'Download Results',
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
            downloadButton(outputId = ns("downloadScoreCSV"), 
                           "Download Table as CSV file"),
            br(), br(),
            downloadButton(outputId = ns("downloadScoreRDS"), 
                           "Download Table as RDS file"),
          )
        )
      })
      ## Download CSV File ####
      output$downloadScoreCSV <- downloadHandler(

        filename = function() {
          paste0("CaDrA-Observed-Input-Scores.csv")
        },

        content = function(file) {

          Scores <- input_score_data() %>% signif(., digits = 4)

          score_table <- data.frame(
            Samples = names(Scores),
            Scores = Scores,
            stringsAsFactors = FALSE
          )
          write.csv(score_table, file, row.names=FALSE)
        }
      )
      ## Download RDS File ####
      output$downloadScoreRDS <- downloadHandler(

        filename = function() {
          paste0("CaDrA-Observed-Input-Scores.rds")
        },
        content = function(file) {

          Scores <- input_score_data() %>% signif(., digits = 4)

          score_table <- data.frame(
            Samples = names(Scores),
            Scores = Scores,
            stringsAsFactors = FALSE
          )
          saveRDS(score_table, file)
        }
      )
      output$meta_plot_title <- renderUI({

        req(rVal$candidate_search_result)

        h2("Best Meta-Feature Plot")
      })
      output$meta_plot <- renderPlot({

        req(rVal$candidate_search_result)

        topn_res <- rVal$candidate_search_result
        topn_best_meta <- CaDrA::topn_best(topn_res)

        CaDrA::meta_plot(topn_best_list = topn_best_meta)
      })
      output$topn_plot_title <- renderUI({

        req(rVal$candidate_search_result)

        if(length(rVal$candidate_search_result) == 1){

          div(
            h2("Top N Overlapping Heatmap"),
            br(),
            h4(style="color: red; font-weight: bold;", 
               "NOTE: Cannot plot overlap matrix with provided top N seed = 1 or the number of provided feature names = 1.")
          )
        } else {
          h2("Top N Overlapping Heatmap")
        }
      })
      output$topn_plot <- renderPlot({

        req(rVal$candidate_search_result)

        topn_res <- rVal$candidate_search_result

        CaDrA::topn_plot(topn_res)
      })
      output$permutation_plot_title <- renderUI({

        req(rVal$cadra_permutation_result)

        h2("Permutation-Based Testing")
      })
      output$permutation_plot <- renderPlot({

        req(rVal$cadra_permutation_result)

        perm_res <- rVal$cadra_permutation_result
        permutation_plot(perm_res)
      })
    }
  )
}
#' Run both Shiny UI and Server Modules 
#' 
#' @return Shiny application 
#'
#' @examples
#' 
#' # Load R library
#' library(shiny)
#'
#' # Create a Shiny app
#' app <- CaDrA::CaDrA_App()
#' 
#' # Launch and deploy Shiny app (NOT RUN)
#' # shiny::runApp(app, host='0.0.0.0', port=3838)
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
    
    ## Reconnecting to new sessions after grey out
    session$allowReconnect("force")
    
    ## Close app when browser closing
    session$onSessionEnded(function() {
     #print("Closed Shiny Session.")
     stopApp()
    })
    
  }
  shinyApp(ui=ui, server=server)
}







