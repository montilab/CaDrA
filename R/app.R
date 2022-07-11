
# function to hover columns in data table
create_hover_txt <- function(table){
  
  column_names <- colnames(table)
  
  th_tr <- lapply(seq_along(column_names), function(l){ 
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

dataset_choices <- list(
  "TTGA BRCA Mutation Signatures (BRCA_GISTIC_MUT_SIG)" = "BRCA_GISTIC_MUT_SIG",
  "CCLE SCNA Mutations in Cancers (CCLE_MUT_SCNA)" = "CCLE_MUT_SCNA", 
  "Simulated Data (sim.ES)" = "sim.ES"
)

score_choices <- list(
  "YAP/TAZ Activity in Human Breast Cancer (TAZYAP_BRCA_ACTIVITY)" = "TAZYAP_BRCA_ACTIVITY",
  "Activation of B-catenin in Cancers (CTNBB1_reporter)" = "CTNBB1_reporter",
  "Random simulated input scores (sim.Scores)" = "sim.Scores"
)

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
#' 
#' server <- function(input, output, session) {
#'    CaDrA::CaDrA_Server(id = id)
#' }
#' 
#' app <-  shinyApp(ui=ui, server=server)
#' 
#' # Launch Shiny app (NOT RUN)
#' # shiny::runApp(app, host='0.0.0.0', port=3838)
#' 
#' @import shiny markdown
#' @importFrom htmltools includeMarkdown
#' 
#' @export 
CaDrA_UI <- function(id){
  
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
          width: 120px;
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
          "Shiny.addCustomMessageHandler('ToggleOperation', function(message) {",
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
          label = "Choose a BINARY feature set:", 
          choices = c(dataset_choices, "Import Data"),
          selected = dataset_choices[1], 
          width = "100%"
        ),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == '%s'", ns("dataset"), dataset_choices[1]),
          selectInput(
            inputId = ns(paste0(dataset_choices[1], "_scores")), 
            label = "Choose an input_score:", 
            choices = c(score_choices[1], "Import Data"), 
            selected = score_choices[1], 
            width = "100%"
          )
        ),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == '%s'", ns("dataset"), dataset_choices[2]),
          selectInput(
            inputId = ns(paste0(dataset_choices[2], "_scores")), 
            label = "Choose an input_score:", 
            choices = c(score_choices[2], "Import Data"), 
            selected = score_choices[2], 
            width = "100%"
          )
        ),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == '%s'", ns("dataset"), dataset_choices[3]),
          selectInput(
            inputId = ns(paste0(dataset_choices[3], "_scores")), 
            label = "Choose an input_score:", 
            choices = c(score_choices[3], "Import Data"), 
            selected = score_choices[3], 
            width = "100%"
          )
        ),
          
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Import Data'", ns("dataset")),
          fileInput(
            inputId = ns("ES_file"), 
            label = strong(span(style = "color: red;", "*"), "Choose a binary feature file:"), 
            width = "100%"
          ),
          radioButtons(
            inputId = ns("ES_file_type"), 
            label = HTML("File type", paste0('<a class="tooltip-txt" data-html="true" data-tooltip-toggle="tooltip" data-placement="top" title=\"NOTE: If file is in csv format, the \'BINARY feature set\' must be a data frame including a \'Features\' column name that contains unique names or labels to search for best features. Otherwise, Feature Set must be an object of class ExpressionSet from BioBase package.\">?</a>')), 
            choices = c(".csv", ".rds"), 
            selected = ".csv", 
            inline = TRUE
          )
        ),

        conditionalPanel(
          condition = sprintf("input['%s'] == 'Import Data' | (input['%s'] == '%s' & input['%s'] == 'Import Data') | (input['%s'] == '%s' & input['%s'] == 'Import Data') | (input['%s'] == '%s' & input['%s'] == 'Import Data')", ns("dataset"), ns("dataset"), dataset_choices[1], ns(paste0(dataset_choices[1], "_scores")), ns("dataset"), dataset_choices[2], ns(paste0(dataset_choices[2], "_scores")), ns("dataset"), dataset_choices[3], ns(paste0(dataset_choices[3], "_scores"))),
          fileInput(
            inputId = ns("input_score_file"), 
            label = strong(span(style = "color: red;", "*"), "Choose an input score file:"), 
            width = "100%"
          ),
          radioButtons(
            inputId = ns("input_score_file_type"), 
            label = HTML("File type", paste0('<a class="tooltip-txt" data-html="true" data-tooltip-toggle="tooltip" data-placement="top" title=\"NOTE: The input score file must be a data frame with two columns (Samples and Scores) and the Samples must match the colnames of the binary feature set.\">?</a>')), 
            choices = c(".csv", ".rds"), 
            selected = ".csv", 
            inline = TRUE
          )
        ),
        
        uiOutput(outputId = ns("min_cutoff_tooltip")),
        
        numericInput(
          inputId = ns("min_cutoff"), 
          label = NULL,
          value = 30,
          min = 5, 
          max = Inf, 
          step = 1,
          width = "100%"
        ),
        
        uiOutput(outputId = ns("max_cutoff_tooltip")),
        
        numericInput(
          inputId = ns("max_cutoff"), 
          label = NULL,
          value = 60, 
          min = 0, 
          max = 100, 
          step = 1, 
          width = "100%"
        ),
        
        radioButtons(inputId = ns("method"), label = strong(span(style="color:red;", "*"), "Scoring method:"), choices = c("ks", "wilcox", "revealer"), selected = "ks", inline = TRUE),
        
        conditionalPanel(
          condition = sprintf("input['%s'] == 'ks'", ns("method")),
          checkboxInput(inputId = ns("weighted_ks"), label = "Compute the weighted ks?", value = FALSE), 
          conditionalPanel(
            condition = sprintf("input['%s'] == true", ns("weighted_ks")),
            fileInput(
              inputId = ns("weights_file"), 
              label = strong(span(style = "color: red;", "*"), "Choose a weight file:"), 
              width = "100%"
            ),
            radioButtons(
              inputId = ns("weights_file_type"), 
              label = HTML("File type", paste0('<a class="tooltip-txt" data-html="true" data-tooltip-toggle="tooltip" data-placement="top" title=\"NOTE: The weights file must be a data frame with two columns (Samples and Weights) and the Samples must match the colnames of the binary feature set.\">?</a>')), 
              choices=c(".csv", ".rds"), 
              selected = ".csv", 
              inline = TRUE
            )     
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
          textAreaInput(inputId = ns("search_start"), label = strong(span(style = "color:red;", "*"), "Enter a list of character strings (separated by commas) which specifies feature names within the expression set object to start the search with"), value="", width="100%")
        ),
        
        checkboxInput(inputId = ns("permutation_test"), label = strong("Perform permutation testing?"), value = FALSE), 
        
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("permutation_test")),
          numericInput(inputId = ns("n_perm"), label = strong(span(style="color:red;", "*"), "Number of permutations to perform"), min = 1, max = Inf, step = 1, value = 100),
          numericInput(inputId = ns("ncores"), label = strong(span(style="color:red;", "*"), "Number of cores to perform parallelization for permutation testing"), min = 1, max = Inf, step = 1, value = 1)
        ),
        
        br(),
        
        uiOutput(outputId = ns("error_message")),
        
        actionButton(inputId = ns("run_cadra"), label = strong("RUN"), style="background: blue; color: white;"),
        actionButton(inputId = ns("stop_cadra"), label = strong("STOP"), style="background: blue; color: white;"),
        
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
            id = ns("loading_icon"), class = "loading_div", style="display: none;",
            span(
              div(class = "loader"),
              br(),
              p(class = "loading_text", "Loading...")
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
          p("Kartha VK, Kern JG, Sebastiani P, Zhang L, Varelas X, Monti S (2019) CaDrA: A computational framework for performing candidate driver analyses using binary genomic features.", a(href="https://www.frontiersin.org/articles/10.3389/fgene.2019.00121/full", "{Frontiers in Genetics}")),
          
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
#' @import shiny htmltools Biobase
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom dplyr mutate_all
#' @importFrom stats rnorm 
#' @importFrom utils read.csv write.csv data
#' @importFrom methods new
#' 
#' @export 
CaDrA_Server <- function(id){
  
  moduleServer(
    id,
    function(input, output, session) {
      
      # Create reative values
      stop_process <- reactiveVal(FALSE)
      feature_set_description <- reactiveVal()
      feature_set_data <- reactiveVal()
      input_score_data <- reactiveVal()
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
      
      output$min_cutoff_tooltip <- renderUI({
        
        if(input$dataset == "BRCA_GISTIC_MUT_SIG"){
          HTML('<strong>Minimum Samples Frequency</strong>', '<a class="tooltip-txt" data-html="true" data-tooltip-toggle="tooltip" data-placement="top" title=\"The \'Minimum Samples Frequency\' means each feature in \'Feature Set\' must have at least 30 or more samples with presence of \'omic feature\' (i.e., any feature occurs in less than 30 samples will be automatically removed).\n\nNOTE: \'Minimum Samples Frequency\' must be >= 5.\">?</a>')
        }else{
          HTML('<strong>Minimum Samples Frequency</strong>', '<a class="tooltip-txt" data-html="true" data-tooltip-toggle="tooltip" data-placement="top" title=\"The \'Minimum Samples Frequency\' means each feature in \'Feature Set\' must have at least 5 or more samples with presence of \'omic feature\' (i.e., any feature occurs in less than 5 samples will be automatically removed).\n\nNOTE: \'Minimum Samples Frequency\' must be >= 5.\">?</a>')
        }
        
      })
      
      output$max_cutoff_tooltip <- renderUI({
        
        HTML(paste0('<strong>Maximum Percent Cutoff</strong> <a class="tooltip-txt" data-html="true" data-tooltip-toggle="tooltip" data-placement="top" title=\"The \'Maximum Percent Cutoff\' means each feature in \'Feature Set\' must have at least 60% or less of the samples with presence of \'omic feature\' (i.e., any feature occurs in > 60% of the samples will be removed).\">?</a>'))
        
      })
      
      observeEvent(input$dataset, {
        
        if(input$dataset == "BRCA_GISTIC_MUT_SIG"){
          updateNumericInput(session, inputId = "min_cutoff", value = 30)
        }else {
          updateNumericInput(session, inputId = "min_cutoff", value = 5)
        }
        
      })
      
      observeEvent(input$stop_cadra, {
        
        ns <- session$ns
        
        shiny::invalidateLater(300, session)
        stop_process(TRUE)
        error_message("Your process has been interrupted")
        
        ## Whether to display loading icon ####
        session$sendCustomMessage(type = "ToggleOperation", message = list(id=ns("loading_icon"), display="no"))
        
      })
        
      observeEvent(input$run_cadra, {
        
        ns <- session$ns
        
        shiny::invalidateLater(300, session)
        
        if(stop_process()){         
          return(NULL); 
        }
        
        stop_process(FALSE)
        feature_set_description(NULL)
        feature_set_data(NULL)
        input_score_data(NULL)
        candidate_search_result(NULL)
        permutation_result(NULL) 
        instructions_message(FALSE)
        error_message(NULL)
        
        ## Update connectivity option ####
        session$sendCustomMessage(type = "ToggleOperation", message = list(id=ns("loading_icon"), display="yes"))
        
        if(input$dataset == "BRCA_GISTIC_MUT_SIG"){
          
          ## Read in BRCA GISTIC+Mutation ESet object
          utils::data("BRCA_GISTIC_MUT_SIG", envir = environment())
          eset_mut_scna <- get("BRCA_GISTIC_MUT_SIG", envir = environment())
          
          if(input$BRCA_GISTIC_MUT_SIG_scores == "TAZYAP_BRCA_ACTIVITY"){
            
            ## Read in input score
            utils::data("TAZYAP_BRCA_ACTIVITY", envir = environment())
            input_score <- get("TAZYAP_BRCA_ACTIVITY", envir = environment())
              
            ## Samples to keep based on the overlap between the two inputs
            overlap <- intersect(names(input_score), Biobase::sampleNames(eset_mut_scna))
            eset_mut_scna <- eset_mut_scna[,overlap]
            input_score <- input_score[overlap]
            
          }
          
          ## Binarize ES to only have 0's and 1's
          exprs(eset_mut_scna)[exprs(eset_mut_scna) >= 1] <- 1.0
          
          ES <- eset_mut_scna
            
        }
        
        if(input$dataset == "CCLE_MUT_SCNA"){
          
          ## Read in SCNA ESet object
          utils::data("CCLE_MUT_SCNA", envir = environment())
          ES <- get("CCLE_MUT_SCNA", envir = environment())
            
          if(input$CCLE_MUT_SCNA_scores == "CTNBB1_reporter"){
            utils::data("CTNBB1_reporter", envir = environment())
            input_score <- get("CTNBB1_reporter", envir = environment())
          }
  
        }
        
        if(input$dataset == "sim.ES"){
          
          ## Read in simulated eset object
          utils::data("sim.ES", envir = environment())
          ES <- get("sim.ES", envir = environment())
            
          if(input$sim.ES_scores == "sim.Scores"){
            utils::data("sim.Scores", envir = environment())
            input_score = get("sim.Scores", envir = environment())
          }
          
        }
        
        if(input$dataset == "Import Data"){
          
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
            Eset <- read.csv(inputfile$datapath, header=TRUE) %>% tibble::column_to_rownames(var="Features") %>% dplyr::mutate_all(as.integer)
            
            # convert Eset to matrix
            Eset <- as.matrix(Eset, nrow=nrow(Eset), ncol=ncol(Eset), byrow=TRUE, dimnames=list(rownames(Eset), colnames(Eset)))
            
            #create phenotypic data
            pData <- data.frame(Samples = colnames(Eset), stringsAsFactors=TRUE)
            rownames(pData) <- pData$Samples
            phenoData <- methods::new("AnnotatedDataFrame", data=pData)
            
            #create feature data
            fData <- data.frame(Features = rownames(Eset), stringsAsFactors = TRUE)
            rownames(fData) <- fData$Features
            featureData <-  methods::new("AnnotatedDataFrame", data=fData)
            
            #create expression set
            ES <- Biobase::ExpressionSet(assayData=Eset, phenoData=phenoData, featureData=featureData)
            
          }else if(inputtype %in% ".rds" & length(rds_ext) > 0){
            
            ES <- readRDS(inputfile$datapath)
            
          }else{
            
            error_message("Incorrect file format. Please check your file again.")
            return(NULL)
            
          }
          
        }
        
        if(input$dataset == "Import Data" | (input$dataset == "BRCA_GISTIC_MUT_SIG" & input$BRCA_GISTIC_MUT_SIG_scores == "Import Data") | (input$dataset == "CCLE_MUT_SCNA" & input$CCLE_MUT_SCNA_scores == "Import Data") | (input$dataset == "sim.ES" & input$sim.ES_scores == "Import Data")){
          
          inputfile <- input$input_score_file;
          inputtype <- input$input_score_file_type;
          
          if(is.null(inputfile)){
            error_message("Please choose a file to import.")
            return(NULL)
          }
          
          csv_ext <-  grep(toupper(".csv"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
          rds_ext <-  grep(toupper(".rds"), toupper(substr(inputfile$datapath, nchar(inputfile$datapath)-4, nchar(inputfile$datapath))), fixed = TRUE)
          
          if(inputtype %in% ".csv" & length(csv_ext) > 0){
            
            dat <- read.csv(inputfile$datapath, header = TRUE) %>% tibble::column_to_rownames(var="Samples") %>% dplyr::mutate_all(as.numeric)
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
        
        # Obtain the pre-filter parameters (min_cutoff)
        min_cutoff <- as.integer(input$min_cutoff)
        
        #print(sprintf("minimum cutoff: %s", min_cutoff))
        
        if(is.na(min_cutoff) || length(min_cutoff)==0 || min_cutoff < 5){
          
          error_message("Please specify an integer value for Minimum Samples Cutoff >= 5\n")
          return(NULL)
          
        }else{
          
          if(ncol(ES) < 5){
            error_message(sprintf("The # of samples in Feature Set must be greater or equal to the Minimum Samples Cutoff = %s \n", min_cutoff))
            return(NULL)
          }
          
          percent_min_cutoff = round(min_cutoff/ncol(ES), 2)
          
        }
        
        # Obtain the pre-filter parameters (max_cutoff)
        max_cutoff <- as.numeric(input$max_cutoff)
        
        #print(sprintf("maximum cutoff: %s", max_cutoff))
        
        if(is.na(max_cutoff) || length(max_cutoff)==0 || max_cutoff < 0 || max_cutoff > 100){
          
          error_message("Please specify a value for Maximum Percent Cutoff between 0 to 100\n")
          return(NULL)
          
        }else{
          
          max_cutoff <- max_cutoff/100
          
        }
        
        #print(sprintf("number of samples: %s", ncol(ES)))
        #print(sprintf("percent minimum cutoff: %s", percent_min_cutoff))
        #print(sprintf("percent maximum cutoff: %s", max_cutoff))
        
        ## keep a record of the number of features in original ES
        n_orig_features <- nrow(ES)
        
        ## Pre-filter ESet based on occurrence frequency
        ES <- CaDrA::prefilter_data(
          ES = ES,
          max.cutoff = max_cutoff,
          min.cutoff = percent_min_cutoff 
        ) 
        
        #print(sprintf("number of samples after filtered: %s", ncol(ES)))
        
        # Make sure matrix is not empty after removing uninformative features
        if(nrow(exprs(ES)) == 0){
          error_message("After removing features based on the specified occurrence parameters. There are no more features remained for downsteam computation.\n")
          return(NULL)
        }
        
        # Check if the ES is provided and is a BioBase ExpressionSet object
        if(length(ES) == 0 || !is(ES, "ExpressionSet")){
          error_message("'ES' must be an ExpressionSet class argument (required).")
          return(NULL)
        }
        
        # Check if the dataset has only binary 0 or 1 values 
        if(!all(exprs(ES) %in% c(0,1))){
          error_message("The expression matrix (ES) must contain only binary values (0/1) with no NAs.\n")
          return(NULL)
        }
        
        # Make sure the input ES has rownames for features tracking
        if(is.null(rownames(ES))){
          error_message("The ES object does not have rownames or featureData to track the features by. Please provide unique features or rownames for the expression matrix.\n")
          return(NULL)
        }
        
        # Check input_score is provided and are continuous values with no NAs
        if(length(input_score) == 0 || any(!is.numeric(input_score)) || any(is.na(input_score))){
          error_message("input_score must be a vector of continous values (with no NAs) where the vector names match the colnames of the expression matrix.\n")
          return(NULL)
        }
        
        # Make sure the input_score has names or labels that are the same as the colnames of ES
        if(is.null(names(input_score))){
          error_message("The input_score object must have names or labels to track the samples by. Please provide unique sample names or labels that matches the colnames of the expression matrix.\n")
          return(NULL)
        }
        
        # Make sure the input_score has the same length as number of samples in ES
        if(length(input_score) != ncol(ES)){
          error_message("The input_score must have the same length as the number of columns in ES.\n")
          return(NULL)
        }else{
          if(any(!names(input_score) %in% colnames(ES))){
            error_message("The input_score object must have names or labels that match the colnames of the expression matrix.\n")
            return(NULL)
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
              
              dat <- read.csv(inputfile$datapath, header=TRUE) %>% tibble::column_to_rownames(var="Samples") %>% mutate_all(as.numeric)
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
            if(length(weights) == 0 || any(!is.numeric(weights)) || any(is.na(weights))){
              error_message("weights must be a vector of continous values (with no NAs) where the vector names match the colnames of the expression matrix.\n")
              return(NULL)
            }
            
            # Make sure the weights has names or labels that are the same as the colnames of ES
            if(is.null(names(weights))){
              error_message("The weights object must have names or labels to track the samples by. Please provide unique sample names or labels that matches the colnames of the expression matrix.\n")
              return(NULL)
            }
            
            # Make sure the weights has the same length as number of samples in ES
            if(length(weights) != ncol(ES)){
              
              error_message("The weights must have the same length as the number of columns in ES.\n")
              return(NULL)
              
            }else{
              
              if(any(!names(weights) %in% colnames(ES))){
                error_message("The weights object must have names or labels that match the colnames of the expression matrix.\n")
                return(NULL)
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
          top_N = as.integer(input$top_N)
          
          if(is.na(top_N) || length(top_N)==0){
            error_message("Please specify an INTEGER top_N value to evaluate over top N features.\n")
            return(NULL)
          }
          
        }else{
          
          search_start = strsplit(as.character(input$search_start), ",", fixed=TRUE) %>% unlist() %>% trimws()
          top_N = NULL
          
          if(!(search_start %in% rownames(ES))){ 
            error_message("Provided starting feature(s) does not exist among ES's rownames.\n\n")
            return(NULL)
          }
          
        }
        
        permutation = input$permutation_test
        
        if(permutation == TRUE){
          
          n_perm = as.integer(input$n_perm)
          
          if(is.na(n_perm) || length(n_perm)==0){
            error_message("Please specify an INTEGER number of permutations to perform for permutation testings.\n")
            return(NULL)
          }
          
          ncores = as.integer(input$ncores)
          
          if(is.na(ncores) || length(ncores)==0){
            error_message("Please specify the number of cores to perform parallelization for permutation testings.\n")
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
              smooth = TRUE,
              obs_best_score = NULL,
              plot = FALSE, 
              ncores = ncores,
              cache_path = NULL,
              verbose = FALSE
            )
          )
          
        }
        
        # Export the ES and input_score to reactive datase
        feature_set_description(sprintf("After filtering features with Minimum Samples Frequency = %s (or having < %s%% prevalence across all samples) and Maximum Percent Cutoff = %s (or having > %s%% prevalance across all samples), the Binary Feature Set retained %s genomic features out of %s supplied features across %s samples.", min_cutoff, percent_min_cutoff*100, max_cutoff*100, max_cutoff*100, format(nrow(ES), big.mark = ","), format(n_orig_features, big.mark = ","), format(ncol(ES), big.mark =",")))
        feature_set_data(exprs(ES))
        input_score_data(input_score)
        
        # Compute the candidate search algorithm
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
        
        error_message("NONE")
        
      })
      
      output$error_message <- renderUI({
        
        req(error_message())
        
        ns <- session$ns
        
        ## Update loading icon ####
        session$sendCustomMessage(type = "ToggleOperation", message = list(id=ns("loading_icon"), display="no"))
        
        if(error_message() != "NONE"){
          p(style="color: red; font-weight: bold;", error_message())
        }
        
      })
      
      output$featureData_title <- renderUI({
        
        req(feature_set_description())
        
        ns <- session$ns
        
        description <- feature_set_description()
        
        dataset <- isolate({ input$dataset })
        
        if(dataset == "BRCA_GISTIC_MUT_SIG"){
          title <- "Dataset: TTGA BRCA Mutation Signatures (BRCA_GISTIC_MUT_SIG)"
        }else if(dataset == "CCLE_MUT_SCNA"){
          title <- "Dataset: CCLE Mutation in Cancers in Cancers (CCLE_MUT_SCNA)"
        }else if(dataset == "sim.ES"){
          title <- "Dataset: Simulated Data (sim.ES)"
        }else if(dataset == "Import Data"){
          title <- "Dataset: Imported Data"
        }
        
        div(
          h2(title),
          br(),
          p(description),
          downloadButton(outputId = ns("download_featureset"), label="Download Filtered Features Set")
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
        
        req(candidate_search_result())
        
        h2("Best Meta-Features ESet")
      
      })
      
      output$bestFeatureData <- DT::renderDataTable({
        
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
          
          write.csv(Eset_table, file, row.names=FALSE)
          
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
        
        req(input_score_data())
        
        dataset <- isolate({ input$dataset })
        
        if(dataset == "BRCA_GISTIC_MUT_SIG"){
          
          scores <- isolate({ input$BRCA_GISTIC_MUT_SIG_scores })
          
          if(scores == "TAZYAP_BRCA_ACTIVITY"){
            title <- "YAP/TAZ Activity in BRCA (TAZYAP_BRCA_ACTIVITY)"
          }else{
            title <- "Imported Data"
          }
          
        }else if(dataset == "CCLE_MUT_SCNA"){

          scores <- isolate({ input$CCLE_MUT_SCNA_scores })
          
          if(scores == "CTNBB1_reporter"){
            title <- "Activation of B-catenin in Cancers (CTNBB1_reporter)"
          }else{
            title <- "Imported Data"
          }
          
        }else if(dataset == "sim.ES"){
          
          scores <- isolate({ input$sim.ES_scores })
          
          if(scores == "sim.Scores"){
            title <- "Random simulated scores (sim.Scores)"
          }else{
            title <- "Imported Data"
          }
          
        }else if(dataset == "Import Data"){

          title <- "Imported Data"
          
        }
        
        h2("Observed Input Scores:", title)
        
      })
      
      output$inputScoreData <- DT::renderDataTable({
        
        req(input_score_data())
        
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
    
  }
  
  shinyApp(ui=ui, server=server)
  
}







