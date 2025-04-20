# args <- commandArgs(trailingOnly = TRUE)
# args[1] <- tempFolder
# args[2] <- MetaFilePath
# args[3] <- GenoFilePath

library(shiny)
library(readxl)
library(DT)
library(shinyjs)
library(webshot)
library(tidyverse)
library(plotly)
library(shinycssloaders)
require(tidyr)
require(magrittr)
require(dplyr)
require(gt)


tempFolder <- tempfile("")
# userFolder <- basename(tempFolder)

# folder_path <- paste(tempFolder, "Output", sep="\\")
# 
# if (!file.exists(folder_path)) {
#   dir.create(folder_path)
# }

dir.create(tempFolder, recursive = TRUE)

#tempFolder <- "C:\\Users\\RDAS\\OneDrive - CIMMYT\\Desktop\\QAQC"

### outpus files 
tempFolderOut <- paste(tempFolder, "/Out/", sep = "")
MetaFilePath <- tempfile("META_", tmpdir = tempFolder, fileext = ".txt")
GenoFilePath <- tempfile("Geno_", tmpdir = tempFolder, fileext = ".csv")

ms_file <- paste(tempFolderOut, "F1_marker_summary.txt", sep="")

# F1_F1_Summary <- read.delim(paste(tempFolderOut, "F1_F1_Summary.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE)
# Automatically import the dataset
# purity_smry <- read.delim(purity_smry_dat)

# pos_val <- unlist(gregexpr(":", purity_smry$SAMPLE_NAME))
# Geno_name <- substr(purity_smry$SAMPLE_NAME,1,pos_val-1)
# Rep_name <- substr(purity_smry$SAMPLE_NAME,pos_val+1, nchar(purity_smry$SAMPLE_NAME))
# 
# purity_smry <- data.frame(Entry = Geno_name,
#                           Rep = as.numeric(Rep_name),
#                           purity_smry[,-1])
# 
# Mean_Purity <- purity_smry %>%
#   group_by(Entry) %>% 
#   summarise(Reps = length(Rep),
#             Mean_Purity = round(mean(PurityScore,na.rm = T),2)) %>%
#   arrange(Mean_Purity)
# 
# # Display the DataTable with checkboxes and "Select All" option
# Mean_Purity[,1:2] <- lapply( Mean_Purity[,1:2], factor)
# # purity_smry[,1:2] <- lapply( purity_smry[,1:2], factor)
# # Reactive expression for the dataset
# 
# 
# 
# tags$head(
#   tags$style(HTML("
#       /* Custom 3D style for a table with ID #consensus */
#       #consensus .dataTable td {
#         box-shadow: inset 0 0 5px #666;
#         border: 1px solid #ddd;
#       }
#     "))
# )
# 


ui = fluidPage(theme = shinythemes::shinytheme("sandstone"),
               useShinyjs(),  # Initialize shinyjs
               tags$head(
                 tags$style(HTML("
      /* Style the table data elements (cells) */
      .dataTable td {
        border-right: 1px solid #dddddd; /* Add right border to all cells */
      }
      
      /* Optional: Style the last cell in each row to remove the right border if desired */
      .dataTable tr td:last-child {
        border-right: none;
      }
    "))),         
               tags$head(
                 tags$style(HTML("
    .custom-well-panel {
      min-height: 220px; /* Adjust as necessary */
    }
  "))
               ),
               
               navbarPage(
                 title = div(
                   # Include your logo image; adjust the src path to point to your image file
                   img(src = "logoW.png", height = "30px", style = "padding-top: 2px;"),
                   # Optional: Add text or additional elements alongside your logo
                   "Dryland Crops Program QA/QC"
                 ),
                 tabPanel(
                   "Home", 
                   icon = icon("home", lib = "glyphicon"), 
                   
                   div(
                     class = "container-fluid",
                     
                     div(
                       class = "jumbotron",
                       h1("Welcome to the Dryland Crops Program QA/QC"),
                       p("This online web app is developed with the intention of empowering researchers in the use of molecular markers for QA/QC molecular based tests. Perform downstream Bioinformatics analysis for Genetic purity of elite inbred/parental lines and validation of crosses using F1 pedigree verification.", 
                         class = "lead"),
                       hr(),
                       p("This app is packed with 5 modules: Data Upload and Marker Summary, Parent Consensus, Parent QC Selection and F1 verification. These modules are designed to assist researchers in designing hypotheses or answering research questions with little or no expertise in Bioinformatics.")
                     ),
                     
                     
                     div(
                       class = "row",
                       div(
                         class = "col-md-4",
                         title = "Parent Purity Testing", 
                         status = "primary", 
                         solidHeader = TRUE, 
                         background = "teal",  # Custom background color
                         p("This module helps in testing the genetic purity of parent crops.")
                         
                       ),
                       div(
                         class = "col-md-5",
                         title = "F1 Pedigree Verification", 
                         status = "warning", 
                         solidHeader = TRUE, 
                         background = "yellow",  # Custom background color
                         p("Verify the pedigree of F1 crops to ensure correct parentage.")
                         
                       ),
                       div(
                         class = "col-md-3",
                         
                         title = "Data Analysis", 
                         status = "info", 
                         solidHeader = TRUE, 
                         background = "light-blue",  # Custom background color
                         p("Perform comprehensive data analysis to derive meaningful insights.")
                         
                       ))),
                   column(12,imageOutput("image3", width = "300px", height = 233)),
                   # fluidRow(
                   #  column(3,
                   #         imageOutput("image1", width = 86, height = 233),
                   #        ),
                   #  column(6,
                   #         imageOutput("image2"),
                   #        ))
                 ),
                 
                 
                 tabPanel("Data Import",
                          fluidPage(
                            column(5,
                                   wellPanel(class = "custom-well-panel",
                                             fileInput('file1', label = 'Upload sample metadata. Choose a .TXT file',
                                                       accept = c(".txt"), multiple = TRUE, placeholder = "No file selected"),
                                             div(dataTableOutput('contents_file1') ))),
                            
                            fluidPage(
                              column(7,
                                     wellPanel(class = "custom-well-panel",
                                               
                              sidebarLayout(
                                sidebarPanel(width = 6,
                                  fileInput("file2", "Upload intertek snp data. Choose a .CSV file", accept = c(".csv", "text/csv", "text/comma-separated-values,text/plain"))
                                ),
                                  tabsetPanel(
                                    tabPanel("Project Details", gt_output("first12")),
                                    tabPanel("Statistics", DTOutput("remaining"))
                                  )
                                
                              )
                            )))),
                          fluidRow(column(12,
                            # uiOutput("analyze_button_ui"),
                            uiOutput("Visualize_button_ui"),
                            div(id = "processing", style = "display: none;", "Processing, please wait..."),
                            plotlyOutput("plot1",  width = "100%", height = "1200px"),
                            DTOutput("result_ms_table"),
                            verbatimTextOutput("selectedMarkersText"),
                            uiOutput("analyze_button2_ui")         
                            
                            
                            
                          ))
                          
                 ),


                 navbarMenu("Result",
                            tabPanel("Parent Consensus",
                                     fluidRow(
                                       column(12,
                                              DTOutput("consensus")
                                       )
                                     ),
                            ),
                            tabPanel("Parent QC Selection",
                                     fluidRow(
                                       column(12,
                                              plotOutput("display_plot", height = "500px")
                                       ),
                                       column(4,
                                     
                                       tags$div(
                                         tags$h4("Parent Purity Selection Decision:"),
                                         tags$table(
                                           class = "fancy-table table table-bordered", 
                                           tags$thead(
                                             tags$tr(
                                               tags$th("Purity Score"),
                                               tags$th("Comment")
                                             )
                                           ),
                                           tags$tbody(
                                             tags$tr(
                                               tags$td("PurityScore > 80"),
                                               tags$td(style = "background-color: #4CAF50; color: white;", "Pure Line") # Green
                                             ),
                                             tags$tr(
                                               tags$td("PurityScore <= 80"),
                                               tags$td(style = "background-color:  #FFA500; color: white;", "Need to Purify Parent") # Purple
                                             ),

                                           )
                                         )
                                       )),
                                       column(12,
                                              DTOutput("display_table"),
                                       )
                                     )

                            ),
                            tabPanel("F1 Verification",
                                     fluidRow(
                                       column(6,
                                              tags$head(
                                                tags$style(HTML("
      .fancy-table th {
        background-color: #0072B2; /* Header background color */
        color: white; /* Header text color */
      }
    "))
                                              ),
                                              
                                              tags$div(
                                                tags$h4("Plant Selection Decision:"),
                                                tags$table(
                                                  class = "fancy-table table table-bordered", 
                                                  tags$thead(
                                                    tags$tr(
                                                      tags$th("Heterozyous Percentage"),
                                                      tags$th("Parent Score (A & B)"),
                                                      tags$th("Comment")
                                                    )
                                                  ),
                                                  tags$tbody(
                                                    tags$tr(
                                                      tags$td("Percentage_He == 100"),
                                                      tags$td("Both > 80"),
                                                      tags$td(style = "background-color: #4CAF50; color: white;", "successfulF1") # Green
                                                    ),
                                                    tags$tr(
                                                      tags$td("Percentage_He == 100"),
                                                      tags$td("Either <= 80"),
                                                      tags$td(style = "background-color:  #FFA500; color: white;", "PQF") # Purple
                                                    ),
                                                    tags$tr(
                                                      tags$td("Percentage_He >= 60 and Percentage_He < 100"),
                                                      tags$td("Both > 80"),
                                                      tags$td(style = "background-color: #8FD744; color: white;", "possibleF1") # Orange
                                                    ),
                                                    tags$tr(
                                                      tags$td("Percentage_He >= 60 and Percentage_He < 100"),
                                                      tags$td("Either <= 80"),
                                                      tags$td(style = "background-color:  #FFA500; color: white;", "PQF") # Purple
                                                    ),
                                                    tags$tr(
                                                      tags$td("Percentage_He < 60"),
                                                      tags$td("Any"),
                                                      tags$td(style = "background-color: #FF8080; color: white;", "failedF1") # Red
                                                    )
                                                  )
                                                )
                                              ),
                                              
                                       ),
                                       
                                       column(12,
                                              
                                              
                                              DTOutput("F1_verification")
                                       ),
                                     )),

                 ),
                 
                 navbarMenu("More", 
                            tabPanel("Help"
                            ),
                            tabPanel("About Us"
                            ),
                 )
               ))
server = function(input, output){
  
  # MetaFilePath1 <- tempfile("META_", tmpdir = tempFolder, fileext = ".txt")
  output$contents_file1 <- renderDataTable({
    req(input$file1)
    inFile1 <- input$file1
    
    if(is.null(inFile1))
      return(NULL)
    
    # fileext = paste(".",tools::file_ext(inFile1$name), sep = "")
    # MetaFilePath <- tempfile("META_", tmpdir = tempFolder, fileext = fileext)
    file.copy(inFile1$datapath,MetaFilePath)
    dd <- read.delim(inFile1$datapath, header = TRUE, sep = "\t") # or sep = "," etc., depending on your file
    names <- c('SAMPLE_UID','DESIGNATION','SAMPLE_NAME','PARENT_A','PARENT_B')
    dd[,names] <- lapply(dd[,names] , factor)
    datatable(dd, filter = "top", escape = FALSE,
              options = list(dom = 'Blfrtip', autoWidth = TRUE, scrollX = TRUE, pageLength = 11
                             ))
    
  })
  
  
  data1 <- reactive({
    inFile <- input$file2
    if (is.null(inFile)) {
      return(NULL)
    }
    # fileext = paste(".",tools::file_ext(inFile$name), sep = "")
    # GenoFilePath <- tempfile("GENO_", tmpdir = tempFolder, fileext = fileext)
    #file.copy(inFile$datapath, GenoFilePath)
    #read.csv(inFile$datapath, nrows = 12, header = FALSE)
    file.copy(inFile$datapath, GenoFilePath)
    readdata <- read.csv(GenoFilePath, nrows = 12, header = FALSE)
    ##############################################################################################################################
    ## Generate Summary files ##
    ##############################################################################################################################
    #rn2 <- paste("./run_pipe_1.sh", shQuote(tempFolder), shQuote(GenoFilePath), shQuote(MetaFilePath))
    #system(rn2)
    rn2 <- paste("./run_pipe_1.sh", shQuote(tempFolder), shQuote(GenoFilePath), shQuote(MetaFilePath))
    system(rn2)
    readdata
    # ##############################################################################################################################
    # 
  })
  
#  data_temp <- reactive({
#  req(data1())
#  rn2 <- paste("./run_pipe_1.sh", shQuote(tempFolder), shQuote(GenoFilePath), shQuote(MetaFilePath))
#  system(rn2)
#  })
  
  data2 <- reactive({
    inFile <- input$file2
    if (is.null(inFile)) {
      return(NULL)
    }
    
    read.csv(inFile$datapath, skip = 14)
  })
  
  # Output the first 12 rows
  output$first12 <- render_gt({
    df <- data1()
    if (is.null(df)) {
      return()
    }
    gt(df)
    
  })
  
  # Output the remaining rows
  output$remaining <- renderDT({
    df <- data2()
    
    datatable(df, options = list(pageLength = 11))
  })
  
  # GenoFilePath1 <- tempfile("GENO_", tmpdir = tempFolder, fileext = ".csv")
  # output$contents_file2 <- render_gt({
  #   req(input$file2)
  #   
  #   inFile2 <- input$file2
  #   
  #   if(is.null(inFile2))
  #     return(NULL)
  #   
  #   fileext = paste(".",tools::file_ext(inFile2$name), sep = "")
  #   # print(paste("File extension for GENO:", fileext))
  #   
  #   GenoFilePath <- tempfile("GENO_", tmpdir = tempFolder, fileext = fileext)
  #   file.copy(inFile2$datapath, GenoFilePath)
  # 
  #   ###############################################################################################################
  #   ################################## Calling the LGCTool Pipeline ###############################################
    # rn2 <- paste("./run_pipeline.sh", shQuote(tempFolder), shQuote(GenoFilePath1), shQuote(MetaFilePath1))
    # system(rn2)
  #   ################################################################################################################
  #   ################################################################################################################
  #   dd2 <- read.csv(inFile2$datapath, header = FALSE, sep = ",", stringsAsFactors = FALSE, nrows = 12)
  #   
  #   gt(dd2)
  # })
  # 
  # 
  
  output$image1 <- renderImage({
    filename <- normalizePath(file.path('.',
                                        paste('crops', '.png', sep='')))
    
    # Return a list containing the filename
    list(src = filename)
  }, deleteFile = FALSE)
  
  output$image2 <- renderImage({
    filename <- normalizePath(file.path('.',
                                        paste('Linking', '.png', sep='')))
    
    # Return a list containing the filename
    list(src = filename)
  }, deleteFile = FALSE)
  output$image5 <- renderImage({
    filename <- normalizePath(file.path('.',
                                        paste('leaf', '.png', sep='')))
    
    # Return a list containing the filename
    list(src = filename, style = 'position: absolute; opacity: 0.8;')
  }, deleteFile = FALSE)
  
  
  
  file_statuses <- reactiveValues(file1_uploaded = FALSE, file2_uploaded = FALSE)
  
  observeEvent(input$file1, {
    file_statuses$file1_uploaded <- TRUE
  })
  
  observeEvent(input$file2, {
    file_statuses$file2_uploaded <- TRUE
  })
  # Generate the Analyze button UI based on the status of file uploads
  output$analyze_button_ui <- renderUI({
    if(file_statuses$file1_uploaded && file_statuses$file2_uploaded) {
      actionButton("analyze_button", "Run Marker Summary")
    }
  })
  
  output$Visualize_button_ui <- renderUI({
    if(file_statuses$file1_uploaded && file_statuses$file2_uploaded) {
      actionButton("Visualize_button", "Get SNP Data Visualization")
      
    }
  })
  
  
  isTableDataAvailable <- reactive({
    # Check if the data frame 'dd' is not empty
    !is.null(dd) && nrow(dd) > 0
  })
  
  observeEvent(input$analyze_button2, {
    shinyjs::runjs('document.getElementById("analyze_button2").style.backgroundColor = "green";')
  })
  
  # Conditional UI for the update action button
  output$analyze_button2_ui <- renderUI({
    req(tableRendered())  # Assuming tableRendered is defined elsewhere in your server logic
    if(isTableDataAvailable()) {  # Assuming isTableDataAvailable is defined elsewhere in your server logic
      actionButton("analyze_button2", "Analyze")
    }
  })
  
  # Example placeholders for tableRendered() and isTableDataAvailable()
  tableRendered <- reactiveVal(FALSE)
  isTableDataAvailable <- reactive({ TRUE })
  # ... [rest of the server code, including rendering data tables]
  
  observeEvent(input$analyze_button, {
    # Code to perform when analyze_button is clicked
  })
  
  
  
  output$image3 <- renderImage({
    filename <- normalizePath(file.path('.',
                                        paste('DNA', '.png', sep='')))
    
    # Return a list containing the filename
    list(src = filename,  style = 'position: absolute; opacity: 0.9;')
  }, deleteFile = FALSE)
  tableRendered <- reactiveVal(FALSE)
  
  
  jsCallback <- JS(
    "$(document).on('change', 'input.excludeMarkerCheckbox', function() {
    var selectedMarkers = $('input.excludeMarkerCheckbox:checked').map(function() {
      return this.value;
    }).get();
    Shiny.setInputValue('selectedMarkers', selectedMarkers);
  });"
  )
  
  dd <- reactiveVal()
  
  observeEvent(input$Visualize_button, {
    # Assuming ms_file is defined somewhere in your app and accessible here
    dd(read.table(ms_file, header = TRUE, stringsAsFactors = FALSE))
  })
  
  observeEvent(input$Visualize_button, {
    req(dd())
    
    # Correctly updating dd with new values
    updated_dd <- dd() %>%
      mutate(exclude_marker = '') %>%
      select(marker_name, exclude_marker, everything())
    
    dd(updated_dd)  # Update the reactiveVal with the new data frame
    output$result_ms_table <- renderDT({
      datatable(dd(), escape = FALSE,  selection = 'none', extensions = 'Buttons',filter  = "top",options = list( pageLength = 10, dom = 'Blfrtip', # Set initial page length
                                                                                                                lengthMenu = list(c(5, 10, 15, -1), c('5 rows', '10 rows', '15 rows', 'All')) , # Customize page length options
                                                                                                                scrollX = TRUE,buttons = c("selectAll", "copy", "excel", "print"),
                                                                                                                
                                                                                                                columnDefs = list(
                                                                                                                  list(
                                                                                                                    targets = 2,  # Adjust based on your dataframe's structure
                                                                                                                    defaultContent = '',
                                                                                                                    orderable = FALSE,
                                                                                                                    render = JS(
                                                                                                                      "function(data, type, full, meta) {
            return '<input type=\"checkbox\" class=\"excludeMarkerCheckbox\" value=\"' + full[0] + '\">'; 
          }"
                                                                                                                    )   ,   initComplete = JS(jsCallback)
                                                                                                                    
                                                                                                                  )
                                                                                                                ),
                                                                                                                select = list(style = 'os', selector = 'td:nth-child(3)'),
                                                                                                                order = list(list(1, 'asc'))
                                                                                                                
      )) %>%
        
        formatStyle('missing_percentage', 
                    backgroundColor = styleInterval(c(0, 0.2, 0.4, 0.6, 0.8), 
                                                    c('#fcbad3', '#fdb7c2', '#fda2b1', '#fd8da0', '#fd788f', '#ff6b6b'))) %>%
        formatStyle('missing_percentage', 
                    backgroundColor = styleInterval(c(0, 0.2, 0.4, 0.6, 0.8), 
                                                    c('#fcbad3', '#fdb7c2', '#fda2b1', '#fd8da0', '#fd788f', '#ff6b6b'))) %>%
        formatStyle('major_allele_freq', 
                    backgroundColor = styleInterval(c(0, 0.2, 0.4, 0.6, 0.8), 
                                                    c('#ffffba', '#fff7a8', '#fff096', '#ffe984', '#ffe272', '#ff9f1c'))) %>%
        formatStyle('minor_allele_freq', 
                    backgroundColor = styleInterval(c(0, 0.2, 0.4, 0.6, 0.8), 
                                                    c('#baffc9', '#a4f6be', '#8ef2b3', '#78eea8', '#62ea9d', '#17c3b2'))) %>%
        formatStyle('het_proportion', 
                    backgroundColor = styleInterval(c(0, 0.2, 0.4, 0.6, 0.8), 
                                                    c('#e2caff', '#d8c1d1', '#ceb8d3', '#c4afd5', '#baa6d7', '#c3aed6'))) %>%
        formatStyle('PIC', 
                    backgroundColor = styleInterval(c(0, 0.2, 0.4, 0.6, 0.8), 
                                                    c('#e0e0e0', '#c4c4c4', '#a8a8a8', '#8c8c8c', '#707070', '#484848')))
      
    })
    tableRendered(TRUE)
  })
  showProgress <- reactiveVal(FALSE)
  
  output$showProgress <- reactive({
    showProgress()
  })
  outputOptions(output, "showProgress", suspendWhenHidden = FALSE)
  
  
  
  plotdd <- reactiveVal()
  observeEvent(input$Visualize_button, {
  plotdd(read.csv(paste(tempFolderOut, "plot_file.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE))
  })
  
  
  observeEvent(input$Visualize_button, {
    shinyjs::show("processing")  # Show the processing message
    req(plotdd())
    output$plot1 <- renderPlotly({
      Sys.sleep(5)  # Assume it takes 5 seconds to generate the plot, adjust based on your needs
      
      
      # Generate the plot
      p1 <- ggplot(plotdd(), aes(X, Y, col = Call, text = paste("MasterPlate:", MasterPlate,
                                                              "<br>MasterWell:", MasterWell,
                                                              "<br>SubjectID:", SubjectID
                                                              )))  + 
        geom_point(alpha = 0.4) + 
        theme_bw() +
        facet_wrap(~SNPID, ncol = 4)

      
      ggplotly(p1)
    })
    
    shinyjs::delay(100, shinyjs::hide("processing"))
  })
  
  
 
  puritySummary <- reactiveVal()
  observeEvent(input$analyze_button2, {
    puritySummary(read.delim(paste(tempFolderOut, "F1_PuritySummary.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE))
  })
  
  
  # Display the heatmap of Summary dataset
  # Display the heatmap of Summary dataset
  output$display_plot <- renderPlot({
    puritysmry <- req(puritySummary())
    pos_val <- unlist(gregexpr(":", puritysmry$SAMPLE_NAME))
    Geno_name <- substr(puritysmry$SAMPLE_NAME,1,pos_val-1)
    Rep_name <- substr(puritysmry$SAMPLE_NAME,pos_val+1, nchar(puritysmry$SAMPLE_NAME))
    
    puritysmry <- data.frame(Entry = Geno_name,
                             Rep = as.numeric(Rep_name),
                             puritysmry[,-1])

    
    ptbl <- ggplot(puritysmry, aes(x = Entry, y = Rep, group = Entry)) +
      geom_tile(color = "black", aes(fill = PurityScore)) +
      scale_fill_gradient(low = "#F8766D", high = "#66C2A5") +  # Pastel blue to pastel pink
      theme(panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.ontop = TRUE,
            panel.background = element_rect(fill = "transparent"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text.y = element_blank(), legend.position = "top")
    
    ptbl
  })



  Parent_consensus <- reactiveVal()
  observeEvent(input$analyze_button2, {

    Parent_consensus(read.csv(paste(tempFolderOut, "F1_consensus_grid.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE))

  })

  # Parent_consensus <- read.csv(paste(tempFolderOut, "F1_consensus_grid.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  output$consensus <- renderDT({
    parentConsensus <- req(Parent_consensus())
    genotypes_vector <- as.vector(as.matrix(parentConsensus[, 2:ncol(parentConsensus)]))
    unique_genotypes <- unique(genotypes_vector)
    generate_spring_colors <- function(n) {
      # Define the hues for spring colors with a wider range
      hues <- c(seq(0, 360, length.out = n)) # Use a full range of hues
      
      # Generate colors with higher chroma for vibrancy
      colors <- hcl(h = hues, l = 65, c = 50) # Adjust l and c for desired brightness and saturation
      
      return(colors)
    }
    
    # Assuming you've already defined 'genotypes_vector'
    colors <- generate_spring_colors(length(unique_genotypes))
    names(colors) <- unique_genotypes
    parentConsensus$Type <- ifelse(grepl(":", parentConsensus$DNA...Assay), "Sample", "Consensus")
    parentConsensus <- relocate(parentConsensus, Type, .after = DNA...Assay)
    
    parentConsensus[, 1:ncol(parentConsensus)] <- lapply(parentConsensus[, 1:ncol(parentConsensus)], factor)
    datatable(parentConsensus, extensions = 'Buttons', filter = "top",options = list(
      dom = 'Blfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf'),
      scrollX = TRUE,
      lengthMenu = list(c(15, 30, 50, -1), c(15, 30, 50, "All")),
      columnDefs = list(
        list(
          targets = 3:ncol(parentConsensus), 
          render = JS(
            "function(data, type, row, meta) {",
            "  if(type === 'display' && data !== null && data !== '') {",
            "    var color = '#FFC425';", 
            "    switch(data) {",
            "      case 'A:A': color = '#7CAE00'; break;",
            "      case 'T:T': color = '#C77CFF'; break;",
            "      case 'G:G': color = '#F8766D'; break;",
            "      case 'C:C': color = '#00BFC4'; break;",
            "      case '?': color = 'grey'; break;",
            "      case 'N:N': color = 'grey'; break;",
            "      case 'Uncallable': color = 'grey'; break;",
            "      case 'NTC': color = 'grey'; break;",
            "    }",
            "    return '<span style=\"background-color:' + color + '; color: white; padding: 0 4px;\">' + data + '</span>';",
            "  }",
            "  return data;",
            "}"
          )
        )
      )
    )) %>%
      formatStyle(
        'Type',
        backgroundColor = styleEqual(c("Consensus", "Sample"), c("#b2dfee", "#dedede")))
  })
  
  
 
  rv <- reactiveValues(F1SummaryData = NULL)
  
  observeEvent(input$analyze_button2, {
    # Simulate reading data
    tempData <- read.delim(paste(tempFolderOut, "F1_QAQC.txt", sep = ""), sep = ",", header = TRUE, stringsAsFactors = FALSE, na.strings = "na")
    tempData$`Plant Selection` <- ifelse(tempData$Comment %in% c("possibleF1", "successfulF1"), "Select", "No Selection")
    rv$F1SummaryData <- tempData
  })
  
  output$F1_verification <- renderDT({
    req(rv$F1SummaryData)
    # names(tempData)[1] <- "F1"
    columns_to_exclude <- c("Parent_A_Score", "Parent_B_Score", "TotalMarkers" , "Polymorphic", "Call_Hetero" , "Call_Parent_A","Call_Parent_B" , "Call_Missing","Percentage_He")
    rv$F1SummaryData <- rv$F1SummaryData %>%
      mutate(across(-all_of(columns_to_exclude), factor))
    
    rv$F1SummaryData$`Plant Selection` <- as.factor(rv$F1SummaryData$`Plant Selection`)
    plantSelectionColIndex <- which(names(rv$F1SummaryData) == "Plant Selection")
    snp_columns_indices <- grep("^snp", names(rv$F1SummaryData))
    
    # Find the maximum index (the last 'snp' column)
    last_snp_column_index <- max(snp_columns_indices)
    first_snp_column_index <- min(snp_columns_indices)
    
    datatable(rv$F1SummaryData, 
              filter = "top", 
              extensions = "Buttons",
              options = list(dom = 'Blfrtip',  scrollX = TRUE, buttons = c('copy', 'excel', 'print'), lengthMenu = list(c(15, 30, 50, -1), c(15, 30, 50, "All")),
                             columnDefs = list(list(targets = plantSelectionColIndex - 1, className = 'editable'),
                                               list(
                                                 targets = first_snp_column_index:last_snp_column_index, 
                                                 render = JS(
                                                   "function(data, type, row, meta) {",
                                                   "  if(type === 'display' && data !== null && data !== '') {",
                                                   "    var color = '#FFC425';", 
                                                   "    switch(data) {",
                                                   "      case 'A:A': color = '#7CAE00'; break;",
                                                   "      case 'T:T': color = '#C77CFF'; break;",
                                                   "      case 'G:G': color = '#F8766D'; break;",
                                                   "      case 'C:C': color = '#00BFC4'; break;",
                                                   "      case '?': color = 'grey'; break;",
                                                   "      case 'N:N': color = 'grey'; break;",
                                                   "      case 'Uncallable': color = 'grey'; break;",
                                                   "      case 'NTC': color = 'grey'; break;",
                                                   "    }",
                                                   "    return '<span style=\"background-color:' + color + '; color: white; padding: 0 4px;\">' + data + '</span>';",
                                                   "  }",
                                                   "  return data;",
                                                   "}"
                                                 )
                                               )
                             )),
              editable = list(target = 'cell', disable = list(columns = c(0:ncol(rv$F1SummaryData)-1))))%>%
      formatStyle('Comment', target = 'cell', backgroundColor = styleEqual(c("successfulF1", "possibleF1", "failedF1", "PQF"), c("#4CAF50", "#8FD744", "#FF8080", "#FFA500")), color = "white")%>%
      formatStyle(
        columns = c('PLATE_ID', 'WELL', 'SUBJECT_ID', 'BMS_ID', 'DESIGNATION', "Comment", "Plant Selection", "PARENT_A", "PARENT_B"),
        `white-space` = 'nowrap'
      )
  })
  observeEvent(input$F1_verification_cell_edit, {
    info <- input$F1_verification_cell_edit
    str(info) # For debugging
    # Update the reactiveValues data
    rv$F1SummaryData[info$row, info$col] <- DT::coerceValue(info$value, rv$F1SummaryData[info$row, info$col])
  })


  # F1_F1_Summary <- reactiveVal()
  # observeEvent(input$analyze_button2, {
  #   F1_F1_Summary(read.delim(paste(tempFolderOut, "F1_F1_Summary.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE))
  #   
  # })
  # 
  # 
  # output$F1_verification <- renderDT({
  #   F1SummaryData <- req(F1_F1_Summary())
  #   names(F1SummaryData)[1] <- "F1"
  #   F1SummaryData$F1 <- as.factor(F1SummaryData$F1)
  #   F1SummaryData$Parent_A <- as.factor(F1SummaryData$Parent_A)
  #   F1SummaryData$Parent_B <- as.factor(F1SummaryData$Parent_B)
  #   F1SummaryData$`Plant Selection` <- ifelse(F1SummaryData$Comment %in% c("possibleF1", "successfulF1"), "Select", "No Selection")
  #   F1SummaryData$`Plant Selection` <- as.factor(F1SummaryData$`Plant Selection`)
  #   
  #   datatable(F1SummaryData, 
  #             filter = "top", 
  #             extensions = "Buttons", 
  #             options = list(dom = 'Blfrtip', 
  #                            buttons = c('copy', 'excel', 'print'), 
  #                            lengthMenu = list(c(15, 30, 50, -1), c(15, 30, 50, "All")),
  #                            columnDefs = list(list(targets = ncol(F1SummaryData), # Target the last column (Plant Selection)
  #                                                   className = 'editable'))),
  #             editable = TRUE) %>% # Make the DataTable editable
  #     formatStyle(
  #       'Comment',
  #       target = 'cell',
  #       backgroundColor = styleEqual(
  #         c("successfulF1", "possibleF1", "failedF1", "PQF"),
  #         c("#4CAF50", "#8FD744", "#FF8080", "#FFA500")
  #       ),
  #       color = "white"
  #     )
  # })
  # observeEvent(input$F1_verification_cell_edit, {
  #   info <- input$F1_verification_cell_edit
  #   str(info) # For debugging, to see the structure of the edit event
  #   # Update the data frame with the new value
  #   F1SummaryData[info$row, info$col] <<- DT::coerceValue(info$value, F1SummaryData[info$row, info$col])
  # })
  # 
  # 

  
  
  
  rv1 <- reactiveValues(ParentPuritySummaryData = NULL)
  
  observeEvent(input$analyze_button2, {
    # Simulate reading data
    tempData <- read.delim(paste(tempFolderOut, "Parent_QAQC.txt", sep = ""), sep = ",", header = TRUE, stringsAsFactors = FALSE)
    # tempData <- tempData %>%
    #   mutate(Parent = sub(":.*", "", SAMPLE_NAME))  
    # tempData <- relocate(tempData, Parent, .before = SAMPLE_NAME)
    tempData$`Plant Selection` <- ifelse(tempData$Comment %in% c("Pure Line"), "Select", "No Selection")
    
    rv$ParentPuritySummaryData <- tempData
    
    
  })
  
  output$display_table <- renderDT({
    req(rv$ParentPuritySummaryData)
    columns_to_exclude <- c("TotalMakers" , "Match" , "MisMatch","BothMissing","Consensus_Missing", "SampleMissing","PurityScore")
    rv$ParentPuritySummaryData <- rv$ParentPuritySummaryData %>%
      mutate(across(-all_of(columns_to_exclude), factor))
    rv$ParentPuritySummaryData$`Plant Selection` <- as.factor(rv$ParentPuritySummaryData$`Plant Selection`)
    
    unique_designations <- unique( rv$ParentPuritySummaryData$DESIGNATION)
    colors <- grDevices::rainbow(length(unique_designations), s = 0.6, v = 0.85, alpha = 0.5) # Adjust 'alpha' for transparency
    
    # Convert colors to RGBA
    colors_rgba <- sapply(colors, function(col) {
      col_rgb <- col2rgb(col, alpha = TRUE)
      sprintf("rgba(%d,%d,%d,%f)", col_rgb[1,], col_rgb[2,], col_rgb[3,], col_rgb[4,]/255)
    })
    colors_map <- setNames(colors_rgba, unique_designations)
    
    snp_columns_indices <- grep("^snp", names(rv$ParentPuritySummaryData))
    plantSelectionColIndex <- which(names(rv$ParentPuritySummaryData) == "Plant Selection")
    
    # Find the maximum index (the last 'snp' column)
    last_snp_column_index <- max(snp_columns_indices)
    first_snp_column_index <- min(snp_columns_indices)
    
    # Map the generated colors to the unique groups
    datatable(rv$ParentPuritySummaryData, 
              filter = "top", 
              extensions = "Buttons", 
              options = list(dom = 'Blfrtip', scrollX = TRUE,buttons = c('copy', 'excel', 'print'), lengthMenu = list(c(15, 30, 50, -1), c(15, 30, 50, "All")),
                             columnDefs = list(list(targets = plantSelectionColIndex - 1, className = 'editable'),
                                               list(
                                                 targets = first_snp_column_index:last_snp_column_index, 
                                                 render = JS(
                                                   "function(data, type, row, meta) {",
                                                   "  if(type === 'display' && data !== null && data !== '') {",
                                                   "    var color = '#FFC425';", 
                                                   "    switch(data) {",
                                                   "      case 'A:A': color = '#7CAE00'; break;",
                                                   "      case 'T:T': color = '#C77CFF'; break;",
                                                   "      case 'G:G': color = '#F8766D'; break;",
                                                   "      case 'C:C': color = '#00BFC4'; break;",
                                                   "      case '?': color = 'grey'; break;",
                                                   "      case 'N:N': color = 'grey'; break;",
                                                   "      case 'Uncallable': color = 'grey'; break;",
                                                   "      case 'NTC': color = 'grey'; break;",
                                                   "    }",
                                                   "    return '<span style=\"background-color:' + color + '; color: white; padding: 0 4px;\">' + data + '</span>';",
                                                   "  }",
                                                   "  return data;",
                                                   "}"
                                                 )
                                               )
                             )),
              editable = list(target = 'cell', disable = list(columns = c(0:ncol(rv$ParentPuritySummaryData)-1))))%>%
      formatStyle('DESIGNATION', 
                  backgroundColor = styleEqual(unique_designations, colors)) %>%
      formatStyle('Comment', target = 'cell', backgroundColor = styleEqual(c("Pure Line", "Need to Purify Parent"), c("#4CAF50", "#FFA500")), color = "white")%>%
      
      formatStyle(
        columns = c('PLATE_ID', 'WELL', 'SUBJECT_ID', 'BMS_ID', 'DESIGNATION', "Comment", "Plant Selection"),
        `white-space` = 'nowrap')
  })
  
  observeEvent(input$display_table_cell_edit, {
    info <- input$display_table_cell_edit
    str(info) # For debugging
    # Update the reactiveValues data
    rv$ParentPuritySummaryData[info$row, info$col] <- DT::coerceValue(info$value, rv$ParentPuritySummaryData[info$row, info$col])
  })
  
  
  
  # Reactive expression to fetch selected marker names based on indices
  # Reactive expression to fetch selected marker names based on indices
  selected_marker_names <- reactive({
    req(dd())
    # Ensure input$selectedMarkers is not NULL and has at least one selection
    if (!is.null(input$selectedMarkers) && length(input$selectedMarkers) > 0) {
      selected_indices <- as.numeric(input$selectedMarkers)  # Convert indices to numeric
      marker_names <- dd()[selected_indices, "marker_name", drop = TRUE]  # Extract marker names
      return(marker_names)
    } else {
      return(character(0))  # Return an empty character vector if no selection
    }
  })
  selectionMade <- reactiveVal(FALSE)
  observeEvent(input$selectedMarkers, {
    # User has made a selection, update the reactive value
    selectionMade(TRUE)
  })
  output$selectedMarkersText <- renderText({
    if (!selectionMade()) {
      ""  # Don't display any message by default
    } else {
      # Fetch the reactive value of selected marker names
      marker_names <- selected_marker_names()
      
      if (length(marker_names) == 0) {
        "No markers selected."
      } else {
        paste("Selected Marker Names:", paste(marker_names, collapse = ", "))
      }
    }
  })
  
  observe({
    # Fetch the reactive value of selected marker names
    marker_names <- selected_marker_names()
    
    if (length(marker_names) > 0) {
      # Specify the path to the text file where you want to save the names
      save_path <-paste(tempFolderOut,"excluded_markers.txt", sep = "") # Replace with your directory path
      
      # Save the selected marker names to the file
      writeLines(marker_names, con = save_path)
      
    }
    
    else {
      excluded_markers <- "No markers selected"
      save_path <-paste(tempFolderOut,"excluded_markers.txt", sep = "") # Replace with your directory path
      writeLines(excluded_markers, save_path)
    }
    
    ##############################################################################################################################
    #### Generate F1 summary
    ##############################################################################################################################
    rn3 <- paste("./run_pipe_2.sh", shQuote(tempFolder))
    system(rn3)
    ##############################################################################################################################
  })
  

}

shinyApp(ui = ui, server = server) 
