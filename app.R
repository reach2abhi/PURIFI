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
#############################
######################################################################################################################################
###Read metadata and intertek fle
#setwd(tempdir())
read_datasets <- function(file_path1, file_path2) {
  # Check if the file exists
  if (!file.exists(file_path1)) {
    stop("File does not exist.")
  }
  if (!file.exists(file_path2)) {
    stop("File does not exist.")
  }
  meta_data<-read.table(file_path1, header = TRUE, sep = "\t")
  data <- readLines(file_path2)
  skip_lines<-grep("Data", data)
  skip_lines
  #Creating plot_file
  Geno_data <- read_csv(file_path2, skip = skip_lines)
  start_index <- which(grepl("Statistics", data)) + 1
  end_index <- which(grepl("DNA", data)) - 3
  alleles_data<-read.csv(file_path2, skip = start_index - 1, nrows = end_index - start_index + 1, header = TRUE)
  alleles_data$alleles <- paste(alleles_data$Allele.Y, alleles_data$Allele.X, sep = "/")
  alleles_d <- alleles_data %>%
    select(SNP, alleles)%>%distinct()
  #Return data that can be used as input for plotting
  return(list(metadata = meta_data, lgc_data = Geno_data,allele=alleles_d))
}

###Create HapMap
create_hapmap<-function(datasets){
  grid_df<-create_grid_file(datasets$lgc_data,datasets$metadata)
  snp_vector <- datasets$allele$SNP
  hap1<-grid_df%>%select(SAMPLE_NAME,all_of(snp_vector))
  rownames(hap1) <- NULL
  hap2 <- t(hap1)
  hap2 <- as.data.frame(hap2)
  colnames(hap2) <- hap2[1, ]
  hap2 <- hap2[-1, ]
  hap2$`SNP` <- rownames(hap2)
  hap2 <- hap2[, c(ncol(hap2), 1:(ncol(hap2)-1))]
  rownames(hap2) <- NULL
  hap3<-merge(hap2, datasets$allele, by = "SNP")
  hap3<-hap3 %>%
    relocate(alleles, .after = SNP)

  hap4 <- hap3 %>%
    mutate(
      chrom = 0,
      pos = 0,
      strand = "+",
      `assembly#` = NA,
      center = NA,
      protLSID = NA,
      assayLSID = NA,
      panelLSID = NA,
      QCcode = NA
    )

    hapmap <- hap4 %>%
    relocate(chrom, pos, strand, `assembly#`, center, protLSID, assayLSID, panelLSID, QCcode, .after = alleles)
    return(hapmap)

}


#################################################################################################
#Create a grid file with metadata
##################################################################################################
create_grid_file <- function(geno_data,meta_data) {
  MISSING_CALLS<-c('N:N', '?:?', "?", "Uncallable", "Unused", "missing", "Empty", "NTC")
  BLANK_SAMPLES <-c('NTC', 'Empty', "")
  #replace missing values
  geno_data$Call[geno_data$Call %in% MISSING_CALLS] <- "?:?"
  geno_data$SubjectID[geno_data$SubjectID %in% BLANK_SAMPLES] <- "NA"
  grid_draft <- geno_data %>%
    select("SubjectID", "SNPID", "Call") %>%
    spread(key = "SNPID", value = "Call")%>%rename("DNA /\ Assay" = "SubjectID")
  grid2<-grid_draft%>%rename("SAMPLE_UID"="DNA /\ Assay")
  grid_df <- merge(meta_data, grid2, by = "SAMPLE_UID")
  
  #write.csv(grid_draft, "Test2_Grid_file.csv", row.names = FALSE, quote = FALSE)
  return (grid_df)
}
#grid_df<-create_grid_file(datasets$lgc_data,datasets$metadata)
#####################################################################################################################
## Calculate SNP summary
calculate_snp_info <- function(alleles, total_samples) {
  missing_count <- sum(alleles == "?:?")
  non_missing_alleles <- alleles[alleles != "?:?"]
  
  # Flatten the alleles and count individual allele occurrences
  allele_list <- unlist(strsplit(non_missing_alleles, ":"))
  allele_counts <- table(allele_list)
  total_alleles <- sum(allele_counts)
  
  # Major and Minor Alleles (most frequent and second most frequent)
  sorted_alleles <- sort(allele_counts, decreasing = TRUE)
  major_allele <- names(sorted_alleles)[1]
  minor_allele <- names(sorted_alleles)[length(sorted_alleles)]
  
  # Frequencies of major and minor alleles
  major_allele_count <- sorted_alleles[1]
  major_allele_freq <- round(major_allele_count / total_alleles, 2)
  minor_allele_count <- sorted_alleles[length(sorted_alleles)]
  minor_allele_freq <- round(minor_allele_count / total_alleles, 2)
  
  # Heterozygote count: Check if the allele pair consists of two different alleles
  het_count <- sum(sapply(non_missing_alleles, function(x) {
    alleles_pair <- strsplit(x, ":")[[1]]
    length(alleles_pair) == 2 && alleles_pair[1] != alleles_pair[2]
  }))
  
  het_proportion <- round(het_count / total_samples, 2)
  
  # PIC Calculation: Based on all alleles present
  #pic <- round(1 - (major_allele_freq^2 + minor_allele_freq^2),2)
  pic <- round(1 - ((major_allele_freq^2) + (minor_allele_freq^2)) - (2 * (major_allele_freq^2) * (minor_allele_freq^2)),2)
  pic <- max(pic, 0)
  allele_frequencies <- allele_counts / total_alleles
  #pic <- round(1 - sum(allele_frequencies^2), 2)
  
  return(data.frame(
    total_samples = total_samples,
    missing_count = missing_count,
    missing_percentage = round(missing_count / total_samples * 100, 2),
    major_allele = major_allele,
    major_allele_count = major_allele_count,
    major_allele_freq = major_allele_freq,
    minor_allele = minor_allele,
    minor_allele_count = minor_allele_count,
    minor_allele_freq = minor_allele_freq,
    het_count = het_count,
    het_proportion = het_proportion,
    PIC = pic
  ))
}
## Calculate Sample summary

compute_sample_summary <- function(df) {
  results <- list()
  for (i in 1:nrow(df)) {
    sample_name <- df[i, 2] 
    #print(sample_name)
    snps <- df[i, -c(1:6)]  
    #print(snps)
    total_snps <- length(snps)
    missing_count <- sum(snps == "?:?")
    missing_percentage <- round((missing_count / total_snps) * 100,3)
    alleles <- unlist(strsplit(as.character(snps), ":"))
    alleles <- alleles[alleles != "?"]  
    allele_freq <- table(alleles)
    
    if (length(allele_freq) > 0) {
      major_allele <- names(which.max(allele_freq))
      major_allele_count <- allele_freq[major_allele]
      minor_allele <- names(which.min(allele_freq))
      minor_allele_count <- allele_freq[minor_allele]
      major_allele_freq <- round(major_allele_count / sum(allele_freq),3)
      minor_allele_freq <- round(minor_allele_count / sum(allele_freq),3)
    } else {
      major_allele <- NA
      major_allele_count <- 0
      minor_allele <- NA
      minor_allele_count <- 0
      major_allele_freq <- 0
      minor_allele_freq <- 0
    }
    het_count <- sum(sapply(snps, function(x) {
      if (x == "?:?") {
        return(FALSE) 
      }
      alleles <- unlist(strsplit(x, ":"))
      return(alleles[1] != alleles[2]) 
    }))
    het_proportion <- het_count / total_snps
    results[[i]] <- data.frame(
      samplename = sample_name,
      total_snps = total_snps,
      missing_count = missing_count,
      missing_percentage = round(missing_percentage,2),
      major_allele = major_allele,
      major_allele_count = major_allele_count,
      minor_allele = minor_allele,
      minor_allele_count = minor_allele_count,
      major_allele_frequency = round(major_allele_freq,2),
      minor_allele_frequency = round(minor_allele_freq,2),
      het_count = het_count,
      het_proportion = round(het_proportion,2)
    )
  }
  summary_df <- do.call(rbind, results)
  return(summary_df)
}
##############################################################################################
##Create parents consensus
########################################################################################
create_parents_consensus<-function(grid_df){
  parents_df<-grid_df%>%filter(GENERATION=="Parent")
  header<-names(parents_df)
  snp_s<-header[7:length(header)]
  k<-c('DESIGNATION','SAMPLE_NAME',snp_s)
  parents_df2<-parents_df%>%select(all_of(k))
  #head(parents_df2)
  # Create a grouping by designation
  grouped <- parents_df2 %>%
    group_by(DESIGNATION)
  consensus_df <- data.frame()
  # Iterate through each designation
  for (group_name in unique(parents_df2$DESIGNATION)) {
    group_data <- parents_df2 %>% filter(DESIGNATION == group_name)
    
    # Calculate the mode for each column (ignoring Designation and SAMPLE_NAME)
    column_frequencies <- list()
    excl <- c("DESIGNATION", "SAMPLE_NAME")
    #excl
    
    for (column in colnames(group_data)) {
      if (!(column %in% excl)) {
        freq_table <- table(group_data[[column]])
        max_freq <- max(freq_table)
        most_common_entries <- names(freq_table[freq_table == max_freq])
        #print(length(most_common_entries))
        
        if (length(most_common_entries) == 1) {
          column_frequencies[[column]] <- most_common_entries[1]
        } else {
          column_frequencies[[column]] <- "Inconclusive"
        }
      }
    }
    #print(column_frequencies)
    # Create a dataframe from column frequencies
    result_df <- data.frame(t(column_frequencies))
    # Add the new row for "Prefix" and "DNA \ Assay" with "Consensus"
    new_row <- data.frame(DESIGNATION = group_name, "SAMPLE_NAME" = "Consensus")
    df1 <- bind_cols(new_row, result_df)
    consensus_df <- bind_rows(consensus_df, df1)
  }
  parents_consensus_df<-rbind(parents_df2,consensus_df)
  #sorted_final_df <- final_df[order(final_df$DESIGNATION), ]
  sorted_parent_consensus_df <- parents_consensus_df[order(parents_consensus_df$DESIGNATION,
                                                           parents_consensus_df$SAMPLE_NAME), ]
  #return(sorted_parent_consensus_df)
  return(list(parent_cons_df = sorted_parent_consensus_df, cons_df = consensus_df))
  
}

########################################################################################################
#Compute parental Purity
##################################################################################################
parental_purity<-function(parent_consensus){
  # Group data by the 'DESIGNATION' column
  group_by_desig <- parent_consensus %>% group_by(DESIGNATION)
  parents_purity <- file(file.path(tempdir(),"Sample_Purity_Summary.txt"), "wt")
  purity_header <- c("SAMPLE_NAME", "TotalMarkers", "Match", "Mismatch", "BothMissing", "ConsensusMissing", "SampleMissing", "PurityScore")
  writeLines(paste(purity_header, collapse = ","), parents_purity)
  #Loop over each group
  group_by_desig %>%mutate(across(everything(), as.character))%>%
    group_walk(function(group_df, group_name) {
      #print(as.character(group_df[1, 1:ncol(group_df)]))
      first_row <- as.character(group_df[1, 2:ncol(group_df)])
      #ww<-as.character(group_df[1, 1:ncol(group_df)])
      for (i in 2:nrow(group_df)) {
        row <- group_df[i, ]
        row_values <- as.character(row[2:ncol(row)])
        matches <- 0
        mismatches <- 0
        missing <- 0
        both_mis <- 0
        cons_miss <- 0
        mark_miss <- 0
        # Compare the first row with the current row
        for (j in seq_along(first_row)) {
          ref_value <- first_row[j]
          current_value <- row_values[j]
          #print(c(j,ref_value,current_value))
          
          # Handle missing values as "?" (missing values)
          if (ref_value == "?" && current_value == "?") {
            both_mis <- both_mis + 1
          } else if (ref_value == current_value) {
            matches <- matches + 1
          } else if (ref_value != current_value) {
            mismatches <- mismatches + 1
            if (ref_value == "Inconclusive") {
              cons_miss <- cons_miss + 1
            } else if (current_value == "?") {
              mark_miss <- mark_miss + 1
            }
          }
        }
        
        # Calculate total number of markers
        total_rows <- length(row_values)
        percent_match <- (matches / total_rows) * 100
        
        # Prepare final data
        purity_data <- c(
          as.character(row[1]), 
          as.character(total_rows),
          as.character(matches),
          as.character(mismatches),
          as.character(both_mis),
          as.character(cons_miss),
          as.character(mark_miss),
          sprintf("%.2f", percent_match)
        )
        #print(purity_data)
        writeLines(paste(purity_data, collapse = ","), parents_purity)
      }
    })
  
  # Close the output file
  close(parents_purity)
  Purity_df<- read.csv(file.path(tempdir(),"Sample_Purity_Summary.txt"))
  #"Sample_Purity_Summary.txt")
  return (Purity_df)
  
}
#Purity_df<-parental_purity(parent_consensus)


###########################################################################################################
##Create purity Report
##########################################################################################################
create_purity_report<-function(Purity_df,consensus_df){
  purity_report <- Purity_df %>%
    mutate(DESIGNATION = sub("(:.*)", "", SAMPLE_NAME)) %>% # Extract the prefix before ':'
    group_by(DESIGNATION) %>%
    summarise(
      NoReps = n(), 
      PurityScore = mean(PurityScore)
    ) %>%ungroup()
  #head(purity_report)
  consensus_report <- merge(purity_report,consensus_df, by = "DESIGNATION")
  consensus_report[] <- lapply(consensus_report, function(x) {
    if (is.list(x)) {
      return(sapply(x, function(y) paste(unlist(y), collapse = ", ")))  # Flatten the list column
    }
    return(x)  # Return non-list columns as is
  })
  write.csv(consensus_report, file.path(tempdir(),"consensus_report_grid2.csv"),row.names = FALSE, quote = FALSE)
  return(consensus_report)
  
}
#cons_report<-create_purity_report(Purity_df,consensus_df)
###########################################################################################
######################### Find poly SNPs
###########################################################################################
find_homo_poly_snps <- function(parent_a, parent_b, call) {
  poly_count <- 0
  het_poly_count <- 0
  call_parA <- 0
  call_parB <- 0
  score_A <- parent_a[1]
  score_B <- parent_b[1]
  #print(call)
  total_markers <- length(parent_a) - 1 #Exclude the score in the length
  call_missing <- 0  # Count of missing data in call
  #print(total_markers)
  for (i in 2:length(parent_a)) {
    snp1 <- parent_a[i]
    snp2 <- parent_b[i]
    snp_call <- call[i - 1] #i-i since parent have score in the first index, hence begin at index2 and f1 begin at index1
    #print(snp_call)
    
    if (snp1 == 'Inconclusive' || snp2 == 'Inconclusive') {
      next
    }
    
    if (snp_call %in% c('Uncalled', 'N:N', '?')) {
      call_missing <- call_missing + 1
      next
    }
    
    a1_1 <- unlist(strsplit(snp1, ":"))[1]
    a1_2 <- unlist(strsplit(snp1, ":"))[2]
    a2_1 <- unlist(strsplit(snp2, ":"))[1]
    a2_2 <- unlist(strsplit(snp2, ":"))[2]
    
    if (grepl(":", snp_call)) {
      c1_1 <- unlist(strsplit(snp_call, ":"))[1]
      c1_2 <- unlist(strsplit(snp_call, ":"))[2]
    } else {
      c1_1 <- c1_2 <- NA
    }
    
    if (snp1 != snp2 && a1_1 == a1_2 && a2_1 == a2_2) {
      poly_count <- poly_count + 1
      if (c1_1 != c1_2) {
        het_poly_count <- het_poly_count + 1
      } else {
        if (c1_1 == a1_1 && c1_2 == a1_2) {
          call_parA <- call_parA + 1
        }
        if (c1_1 == a2_1 && c1_2 == a2_2) {
          call_parB <- call_parB + 1
        }
      }
    }
  }
  
  perc_het <- ifelse(poly_count > 0, round(100 * het_poly_count / poly_count, 2), 0.00)
  #print(c(score_A, score_B, total_markers, poly_count, het_poly_count, call_parA, call_parB, call_missing, perc_het))
  
  return(c(score_A, score_B, total_markers, poly_count, het_poly_count, call_parA, call_parB, call_missing, perc_het))
}

###########################################################
##Function to create F1 Report
#################################################
create_f1<-function(grid_df){
  f1_grid<-grid_df%>%filter(GENERATION=="F1")
  #write.csv(f1_grid, "f11_grid.csv",row.names = FALSE, quote = FALSE)
  return(f1_grid)
}

create_f1_report<-function(cons_report,f1_grid,het_threshold, purity_threshold){
  parents <- list()
  # Store the parents data in a dictionary
  for (i in 1:nrow(cons_report)) {
    dat <- unname(unlist(cons_report[i, ]))  
    #print(dat)
    parents[[dat[1]]] <- c(dat[3], dat[5:length(dat)])  
  }

  ##Create f1 report dataframe
  f_header<-c("F1_one_Name","Parent_A","Parent_A_Score","Parent_B", "Parent_B_Score", "TotalMarkers", "Polymorphic", "Call_Hetero", "Call_Parent_A", "Call_Parent_B", "Call_Missing", "Percentage_Het")
  f1_report_out <- data.frame(matrix(ncol = length(f_header), nrow = 0))
  for (i in 1:nrow(f1_grid)) {
    data<-unname(unlist(f1_grid[i,]))
    #  #names are keys in R
    if (data[4] %in% names(parents) && data[5] %in% names(parents)){
      #print(data[4])
      #print(data[5])
      PA<-as.character(unlist(parents[data[4]]))
      PB<-as.character(unlist(parents[data[5]]))
      f1_out<-find_homo_poly_snps(PA,PB,data[7:length(data)])
      #print(c(data[2],data[4],data[5]))
      f1_out2<-c(data[2],data[4],f1_out[1],data[5],f1_out[2], f1_out[3:length(f1_out)])
      #print(f1_out2)
      f1_report_out <- rbind(f1_report_out, f1_out2)
    }
    
  }
  
  colnames(f1_report_out) <- f_header
  #print(names(f1_report_out))
  
  f1_report_out <- f1_report_out %>%
    mutate(comment = ifelse(Percentage_Het == 100, ifelse(Parent_A_Score >= purity_threshold & Parent_B_Score >= purity_threshold,"successfulF1","PQF"), 
                            ifelse(Percentage_Het >= het_threshold, 
                                   ifelse(Parent_A_Score >= purity_threshold & Parent_B_Score >= purity_threshold,"PossibleF1","PQF"),
                                   "FailedF1")
    )
    )

  return(f1_report_out)
}

##################################################
ui = fluidPage(theme = shinythemes::shinytheme("flatly"),
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
      
      /* Style for the checkbox label in red */
      .checkbox-label-red {
        color: red;
        font-weight: bold;
      }
      
    "))),         
               tags$head(
                 tags$style(HTML("
    .custom-well-panel {
      min-height: 250px; /* Adjust as necessary */
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
                     
                     column(8,
                            class = "jumbotron",
                            h1("Dryland Crops PURIFI program"),
                            h2("This online web app is developed with the intention of empowering researchers in the use of molecular markers for QA/QC molecular based tests. Perform downstream Bioinformatics analysis for genetic purity of elite inbred/parental lines and validation of crosses using F1 pedigree verification.", 
                               class = "lead"),
                            hr(),
                            p("
                         Featuring four integral components: Data Upload and Marker Summary, Parent Consensus, Parent QC Selection, and F1 Verification, the application is tailored to support researchers in formulating hypotheses or addressing scientific inquiries, even those with minimal bioinformatics background. The Parent QC Selection tool specifically enables the detailed examination of parent plant genetic consistency, crucial for the propagation of desired traits. Meanwhile, the F1 Verification function offers a critical assessment tool, ensuring the genetic lineage and purity of F1 hybrids, confirming that they are true progenies of their intended parentage.")
                     ),
                     
                     fluidRow(
                       column(4,
                              wellPanel(
                                class = "well-custom",
                                # Use an inline image next to the "Parent Consensus" title
                                h4(img(src = "img1.png", height = "75px", style = "vertical-align:middle; margin-right: 5px;"), "Parent Consensus", class="well-title", style="display: inline;"),
                                p("Parent Consensus analysis aggregates genetic data from parent specimens to identify shared markers, ensuring breeding consistency. It confirms parental line homogeneity, enhancing offspring trait predictability."),
                                style = "background-color: #aca4e2; color: white;" # Pastel Blue
                              )
                       ),
                       column(4,
                              wellPanel(
                                class = "well-custom",
                                h4(img(src = "img2.png", height = "75px", style = "vertical-align:middle; margin-right: 5px;"), "Parent QC Selection", class="well-title", style="display: inline;"),
                                p("Parent Purity testing ensures genetic consistency in parent plants, vital for transmitting desired traits to offspring. Utilizing molecular markers, it underpins breeding predictability and the development of superior hybrids."),
                                style = "background-color: #7db0dd; color: white;" # Pastel Orange
                              )
                       ),
                       column(4,
                              wellPanel(
                                class = "well-custom",
                                h4(img(src = "img3.png", height = "75px", style = "vertical-align:middle; margin-right: 5px;"), "F1 Verification", class="well-title", style="display: inline;"),
                                p("F1 Pedigree Verification uses molecular markers to confirm the genetic lineage of F1 hybrids, ensuring offspring are true products of intended parental crosses. This precise verification process is pivotal in upholding the integrity of hybrid breeding programs."),
                                style = "background-color: #4cb9cc; color: white;" # Pastel Green
                              )
                       ))),
                   tags$head(
                     tags$style(HTML("
      .full-width-image img {
        width: 98.5%;
        height: auto;
      }
    "))
                   ),
                   column(12, div(class = "full-width-image", imageOutput("image3")))
                   # fluidRow(
                   #  column(3,
                   #         imageOutput("image1", width = 86, height = 233),
                   #        ),
                   #  column(6,
                   #         imageOutput("image2"),
                   #        ))
                 ),
                 
                 
                 tabPanel("Data Upload & Marker Visualization",
                          fluidPage(
                            column(5,
                                   wellPanel(class = "custom-well-panel",
                                             fileInput('file1', label = 'Upload sample metadata. Choose a .TXT file',
                                                       accept = c(".txt"), multiple = TRUE, placeholder = "No file selected"),
                                             actionButton("load_example_meta", "Load Sample Meta Example", class = "btn-info"),
                                             
                                             div(dataTableOutput('contents_file1') ))),
                            
                            column(7,
                                   wellPanel(class = "custom-well-panel",
                                             
                                             sidebarLayout(
                                               sidebarPanel(width = 6,
                                                            fileInput("file2", "Upload intertek snp data. Choose a .CSV file", accept = c(".csv", "text/csv", "text/comma-separated-values,text/plain")),
                                                            actionButton("load_example_geno", "Load Intertek Marker Example", class = "btn-info"),
                                                            
                                               ),
                                               tabsetPanel(
                                                 tabPanel("Project Details", gt_output("first12")),
                                                 tabPanel("Statistics", DTOutput("remaining"))
                                               )
                                               
                                             )
                                   ))),
                          
                          tags$head(
                            tags$style(HTML(
                              "
      .visualize-button {
        background-color: #00BFC4; /* Green background */
        color: white; /* White text */
      }
      .download-button {
        background-color: #7CAE00; /* Red background */
        color: white; /* White text */
      }
      .visualize-button:hover {
        background-color: #45a049; /* Darker green on hover */
      }
      .download-button:hover {
        background-color: #e53935; /* Darker red on hover */
      }
      "
                            ))
                          ),
                          fluidRow(column(12,
                                          # uiOutput("analyze_button_ui"),
                                          uiOutput("Visualize_button_ui"),
                                          
                                          div(id = "processing", style = "display: none;", "Processing, please wait..."),
                                          plotlyOutput("plot1",  width = "100%", height = "1200px"),
                                          DTOutput("result_ms_table"),
                                          verbatimTextOutput("selectedMarkersText"),
                                          uiOutput("analyze_button2_ui"),
                                          uiOutput("download_button_ui")
                                          
                                          
                                          
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
                                              wellPanel(
                                                h4("Adjust Purity Score Threshold:"),
                                                sliderInput("purity_threshold", "Select Threshold:", min = 0, max = 100, value = 80),
                                                # checkboxInput("drop_below_threshold", 
                                                #               label = tags$span("Drop parents below threshold", class = "checkbox-label-red"), 
                                                #               value = FALSE),   
                                                style = "height: 200px;" # Set the fixed height here
                                                
                                              )
                                       ),
                                       column(4,
                                              wellPanel(
                                                
                                                uiOutput("purity_table_ui"),
                                                style = "height: 200px;" # Set the same fixed height here
                                                
                                              )
                                       ),
                                       # column(4, 
                                       #        verbatimTextOutput("droppedBMSIDs")
                                       # ),
                                       column(12,
                                              DTOutput("display_table")
                                       )
                                     )
                            ),
                            tabPanel("F1 Verification",
                                     
                                     fluidRow(
                                       column(4,
                                              wellPanel(style = "height: 280px;",
                                                        h4("Adjust F1 Heterozygosity Percentage Threshold:"),
                                                        
                                                        sliderInput("percentage_he_threshold", "Select Threshold:", min = 0, max = 100, value = 60))),
                                       column(8,
                                              wellPanel(style = "height: 280px;",
                                                        
                                                        uiOutput("decision_table"))),
                                       column(12,       
                                              DTOutput("F1_verification"))
                                       
                                     ),
                            ),
                            
                 ),
                 
                 navbarMenu("Help", 
                            tabPanel("About Us",
                                     fluidPage(
                                       titlePanel("About Us"),
                                       # Update styles to remove gray color and use a neutral background
                                       tags$style(HTML("
                        .about-us-container {
                          background-color: #ffffff; /* White background for sections */
                          padding: 20px;
                          border-radius: 10px;
                          box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
                          min-height: 300px; /* Ensure sections have consistent height */
                        }
                        .about-us-header {
                          background-color: #003366; /* Matching navigation panel color */
                          color: white;
                          padding: 15px;
                          border-radius: 10px;
                          text-align: center;
                          box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2);
                        }
                        .footer-info {
                          background-color: #003366; /* Matching navigation panel color */
                          color: white;
                          padding: 15px;
                          border-radius: 10px;
                          text-align: center;
                          margin-top: 20px;
                        }
                        .logo-img {
                          height: 60px;
                          vertical-align: middle;
                        }
                        .team-icon, .contact-icon {
                          font-size: 30px;
                          color: #003366; /* Matching icon color */
                        }
                        .contact-info, .team-info {
                          color: #333;
                        }
                      ")),
                                       fluidRow(
                                         # Header with matching color background
                                         column(12,
                                                div(class = "about-us-header",
                                                    img(src = "logoW.png", class = "logo-img"),
                                                    br(),
                                                )
                                         )
                                       ),
                                       fluidRow(
                                         # Column for Mission and Team
                                         column(6,
                                                div(class = "about-us-container",
                                                    h3(class = "team-icon", icon("flag-checkered"), "Our Mission"),
                                                    p(class = "team-info", "Our mission is to provide researchers with robust tools for molecular marker analysis and quality control in genetic studies. We aim to facilitate accurate genetic evaluations and enhance breeding programs."),
                                                    br(),
                                                    h3(class = "team-icon", icon("users"), "Our Development Team"),
                                                    tags$ul(
                                                      tags$li("Dr. Abhishek Rathore - Principal Scientist, Breeding Data and Informatics Expert"),
                                                      tags$li("Mr. Peter Kimathi - Bioinformatics and Software Developer Specialist"),
                                                      tags$li("Dr. Roma Rani Das - Biometrician")
                                                    )
                                                )
                                         ),
                                         # Column for Contact and Acknowledgements
                                         column(6,
                                                div(class = "about-us-container",
                                                    h3(class = "contact-icon", icon("envelope"), "Contact Us"),
                                                    p(class = "contact-info", "For inquiries, feedback, or support, please contact us at:"),
                                                    tags$ul(
                                                      tags$li(a(href = "mailto:abhishek.rathore@cgiar.org", "Email: abhishek.rathore@cgiar.org")),
                                                      tags$li("Phone: +91 9676 403 413"),
                                                      tags$li("Address: CIMMYT, Nairobi")
                                                    ),
                                                    br(),
                                                    h3(class = "contact-icon", icon("handshake"), "Acknowledgements"),
                                                    p(class = "contact-info", "We thank the CIMMYT and NARS partners for their support and collaboration.")
                                                )
                                         )
                                       ),
                                       # Add some spacing and a footer with matching color
                                       div(class = "footer-info",
                                       )
                                     )
                            )
                 )
                 
                 
                 
                 
               ))
server = function(input, output){
  
  
  #RV to store paths
  rv <- reactiveValues(
    tempFolderOut = NULL,
    MetaFilePath = NULL,
    GenoFilePath = NULL,
    ms_file = NULL,
    tempFolder=NULL
  )
  
  # Creating a session-specific temp folder with an "Out" subfolder
  createtempFolder <- reactive({
    tempFolder <- tempfile("")
    #setwd(tempdir())
    dir.create(tempFolder, recursive = TRUE)
    tempFolderOut <- paste(tempFolder, "/Out/", sep = "")
    dir.create(tempFolderOut, recursive=TRUE)
    #rv$tempFolderOut <- tempFolder
    return(tempFolder)
  })
  Mdata <- reactiveVal()
  Gdata <- reactiveVal()
  observeEvent(list(input$file1, input$file2), {
    req(input$file1, input$file2)
    tempFolder=createtempFolder()
    #Set working directory
    rv$tempFolder<-tempFolder
    setwd(rv$tempFolder)
    rv$tempFolderOut <- paste(tempFolder, "/Out/", sep = "")
    rv$MetaFilePath <- tempfile("META_", tmpdir = tempFolder, fileext = ".txt")
    rv$GenoFilePath <- tempfile("Geno_", tmpdir = tempFolder, fileext = ".csv")
    rv$ms_file <- file.path(tempFolder, "Out", "F1_marker_summary.txt")
    
    inFile1 <- input$file1
    inFile2 <- input$file2
    
    if (!is.null(inFile1)) {
      # Read the uploaded file
      file.copy(inFile1$datapath,rv$MetaFilePath)
      fileData <- read.delim(inFile1$datapath, header = TRUE, sep = "\t") # Adjust the 'sep' as necessary
      Mdata(fileData) # Update the reactive variable
      ################################################
      output$contents_file1 <- renderDataTable({
        dd <- req(Mdata()) # Ensure data is available before rendering
        
        names <- c('SAMPLE_UID','DESIGNATION','SAMPLE_NAME','PARENT_A','PARENT_B')
        dd[,names] <- lapply(dd[,names] , factor)
        datatable(dd, filter = "top", escape = FALSE,
                  options = list(dom = 'Blfrtip', autoWidth = TRUE, scrollX = TRUE, pageLength = 11,
                                 lengthMenu = list(c(11, 15, -1), c('11', '15', 'All'))  # Customize page length options
                                 
                  )) %>%
          formatStyle(
            columns = c('SAMPLE_UID','SAMPLE_NAME','DESIGNATION', 'PARENT_A', 'PARENT_B'),
            `white-space` = 'nowrap')
      })
      ####################################################
    }
    if (!is.null(inFile2)) {
      # Read the uploaded file
      file.copy(inFile2$datapath, rv$GenoFilePath)
      fileData2 <- read.csv(inFile2$datapath, nrows = 12, header = F) # Adjust the 'sep' as necessary
      Gdata(fileData2) # Update the reactive variable
    }
  })
  
  ##################################################################################################################
  #################################################################################################################
  observeEvent(list(input$load_example_meta,input$load_example_geno), {
    req(input$load_example_meta, input$load_example_geno)
    tempFolder=createtempFolder()
    rv$tempFolder<-tempFolder
    #setwd(rv$tempFolder)
    rv$tempFolderOut <- paste(tempFolder, "/Out/", sep = "")
    rv$MetaFilePath <- tempfile("META_", tmpdir = tempFolder, fileext = ".txt")
    rv$GenoFilePath <- tempfile("Geno_", tmpdir = tempFolder, fileext = ".csv")
    #rv$ms_file <- file.path(tempFolder, "Out", "F1_marker_summary.txt")
    
    # Path to the example file in the 'www' folder
    exampleFilePath <- "www/meta.txt" # Adjust the path/name as necessary
    file.copy(exampleFilePath,rv$MetaFilePath)
    # Read the example file
    exampleData <- read.delim(exampleFilePath, header = TRUE, sep = "\t") # Adjust the 'sep' as necessary
    Mdata(exampleData) # Update the reactive variable
    
    
    output$contents_file1 <- renderDataTable({
      dd <- req(Mdata()) # Ensure data is available before rendering
      
      names <- c('SAMPLE_UID','DESIGNATION','SAMPLE_NAME','PARENT_A','PARENT_B')
      dd[,names] <- lapply(dd[,names] , factor)
      datatable(dd, filter = "top", escape = FALSE,
                options = list(dom = 'Blfrtip', autoWidth = TRUE, scrollX = TRUE, pageLength = 11,
                               lengthMenu = list(c(11, 15, -1), c('11',  '15', 'All'))  # Customize page length options
                               
                )) %>%
        formatStyle(
          columns = c('SAMPLE_UID','SAMPLE_NAME','DESIGNATION', 'PARENT_A', 'PARENT_B'),
          `white-space` = 'nowrap')
    })
    # Path to the example file in the 'www' folder
    exampleFilePath <- "www/geno.csv" # Adjust the path/name as necessary
    file.copy(exampleFilePath, rv$GenoFilePath)   
    # Read the example file
    exampleData <- read.csv(exampleFilePath, nrows = 12, header = F) # Adjust the 'sep' as necessary
    Gdata(exampleData) # Update the reactive variable
  })
  ##################################################################################################################
  ###################################################################################################################
  # Reading intertek data 
  Gdata1 <- reactiveVal()
  
  observe({
    req(input$file2)
    inFile <- input$file2
    
    if (!is.null(inFile)) {
      # Read the uploaded file
      fileData <- read.csv(inFile$datapath, skip = 14) # Adjust the 'sep' as necessary
      Gdata1(fileData) # Update the reactive variable
    }
  })
  
  
  observeEvent(input$load_example_geno, {
    # Path to the example file in the 'www' folder
    exampleFilePath <- "www/geno.csv" # Adjust the path/name as necessary
    
    # Read the example file
    exampleData <- read.csv(exampleFilePath, skip = 14) # Adjust the 'sep' as necessary
    Gdata1(exampleData) # Update the reactive variable
  })
  
  
  
  # Output the first 12 rows
  output$first12 <- render_gt({
    df <- Gdata()
    if (is.null(df)) {
      return()
    }
    gt(df)
    
  })
  
  # Output the remaining rows
  output$remaining <- renderDT({
    df <- Gdata1()
    
    if("Plate" %in% colnames(df)) {
      datatable(df, options = list(pageLength = 11, scrollX = TRUE)) %>%
        formatStyle(
          columns = c('Plate'),
          `white-space` = 'nowrap'
        )
    } else {
      # Optionally, return a message or an empty data table if the column is not found
      datatable(data.frame(), options = list(pageLength = 11, scrollX = TRUE))
    }
  })
  
  
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
  
  
  
  
  file_statuses <- reactiveValues(file1_uploaded = FALSE, file2_uploaded = FALSE, example1_uploaded = FALSE, example2_uploaded =FALSE)
  
  observeEvent(input$file1, {
    file_statuses$file1_uploaded <- TRUE
  })
  
  observeEvent(input$file2, {
    file_statuses$file2_uploaded <- TRUE
  })
  
  
  observeEvent(input$load_example_meta, {
    file_statuses$example1_uploaded <- TRUE
  })
  
  observeEvent(input$load_example_geno, {
    file_statuses$example2_uploaded <- TRUE
  })
  # Generate the Analyze button UI based on the status of file uploads
  output$analyze_button_ui <- renderUI({
    if(file_statuses$file1_uploaded && file_statuses$file2_uploaded) {
      actionButton("analyze_button", "Run Marker Summary")
    } else
      if(file_statuses$example1_uploaded && file_statuses$example2_uploaded) {
        actionButton("analyze_button", "Run Marker Summary")
      }
  })
  
  output$Visualize_button_ui <- renderUI({
    if(file_statuses$file1_uploaded && file_statuses$file2_uploaded) {
      tagList(
        actionButton("Visualize_button", "Get SNP Data Visualization",class = "visualize-button"),
      )
      
    } else
      
      if(file_statuses$example1_uploaded && file_statuses$example2_uploaded) {
        tagList(
          actionButton("Visualize_button", "Get SNP Data Visualization", class = "visualize-button"),
        )
        
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
  
  observeEvent(input$analyze_button2, {
    # Render the download button UI after the analyze button is clicked
    output$download_button_ui <- renderUI({
      tagList(
        div(
          style = "display: flex; align-items: center;", # Flexbox to align items in a row
          downloadButton("download_snp_data", "Download Raw Data in HAPMAP Format", class = "download-button"),
          div(
            style = "margin-left: 10px; font-weight: bold; color: blue;", # Spacing and styling for the message
            "Please go to the Result Tab to view analysed outputs."
          )
        )
      )
    })
  })
  
  
  output$analysis_message_ui <- renderUI({
    div(
      style = "margin-top: 10px; font-weight: bold; color: blue;",
      "Please go to the Result Tab to view your analysis."
    )
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

  #print(plotdd())
  plotdd <- reactive({
  #req(input$file1, input$file2)
  req( (!is.null(input$file1) && !is.null(input$file2)) ||(input$load_example_meta > 0 && input$load_example_geno > 0))

  if(file_statuses$file1_uploaded && file_statuses$file2_uploaded) {
  tempData <- read_datasets(input$file1$datapath, input$file2$datapath)$lgc_data
  #print(input$file1$datapath)
  #print(input$file2$datapath)
  }
  if(file_statuses$example1_uploaded && file_statuses$example2_uploaded){
  #print(rv$MetaFilePath)
  #print(rv$GenoFilePath)
  tempData <- read_datasets(rv$MetaFilePath,rv$GenoFilePath)$lgc_data
  }
  tm <- tempData
  print(rv$tempFolderOut)
  write.table(tm, paste(rv$tempFolderOut, 'plot_file.csv', sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
  tempData <- read_csv(paste(rv$tempFolderOut, "plot_file.csv", sep = ""), show_col_types = FALSE)
  return(tempData)  # Return the dataframe from the reactive expression
   })

  observeEvent(input$Visualize_button, {
    shinyjs::show("processing")  # Show the processing message
    req(plotdd())
    output$plot1 <- renderPlotly({
      Sys.sleep(5)  # Assume it takes 5 seconds to generate the plot, adjust based on your needs

      p1 <- ggplot(plotdd(), aes(X, Y, col = Call, text = paste("MasterPlate:", MasterPlate,
                                                                "<br>MasterWell:", MasterWell,
                                                                "<br>SubjectID:", SubjectID
      )))  + 
        geom_point(alpha = 0.4) + 
        theme_bw() +
        facet_wrap(~SNPID, ncol = 4, scales = "fixed")


      ggplotly(p1)  

    })
    
    shinyjs::delay(100, shinyjs::hide("processing"))
  })


  ##########################################################
  grid_df <- reactive({
  req( (!is.null(input$file1) && !is.null(input$file2)) ||(input$load_example_meta > 0 && input$load_example_geno > 0))

  if(file_statuses$file1_uploaded && file_statuses$file2_uploaded) {
  datasets <- read_datasets(input$file1$datapath, input$file2$datapath)
  }
  if(file_statuses$example1_uploaded && file_statuses$example2_uploaded){
  datasets <- read_datasets(rv$MetaFilePath,rv$GenoFilePath)
  }
  #req(input$file1, input$file2)
  #datasets<-read_datasets(input$file1$datapath,input$file2$datapath)
  tmpdata<-create_grid_file(datasets$lgc_data,datasets$metadata)
  #print(tmpdata)
  return(tmpdata)  # Return the dataframe from the reactive expression
   })


  ms_file <- reactive({
    req(grid_df())
    ##Get SNP summary
    snp_summary <- lapply(grid_df()[, 7:ncol(grid_df())], function(col) calculate_snp_info(col, nrow(grid_df())))
    #snp_summary <- lapply(grid_df[, 7:ncol(grid_df)], function(col) round(calculate_snp_info(col, nrow(grid_df)), 2))
    # Combine results into a single dataframe
    snp_summary_df <- do.call(rbind, snp_summary)
    snp_summary_df$marker_name <- rownames(snp_summary_df)
    #Now drop the index column
    rownames(snp_summary_df) <- NULL 
    ##Re-order to make "marker_name" to appear in first column
    snp_summary_df <- snp_summary_df[, c("marker_name", setdiff(names(snp_summary_df), "marker_name"))]
    write.table(snp_summary_df, paste0(rv$tempFolderOut,"/F1_marker_summary.txt"), sep = "\t", row.names = FALSE,col.names = TRUE,quote = FALSE)
    ms_file <- paste(rv$tempFolderOut, "/F1_marker_summary.txt", sep="")
    return(ms_file)
    #dd(read.table(rv$ms_file, header = TRUE, stringsAsFactors = FALSE))

  })

  observeEvent(input$Visualize_button, {
    req(ms_file())
    #print(ms_file)
    dd<-read.table(ms_file(), header = TRUE, stringsAsFactors = FALSE)
    #print("YEs")
    #print(dd)
    # Correctly updating dd with new values
    updated_dd <- dd %>%
      mutate(exclude_marker = '') %>%
      select(marker_name, exclude_marker, everything())
    
    dd <- updated_dd
     #dd(updated_dd)  # Update the reactiveVal with the new data frame
     output$result_ms_table <- renderDT({
       datatable(dd, escape = FALSE,  selection = 'none', extensions = 'Buttons',filter  = "top",options = list( pageLength = 10, dom = 'Blfrtip', # Set initial page length
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
  #   showProgress()
  })
  outputOptions(output, "showProgress", suspendWhenHidden = FALSE)

 

  puritySummary <- reactiveVal()
  observeEvent(input$analyze_button2, {
    req(grid_df())
    #req(input$file1, input$file2)
    # Assuming ms_file is defined somewhere in your app and accessible here
    parent_consensus<-create_parents_consensus(grid_df())$parent_cons_df
    purity_df<-parental_purity(parent_consensus)
    #print("Yes")
    #print(head(purity_df))
    write.table(purity_df, paste0(rv$tempFolderOut,"/F1_PuritySummary.txt"), sep = "\t", row.names = FALSE,col.names = TRUE,quote = FALSE)
    puritySummary(read.delim(paste(rv$tempFolderOut, "F1_PuritySummary.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE))
  })
  

  # Display the heatmap of Summary dataset
  output$display_plot <- renderPlot({
    puritysmry <- req(puritySummary())
    pos_val <- unlist(gregexpr(":", puritysmry$SAMPLE_NAME))
    Geno_name <- substr(puritysmry$SAMPLE_NAME,1,pos_val-1)
    Rep_name <- substr(puritysmry$SAMPLE_NAME,pos_val+1, nchar(puritysmry$SAMPLE_NAME))
    #print(head(Geno_name))
    #print(Rep_name)
    puritysmry <- data.frame(Entry = Geno_name,
                             Rep = as.numeric(Rep_name),
                             puritysmry[,-1])
    
    #print(head(puritysmry))
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
    print(ptbl)
  })

  

  #Parent_consensus <- reactiveVal()
  Parent_consensus<-reactive({
    req( (!is.null(input$file1) && !is.null(input$file2)) ||(input$load_example_meta > 0 && input$load_example_geno > 0))
    req(grid_df())
    if(file_statuses$file1_uploaded && file_statuses$file2_uploaded) {
    datasets <- read_datasets(input$file1$datapath, input$file2$datapath)
    }
    if(file_statuses$example1_uploaded && file_statuses$example2_uploaded){
    datasets <- read_datasets(rv$MetaFilePath,rv$GenoFilePath)
    }


    #req(input$file1, input$file2)
    
    excluded <- read.table(paste(rv$tempFolderOut,"excluded_markers.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE)
    exc_markers<-excluded[[1]]
    #grid_df() <- grid_df()[, setdiff(names(grid_df()), exc_markers)]
    grid_df<-grid_df()
    grid_df <- grid_df[, setdiff(names(grid_df), exc_markers)]
    print(names(grid_df))
    parent_consensus<-create_parents_consensus(grid_df)$parent_cons_df
    consensus_df<-create_parents_consensus(grid_df)$cons_df
    Purity_df<-parental_purity(parent_consensus)
    #Replace the consensus designations with their actual names
    parent_consensus$SAMPLE_NAME[parent_consensus$SAMPLE_NAME == "Consensus"] <- parent_consensus$DESIGNATION[parent_consensus$SAMPLE_NAME == "Consensus"]
    #print(Purity_df)
    temp1<-as.data.frame(parent_consensus[,-1])
    #print(typeof(temp1))
    #Rename Sample_Name to "DNA\Assay"
    colnames(temp1)[1] <- "DNA...Assay"
    #unlist columns
    unlist_columns <- function(df) {
    df[] <- lapply(df, function(x) if (is.list(x)) unlist(x) else x)
    return(df)
    }

    # Unlist the columns if necessary
    temp1 <- unlist_columns(temp1)
    pr_df<-unlist_columns(Purity_df)
    cn_df<-unlist_columns(consensus_df)
    #datasets<-read_datasets(input$file1$datapath,input$file2$datapath)
    cons_report<-create_purity_report(pr_df,cn_df)
    het_threshold<-60
    purity_threshold<-80
    f1_grid<-create_f1(grid_df)
    f1_report_output<-create_f1_report(cons_report,f1_grid, het_threshold, purity_threshold)
    #print(f1_report_output)
    write.table(temp1, paste0(rv$tempFolderOut,"/F1_consensus_grid.csv", sep = ""), row.names = FALSE,col.names = TRUE,quote = FALSE)
    write.table(pr_df, paste0(rv$tempFolderOut,"/F1_PuritySummary.txt",sep=""), sep = "\t", row.names = FALSE,col.names = TRUE,quote = FALSE)
    write.table(cons_report, paste0(rv$tempFolderOut,"/F1_PurityReport.txt",sep=""), sep = "\t", row.names = FALSE,col.names = TRUE,quote = FALSE)
    write.table(f1_report_output, paste0(rv$tempFolderOut,"/F1_F1_summary.txt",sep=""), sep = "\t", row.names = FALSE,col.names = TRUE,quote = FALSE)

    #datasets<-read_datasets(file_in1,file_in2) 
    #subset lgc_data
    d1 <- datasets$lgc_data[, c("MasterPlate", "MasterWell", "SubjectID")]
    #Rename columns
    colnames(d1) <- c("PLATE_ID", "WELL", "SAMPLE_UID")
  
    d2<-datasets$metadata
    colnames(d2)<-c("SUBJECT_ID","SAMPLE_NAME","DESIGNATION","PARENT_A","PARENT_B","FN" )
    final_dataset <- merge(d1, grid_df, by = "SAMPLE_UID", all = TRUE)
    colnames(final_dataset)[colnames(final_dataset) == "SAMPLE_UID"] <- "SUBJECT_ID"
    colnames(final_dataset)[colnames(final_dataset) == "GENERATION"] <- "FN"
    
    final_dataset1 <- final_dataset %>%relocate(PLATE_ID, WELL, SUBJECT_ID, SAMPLE_NAME, DESIGNATION, PARENT_A, PARENT_B, FN)
 
    
    
    ##F1_QAQC
    final_dataset_F1 <- final_dataset1[final_dataset1$FN == "F1", ]
    f1_report_file <- f1_report_output[, c("F1_one_Name", "Parent_A_Score", "Parent_B_Score", 
                    "Call_Hetero", "Call_Parent_A", "Call_Parent_B", 
                    "Call_Missing", "Percentage_Het", "comment")]
    colnames(f1_report_file)[colnames(f1_report_file) == "F1_one_Name"] <- "SAMPLE_NAME"
    final_f1_report <- merge(final_dataset_F1, f1_report_file, by = "SAMPLE_NAME", all = TRUE)
    #print(names(final_f1_report))
    write.table(final_f1_report, paste0(rv$tempFolderOut,"/F1_QAQC.txt",sep=""), sep = ",", row.names = FALSE,col.names = TRUE,quote = FALSE)

    ##Parent QAQC
    final_dataset_Parent <- final_dataset1[final_dataset1$FN == "Parent", ]
    final_parent_purity_report <- merge(final_dataset_Parent, Purity_df, by = "SAMPLE_NAME", all = TRUE)
    #print(names(final_parent_purity_report))
    write.table(final_parent_purity_report, paste0(rv$tempFolderOut,"/Parent_QAQC.txt",sep=""), sep = ",", row.names = FALSE,col.names = TRUE,quote = FALSE)

    return(temp1)
    #Parent_consensus()(read.csv(paste(rv$tempFolderOut, "F1_consensus_grid.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE))

  })
  hapmap<-reactive({
    #req(input$file1, input$file2)
    req( (!is.null(input$file1) && !is.null(input$file2)) ||(input$load_example_meta > 0 && input$load_example_geno > 0))
    if(file_statuses$file1_uploaded && file_statuses$file2_uploaded) {
    hapmap <- create_hapmap(read_datasets(input$file1$datapath, input$file2$datapath))
    }
    if(file_statuses$example1_uploaded && file_statuses$example2_uploaded){
    hapmap <- create_hapmap(read_datasets(rv$MetaFilePath,rv$GenoFilePath))
    }
    #hapmap<-create_hapmap(read_datasets(input$file1$datapath, input$file2$datapath))
    write.table(hapmap, paste0(rv$tempFolderOut,"/hapmap.hmp.txt",sep=""), sep = "\t", row.names = FALSE,col.names = TRUE,quote = FALSE)
    })

  observeEvent(input$analyze_button2, {
    req(Parent_consensus())
    req(hapmap())
    #print(Parent_consensus())
    #Parent_consensus(read.csv(paste(tempFolderOut, "F1_consensus_grid.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE))
    read.csv(paste(rv$tempFolderOut, "F1_consensus_grid.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    
  })
  

  # Parent_consensus <- read.csv(paste(tempFolderOut, "F1_consensus_grid.csv", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  output$consensus <- renderDT({
    req(Parent_consensus())
    parentConsensus <- Parent_consensus()
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
    datatable(parentConsensus, 
              extensions = 'FixedColumns', 
              filter = "top",
              options = list(
                fixedColumns = list(leftColumns = 3),
                dom = 'Blfrtip',
                buttons = list(
                  list(extend = 'copy', filename = 'parent consensus'),
                  list(extend = 'csv', filename = 'parent consensus'),
                  list(extend = 'excel', filename = 'parent consensus'),
                  list(extend = 'pdf', filename = 'parent consensus')
                ),  
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
              )
    )%>%
      formatStyle(
        'Type',
        backgroundColor = styleEqual(c("Consensus", "Sample"), c("#b2dfee", "#dedede")))
  })


  #fs <- reactiveValues(F1SummaryData = NULL)
  fs <- reactive({
  req( (!is.null(input$file1) && !is.null(input$file2)) ||(input$load_example_meta > 0 && input$load_example_geno > 0))
  req(grid_df())
  #if(file_statuses$file1_uploaded && file_statuses$file2_uploaded) {
  # datasets <- read_datasets(input$file1$datapath, input$file2$datapath)
  # }
  #if(file_statuses$example1_uploaded && file_statuses$example2_uploaded){
  # datasets <- read_datasets(rv$MetaFilePath,rv$GenoFilePath)
  # }
  #req(input$file1, input$file2)
  #req(grid_df())
  #req(hapmap())
  req(Parent_consensus())
  req(paste(rv$tempFolderOut, "/F1_QAQC.txt", sep = ""))
  #tempData <- read.delim(paste(rv$tempFolderOut, "/F1_QAQC.txt", sep = ""), sep = ",", header = TRUE, stringsAsFactors = FALSE, na.strings = "na")
  tempData <-read.csv(paste(rv$tempFolderOut, "F1_QAQC.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  tempData <- tempData %>%relocate(PLATE_ID, WELL, SUBJECT_ID, SAMPLE_NAME, DESIGNATION, PARENT_A, PARENT_B, FN)%>%
                  rename(BMS_ID = SAMPLE_NAME)
  return(tempData)
  })
  parent_qc_threshold <- reactive({
    input$purity_threshold
  })
  
  percentage_he_threshold <- reactive({
    input$percentage_he_threshold
  })
  fs3 <- reactiveValues(F1SummaryData = NULL)

  fs2 <- reactive({
  req(fs())
  tempData <-fs()
  #print(names(tempData))
  # Update the Comment field based on the selected thresholds
  tempData <- tempData %>%
     mutate(
       Status = case_when(
         is.na(Parent_A_Score) | is.na(Parent_B_Score) ~ "Failed F1",
         Percentage_Het == 100 & Parent_A_Score > parent_qc_threshold() & Parent_B_Score > parent_qc_threshold() ~ "Successful F1",
         (Parent_A_Score <= parent_qc_threshold() | Parent_B_Score <= parent_qc_threshold()) ~ "Parent Quality Failed",
         Percentage_Het >= percentage_he_threshold() & Percentage_Het < 100 & Parent_A_Score > parent_qc_threshold() & Parent_B_Score > parent_qc_threshold() ~ "Possible F1",
         Percentage_Het < percentage_he_threshold() ~ "Failed F1",
        TRUE ~ "Unknown"
       )
     )
    
  tempData <- tempData[, !(names(tempData) %in% c("Comment"))]
  tempData$`F1 Selection` <- ifelse(tempData$Status %in% c("Possible F1", "Successful F1"), "Select", "No Selection") 
  return (tempData)
  })

  observeEvent(input$analyze_button2, {
    req(fs2())
    fs3$F1SummaryData <- fs2()
  })
 

  output$F1_verification <- renderDT({
    req(fs())
    #req(fs2())
    req(fs3$F1SummaryData)
    #req(hapmap())
    char_columns <- sapply(fs3$F1SummaryData, is.character)
    #print(char_columns)
    fs3$F1SummaryData[, char_columns] <- lapply(fs3$F1SummaryData[, char_columns], as.factor)
    fs3$F1SummaryData$`F1 Selection` <- as.factor(fs3$F1SummaryData$`F1 Selection`)
    plantSelectionColIndex <- which(names(fs3$F1SummaryData) == "F1 Selection")
    fs3$F1SummaryData <- fs3$F1SummaryData[, !(names(fs3$F1SummaryData) %in% c("Comment"))]
    #print(names(fs3$F1SummaryData))
    
    # Find the SNP columns
    snp_columns_indices <- grep("^snp", names(fs3$F1SummaryData))
    first_snp_column_index <- min(snp_columns_indices)
    last_snp_column_index <- max(snp_columns_indices)
    
    # Render the DataTable with the updated comments
    datatable(fs3$F1SummaryData, 
              filter = "top", 
              extensions = "FixedColumns",
              options = list(dom = 'Blfrtip',  
                             fixedColumns = list(leftColumns = 9), 
                             scrollX = TRUE, 
                             buttons = list(
                               list(extend = 'copy', filename = 'F1 verification'),
                               list(extend = 'csv', filename = 'F1 verification'),
                               list(extend = 'excel', filename = 'F1 verification'),
                               list(extend = 'pdf', filename = 'F1 verification')
                             ),                             
                             lengthMenu = list(c(15, 30, 50, -1), c(15, 30, 50, "All")),
                             columnDefs = list(
                               list(targets = plantSelectionColIndex - 1, className = 'editable'),
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
              editable = list(target = 'cell', disable = list(columns = c(0:ncol(fs3$F1SummaryData)-1)))) %>%
      formatStyle('Status', 
                  target = 'cell', 
                  backgroundColor = styleEqual(c("Successful F1", "Possible F1", "Failed F1", "Parent Quality Failed"), c("#4CAF50", "#8FD744", "#FF8080", "#FFA500")), 
                  color = "white") %>%
      formatStyle(
        columns = c('PLATE_ID', 'WELL', 'SUBJECT_ID', 'BMS_ID', 'DESIGNATION',  "F1 Selection", "PARENT_A", "PARENT_B"),
        `white-space` = 'nowrap'
      )
      #}
  })

  observeEvent(input$F1_verification_cell_edit, {
    info <- input$F1_verification_cell_edit
    str(info) # For debugging
    # Update the reactiveValues data
    fs3$F1SummaryData[info$row, info$col] <- DT::coerceValue(info$value, fs3$F1SummaryData[info$row, info$col])
  })
  #}

  ParentPuritySummaryData <- reactive({
    req(input$purity_threshold) # Ensure the slider input is available
    # Read the data
    tempData <-read.csv(paste(rv$tempFolderOut, "Parent_QAQC.txt", sep = ""), header = TRUE, stringsAsFactors = FALSE)
    tempData$Status <- ifelse(tempData$PurityScore > input$purity_threshold, "Pure Line", "Need to Purify Parent")
    tempData$`Parent Selection` <- ifelse(tempData$Status %in% c("Pure Line"), "Select", "No Selection")
    tempData <- tempData %>%relocate(PLATE_ID, WELL, SUBJECT_ID, SAMPLE_NAME, DESIGNATION, PARENT_A, PARENT_B, FN)%>%
                  rename(BMS_ID = SAMPLE_NAME)%>%select(-PARENT_A, -PARENT_B)
    
    return(tempData) # Return the updated data
  })
    

  
  output$display_table <- renderDT({
    req(ParentPuritySummaryData())
    data <- ParentPuritySummaryData()
    #print(head(data))
    
    columns_to_exclude <- c("TotalMarkers" , "Match" , "Mismatch","BothMissing","ConsensusMissing", "SampleMissing","PurityScore")
    data <- data %>% mutate(across(-all_of(columns_to_exclude), factor))
    
    unique_designations <- unique(data$DESIGNATION)
    
    # Generate colors
    colors <- grDevices::rainbow(length(unique_designations), s = 0.6, v = 0.85)
    
    # Convert colors to RGBA
    colors_rgba <- sapply(colors, function(col) {
      col_rgb <- col2rgb(col, alpha = TRUE)
      sprintf("rgba(%d,%d,%d,%f)", col_rgb[1,], col_rgb[2,], col_rgb[3,], col_rgb[4,]/255)
    })
    colors_map <- setNames(colors_rgba, unique_designations)
    
    snp_columns_indices <- grep("^snp", names(data))
    plantSelectionColIndex <- which(names(data) == "Parent Selection")
    data <- data[, !(names(data) %in% c("Comment"))]
    # Render the DataTable
    datatable(data, 
              filter = "top", extensions = 'FixedColumns',
              options = list(dom = 'Blfrtip', scrollX = TRUE, fixedColumns = list(leftColumns = 7), 
                             buttons = list(
                               list(extend = 'copy', filename = 'parent purity'),
                               list(extend = 'csv', filename = 'parent purity'),
                               list(extend = 'excel', filename = 'parent purity'),
                               list(extend = 'pdf', filename = 'parent purity')
                             ),                             lengthMenu = list(c(15, 30, 50, -1), c(15, 30, 50, "All")),
                             columnDefs = list(list(targets = plantSelectionColIndex - 1, className = 'editable'),
                                               list(
                                                 targets = snp_columns_indices, 
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
              editable = list(target = 'cell', disable = list(columns = c(0:ncol(data)-1)))) %>%
      formatStyle('DESIGNATION', 
                  backgroundColor = styleEqual(unique_designations, colors)) %>%
      formatStyle(
        columns = c('PLATE_ID', 'WELL', 'SUBJECT_ID', 'BMS_ID', 'DESIGNATION',  "Parent Selection"),
        `white-space` = 'nowrap')%>%
      formatStyle(
        'Status',
        backgroundColor = styleEqual(
          c("Pure Line", "Need to Purify Parent"), 
          c("#4CAF50", "#FFA500")  # Green for "Pure Line", Orange for "Need to Purify Parent"
        ),
        color = "white"
      )
  })
  
  
  
  # Handle cell edits in the Parent Purity DataTable
  observeEvent(input$display_table_cell_edit, {
    info <- input$display_table_cell_edit
    str(info) # For debugging
    # Update the reactiveValues data
    rv$ParentPuritySummaryData[info$row, info$col] <- DT::coerceValue(info$value, rv$ParentPuritySummaryData[info$row, info$col])
  })
  
  
  # Reactive expression to fetch selected marker names based on indices
  selected_marker_names <- reactive({
    #req(dd())
    req(ms_file())
    unlist_columns <- function(df) {
    df[] <- lapply(df, function(x) if (is.list(x)) unlist(x) else x)
    return(df)
    }
    d1<-read.table(ms_file(), header = TRUE, stringsAsFactors = FALSE)
    dd<-unlist_columns(d1)
    #print(head(dd))
    # Ensure input$selectedMarkers is not NULL and has at least one selection
    if (!is.null(input$selectedMarkers) && length(input$selectedMarkers) > 0) {
      selected_indices <- as.numeric(input$selectedMarkers)  # Convert indices to numeric
      marker_names <- dd[selected_indices, "marker_name", drop = TRUE]  # Extract marker names
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
        paste("List of Markers Excluded for Further Analysis:", paste(marker_names, collapse = ", "))
      }
    }
  })
  
  output$purity_table_ui <- renderUI({
    threshold <- input$purity_threshold
    
    tags$div(
      tags$h4("Parent Purity Selection Decision:"),
      tags$table(
        class = "fancy-table table table-bordered", 
        tags$thead(
          tags$tr(
            tags$th("Status"),
            tags$th("Criteria")
            
          )
        ),
        tags$tbody(
          tags$tr(
            tags$td(style = "background-color: #4CAF50; color: white;", "Pure Line"),
            tags$td(paste("Parent Purity Score >", threshold))
            # Green
          ),
          tags$tr(
            tags$td(style = "background-color:  #FFA500; color: white;", "Need to Purify Parent"), # Orange
            tags$td(paste("Parent Purity Score <=", threshold))
            
          )
        )
      )
    )
  })
  
  output$decision_table <- renderUI({
    tags$div(
      tags$h4("F1 Cross Selection Decision:"),
      tags$table(
        class = "fancy-table table table-bordered", 
        style = "width: 60%;",  # Adjust the percentage or use px (pixels) to fix the width
        
        tags$thead(
          tags$tr(
            # tags$th("Parent Purity Score"),
            tags$th("Status"),
            tags$th("Criteria")
          )
        ),
        tags$tbody(
          tags$tr(
            # tags$td(paste("Both Parents >", parent_qc_threshold())),
            tags$td(style = "background-color: #4CAF50; color: white;", "Successful F1"), # Green
            tags$td(paste("F1 Heterozygosity = 100 %"))
            
          ),
          tags$tr(
            # tags$td(paste("Both Parents >", parent_qc_threshold())),
            tags$td(style = "background-color: #8FD744; color: white;", "Possible F1"),
            tags$td(paste("F1 Heterozygosity >=", percentage_he_threshold(),"%", "and F1 Heterozygosity < 100 %"))
            # Orange
          ),
          
          tags$tr(
            # tags$td("Either or Both Parents"),
            tags$td(style = "background-color: #FF8080; color: white;", "Failed F1"),
            tags$td(paste("Either Parent Purity Score is missing or F1 Heterozygosity <", percentage_he_threshold() ,"%"))
            # Red
          ),
          tags$tr(
            # tags$td(paste("Either Parent <=", parent_qc_threshold())),
            # tags$td(paste("F1 Heterozygosity >=", percentage_he_threshold(),"%", "and F1 Heterozygosity <= 100 %")),
            
            tags$td(style = "background-color:  #FFA500; color: white;", "Parent Quality Failed"),
            tags$td(paste("Either Parent Purity Score <=", parent_qc_threshold(), "%"))
            # Purple
          )
        )
      )
    )
  })
  #req(hapmap())
  output$download_snp_data <- downloadHandler(
    filename = function() {
      paste("hapmap-", Sys.Date(), ".hmp.txt", sep = "")  # Generates the filename dynamically based on the current date
    },
    content = function(file) {
      file.copy(file.path(rv$tempFolderOut, "hapmap.hmp.txt"), file)  # Copy the file from the output directory to the download path
    },
    contentType = "text/plain"  # This line makes sure that the browser handles the file as a plain text file
  )
 
 observe({
    # Fetch the reactive value of selected marker names
    marker_names <- selected_marker_names()
    
    if (length(marker_names) > 0) {
      # Specify the path to the text file where you want to save the names
      save_path <-paste(rv$tempFolderOut,"excluded_markers.txt", sep = "") # Replace with your directory path
      
      # Save the selected marker names to the file
      writeLines(marker_names, con = save_path)
      
    }
    else {
      excluded_markers <- "No markers selected"
      save_path <-paste(rv$tempFolderOut,"excluded_markers.txt", sep = "") # Replace with your directory path
      writeLines(excluded_markers, save_path)
    }
    
  })



#}
  
   
}

shinyApp(ui = ui, server = server) 
