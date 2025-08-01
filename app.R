#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# app.R
library(shiny)
library(tidyverse)

# Increase upload limit
options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB

ui <- fluidPage(
  titlePanel("Phospho Summary Stats"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload pr_matrix.tsv or phosphosites_90.tsv File (max. 100MB)", accept = ".tsv"),
      downloadButton("downloadData", "Download Results")
    ),
    
    mainPanel(
      tableOutput("resultTable")
    )
  )
)

server <- function(input, output) {
  
  data_processed <- reactive({
    req(input$file)
    
    df <- read_tsv(input$file$datapath)
    fname <- input$file$name
    
    # Determine mode
    if (grepl("pr_matrix", fname)) {
      # ==== Original logic ====
      n <- ncol(df)
      results_table <- data.frame(FileName = character(), Total_Protein = numeric(), Total_Peptide = numeric(), stringsAsFactors = FALSE)
      
      for (i in 11:n) {
        sample <- df[, c(1, 8, i)]
        sample_no_na <- na.omit(sample)
        sample_protein <- length(unique(sample_no_na$Protein.Group))
        sample_peptide <- length(unique(sample_no_na$Modified.Sequence))
        sample_name <- colnames(sample_no_na[, 3])
        
        results_table <- rbind(results_table, data.frame(
          FileName = sample_name,
          Total_Protein = sample_protein,
          Total_Peptide = sample_peptide
        ))
      }
      
      sty_df <- df[grep("(UniMod:21)", df$Modified.Sequence), ]
      phos_results_table <- data.frame(FileName = character(), Phos_Protein = numeric(), Phos_Peptide = numeric(), Count_S = numeric(), Count_T = numeric(), Count_Y = numeric(), stringsAsFactors = FALSE)
      
      for (i in 11:n) {
        sample <- sty_df[, c(1, 8, i)]
        sample_no_na <- na.omit(sample)
        sample_protein <- length(unique(sample_no_na$Protein.Group))
        sample_peptide <- length(unique(sample_no_na$Modified.Sequence))
        
        count_S <- sum(grepl("S\\(UniMod:21\\)", sample_no_na$Modified.Sequence))
        count_T <- sum(grepl("T\\(UniMod:21\\)", sample_no_na$Modified.Sequence))
        count_Y <- sum(grepl("Y\\(UniMod:21\\)", sample_no_na$Modified.Sequence))
        sample_name <- colnames(sample_no_na[, 3])
        
        phos_results_table <- rbind(phos_results_table, data.frame(
          FileName = sample_name,
          Phos_Protein = sample_protein,
          Phos_Peptide = sample_peptide,
          Count_S = count_S,
          Count_T = count_T,
          Count_Y = count_Y
        ))
      }
      
      results_all <- merge(results_table, phos_results_table, by = "FileName")
      results_all$Phos_Peptide_Ratio <- results_all$Phos_Peptide / results_all$Total_Peptide
      
      results_all_final <- select(results_all, FileName, Total_Protein, Phos_Protein, Total_Peptide, Phos_Peptide, Phos_Peptide_Ratio, Count_S, Count_T, Count_Y)
      return(results_all_final)
      
    } else if (grepl("phosphosites_90", fname)) {
      # ==== 90% phosphosite filtering logic ====
      n <- ncol(df)
      results_table <- data.frame(FileName = character(), Phos_Protein = numeric(), Phos_Peptide = numeric(), Count_S = numeric(), Count_T = numeric(), Count_Y = numeric(), stringsAsFactors = FALSE)
      
      for (i in 7:n) {
        sample <- df[, c(1, 4, 6, i)]
        sample_no_na <- subset(sample, sample[, 4] != 0)
        sample_protein <- length(unique(sample_no_na$Protein))
        sample_peptide <- length(unique(sample_no_na$Sequence))
        count_S <- sum(grepl("S", sample_no_na$Residue))
        count_T <- sum(grepl("T", sample_no_na$Residue))
        count_Y <- sum(grepl("Y", sample_no_na$Residue))
        
        sample_name <- colnames(sample_no_na[, 4])
        
        results_table <- rbind(results_table, data.frame(
          FileName = sample_name,
          Phos_Protein = sample_protein,
          Phos_Peptide = sample_peptide,
          Count_S = count_S,
          Count_T = count_T,
          Count_Y = count_Y
        ))
      }
      
      return(results_table)
      
    } else {
      # Unrecognized file format
      showNotification("Unrecognized file type. File must include 'pr_matrix' or 'phosphosites_90' in name.", type = "error")
      return(NULL)
    }
  })
  
  output$resultTable <- renderTable({
    data_processed()
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      if (grepl("pr_matrix", input$file$name)) {
        "Phospho_count_results.csv"
      } else {
        "Phospho_90_count_results.csv"
      }
    },
    content = function(file) {
      write.csv(data_processed(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)

