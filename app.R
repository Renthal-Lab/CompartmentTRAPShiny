library(shiny)
library(tidyverse)
library(DT)
library(ggplot2)
library(ggrepel)

# Set working directory
#setwd("~/Partners HealthCare Dropbox/Renthal Lab/Team/Shams Bhuiyan/Liz_shiny/")
#EVERYTHING BUT THE ASTERIKS
# Load Data
PA_counts <- read_csv("PA_counts.csv")
CA_counts <- read_csv("CA_counts.csv")
DRG_counts <- read_csv("DRG_counts.csv")

PA_DE <- read_csv("PA_DE.csv")
CA_DE <- read_csv("CA_DE.csv")
DRG_DE <- read_csv("DRG_DE.csv")

GO_PA <- read_csv("GO_PA.csv")
GO_CA <- read_csv("GO_CA.csv")

# Function to transform counts to long format for plotting
convert_counts_long <- function(df, condition) {
  df_long <- df %>%
    pivot_longer(cols = starts_with(condition) | starts_with("X"), 
                 names_to = "Sample", 
                 values_to = "TPM") %>%
    mutate(
      Condition = condition,
      GFP_Status = factor(ifelse(grepl("^X", Sample), "GFP-", "GFP+"), levels = c("GFP+", "GFP-"))
    )
  return(df_long)
}

# Convert count tables to long format
PA_long <- convert_counts_long(PA_counts, "PA")
CA_long <- convert_counts_long(CA_counts, "CA")
DRG_long <- convert_counts_long(DRG_counts, "DRG")

# Merge all conditions
all_counts_long <- bind_rows(PA_long, CA_long, DRG_long)

# UI
ui <- fluidPage(
  titlePanel("TRAP-seq Data Explorer"),
  
  tabsetPanel(
    tabPanel("Home",
             fluidRow(
               column(12,
                      h2("Graphical Abstract"),
                      img(src = "graphical_abstract.png", width = "50%") 
               )
             )
    ),
    
    tabPanel("Gene Expression",
             sidebarLayout(
               sidebarPanel(
                 textInput("gene_search", "Search for a Gene:", ""),
                 actionButton("search_gene", "Search")
               ),
               mainPanel(
                 plotOutput("expression_plot"),
                 fluidRow(

                   column(4, plotOutput("volcano_ca")),
                   column(4, plotOutput("volcano_drg")),
                   column(4, plotOutput("volcano_pa"))
                 ),
                 textOutput("no_gene_found")
               )
             )
    ),
    
    tabPanel("GO Term Analysis",
             sidebarLayout(
               sidebarPanel(
                 selectInput("go_term", "Select GO Term:", 
                             choices = unique(c(GO_PA$Description, GO_CA$Description))),
                 actionButton("search_go", "Find GO Genes")
               ),
               mainPanel(
                 DT::dataTableOutput("go_gene_table"),
                 textOutput("no_go_genes_found")
               )
             )
    )
  )
)

# Server
server <- function(input, output) {
  
  # Reactive function for gene selection
  selected_gene_data <- reactive({
    req(input$gene_search)
    
    gene_name <- input$gene_search
    
    filtered_data <- all_counts_long %>% filter(gene_name == !!gene_name)
    
    if (nrow(filtered_data) == 0) {
      return(NULL)
    }
    
    return(filtered_data)
  })
  
  # Barplot for gene expression (fixes missing conditions issue)
  output$expression_plot <- renderPlot({
    data <- selected_gene_data()
    
    if (is.null(data)) {
      return(NULL)
    }
    
    required_conditions <- c("PA", "DRG", "CA")
    missing_conditions <- setdiff(required_conditions, unique(data$Condition))
    
    if (length(missing_conditions) > 0) {
      missing_data <- tibble(
        Condition = missing_conditions,
        Sample = NA,
        TPM = 0,  
        GFP_Status = factor("GFP+", levels = c("GFP+", "GFP-"))  
      )
      data <- bind_rows(data, missing_data)
    }
    data <- data %>% mutate(Condition = case_when(
      Condition == "PA" ~ "Peripheral Axon",
      Condition == "CA" ~ "Central Axon",
      TRUE ~ Condition
    ))
    print(data)
    gene_name<-data[1,3] %>% pull()
    print(gene_name)
    
    PA_sig = "n.s"
    
    if(gene_name %in% PA_DE$symbol & PA_DE$padj<0.05){
      PA_sig = "*"
    }
    
    CA_sig = "n.s"
    
    if(gene_name %in% CA_DE$symbol & CA_DE$padj<0.05){
      CA_sig = "*"
    }
    
    DRG_sig = "n.s"
    
    if(gene_name %in% DRG_DE$symbol & DRG_DE$padj<0.05){
      DRG_sig = "*"
    }
    
    ggplot(data, aes(x = Condition, y = TPM, fill = GFP_Status)) +
      geom_bar(stat = "summary", fun = "mean", position = position_dodge(), color = "black", width = 0.5) +
      geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.1, position = position_dodge(0.5)) +
      theme_minimal() +
      labs(
        title = paste("Expression of", input$gene_search),
        y = "Normalized counts",
        x = "Compartment",
        caption = "Error bars represent standard error. Missing conditions are set to 0."
      ) +
      scale_fill_manual(values = c("GFP+" = "#A3D9A5", "GFP-" = "#F4B6C2")) +
      theme(
        text = element_text(size = 15),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.9, 0.9)
      ) +
      # Adding text annotation for "Peripheral Axon"
      geom_text(
        data = data %>% filter(Condition == "Peripheral Axon" & GFP_Status == "GFP+"), 
        aes(label = PA_sig, y = max(TPM) * 1.01), 
        position = position_dodge(width = 0.5),
        size = 5
      ) +
      # Adding text annotation for "Peripheral Axon"
      geom_text(
        data = data %>% filter(Condition == "Central Axon"& GFP_Status == "GFP+"), 
        aes(label = CA_sig, y = max(TPM) * 1.01), 
        position = position_dodge(width = 0.5),
        size = 5
      ) +
      geom_text(
        data = data %>% filter(Condition == "DRG"& GFP_Status == "GFP+"), 
        aes(label = DRG_sig, y = max(TPM) * 1.01), 
        position = position_dodge(width = 0.5),
        size = 5
      )
    
  }) 
  
  # Function to generate Volcano plots
  generate_volcano_plot <- function(de_data, condition) {
    req(input$gene_search)
    gene_name <- input$gene_search
    
    de_data <- de_data %>% mutate(condition = case_when(
      condition == "PA" ~ "Peripheral Axon",
      condition == "CA" ~ "Central Axon",
      TRUE ~ condition
    ))
    print(condition)
    if(condition == "CA"){
      condition = "Central axon"
    }
    else if(condition == "PA"){
      condition = "Peripheral axon"
    }
    de_data <- de_data %>%
      mutate(highlight = ifelse(symbol == gene_name, "Highlighted", "Normal"))
    print((de_data %>% filter(symbol=="Calca")))
    ggplot(de_data, aes(x = log2FoldChange, y = -log10(padj), color = highlight)) +
      geom_point(alpha = 0.6) +
      geom_text_repel(data = filter(de_data, highlight == "Highlighted"), 
                      aes(label = symbol), size = 5, box.padding = 0.3) +
      scale_color_manual(values = c("Highlighted" = "red", "Normal" = "black")) +
      theme_minimal() +
      labs(
        title = paste(condition, "GFP+ vs GFP-"),
        x = "Log2 Fold Change",
        y = "-log10(padj)"
      ) +
      theme(legend.position = "none", text = element_text(size = 15),plot.title = element_text(hjust = 0.5)) 
  }
  
  # Volcano plots

  output$volcano_ca <- renderPlot({
    generate_volcano_plot(CA_DE, "CA")
  })
  
  output$volcano_drg <- renderPlot({
    generate_volcano_plot(DRG_DE, "DRG")
  })
  
  output$volcano_pa <- renderPlot({
    generate_volcano_plot(PA_DE, "PA")
  })
  
  # Reactive function for GO term selection
  selected_go_data <- reactive({
    req(input$go_term)
    
    go_pa_genes <- GO_PA %>% filter(Description == input$go_term) %>% pull(geneID)
    go_ca_genes <- GO_CA %>% filter(Description == input$go_term) %>% pull(geneID)
    
    pa_genes_list <- unlist(strsplit(go_pa_genes, "/"))
    ca_genes_list <- unlist(strsplit(go_ca_genes, "/"))
    
    all_genes <- unique(c(pa_genes_list, ca_genes_list))
    
    go_table <- tibble(
      Gene = all_genes,
      PA_Enrichment = ifelse(Gene %in% pa_genes_list, "Enriched", "Not Enriched"),
      CA_Enrichment = ifelse(Gene %in% ca_genes_list, "Enriched", "Not Enriched")
    )
    
    return(go_table)
  })
  
  # Render GO gene table
  output$go_gene_table <- DT::renderDataTable({
    req(selected_go_data())
    
    datatable(selected_go_data(), options = list(pageLength = 10)) %>%
      formatStyle(
        "PA_Enrichment",
        backgroundColor = styleEqual(c("Enriched", "Not Enriched"), c("#89A8B2", "#E4E0E1"))
      ) %>%
      formatStyle(
        "CA_Enrichment",
        backgroundColor = styleEqual(c("Enriched", "Not Enriched"), c("#89A8B2", "#E4E0E1"))
      )
  })
}

# Run the app
shinyApp(ui = ui, server = server)
