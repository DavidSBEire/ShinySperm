# Load required libraries for the app functionality
library(shiny)
library(DT)
library(plotly)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(shinythemes)

# Define the user interface for the app
ui <- fluidPage(
  theme = shinytheme("flatly"), # Using a predefined theme from shinythemes - optional
  titlePanel(img(src = "Shiny_Sperm_Logo.png", height = "100px"), "Epididymal Sperm Proteome"),
  
  # Define navigation tabs
  tabsetPanel(
    tabPanel("Data & Charts",# Next adding guide pieces of information for users
             p(em("Hints:"),"Interactive with the panel or table for proteins of interest (e.g. kinases); Each tab will adjust to your filtering",
               style = "color: #55B69A;"),
             p("Hover of the pie slices for more info",
               style = "color: #55B69A;"),
             p(em("Legend:"), "Click on types and locations to remove",
               style = "color: #55B69A;"),
             p(a("Click here for the full publication", href = "https://www.cell.com/cell-reports/fulltext/S2211-1247%2822%2901526-1",
                 style = "color: #0571b0;")),
             sidebarLayout(
               sidebarPanel(width = 3,
                # Input controls for filtering data
                 textInput("gene_search", "Search by Gene Symbol:", placeholder = "'Enter gene symbol..."),
                 selectInput("protein_type_filter", "Filter by Protein Type:",
                             choices = ""),
                 selectInput("location_filter", "Filter by Location:",
                             choices = ""),
                 selectInput("sperm_location_filter", "Filter by Sperm Location:",
                             choices = ""),
                 # Filter by Caput
                 checkboxInput("caput_filter", "Filter by Caput"),
                 # Filter by NC Cauda
                 checkboxInput("nc_cauda_filter", "Filter by NC Cauda"),
                 # Filter by Core
                 checkboxInput("core_filter", "Filter by Core"),
                 br(),
                 actionButton("reset_filters", "Reset Filters"),
                 br(),
                 br(),
                 # Feedback button
                 actionButton("feedback_button", "Feedback"),
                 br(),
                 br(),
                 # Export button
                 downloadButton("export_data", "Export Filtered Data")
               ),
               
               mainPanel(
                 # Displays and charts for data visualization
                 DTOutput("data_table"),
                 fluidRow(
                   column(width = 6,
                          # Output: Protein Type Pie Chart
                          plotlyOutput("protein_type_pie")),
                   column(width = 6,
                          # Output: Sperm Location Pie Chart
                          plotlyOutput("sperm_location_pie"))
                 )
               )
               ),
    ),
    # Defining additional tabs
    # First tab - Volcano Plot
    tabPanel("Volcano Plot",
             p(em("Hints:"), "Hover over the dots to see gene symbols",
              style = "color: #55B69A;"),
             p("You can export the adjusted plot",
               style = "color: #55B69A;"),
             fluidRow(
               column(width = 12,
               plotlyOutput("volcano_plot")),
               
             )
    ),
    
    # Next tab for Phenotypes
    tabPanel("Phenotypes",
             p(em("Hints:"),"Hover of the pie slices for more info",
               style = "color: #55B69A;"),
             p(em("Legend:"), "Click on phenotypes to remove",
              style = "color: #55B69A;"),
             p("Hou can export the associated phenotypes below",
              style = "color: #55B69A;"),
             p(a("Click here for all mouse phenotypes", href = "https://www.mousephenotype.org/",
                 style = "color: #0571b0;")),
             fluidRow(
               column(width = 6,
                      plotlyOutput("phenotypes_pie")),
               column(width = 6,
                      plotlyOutput("phenotypes_legend")),
               downloadButton("export_phenotypes_data", "Export Phenotypes Data")
             )
    ),
    
    # Next tab for Heatmap
    tabPanel("Functions",
             p(em("Hints:"),"Export the associated repro functions with the filtered table",
               style = "color: #55B69A;"),
             fluidRow(
               plotlyOutput("heatmap_plot")
             ),
             div(style = "text-align: center; margin-top: 20px;",  # Center-align the button and add some margin at the top
                 downloadButton("export_heatmap_data", "Export Heatmap Data", style = "width: 200px;")  # Set a fixed width for the button
             )
    )
  )
)

# Define server logic required to render UI
server <- function(input, output, session) {
  # Load the datasets
  # Data reactive expressions and output definitions
  # Example: Filtering data based on inputs
  dataset <- reactive({
    read.csv("Epididymal_Sperm_Proteome.csv")
  })
  
  phenotypes_data <- read.csv("Phenotypes.csv")
  
  epididymal_ipa_data <- read.csv("Epididymal_IPA.csv")
  
  # More filter choices
  observe({
    updateSelectInput(session, "protein_type_filter", choices = c("", unique(dataset()$Protein_Type)))
    updateSelectInput(session, "location_filter", choices = c("", unique(dataset()$Location)))
    updateSelectInput(session, "sperm_location_filter", choices = c("", unique(unlist(strsplit(as.character(dataset()$Sperm_Location), "; ")))))
  })
  
  # Reactive filtered dataset
  filtered_dataset <- reactive({
    filtered <- dataset()
    
    if (!is.null(input$gene_search) && input$gene_search != "") {
      filtered <- filtered[grep(input$gene_search, filtered$Gene_symbol), ]
    }
    if (input$protein_type_filter != "") {
      filtered <- filtered[filtered$Protein_Type == input$protein_type_filter, ]
    }
    if (input$location_filter != "") {
      filtered <- filtered[filtered$Location == input$location_filter, ]
    }
    if (input$sperm_location_filter != "") {
      filtered <- filtered[grepl(input$sperm_location_filter, filtered$Sperm_Location), ]
    }
    if (input$caput_filter) {
      filtered <- filtered[filtered$Caput == "Yes" | is.na(filtered$Caput), ]
    }
    if (input$nc_cauda_filter) {
      filtered <- filtered[filtered$NC_Cauda == "Yes" | is.na(filtered$NC_Cauda), ]
    }
    if (input$core_filter) {
      filtered <- filtered[filtered$CORE == "Yes" | is.na(filtered$CORE), ]
    }
    
    # Ensures that log10p_value column is retained
    filtered$log10p_value <- dataset()$log10p_value
    
    filtered
  })
  
  
  # Merge the Phenotypes data with the filtered dataset from tab1 based on the common columns
  merged_phenotypes <- reactive({
    filtered_data <- filtered_dataset()
    merge(filtered_data, phenotypes_data, by = "Accession", all.x = TRUE)
  })
  
  # Reactive expression for merged heatmap data
  merged_heatmap_data <- reactive({
    filtered_data <- filtered_dataset()
    # Merge with the epididymal IPA data based on the "Accession" column
    merge(filtered_data, epididymal_ipa_data, by = "Accession", all.x = TRUE)
  })
  
  # Render the interactive DataTables and plots
  output$data_table <- renderDT({
    datatable(
      filtered_dataset(), 
      filter = "top",
      selection = 'single', 
      rownames = FALSE,
      options = list(
        pageLength = 5, 
        lengthMenu = c(5, 10, 15),
        responsive = TRUE, # Enable responsive option for DataTable
        columnDefs = list(
          list(visible = FALSE, targets = c(5, 6, 7))  # Assuming "Caput", "NC_Cauda", "CORE" are the 6th, 7th, and 8th columns (0-indexed)
        )
      )
    )
  })
  
  # Render the protein type pie chart
  output$protein_type_pie <- renderPlotly({
    protein_type_counts <- table(filtered_dataset()$Protein_Type)
    pie_data <- data.frame(Protein_Type = names(protein_type_counts), Count = as.numeric(protein_type_counts))
    
    # Remove rows with zero count or empty Protein_Type
    pie_data <- pie_data[pie_data$Count > 0 & pie_data$Protein_Type != "", ]
    
    # Define pastel colors from RColorBrewer
    pastel_colors <- brewer.pal(length(pie_data$Count), "Pastel1")
    
    pie_chart <- plot_ly(pie_data, labels = ~Protein_Type, values = ~Count, type = "pie", 
                         title = "Protein Type",
                         marker = list(colors = pastel_colors),
                         textinfo = "none") %>%
      layout(showlegend = TRUE) %>%
      config(displayModeBar = FALSE)
    
    pie_chart
  })
  
  # Render the sperm location pie chart
  output$sperm_location_pie <- renderPlotly({
    # Split and unlist the locations
    locations <- unlist(strsplit(as.character(filtered_dataset()$Sperm_Location), "; "))
    
    # Count the occurrences of each location
    location_counts <- table(locations)
    pie_data <- data.frame(Location = names(location_counts), Count = as.numeric(location_counts))
    
    # Remove rows with zero count or empty Location
    pie_data <- pie_data[pie_data$Count > 0 & pie_data$Location != "", ]
    
    # Define pastel colors from RColorBrewer
    pastel_colors <- brewer.pal(length(pie_data$Count), "Pastel2")
    
    pie_chart <- plot_ly(pie_data, labels = ~Location, values = ~Count, type = "pie", 
                         title = "Sperm Location",
                         marker = list(colors = pastel_colors),
                         textinfo = "none") %>%
      layout(showlegend = TRUE) %>%
      config(displayModeBar = FALSE)
    
    pie_chart
  })
  
  # Render the volcano plot
  output$volcano_plot <- renderPlotly({
    volcano_data <- filtered_dataset()
    
    # Create ggplot object
    gg_volcano <- ggplot(volcano_data, aes(x = FC, y = NegLog10_p_value, text = Gene_symbol, fill = factor(ifelse(FC > 0.585 & NegLog10_p_value > 1.3, "Caput",
                                                                                                         ifelse(FC < -0.585 & NegLog10_p_value > 1.3, "NC_Cauda", "Other"))))) +
      geom_point(shape = 21, size = 3, color = "black") +
      scale_fill_manual(values = c("Caput" = "lightgreen", "NC_Cauda" = "lightblue", "Other" = "grey"),
                        labels = c("Caput", "NC_Cauda", "Other")) +
      labs(title = "Volcano Plot", x = "Log2 Fold Change (NC Cauda / Caput)", y = "-log10 (p-value)") +
      theme_minimal() +
      theme(axis.text = element_text(color = "black"),
            axis.line = element_line(color = "black", size = 1),
            plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white"),
            plot.title = element_text(color = "black")) +
      ylim(0, NA) +  # Limit y-axis minimum to be at zero and no lower
      geom_vline(xintercept = c(0.585, -0.585), linetype = "dotted", size = 0.3) +  # Add dotted lines at specified positions
      geom_hline(yintercept = 1.3, linetype = "dotted", size = 0.3)  # Add dotted lines at specified positions
    
    # Convert ggplot object to Plotly plot
    volcano_plot <- ggplotly(gg_volcano, tooltip = "text") %>%
      layout(height = 600, width = 1000,
        showlegend = FALSE)
    # Remove the legend
    
    volcano_plot
  })
  
  # Render the pie chart for Phenotypes
  output$phenotypes_pie <- renderPlotly({
    phenotypes_data <- merged_phenotypes()
    
    # Filter out rows where phenotype is null
    phenotypes_data <- phenotypes_data[!is.na(phenotypes_data$Phenotype), ]
    
    
    # Create the pie chart
    pie_chart <- plot_ly(phenotypes_data, labels = ~Phenotype, type = "pie",
                         textinfo = "none")  # Remove numbers from the pie chart
    
    # Customize layout
    pie_chart <- pie_chart %>%
      layout(title = "Associated Phenotypes")
    
    # Return the plotly object
    pie_chart
  })
  
  # Render the heatmap plot
  output$heatmap_plot <- renderPlotly({
    # Get the merged data for the heatmap
    heatmap_data <- merged_heatmap_data() %>%
      # Group by Function_Annotation and summarize the data, in case there are multiple entries per function
      group_by(Function_Annotation) %>%
      summarize(log10p_value = mean(log10p_value, na.rm = TRUE)) %>%
      # Take the top entries based on log10p_value
      top_n(20, log10p_value) %>%
      ungroup() %>%
      # Arrange might not be necessary here as Plotly will handle ordering, but included for clarity
      arrange(desc(log10p_value))
    
    # Create the heatmap plot
    heatmap_plot <- ggplot(heatmap_data, aes(x = factor(1), y = reorder(Function_Annotation, log10p_value), fill = log10p_value)) +
      geom_tile(color = "black", size = 0.15) +
      scale_fill_gradient(low = "#bdc9e1", high = "#016450", limits = c(1.3, 8), na.value = "#016450") +
      labs(x = NULL, y = NULL, fill = "-log10 (p-value)") +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.background = element_blank())
    
    # Convert ggplot to Plotly
    heatmap_plot <- ggplotly(heatmap_plot, tooltip = "y+fill")
    
    # Adjust layout options for height and width
    heatmap_plot <- layout(heatmap_plot, height = 600, width = 600)  # Set desired height and width
    
    heatmap_plot
  })
  # Export data from the heatmap
  output$export_heatmap_data <- downloadHandler(
    filename = function() {
      paste("heatmap_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(merged_heatmap_data(), file, row.names = FALSE)
    }
  )
  
  # Export data from the pie chart
  output$export_phenotypes_data <- downloadHandler(
    filename = function() {
      "phenotypes_data.csv"
    },
    content = function(file) {
      phenotypes_data <- merged_phenotypes()
      phenotypes_data <- phenotypes_data[!is.na(phenotypes_data$Phenotype), ]
      write.csv(phenotypes_data, file, row.names = FALSE)
    }
  )
  
  
  # Submit feedback
  observeEvent(input$feedback_button, {
    # Show modal dialog for feedback
    showModal(modalDialog(
      title = "Feedback",
      p("Please email your feedback to ", a(href = "mailto:David.Skerrett-Byrne@helmholtz-munich.de", "David.Skerrett-Byrne@helmholtz-munich.de"), ".")
    ))
  })
  
  # Submit feedback
  observeEvent(input$submit_feedback, {
    feedback <- isolate(input$feedback_text)
    
    # Send email
    
    # Show confirmation message
    showModal(modalDialog(
      title = "Thank you for your feedback!",
      footer = modalButton("Close")
    ))
  })
  
  # Export filtered data
  output$export_data <- downloadHandler(
    filename = function() {
      paste("filtered_data", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_dataset(), file, row.names = FALSE)
    }
  )
  
  # Reset filters
  observeEvent(input$reset_filters, {
    updateSelectInput(session, "protein_type_filter", selected = "")
    updateSelectInput(session, "location_filter", selected = "")
    updateSelectInput(session, "sperm_location_filter", selected = "")
    updateCheckboxInput(session, "caput_filter", value = FALSE)
    updateCheckboxInput(session, "nc_cauda_filter", value = FALSE)
    updateCheckboxInput(session, "core_filter", value = FALSE)
    updateTextInput(session, "gene_search", value = "")
  })
}

# Run the application
shinyApp(ui = ui, server = server)
