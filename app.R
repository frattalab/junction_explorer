# Load required libraries
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(tidyverse)
library(arrow)
library(pROC)
library(ggrepel)
library(data.table)
library(glue)

# --- Plotting Functions ---

# NYGC plot function
plot_junction_nygc <- function(junc, plotin_table, measure = 'spliced_reads',
                               gene_name = "") {
    
    vals <- c(`ALS-TDP` = "#E1BE6A", Control = "#40B0A6", `FTLD-TDP` = "#E1BE6A", 
              `ALS-non-TDP` = "#408A3E", `FTLD-non-TDP` = "#408A3E")
    
    plotin_table <- plotin_table %>% 
        mutate(tdp_proteinopathy = ifelse(tdp_proteinopathy %in% c("FTLD-TAU","FTLD-FUS"), "FTLD-non-TDP", tdp_proteinopathy))
    
    if(gene_name != ""){
        plot_title <- glue::glue("{gene_name} - {junc}")
    }else{
        plot_title <- glue::glue("{junc}")
    }
    
    plotin_table =plotin_table %>% 
        filter(rna_tissue_source_simplified %in% c("Spinal_Cord_Thoracic", "Spinal_Cord_Cervical", "Spinal_Cord_Lumbar", 
                                                   "Cortex_Frontal", "Cortex_Motor", "Cortex_Temporal", 
                                                   "Medulla", "Hippocampus", "Choroid", 
                                                   "Cortex_Unspecified", 
                                                   "Spinal_Cord_Unspecified","Cerebellum")) %>% 
        filter(tdp_proteinopathy != "Unknown")
    
    
    if(measure == 'spliced_reads') {
        plt <- plotin_table %>%
            mutate(tdp_proteinopathy = fct_relevel(tdp_proteinopathy, "Control", "ALS-non-TDP", "ALS-TDP", "FTLD-non-TDP", "FTLD-TDP")) %>%
            ggplot(aes(x = tdp_proteinopathy, y = junction_count, fill = tdp_proteinopathy)) +
            geom_boxplot(outlier.colour = NA) +
            geom_jitter(height = 0, alpha = 0.7, pch = 21) +
            facet_wrap(~rna_tissue_source_simplified, scales = "free_y") +
            scale_fill_manual(values = vals) +
            ylab("N spliced reads") +
            xlab("") +
            ggtitle(plot_title) +
            ggpubr::theme_pubr() +
            theme(legend.position = 'none') +
            theme(text = element_text(size = 10)) +
            scale_x_discrete(guide = guide_axis(n.dodge = 2))
    } else {
        plt <- plotin_table %>%
            mutate(tdp_proteinopathy = fct_relevel(tdp_proteinopathy, "Control", "ALS-non-TDP", "ALS-TDP", "FTLD-non-TDP", "FTLD-TDP")) %>%
            ggplot(aes(x = tdp_proteinopathy, y = psi, fill = tdp_proteinopathy)) +
            geom_boxplot(outlier.colour = NA) +
            geom_jitter(height = 0, alpha = 0.7, pch = 21) +
            facet_wrap(~rna_tissue_source_simplified, scales = "free_y") +
            scale_fill_manual(values = vals) +
            ylab("PSI") +
            xlab("") +
            ggtitle(plot_title) +
            ggpubr::theme_pubr() +
            theme(legend.position = 'none') +
            theme(text = element_text(size = 10)) +
            scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
            scale_y_continuous(labels = scales::percent_format())
    }
    
    return(plt)
}

plot_junction_rimod <- function(junc, table, measure = 'spliced_reads',
                                gene_name = "") {
    
    vals <- c(`FTD-GRN` = "#E1BE6A", control = "#40B0A6", `FTD-C9` = "#E1BE6A", 
              `FTD-MAPT` = "#408A3E")
    if(gene_name != ""){
        plot_title <- glue::glue("{gene_name} - {junc}")
    }else{
        plot_title <- glue::glue("{junc}")
    }
    
    if(measure == 'psi'){
        p = table %>% 
            mutate(DiseaseCode = fct_relevel(DiseaseCode,'control','FTD-MAPT')) %>% 
            as.data.table() %>% 
            mutate(psi = junction_count / cluster_count) %>% 
            ggplot(aes(y = psi,x = DiseaseCode,fill = DiseaseCode)) + 
            geom_boxplot() + 
            geom_jitter(height = 0) +
            theme_classic() + 
            ylab("Junction PSI") +
            scale_fill_manual(values = vals) +
            ylab("PSI") +
            xlab("") +
            ggtitle(plot_title) +
            ggpubr::theme_pubr() +
            theme(legend.position = 'none') +
            theme(text = element_text(size = 10)) +
            scale_y_continuous(labels = scales::percent_format())
        
    }else{
        
        p = table %>% 
            mutate(DiseaseCode = fct_relevel(DiseaseCode,'control','FTD-MAPT')) %>% 
            as.data.table() %>% 
            ggplot(aes(y = junction_count,x = DiseaseCode,fill = DiseaseCode)) + 
            geom_boxplot() + 
            geom_jitter(height = 0) +
            theme_classic() + 
            ylab("N spliced reads") +
            scale_fill_manual(values = vals) +
            xlab("") +
            ggtitle(plot_title) +
            ggpubr::theme_pubr() +
            theme(legend.position = 'none') +
            theme(text = element_text(size = 10)) 
        
    }
    return(p)
}

# query the parquets files
query_junction = function(junc, dataset_parquet){
    parts <- unlist(strsplit(junc, "[:\\-]"))
    queried_chrom = parts[[1]]
    queried_start = as.numeric(parts[[2]])
    queried_end = as.numeric(parts[[3]])
    
    dataset_parquet_query = arrow::open_dataset(dataset_parquet)
    
    junction_data <- dataset_parquet_query %>%
        filter(chrom == queried_chrom, 
               start == queried_start, 
               end == queried_end) %>%
        collect()
    return(junction_data)
}

# --- UI ---

ui <- dashboardPage(
    dashboardHeader(title = "Splice Junction Expression Viewer"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Junction Plotter", tabName = "plotter", icon = icon("chart-bar")),
            menuItem("Compendium events", tabName = "input_junctions", icon = icon("table-list")),
            menuItem("About", tabName = "about", icon = icon("info-circle"))
        )
    ),
    
    
    dashboardBody(
        tabItems(
            # Junction Plotter Tab
            tabItem(tabName = "plotter",
                    fluidRow(
                        box(
                            title = "Input Parameters", 
                            status = "primary", 
                            solidHeader = TRUE,
                            width = 4,
                            
                            textInput("junction_input", 
                                      "Enter Junction Coordinates:", 
                                      value = "chr2:241668985-241670726",
                                      placeholder = "e.g., chr1:1000-2000"),
                            
                            radioButtons("measure_type", 
                                         "Measure Type:",
                                         choices = list(
                                             "Spliced Reads" = "spliced_reads",
                                             "PSI" = "psi"
                                         ),
                                         selected = "spliced_reads"),
                            
                            actionButton("plot_button", 
                                         "Generate Plots", 
                                         class = "btn-primary",
                                         style = "width: 100%;"),
                            
                            br(), br(),
                            
                            downloadButton("download_nygc_plot", 
                                           "Download NYGC Plot",
                                           class = "btn-success",
                                           style = "width: 100%; margin-bottom: 5px;"),
                            
                            downloadButton("download_rimod_plot", 
                                           "Download RIMOD Plot",
                                           class = "btn-success",
                                           style = "width: 100%; margin-bottom: 5px;"),
                            
                            downloadButton("download_nygc_data", 
                                           "Download NYGC Data",
                                           class = "btn-info",
                                           style = "width: 100%; margin-bottom: 5px;"),
                            
                            downloadButton("download_rimod_data", 
                                           "Download RIMOD Data",
                                           class = "btn-info",
                                           style = "width: 100%;")
                        ),
                        
                        box(
                            title = "NYGC Junction Expression Plot", 
                            status = "primary", 
                            solidHeader = TRUE,
                            width = 8,
                            
                            conditionalPanel(
                                condition = "output.nygc_plot_ready",
                                plotOutput("nygc_plot", height = "600px")
                            ),
                            
                            conditionalPanel(
                                condition = "!output.nygc_plot_ready",
                                div(
                                    style = "text-align: center; padding: 50px;",
                                    h4("Enter a junction coordinate and click 'Generate Plots' to visualize NYGC expression data"),
                                    icon("chart-bar", style = "font-size: 48px; color: #ccc;")
                                )
                            )
                        ),
                        
                        box(
                            title = "RIMOD Junction Expression Plot", 
                            status = "primary", 
                            solidHeader = TRUE,
                            width = 8,
                            
                            conditionalPanel(
                                condition = "output.rimod_plot_ready",
                                plotOutput("rimod_plot", height = "600px")
                            ),
                            
                            conditionalPanel(
                                condition = "!output.rimod_plot_ready",
                                div(
                                    style = "text-align: center; padding: 50px;",
                                    h4("Enter a junction coordinate and click 'Generate Plots' to visualize RIMOD expression data"),
                                    icon("chart-bar", style = "font-size: 48px; color: #ccc;")
                                )
                            )
                        )
                    )
            ),
            
            # Possible Junctions Tab
            tabItem(tabName = "input_junctions",
                    fluidRow(
                        box(
                            title = "Possible Junctions to Plot",
                            status = "primary",
                            solidHeader = TRUE,
                            width = 12,
                            DTOutput("junction_table"),
                            br(), 
                            downloadButton("download_junction_table", "Download Table", class = "btn-success", style = "width: 100%;")
                        )
                    )
            ),
            
            # About Tab
            tabItem(tabName = "about",
                    fluidRow(
                        box(
                            title = "About This App", 
                            status = "primary", 
                            solidHeader = TRUE,
                            width = 12,
                            
                            h3("Splice Junction Expression Viewer"),
                            p("This application allows you to visualize splice junction expression data from two different datasets, NYGC and RIMOD, across various disease conditions and tissue types."),
                            
                            h4("Features:"),
                            tags$ul(
                                tags$li("Interactive plotting of splice junction expression for both NYGC and RIMOD datasets"),
                                tags$li("Support for both spliced reads count and PSI (Percent Spliced In) metrics"),
                                tags$li("Comparison across different disease proteinopathy conditions"),
                                tags$li("Data exploration and export capabilities for each dataset")
                            ),
                            
                            h4("Usage:"),
                            p("1. Navigate to the 'Junction Plotter' tab"),
                            p("2. Enter a junction coordinate (format: chr:start-end)"),
                            p("3. Select your preferred measure type"),
                            p("4. Click 'Generate Plots' to visualize the data for both datasets"),
                            
                            h4("Data Source:"),
                            p("Data is sourced from parquet files containing splice junction analysis results from neurodegeneration studies.")
                        )
                    )
            )
        )
    )
)

# --- Server ---

server <- function(input, output, session) {
    
    # Reactive values to store data and plots for both datasets
    values <- reactiveValues(
        nygc_data = NULL,
        nygc_plot = NULL,
        rimod_data = NULL,
        rimod_plot = NULL
    )
    
    # Read the possible junctions data
    possible_events <- reactive({
        fread("data/possible_junctions.csv")
    })
    
    # Render the interactive table for possible junctions
    output$junction_table <- renderDT({
        dt_data <- possible_events()
        datatable(
            dt_data,
            filter = 'top',
            options = list(
                pageLength = 10,
                dom = 'lfrtip', 
                columnDefs = list(list(targets = which(names(dt_data) == 'gene_name') - 1, searchable = TRUE))
            ),
            rownames = FALSE
        )
    })
    
    # Generate plots and data when button is clicked
    observeEvent(input$plot_button, {
        req(input$junction_input)
        
        # Show loading message
        showModal(modalDialog(
            title = "Loading...",
            "Generating plots, please wait...",
            footer = NULL,
            easyClose = FALSE
        ))
        
        tryCatch({
            # Load metadata once
            nygc_metadata <- fread("data/nygc/nygc_metadata.csv")
            rimod_metadata <- fread("data/rimod/rimod_meta.csv")
            possible_events_dt <-  fread("data/possible_junctions.csv")
            
            
            # If we've got a gene name for that event in this table
            junction_name = possible_events_dt %>% 
                filter(junctions_coords == input$junction_input) %>% 
                distinct(gene_name) %>% 
                pull(gene_name)
            row_num = possible_events_dt %>% 
                filter(junctions_coords == input$junction_input) %>% 
                distinct(gene_name) %>% nrow()

            if(row_num == 0){
                junction_name = ""
            }
            
            # Query NYGC data
            nygc_junction_data <- query_junction(junc = input$junction_input,
                                                 dataset_parquet = "data/nygc/parquets")
            nygc_junction_data = nygc_junction_data %>% 
                distinct(chrom,sample,cluster_count,start,end,cluster,strand,junction_count,psi)
            
            # Query RIMOD data
            rimod_junction_data <- query_junction(junc = input$junction_input,
                                                  dataset_parquet = "data/rimod/parquets")

            # Process and store NYGC data
            if (nrow(nygc_junction_data) > 0) {
                nygc_data_with_meta <- nygc_junction_data %>% left_join(nygc_metadata, by = "sample")
                values$nygc_data <- nygc_data_with_meta
                values$nygc_plot <- plot_junction_nygc(input$junction_input, 
                                                       nygc_data_with_meta, 
                                                       input$measure_type,
                                                       junction_name)
            } else {
                values$nygc_data <- NULL
                values$nygc_plot <- NULL
                showNotification("No data found for this junction in NYGC dataset.", type = "warning")
            }
            
            # Process and store RIMOD data
            if (nrow(rimod_junction_data) > 0) {
                rimod_data_with_meta <- rimod_junction_data %>% left_join(rimod_metadata, by = "sample")
                values$rimod_data <- rimod_data_with_meta
                values$rimod_plot <- plot_junction_rimod(input$junction_input, 
                                                         rimod_data_with_meta, 
                                                         input$measure_type,
                                                         junction_name)
            } else {
                values$rimod_data <- NULL
                values$rimod_plot <- NULL
                showNotification("No data found for this junction in RIMOD dataset.", type = "warning")
            }
            
        }, error = function(e) {
            showNotification(paste("Error:", e$message), type = "error")
            values$nygc_plot <- NULL
            values$rimod_plot <- NULL
        })
        
        removeModal()
    })
    
    # --- Plot Rendering ---
    
    output$nygc_plot <- renderPlot({
        req(values$nygc_plot)
        values$nygc_plot
    })
    
    output$rimod_plot <- renderPlot({
        req(values$rimod_plot)
        values$rimod_plot
    })
    
    # --- Plot Ready Status ---
    
    output$nygc_plot_ready <- reactive({
        !is.null(values$nygc_plot)
    })
    outputOptions(output, "nygc_plot_ready", suspendWhenHidden = FALSE)
    
    output$rimod_plot_ready <- reactive({
        !is.null(values$rimod_plot)
    })
    outputOptions(output, "rimod_plot_ready", suspendWhenHidden = FALSE)
    
    # --- Download Handlers ---
    
    # Download NYGC plot
    output$download_nygc_plot <- downloadHandler(
        filename = function() {
            junction_str <- str_replace_all(input$junction_input, "[^A-Za-z0-9]", "_")
            paste0("nygc_plot_", junction_str, "_", Sys.Date(), ".png")
        },
        content = function(file) {
            req(values$nygc_plot)
            ggsave(file, values$nygc_plot, width = 12, height = 8, dpi = 300)
        }
    )
    
    # Download RIMOD plot
    output$download_rimod_plot <- downloadHandler(
        filename = function() {
            junction_str <- str_replace_all(input$junction_input, "[^A-Za-z0-9]", "_")
            paste0("rimod_plot_", junction_str, "_", Sys.Date(), ".png")
        },
        content = function(file) {
            req(values$rimod_plot)
            ggsave(file, values$rimod_plot, width = 12, height = 8, dpi = 300)
        }
    )
    
    # Download NYGC data
    output$download_nygc_data <- downloadHandler(
        filename = function() {
            junction_str <- str_replace_all(input$junction_input, "[^A-Za-z0-9]", "_")
            paste0("nygc_data_", junction_str, "_", Sys.Date(), ".csv")
        },
        content = function(file) {
            req(values$nygc_data)
            readr::write_csv(values$nygc_data, file)
        }
    )
    
    # Download RIMOD data
    output$download_rimod_data <- downloadHandler(
        filename = function() {
            junction_str <- str_replace_all(input$junction_input, "[^A-Za-z0-9]", "_")
            paste0("rimod_data_", junction_str, "_", Sys.Date(), ".csv")
        },
        content = function(file) {
            req(values$rimod_data)
            readr::write_csv(values$rimod_data, file)
        }
    )
    # Possible junctions
    
    output$download_junction_table <- downloadHandler(
        filename = function() {
            paste0("possible_junctions_", Sys.Date(), ".csv")
        },
        content = function(file) {
            req(possible_events())
            readr::write_csv(possible_events(), file)
        }
    )
    
    
}

# Run the application
shinyApp(ui = ui, server = server)