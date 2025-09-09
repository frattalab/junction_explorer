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
library(stringr)

# --- Plotting Functions ---
source("plot_gene_expression.R")

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
    dashboardHeader(title = "Cryptic Creeper"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Compendium events", tabName = "input_junctions", icon = icon("table-list")),
            menuItem("Junction Plotter", tabName = "plotter", icon = icon("chart-bar")),
            menuItem("Gene expression", tabName = "gene_expression_tab", icon = icon("dna")),
            menuItem("About", tabName = "about", icon = icon("info-circle"))
        )
    ),
    
    
    dashboardBody(
        tabItems(
            # Junction Plotter Tab
            tabItem(tabName = "plotter",
                    # First Row: Input Parameters
                    fluidRow(
                        box(
                            title = "Input Parameters",
                            status = "primary",
                            solidHeader = TRUE,
                            width = 12, # Full width for the input box
                            fluidRow( # Use a nested fluidRow for better control
                                column(4,
                                       textInput("junction_input",
                                                 "Enter Junction Coordinates:",
                                                 value = "chr2:241668985-241670726",
                                                 placeholder = "e.g., chr1:1000-2000")),
                                column(4,
                                       radioButtons("measure_type",
                                                    "Measure Type:",
                                                    choices = list(
                                                        "Spliced Reads" = "spliced_reads",
                                                        "PSI" = "psi"
                                                    ),
                                                    selected = "spliced_reads")),
                                column(4,
                                       actionButton("plot_button",
                                                    "Generate Plots",
                                                    class = "btn-primary",
                                                    style = "width: 100%; margin-top: 25px;")) # Align with the other inputs
                            )
                        )
                    ),
                    
                    # Second Row: NYGC Plot
                    fluidRow(
                        box(
                            title = "NYGC Junction Expression Plot",
                            status = "primary",
                            solidHeader = TRUE,
                            width = 12, # Full width for the plot
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
                        )
                    ),
                    
                    # Third Row: RIMOD Plot
                    fluidRow(
                        box(
                            title = "RIMOD Junction Expression Plot",
                            status = "primary",
                            solidHeader = TRUE,
                            width = 12, # Full width for the plot
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
                    ),
                    
                    # Fourth Row: Download Buttons
                    fluidRow(
                        box(
                            title = "Download Options",
                            status = "info",
                            solidHeader = TRUE,
                            width = 12, # Full width for the buttons
                            tags$div(
                                style = "display: flex; flex-wrap: wrap; gap: 10px; justify-content: center;",
                                downloadButton("download_nygc_plot", "Download NYGC Plot", class = "btn-success"),
                                downloadButton("download_rimod_plot", "Download RIMOD Plot", class = "btn-success"),
                                downloadButton("download_nygc_data", "Download NYGC Data", class = "btn-info"),
                                downloadButton("download_rimod_data", "Download RIMOD Data", class = "btn-info")
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
            
            tabItem(tabName = "gene_expression_tab",
                    fluidRow(
                        # Input box - full width at top
                        box(
                            title = "Gene Expression Analysis",
                            status = "primary",
                            solidHeader = TRUE,
                            width = 12,
                            textInput("gene_input", 
                                      "Enter Gene Symbol:", 
                                      value = "ELAVL3",
                                      placeholder = "e.g., ELAVL3",
                                      width = "300px"),  # Constrain input width
                            
                            actionButton("generate_plots", 
                                         "Generate Plots", 
                                         class = "btn-primary",
                                         style = "margin-left: 10px;")
                        )
                    ),
                    
                    fluidRow(
                        # Plot output box - full width in middle
                        box(
                            title = "Expression Plots",
                            status = "primary", 
                            solidHeader = TRUE,
                            width = 12,
                            conditionalPanel(
                                condition = "input.generate_plots > 0",
                                plotOutput("gene_expression_plot", height = "800px")
                            ),
                            conditionalPanel(
                                condition = "input.generate_plots == 0",
                                tags$div(
                                    style = "text-align: center; padding: 50px; color: #999;",
                                    icon("chart-line", "fa-3x"),
                                    br(), br(),
                                    h4("Plots will appear here after clicking 'Generate Plots'")
                                )
                            )
                        )
                    ),
                    
                    fluidRow(
                        # Download buttons - full width at bottom
                        box(
                            title = "Download Datasets",
                            status = "info",
                            solidHeader = TRUE,
                            width = 12,
                            tags$div(
                                style = "display: flex; flex-wrap: wrap; gap: 10px; justify-content: center;",
                                downloadButton("download_rna_expression_nmd", "Fractionation DE", class = "btn-info"),
                                downloadButton("download_sh_exp", "SH-SY5Y Dose Curve", class = "btn-info"),
                                downloadButton("download_sk_exp", "SK-N-BE(2) Dose Curve", class = "btn-info"),
                                downloadButton("download_upf1_rna", "UPF1 RNA", class = "btn-info"),
                                downloadButton("download_chx_rna", "CHX RNA", class = "btn-info")
                            )
                        )
                    )
            ),
            
            # About Tab
            tabItem(tabName = "about",
                    fluidRow(
                        box(
                            title = "What can I do this with?", 
                            status = "primary", 
                            solidHeader = TRUE,
                            width = 12,
                            
                            h3("Cryptic Creeper"),
                            p(""),
                            
                            h4("Features:"),
                            tags$ul(
                                tags$li("Interactive plotting of splice junction expression for both NYGC and RIMOD datasets"),
                                tags$li("Support for both spliced reads count and PSI (Percent Spliced In) metrics"),
                                tags$li("Comparison across different disease proteinopathy conditions"),
                                tags$li("Data exploration and export capabilities for each dataset")
                            ),
                            
                            h4("Usage:"),
                            h4("Compendium events:"),
                            p("TDP-43 cryptic splicing events we discovered in across multiple (14) TDP-43 knockdowns as well as NMD-inhibition"),
                            p("Explore the detection, FPR, and TPR of all the events we discovered in postmortem"),
                            p("Click on a row to send the junction and gene to the 'Junction Plotter' and 'Gene expression' tabs"),
                            
                            
                            h4("Junction plotter:"),
                            p("1. Enter a junction coordinate (format: chr:start-end)"),
                            p("2. Select your preferred measure type"),
                            p("3. Click 'Generate Plots' to visualize the data for both datasets"),
                            
                            h4("Gene expression:"),
                            p("Explore effect of dose-dependent TDP-43 knockdown, NMD-inhibition, and the nuclear/cytoplasmic expression of your favorite gene"),
                            
                            
                            h4("Data Source:"),
                            p("Each dataset had its splice junctions extracted from the BAM files with regtools and then was run through LeafCutter's clustering tool to generate PSI values"),
                            
                            h4("Want more?"),
                            p("Due to space limitations, only a subset of all splicing events are available here, but you can download the full files and run the app on your desktop"),
                            p("To run the app on your desktop for better performance and explore all  splice junctions in both datasets,"),
                            tags$a(href="https://github.com/frattalab/junction_explorer/", "Click here!")
                            
                            
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
        rimod_plot = NULL,
        gene_expr_data = NULL,
        gene_expr_plot = NULL
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
                scrollX = TRUE,
                scrollY = "400px",
                autoWidth = TRUE,
                width = 8,
                columnDefs = list(list(width = '100px', targets = "_all")),
                columnDefs = list(list(targets = which(names(dt_data) == 'gene_name') - 1, searchable = TRUE))
            ),
            class = 'cell-border stripe compact',
            rownames = FALSE,
            selection = 'single'  # Enable single row selection
        )
    })
    
    # Observe clicks on the junction_table and update inputs
    observeEvent(input$junction_table_rows_selected, {
        # Check if a row is selected
        row_index <- input$junction_table_rows_selected
        req(row_index)
        
        # Get the data from the selected row
        selected_row <- possible_events()[row_index, ]
        selected_junction <- selected_row$junction_coords
        selected_gene <- selected_row$gene_name
        
        # Update the text inputs on the other tabs
        updateTextInput(session, "junction_input", value = selected_junction)
        updateTextInput(session, "gene_input", value = selected_gene)
        
        # Switch to the Junction Plotter tab
        updateTabsetPanel(session, "sidebar_tabs", selected = "plotter")
        
        # Simulate a click on the "Generate Plots" button to auto-plot
        # This is not a standard Shiny feature. A better approach is to use a reactive expression
        # that depends on a combination of inputs, including the row selection.
        # But if you must use a button, a better way is to move the plotting logic into a reactive
        # that triggers on a combination of events, not just the button.
    })
    
    # The rest of your server code, starting with the plot_button observer
    # The 'plot_button' observer is now the main entry point for plotting logic.
    # The `observeEvent(input$junction_table_rows_selected, ...)` should now
    # call this logic directly or through a reactive value.
    
    # Reactive expression for plotting based on a trigger
    plot_trigger <- reactive({
        list(
            button = input$plot_button,
            junction = input$junction_input,
            measure = input$measure_type
        )
    })
    
    # Now the main plotting logic depends on `plot_trigger`
    observeEvent(plot_trigger(), {
        # Show loading message
        showModal(modalDialog(
            title = "Loading...",
            "Generating plots, please wait...",
            footer = NULL,
            easyClose = FALSE
        ))
        
        req(input$junction_input)
        tryCatch({
            nygc_metadata <- fread("data/nygc/nygc_metadata.csv")
            rimod_metadata <- fread("data/rimod/rimod_meta.csv")
            ma_table =  fread("data/possible_junctions.csv")
            junction_name = ma_table %>%
                filter(junction_coords == input$junction_input) %>%
                distinct(gene_name) %>%
                pull(gene_name)
            
            if(nchar(junction_name) == 0){
                junction_name = ""
            }
            
            nygc_junction_data <- query_junction(junc = input$junction_input,
                                                 dataset_parquet = "data/nygc/parquets")
            nygc_junction_data = nygc_junction_data %>%
                distinct(chrom,sample,cluster_count,start,end,cluster,strand,junction_count,psi)
            
            rimod_junction_data <- query_junction(junc = input$junction_input,
                                                  dataset_parquet = "data/rimod/parquets")
            
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
    
    
    # Observe event for the gene expression tab button
    observeEvent(input$generate_plots, {
        req(input$gene_input)
        
        # Show loading message
        showModal(modalDialog(
            title = "Loading...",
            "Generating gene expression plots, please wait...",
            footer = NULL,
            easyClose = FALSE
        ))
        
        tryCatch({
            # Load all gene expression data at once
            rna_expression_nmd <- arrow::read_parquet("data/gene_expression/fractionation_de.parquet")
            sh_exp <- arrow::read_parquet("data/gene_expression/shsy5y_curve.parquet")
            sk_exp <- arrow::read_parquet("data/gene_expression/sknbe2_curve.parquet")
            upf1_rna <- arrow::read_parquet("data/gene_expression/upf1_rna_full.parquet")
            chx_rna <- arrow::read_parquet("data/gene_expression/chx_rna.parquet")
            
            # Store the raw data frames for download
            values$gene_expr_data$rna_expression_nmd <- rna_expression_nmd
            values$gene_expr_data$sh_exp <- sh_exp
            values$gene_expr_data$sk_exp <- sk_exp
            values$gene_expr_data$upf1_rna <- upf1_rna
            values$gene_expr_data$chx_rna <- chx_rna
            
            # Generate the plot for the specified gene
            gene_name_to_plot <- toupper(input$gene_input)
            
            # Check if the gene is present in at least one dataset before plotting
            if(any(gene_name_to_plot %in% unique(rna_expression_nmd$gene_name),
                   gene_name_to_plot %in% unique(sh_exp$symbol),
                   gene_name_to_plot %in% unique(sk_exp$symbol),
                   gene_name_to_plot %in% unique(upf1_rna$gene_name),
                   gene_name_to_plot %in% unique(chx_rna$gene_name))) {
                
                values$gene_expr_plot <- plot_expression(
                    this_gene  = gene_name_to_plot,
                    rna_expression_nmd = rna_expression_nmd,
                    sh_exp = sh_exp,
                    sk_exp = sk_exp,
                    upf1_rna = upf1_rna,
                    chx_rna = chx_rna
                )
            } else {
                showNotification(paste("Gene", gene_name_to_plot, "not found in any dataset."), type = "warning")
                values$gene_expr_plot <- NULL
            }
            
        }, error = function(e) {
            showNotification(paste("Error:", e$message), type = "error")
            values$gene_expr_plot <- NULL
        })
        
        removeModal()
    })
    
    # Render the new gene expression plot
    output$gene_expression_plot <- renderPlot({
        req(values$gene_expr_plot)
        values$gene_expr_plot
    })
    
    # Download handlers for the five gene expression datasets
    output$download_rna_expression_nmd <- downloadHandler(
        filename = function() { paste0("fractionation_de_", Sys.Date(), ".csv") },
        content = function(file) { req(values$gene_expr_data$rna_expression_nmd); write_csv(values$gene_expr_data$rna_expression_nmd, file) }
    )
    
    output$download_sh_exp <- downloadHandler(
        filename = function() { paste0("shsy5y_curve_", Sys.Date(), ".csv") },
        content = function(file) { req(values$gene_expr_data$sh_exp); write_csv(values$gene_expr_data$sh_exp, file) }
    )
    
    output$download_sk_exp <- downloadHandler(
        filename = function() { paste0("sknbe2_curve_", Sys.Date(), ".csv") },
        content = function(file) { req(values$gene_expr_data$sk_exp); write_csv(values$gene_expr_data$sk_exp, file) }
    )
    
    output$download_upf1_rna <- downloadHandler(
        filename = function() { paste0("upf1_rna_full_", Sys.Date(), ".csv") },
        content = function(file) { req(values$gene_expr_data$upf1_rna); write_csv(values$gene_expr_data$upf1_rna, file) }
    )
    
    output$download_chx_rna <- downloadHandler(
        filename = function() { paste0("chx_rna_", Sys.Date(), ".csv") },
        content = function(file) { req(values$gene_expr_data$chx_rna); write_csv(values$gene_expr_data$chx_rna, file) }
    )
}

# Run the application
shinyApp(ui = ui, server = server)