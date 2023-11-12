library(shiny)
library(DT)
library(ggplot2)
library(ProjectTemplate)
library(glue)
library(dplyr)
library(stringr)
library(tidyr)
library(ggsci)
library(shinyFiles)

project_folder<-'/media/SSD/Bioinformatics/Projects/gene_coverage_analysis_2022/'
sapply(list.files(glue('{project_folder}/lib'),full.names = T),source,.GlobalEnv)
# Assuming you have all the necessary functions like get_all_vcf_stats_with_stratifications, plot_vcfstats_results, collect_happy_results, generate_sample_table, happy_summary_per_group, etc. loaded

ui <- fluidPage(
  titlePanel("Coverage analysis report"),
  mainPanel(
    div(h3('Select coverage output folder'),
        shinyDirButton('folder', 'Select a folder', 'Please select a folder', FALSE)),
    br(),
    tabsetPanel(
      id = "sectionTab",
      tabPanel("Analysis setup", 
               div(h2('Analysis parameters'),DTOutput("analysis_setup"))),
      tabPanel("Coverage results: All samples", 
               div(# UI for selecting columns
                 selectInput("selected_columns_per_sample", "Select Columns",
                             choices = NULL,
                             selected = NULL,
                             multiple = TRUE),
                 DTOutput("all_samples_table"), style = 'width:1800px;')),
      tabPanel("Coverage results: Aggregated per group", 
               div(br(),
                   selectInput("selected_columns_aggregated", "Select Columns",
                               choices = NULL,
                               selected = NULL,
                               multiple = TRUE),
                   br(),
                   DTOutput("aggregated_table"), 
                   style = 'width:1800px;')),
      tabPanel("Coverage results: Plots", 
               div(br(),
                   selectInput("selected_groups_to_plot", "Select groups",
                               choices = NULL,
                               selected = NULL,
                               multiple = TRUE),
                   selectInput("selected_metric_to_plot", "Select metric to plot",
                               choices = NULL,
                               selected = NULL,
                               multiple = FALSE),
                   br(),
                   plotOutput("metric_plot",height='100%'), 
                   style = 'height:800px;width:1200px;'))
    )
  )
)

server <- function(input, output, session) {
  setwd(project_folder)
  input_folder<-NULL
  shinyDirChoose(input, 'folder', roots=c(wd=glue('{project_folder}/output')), filetypes=c('', 'txt'))
  observeEvent(input$folder,{
    input_folder<-parseDirPath(c(wd=glue('{project_folder}/output')),selection = input$folder)
    if (length(input_folder)!=0){
      main_output_folder<-input_folder
      print(main_output_folder)
      # Parse happy results
      message(glue('Parsing coverage analysis results from {main_output_folder}'))
      analysis_setup_file <- grep('*.analysis_setup*',list.files(main_output_folder),ignore.case = T,value = T)
      picard_coverage_analysis_results_file<-grep('picard_coverage_analysis_results*',list.files(main_output_folder),ignore.case = T,value = T)
      picard_coverage_analysis_results<-readr::read_delim(glue('{main_output_folder}/{picard_coverage_analysis_results_file}'), delim = '\t')
      # Analysis params ####
      analysis_setup<-readr::read_delim(glue('{main_output_folder}/{analysis_setup_file}'), delim = '\t')
      output$analysis_setup <- renderDT({
        DT::datatable(analysis_setup,selection ='none', extensions = c('Buttons'),
                      options = list(dom = 'Bfrtip',scrollX = T,buttons = c('copy', 'csv', 'excel')), filter = list(position = 'top', clear = FALSE))
      })
      # All samples ####
      updateSelectInput(inputId = 'selected_columns_per_sample',
                        choices = colnames(picard_coverage_analysis_results),
                        selected = colnames(picard_coverage_analysis_results)[1:3])
      message(glue('Producing all samples output..'))
      output$all_samples_table <- renderDT({
        DT::datatable(picard_coverage_analysis_results%>%select(input$selected_columns_per_sample),
                      extensions = 'Buttons',
                      options=list(scrollX=T,dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel'),pageLength =25),
                      filter = list(position = 'top', clear = FALSE))
      })
      
      # Aggregated data ####
      message(glue('Producing aggregated output..'))
      agg_coverage_df<-picard_coverage_analysis_results%>%
        group_by(group_name,BAIT_SET)%>%
        summarize(across(where(is.numeric),mean))%>%
        mutate(group_name=factor(group_name),BAIT_SET=factor(BAIT_SET))
      updateSelectInput(inputId = 'selected_columns_aggregated',
                        choices = colnames(agg_coverage_df),
                        selected = c('group_name','BAIT_SET','PCT_TARGET_BASES_20X'))
      
      output$aggregated_table <- renderDT({
        DT::datatable(agg_coverage_df%>%select(input$selected_columns_aggregated),
                      extensions = 'Buttons',
                      options=list(scrollX=T,dom = 'Bfrtip',buttons = c('copy', 'csv', 'excel'),pageLength =25),
                      filter = list(position = 'top', clear = FALSE))  
      })
      
      # Plots ####
      updateSelectInput(inputId = 'selected_groups_to_plot',
                        choices = unique(picard_coverage_analysis_results$group_name),
                        selected = unique(picard_coverage_analysis_results$group_name))
      numeric_cols<-picard_coverage_analysis_results%>%select(where(is.numeric))%>%colnames()
      excluded_cols<-c('BAIT_TERRITORY','BAIT_DESIGN_EFFICIENCY','HS_LIBRARY_SIZE','GENOME_SIZE','TARGET_TERRITORY')
      col_options<-setdiff(numeric_cols,excluded_cols)
      updateSelectInput(inputId = 'selected_metric_to_plot',
                        choices = col_options,
                        selected = 'PCT_TARGET_BASES_20X')
      toListen <- reactive({
        list(input$selected_metric_to_plot,input$selected_groups_to_plot)
      })
      observeEvent(toListen(),{
        if (input$selected_metric_to_plot!=''){
          p<-plot_per_target_metric(picard_coverage_analysis_results%>%filter(group_name%in%input$selected_groups_to_plot),input$selected_metric_to_plot)
          output$metric_plot<-renderPlot({p})
        }
      })
    }
  })
  
  
}

shinyApp(ui, server)
