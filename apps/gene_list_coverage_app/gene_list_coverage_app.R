library(shiny)
library(shinydashboard)
library(DT)
library(tidyr)
library(reactable)
library(dplyr)
library(dashboardthemes)
library(xlsx)
library(stringr)
library(glue)
library(plotly)
library(ggsci)
# read references
reference_files<-list.files('./coverage_references/',full.names = T)
# read accessory data
segdup_intersect_data<-readr::read_delim('./accessory_data/refseq_hg19_curated_cds_vs_segdup.csv.gz')
# filter duplications that are less than 10% of the exon
segdup_intersect_data<-segdup_intersect_data%>%filter(overlap_perc>0.1)

# read data
read_reference_data<-function(group_name_by_file_name=F){
    reference_data<-list('per_transcript_reference_coverage_db'=NULL,
                         'per_exon_reference_coverage_db'=NULL)
    #per_transcript_reference_coverage_db<-NULL
    #per_exon_reference_coverage_db<-NULL
    for (file_name in reference_files){
        message(glue('Parsing {file_name}'))
        if(grepl('per_transcript',file_name)){
            group_per_trascript<-readr::read_delim(file_name)
            if (group_name_by_file_name){
                group_name_from_file<-stringr::str_replace(basename(file_name),'_.+','')
                group_per_trascript$group_name<-group_name_from_file
            }
            reference_data[['per_transcript_reference_coverage_db']]<-reference_data[['per_transcript_reference_coverage_db']]%>%
                bind_rows(group_per_trascript)
        }else{
            group_per_exon<-readr::read_delim(file_name)
            if (group_name_by_file_name){
                group_name_from_file<-stringr::str_replace(basename(file_name),'_.+','')
                group_per_exon$group_name<-group_name_from_file
            }
            reference_data[['per_exon_reference_coverage_db']]<-reference_data[['per_exon_reference_coverage_db']]%>%
                bind_rows(group_per_exon)
        }
    }
    return(reference_data)
}
reference_data<-read_reference_data(group_name_by_file_name = T)

per_transcript_reference_coverage_db<-reference_data[['per_transcript_reference_coverage_db']]%>%
    mutate(coverage_cat=forcats::fct_relevel(factor(coverage_cat),'0x','1-20x','20-50x','50-100x','>100x'),
           above_20x=ifelse(coverage_cat%in%c('0x','1-20x'),F,T))
per_exon_reference_coverage_db<-reference_data[['per_exon_reference_coverage_db']]%>%
    mutate(coverage_cat=forcats::fct_relevel(factor(coverage_cat),'0x','1-20x','20-50x','50-100x','>100x'),
           above_20x=ifelse(coverage_cat%in%c('0x','1-20x'),F,T))%>%
    left_join(segdup_intersect_data%>%select(transcript_name,exon_num,chr,dup_perc))%>%
    mutate(dup_level_above_98_perc=ifelse(!is.na(dup_perc) & dup_perc>0.98,T,F))

# Define UI ####
ui <- dashboardPage(dashboardHeader(title = sprintf('Gene Coverage Calculator'),
                                    titleWidth = 290), 
                    dashboardSidebar(width = 290,
                                     sidebarMenu(menuItem("Gene List Upload",
                                                          tabName = "gene_list_upload", 
                                                          icon = icon('search')),
                                                 menuItem("Per Transcript Coverage",
                                                          tabName = "per_transcript_cov", 
                                                          icon = icon('dna')),
                                                 menuItem("Per Transcript Coverage Plots",
                                                          tabName = "per_transcript_cov_plots", 
                                                          icon = icon('dna')),
                                                 menuItem("Per Exon Coverage",
                                                          tabName = "per_exon_cov", 
                                                          icon = icon('dna')),
                                                 menuItem("Disclaimer table",
                                                          tabName = "disclaimer_table", 
                                                          icon = icon('dna'))
                                     )
                    ),
                    
                    dashboardBody(
                        # stlye definitions from https://stackoverflow.com/questions/52198452/how-to-change-the-background-color-of-the-shiny-dashboard-body
                        
                        
                        tags$head(tags$style(HTML('
                                /* logo */
                                .skin-blue .main-header .logo {
                                background-color: #87a86f;
                                font-weight: 800;
                                }
                                /* navbar (rest of the header) */
                                .skin-blue .main-header .navbar {
                                background-color: #87a86f;
                                }
                                /* body */
                                .content-wrapper, .right-side {
                                background-color: #f5fbf8;
                                }
                                /* main sidebar */
                                .main-sidebar {
                                font-size: 20px;
                                }
                                /* sidebar open */
                                .treeview-menu>li>a {
                                font-size: 18px!important;
                                }
                                .dataTables_scrollBody {
                                    transform:rotateX(180deg);
                                }
                                .dataTables_scrollBody table {
                                    transform:rotateX(180deg);
                                }
                                '))),
                        
                        #shinyDashboardThemes(theme = "poor_mans_flatly"),
                        tabItems(
                            # Panel Browser
                            tabItem('gene_list_upload',
                                    fluidRow(
                                        fileInput("gene_list_file", "Select gene list file",
                                                  multiple = FALSE,
                                                  accept = c("text/csv",
                                                             "text/comma-separated-values,text/plain",
                                                             ".csv",
                                                             "text/tsv",
                                                             ".tsv",'xls','xlsx')),
                                        div(DT::dataTableOutput(outputId = 'gene_list_table'), 
                                            style='font-size:125%;'))
                                    ),
                            
                            
                            tabItem('per_transcript_cov',
                                    div(selectizeInput('lib_prep_kit_transcript', label='Select Library Prep Kits', multiple=F,choices=unique(per_transcript_reference_coverage_db$group_name),width = 800,options= list(maxOptions = 5000)), 
                                        style='font-size:200%;'),
                                    div(DT::dataTableOutput(outputId = 'per_transcript_cov'), 
                                        style='font-size:100%;')),
                            tabItem('per_transcript_cov_plots',
                                    fluidPage(
                                        div(numericInput('select_below20x_thresh', label='Plot transcripts with this rate (0.0-1.0) below 20x', value=0.2,min = 0,max=1), 
                                            style='font-size:200%;'),
                                        div(uiOutput(outputId='per_transcript_cov_plot'),height = '400px')    
                                    )),
                                    #reactableOutput(outputId = 'per_transcript_cov')),
                            tabItem('per_exon_cov',
                                    div(selectizeInput('lib_prep_kit_exon', label='Select Library Prep Kits', multiple=F,choices=unique(per_transcript_reference_coverage_db$group_name),width = 800,options= list(maxOptions = 5000)), 
                                        style='font-size:200%;'),
                                    div(DT::dataTableOutput(outputId = 'per_exon_cov'), 
                                        style='font-size:100%;')),
                            tabItem('disclaimer_table',
                                    div(selectizeInput('lib_prep_kit_disclaimer', label='Select Library Prep Kits', multiple=F,choices=unique(per_transcript_reference_coverage_db$group_name),width = 800,options= list(maxOptions = 5000)), 
                                        style='font-size:200%;'),
                                    div(DT::dataTableOutput(outputId = 'disclaimer_table'), 
                                        style='font-size:100%;'))
                        ),
                    )
                    #tags$head(tags$style(HTML('* {font-family: "Kirnberg"};')))
)



# Define server logic required to draw a histogram
server <- function(input, output) {
    # per transcript ####
    observeEvent(input$gene_list_file,{
        if (grepl('xls',input$gene_list_file$datapath)){
            gene_list<-readxl::read_excel(input$gene_list_file$datapath,col_names = F)
            colnames(gene_list)[1]<-'gene_symbol'
        }else{
            gene_list<-readr::read_delim(input$gene_list_file$datapath,col_names = c('gene_symbol'))
        }
        # verify gene names
        gene_list<-gene_list%>%mutate(gene_found_in_reference=gene_symbol%in%per_transcript_reference_coverage_db$gene_name)
        gene_symbols<- gene_list%>%pull(gene_symbol)
        output$gene_list_table<-DT::renderDataTable(gene_list%>%select(gene_symbol,gene_found_in_reference),
                                                    options=list(scrollX=F,pageLength=10),
                                                    selection='single',
                                                    rownames= FALSE,
                                                    filter = list(position = 'top', clear = FALSE))
        observeEvent(input$lib_prep_kit_transcript,{
            per_transcript_cov<-per_transcript_reference_coverage_db%>%
                filter(gene_name%in%gene_list$gene_symbol,
                       group_name==input$lib_prep_kit_transcript)
            
            good_cov_transcripts<-per_transcript_cov%>%
                group_by(gene_name,transcript_name,chr)%>%
                summarize(perc_above_20x=sum(ifelse(above_20x,perc_at_coverage,0)))
            
            output$per_transcript_cov<-DT::renderDataTable(good_cov_transcripts%>%select(-c(chr)),
                                                        options=list(scrollX=F,pageLength=50),
                                                        selection='single',
                                                        rownames= FALSE,
                                                        filter = list(position = 'top', clear = FALSE))
            
            observeEvent(input$select_below20x_thresh,{
                
                to_plot_data<-per_transcript_cov%>%
                    filter(gene_name%in%(good_cov_transcripts%>%
                                             filter((1-perc_above_20x)>input$select_below20x_thresh)%>%
                                             pull(gene_name)))
                
                plot_coverage_data<-function(gene_name_to_plot){
                    
                    cov_plot<-to_plot_data%>%filter(gene_name==gene_name_to_plot)%>%
                        ggplot(aes(x=transcript_name,y=perc_at_coverage,fill=coverage_cat))+
                        geom_col(alpha=0.7)+
                        coord_flip()+
                        scale_fill_manual(values=c('0x'='darkred','1-20x'='orange','20-50x'='darkolivegreen3','50-100x'='chartreuse4','>100x'='darkgreen'))+
                        theme_minimal()+theme(legend.position = 'top')+
                        labs(x='Percent at coverage',y=NULL)
                    box(title = gene_name_to_plot,renderPlot(cov_plot))
                }
                
                output$per_transcript_cov_plot<-renderUI({
                    lapply(to_plot_data%>%pull(gene_name)%>%unique(), plot_coverage_data)
                })
            })
            
            #output$per_transcript_cov_plot<-renderPlotly(cov_plot)
            
        })
        
        # Per exon tab
        observeEvent(input$lib_prep_kit_exon,{
            per_exon_cov<-per_exon_reference_coverage_db%>%
                filter(gene_name%in%gene_list$gene_symbol,
                       group_name==input$lib_prep_kit_exon)

            good_cov_exons<-per_exon_cov%>%
                group_by(gene_name,transcript_name,exon_num,chr,dup_level_above_98_perc)%>%
                summarize(perc_above_20x=sum(ifelse(above_20x,perc_at_coverage,0)))

            output$per_exon_cov<-DT::renderDataTable(good_cov_exons%>%select(-c(chr)),
                                                           options=list(scrollX=F,pageLength=50),
                                                           selection='single',
                                                           rownames= FALSE,
                                                           filter = list(position = 'top', clear = FALSE))

        
        })
        # Disclaimer table tab
        
        observeEvent(input$lib_prep_kit_disclaimer,{
            
            per_transcript_cov<-per_transcript_reference_coverage_db%>%
                filter(gene_name%in%gene_list$gene_symbol,
                       group_name==input$lib_prep_kit_transcript)
            
            good_cov_transcripts<-per_transcript_cov%>%
                group_by(gene_name,transcript_name,chr)%>%
                summarize(perc_above_20x=sum(ifelse(above_20x,perc_at_coverage,0)))
            
            per_exon_cov<-per_exon_reference_coverage_db%>%
                filter(gene_name%in%gene_list$gene_symbol,
                       group_name==input$lib_prep_kit_exon)
            
            bad_cov_exons<-per_exon_cov%>%
                group_by(gene_name,transcript_name,exon_num,chr)%>%
                summarize(perc_above_20x=sum(ifelse(above_20x,perc_at_coverage,0)))%>%
                filter(perc_above_20x<0.9)%>%ungroup()%>%
                group_by(gene_name,transcript_name,chr)%>%
                summarize(exons_with_more_than_10perc_below_20x=paste0(unique(exon_num),collapse=', '))
                
            high_dup_exons<-per_exon_cov%>%
                filter(dup_level_above_98_perc)%>%
                group_by(gene_name,transcript_name,chr)%>%
                summarize(exons_with_dup_level_above_98_perc=paste0(unique(exon_num),collapse=', '))
            
            disclaimer_table<-good_cov_transcripts%>%
                left_join(bad_cov_exons)%>%
                left_join(high_dup_exons)%>%
                mutate(in_disclaimer=factor(ifelse(!(is.na(exons_with_more_than_10perc_below_20x)&is.na(exons_with_dup_level_above_98_perc)),T,F)))%>%
                rename('Add to disclaimer'=in_disclaimer,
                       'Exons with >10% below 20x (zero-based)'=exons_with_more_than_10perc_below_20x,
                       'Exons with >10% with >98% homology (zero-based)'=exons_with_dup_level_above_98_perc,
                       'Percent above 20x'=perc_above_20x,
                       'Gene name'=gene_name,
                       'Transcript name'=transcript_name)
            output$disclaimer_table<-DT::renderDataTable(DT::datatable(disclaimer_table%>%select(-c(chr)),
                                                     options=list(scrollX=F,pageLength=50,
                                                                  buttons = c('csv', 'excel'),
                                                                  dom = 'Bfrtip'),
                                                     extensions = 'Buttons',
                                                     selection='single',
                                                     rownames= FALSE,
                                                     filter = list(position = 'top', clear = FALSE))%>%
                                                         formatPercentage(c('Percent above 20x')))
            
            
        })
        
    })
    
}
#  disable = list(columns = 1:(ncol(final_panel_output)-1)
# Run the application 
shinyApp(ui = ui, server = server)
