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

#gene_list<-readr::read_delim('/media/SSD/Bioinformatics/Projects/gene_coverage_analysis_2022/data/gene_lists/LVNC.csv',col_names = c('gene_symbol','confidence'))
# read refseq names
refseq_to_gene<-readr::read_delim('./accessory_data/refseq_to_gene.txt')
# remove bad chr
refseq_to_gene<-refseq_to_gene%>%filter(!grepl('_',chrom))

# fix gene names to be the approved ones
standard_genenames_table<-readr::read_delim('./accessory_data/standard_genenames_db_2022-04-07.csv')

# read accessory data
# filter duplications that are less than 10% of the exon
segdup_intersect_data<-readr::read_delim('./accessory_data/refseq_hg19_curated_cds_vs_segdup.csv.gz')
segdup_intersect_data<-segdup_intersect_data%>%
    filter(overlap_perc>0.1)%>%
    #mutate(exon_num=exon_num)%>%
    filter(!grepl('_',chr))

refseq_curated<-readr::read_delim('./accessory_data/refseq_hg19_curated_exons_20220410.bed.gz',
                                  col_names = c('chr','start','end','interval_name','zero','strand'))%>%
    filter(!grepl('_',chr))%>%
    mutate(transcript_name=stringr::str_extract(interval_name,'(.+)_(exon|cds)')%>%stringr::str_replace('_(exon|cds)',''),
           exon_num=as.numeric(stringr::str_extract(interval_name,'(exon|cds)_\\d+')%>%stringr::str_extract('\\d+')))
# add total number of exons
refseq_curated<-refseq_curated%>%left_join(refseq_curated%>%group_by(transcript_name,chr)%>%summarize(total_exons=n()))
refseq_with_segdup<-refseq_curated%>%left_join(segdup_intersect_data%>%select(-c(start,end)),by=c('transcript_name','chr','exon_num'))

# Define UI ####
ui <- dashboardPage(dashboardHeader(title = sprintf('Disclaimer Generator'),
                                    titleWidth = 290), 
                    dashboardSidebar(width = 290,
                                     sidebarMenu(menuItem("Gene List Upload",
                                                          tabName = "gene_list_upload", 
                                                          icon = icon('search')),
                                                 menuItem('Segmental dups plots',
                                                          tabName = "segdup_plots", 
                                                          icon = icon('search')),
                                                 menuItem("Disclaimer table",
                                                          tabName = "disclaimer_table", 
                                                          icon = icon('dna')),
                                                 menuItem("Disclaimer text",
                                                          tabName = "disclaimer_text", 
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
                                                  accept = c(".csv",".tsv",'.xls','.xlsx')),
                                        div(DT::dataTableOutput(outputId = 'gene_list_table'), 
                                            style='font-size:125%;'))
                                    ),
                            

                            tabItem('segdup_plots',
                                    fluidPage(
                                        div(selectizeInput('lib_prep_kit_plots', label='Select Refseq Set', multiple=F,choices=c('select','curated'),width = 800,options= list(maxOptions = 5000)), 
                                            style='font-size:200%;'),
                                        div(sliderInput('segdup_plots_thresh', label='Select Minimal Duplication Level',value = 0.98, min = 0.9,max=1,width = 800), 
                                            style='font-size:200%;'),
                                        div(uiOutput(outputId='per_transcript_cov_plot'),height = '400px')
                                    )),
                                    #reactableOutput(outputId = 'per_transcript_cov')),

                            tabItem('disclaimer_table',
                                    div(selectizeInput('lib_prep_kit_disclaimer_table', label='Select Refseq Set', multiple=F,choices=c('select','curated'),width = 800,options= list(maxOptions = 5000)), 
                                        style='font-size:200%;'),
                                    div(sliderInput('disclaimer_segdup_thresh', label='Select Minimal Duplication Level',value = 0.98, min = 0.9,max=1,width = 800), 
                                        style='font-size:200%;'),
                                    div(DT::dataTableOutput(outputId = 'disclaimer_table'), 
                                        style='font-size:100%;')),
                            tabItem('disclaimer_text',
                                    div(textOutput('disclaimer_text'), 
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
        
        # fix gene names
        # first check if any of the provided gene symbols is not found in the refseq gene name table
        gene_list_genes_not_in_refseq<-gene_list%>%
            filter(!gene_symbol%in%refseq_to_gene$gene_symbol)%>%
            select(provided_gene_symbol=gene_symbol)%>%
            left_join(standard_genenames_table,by=c('provided_gene_symbol'='prev_symbol'))%>%
            filter(!is.na(approved_gene_symbol))# remove those that are not found in both refseq and the standard gene table (bad genes)
        # then, for the genes that are found in the standard gene names table, change the gene symbol to the approved name 
        gene_list<-gene_list%>%
            mutate(gene_name_provided=gene_symbol)%>%
            rowwise()%>%
            mutate(gene_symbol=ifelse(gene_symbol%in%gene_list_genes_not_in_refseq$provided_gene_symbol,
                                      gene_list_genes_not_in_refseq%>%filter(provided_gene_symbol==gene_symbol)%>%slice(1)%>%pull(approved_gene_symbol),#if there are more than one approved gene symbol for the previous symbol, pick the first one.
                                      gene_symbol))
        
        #gene_list<-gene_list%>%mutate(gene_found_in_reference=gene_symbol%in%refseq_to_gene$gene_name)
        gene_list_annotated<-gene_list%>%select(gene_symbol,gene_name_provided)%>%left_join(refseq_to_gene)
        #fix input gene names to match refseq
        gene_symbols<- gene_list%>%pull(gene_symbol)
        output$gene_list_table<-DT::renderDataTable(gene_list_annotated,
                                                    options=list(scrollX=F,pageLength=10),
                                                    selection='single',
                                                    rownames= FALSE,
                                                    filter = list(position = 'top', clear = FALSE))
        # add transcript names
        # join gene list with segdup
        genes_seg_data<-gene_list_annotated%>%
            inner_join(segdup_intersect_data)
        
        observeEvent(input$lib_prep_kit_plots,{
            observeEvent(input$segdup_plots_thresh,{
                
                segdup_thresh<-input$segdup_plots_thresh
                transcripts_with_dups<-refseq_with_segdup%>%filter(dup_perc>segdup_thresh)%>%pull(transcript_name)%>%unique()
                to_plot<-refseq_with_segdup%>%
                    mutate(dup_level=cut(dup_perc,breaks=c(0,segdup_thresh,1),labels=c(glue('<{segdup_thresh}'),glue('>{segdup_thresh}'))))%>%
                    mutate(dup_level=factor(ifelse(is.na(dup_level),glue('<{segdup_thresh}'),as.character(dup_level))))%>%
                    filter(transcript_name %in% transcripts_with_dups)%>%
                    left_join(refseq_to_gene)%>%
                    filter(transcript_name%in%(gene_list_annotated$transcript_name))
                if (input$lib_prep_kit_plots=='select'){
                    to_plot<-to_plot%>%filter(cannonical)
                }
                
                plot_coverage_data<-function(gene_name_to_plot){
                    cov_plot<-to_plot%>%filter(gene_symbol==gene_name_to_plot)%>%
                        ggplot(aes(x=transcript_name,xend=transcript_name,y=start,yend=end,col=dup_level))+
                        geom_segment(size=10)+
                        coord_flip()+
                        scale_color_manual(values=c('gray','darkred'),drop=F)+
                        theme_minimal()+
                        labs(x=NULL,y=NULL,color='Duplication level')+
                        theme(legend.position = 'bottom')
                    
                    box(title = gene_name_to_plot,renderPlot(cov_plot))
                }

                output$per_transcript_cov_plot<-renderUI({
                    lapply(to_plot%>%arrange(gene_symbol)%>%pull(gene_symbol)%>%unique(), plot_coverage_data)
                })
            })

            #output$per_transcript_cov_plot<-renderPlotly(cov_plot)

        })
        
        # Disclaimer table tab
        
        observeEvent(input$lib_prep_kit_disclaimer_table,{
          observeEvent(input$disclaimer_segdup_thresh,{
            disclaimer_segdup_thresh<-input$disclaimer_segdup_thresh
            high_dup_exons<-genes_seg_data%>%
                filter(dup_perc>disclaimer_segdup_thresh)%>%
                group_by(gene_symbol,transcript_name)%>%
                summarize(exons_with_dup_level_above_thresh=paste0(unique(exon_num+1),collapse=', '))
            
            disclaimer_table<-gene_list_annotated%>%
                left_join(refseq_curated%>%select(transcript_name,total_exons)%>%distinct())%>%
                left_join(high_dup_exons)%>%
                mutate(in_disclaimer=factor(ifelse(!(is.na(exons_with_dup_level_above_thresh)),T,F)))
            
            if (input$lib_prep_kit_disclaimer_table=='select'){
                disclaimer_table<-disclaimer_table%>%filter(cannonical)
            }
            
            print(disclaimer_table)
            disclaimer_text<-disclaimer_table%>%filter(in_disclaimer==TRUE)%>%
                mutate(as_text=glue('{gene_symbol}:{transcript_name}'))%>%
                mutate(as_text=ifelse(!is.na(exons_with_dup_level_above_thresh),glue('{as_text}:Exons with >10% with >{disclaimer_segdup_thresh*100}% homology:{exons_with_dup_level_above_thresh}'),as_text))
            
            dup_column_name<-glue('Exons with >10% with >{disclaimer_segdup_thresh*100}% homology')
            disclaimer_table<-disclaimer_table%>%
                rename('Add to disclaimer'=in_disclaimer,
                       !!sym(dup_column_name):=exons_with_dup_level_above_thresh,
                       'Gene name'=gene_symbol,
                       'Transcript name'=transcript_name)
            output$disclaimer_table<-DT::renderDataTable(DT::datatable(disclaimer_table,
                                                     options=list(scrollX=F,pageLength=nrow(disclaimer_table),
                                                                  buttons = c('csv', 'excel'),
                                                                  dom = 'Bfrtip'),
                                                     extensions = 'Buttons',
                                                     selection='single',
                                                     rownames= FALSE,
                                                     filter = list(position = 'top', clear = FALSE)))
                
            output$disclaimer_text<-renderText(paste0(disclaimer_text%>%pull(as_text),collapse=' | '))
            
          })
        })
        
    })
    
}
#  disable = list(columns = 1:(ncol(final_panel_output)-1)
# Run the application 
shinyApp(ui = ui, server = server)
