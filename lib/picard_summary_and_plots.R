# plot per target metric
# options include any metric in the picard output -
# MEAN_TARGET_COVERAGE | MAX_TARGET_COVERAGE | MIN_TARGET_COVERAGE | ZERO_CVG_TARGETS_PCT | MEDIAN_TARGET_COVERAGE
plot_per_target_metric<-function(coverage_df,metric='MEAN_TARGET_COVERAGE'){
  g <-
    coverage_df %>% ggplot(aes_string(x = 'group_name', y = metric, fill ='group_name')) +
    geom_boxplot(alpha = 0.6) + 
    coord_flip()+
    facet_wrap(BAIT_SET ~ .) +
    scale_fill_nejm() +
    theme_minimal()+labs(fill=NULL,title=metric,y=NULL,x=NULL)+
    theme(legend.position = 'top',text = element_text(size=20),strip.text = element_text(size=12),axis.text.x=element_text(size=12))
  if (grepl('PCT',metric)){g<-g+scale_y_continuous(labels = scales::percent)}
  g
  return(g)
}

# Plot percent covered above X
plot_percent_covered_above_x <- function(coverage_df,threshs=c(1,2,10,20,30,40,50,100)) {
  threshs_names<-paste0('PCT_TARGET_BASES_',threshs,'X')
  # keep only threshs present in the data
  threshs_names<-intersect(threshs_names,colnames(coverage_df))
  names(threshs_names)<-as.numeric(stringr::str_extract(threshs_names,'\\d+'))
  to_plot<-coverage_df%>%select(group_name,BAIT_SET,threshs_names)%>%
    pivot_longer(-c(group_name,BAIT_SET))%>%
    mutate(name=factor(glue('{name}X')))%>%
    mutate(name=forcats::fct_relevel(name,paste0(names(threshs_names),'X')))#%>%
    #mutate(name=glue('{name}X'))
    # mutate(name=stringr::str_replace(name,'PCT_TARGET_BASES_',''))#%>%
    # mutate(name=forcats::fct_relevel(name,'1X','2X','10X','20X','30X','40X','50X','100X'))
  
  g<-to_plot%>%#filter(BAIT_SET=='refseq_hg19_select_cds_20220307')%>%
    ggplot(aes(x=name,y=value,fill=group_name))+
    geom_boxplot()+facet_wrap(BAIT_SET~.)+
    scale_y_continuous(labels = scales::percent)+
    scale_fill_nejm()+
    labs(x='',y='Percent bases above threshold',fill='')+
    theme_minimal()+
    theme(legend.position = 'top',
          axis.text = element_text(size=15),
          legend.text = element_text(size=15),
          strip.text = element_text(size=15),
          axis.title = element_text(size=15),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          plot.background = element_rect(color='white'))
  g
  return(g)
}

# generate plots for basic params
plot_basic_seq_and_map_metrics <- function(coverage_df,
                                           basic_params=c('TOTAL_READS','PF_UNIQUE_READS','PF_UQ_READS_ALIGNED','PCT_EXC_MAPQ','PCT_EXC_DUPE','PCT_EXC_BASEQ','AT_DROPOUT','GC_DROPOUT')) {
  to_plot<-coverage_df%>%select(group_name,sample_name,BAIT_SET,basic_params)%>%
    pivot_longer(-c(group_name,sample_name,BAIT_SET))%>%
    mutate(sample_name=forcats::fct_reorder(sample_name,as.numeric(factor(group_name))))
  g<-to_plot%>%ggplot(aes(x=sample_name,y=value,fill=group_name))+
    geom_col(alpha=0.6)+facet_wrap(name~.,scales='free')+
    scale_fill_nejm()+coord_flip()+
    labs(x='',fill='')+
    theme_minimal()+
    theme(legend.position = 'top',
          plot.background = element_rect(color='white'),
          strip.text = element_text(size=12),
          axis.text = element_text(size=12),
          legend.text = element_text(size=14))
  g
  return(g)
}
