
collect_per_group_coverage_metrics <- function(per_transcript_data) {
  per_transcipt_cov_perp_summary<-per_transcript_data%>%
    select(group_name,sample_name,transcript_len,all_of(grep('^above_',colnames(per_transcript_data),value=T)))%>%
    pivot_longer(-c(group_name,sample_name,transcript_len))%>%
    mutate(name=as.numeric(stringr::str_replace(name,'above_','')))
  
  per_transcipt_cov_summary<-per_transcipt_cov_perp_summary%>%
    group_by(group_name,name,sample_name)%>%
    summarize(sum_above=sum(value),
              perc_above=sum_above/sum(transcript_len))%>%
    ungroup()%>%group_by(group_name,name)%>%
    summarize(n=n(),
              mean_perc_above=mean(perc_above),
              sd_perc_above=sd(perc_above),
              se_perc_above=sd_perc_above/sqrt(n))
  return(per_transcipt_cov_summary)
}

plot_per_group_cov <- function(per_transcript_data,title=NULL,with_table=T) {
  per_transcipt_cov_summary<-collect_per_group_coverage_metrics(per_transcript_data)
  g<-per_transcipt_cov_summary%>%ggplot(aes(x=factor(name,labels = paste0('>x',unique(per_transcipt_cov_summary$name))),y=mean_perc_above,fill=group_name))+
    geom_col(position='dodge',alpha=0.7)+
    scale_fill_nejm()+
    labs(x='Threshold',y='Percent above threshold',fill='Group',title=title)+
    scale_y_continuous(labels = scales::percent)+
    theme_minimal() +
    theme(legend.position = 'top')
  print(g)
  summary_table<-per_transcipt_cov_summary%>%
    select(group_name,name,n,mean_perc_above)%>%
    mutate(mean_perc_above=glue('{100*round(mean_perc_above,4)}%'))%>%
    pivot_wider(id_cols = name,names_from = group_name,values_from = mean_perc_above)%>%
    rename('Threshold'=name)
  # tg<-per_transcript_data%>%select(group_name,all_of(grep('perc_above_',colnames(per_transcript_data),value=T)))%>%
  #   gtsummary::tbl_summary(by = group_name,
  #                          statistic = gtsummary::all_continuous() ~ "{mean} ({sd})")%>%
  #   gtsummary::as_tibble()
  if (with_table){
    final_plot<-g/gridExtra::tableGrob(summary_table,theme=gridExtra::ttheme_minimal())
  }else{final_plot<-g}
  print(final_plot)
  #g/gridExtra::tableGrob(tg)
  return(final_plot)
}
