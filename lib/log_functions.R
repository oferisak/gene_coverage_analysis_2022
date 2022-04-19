add_to_log<-function(log_file,to_log_title,to_log_data){
  write(glue('## {to_log_title}'),file=log_file,append = T)
  if (class(to_log_data)%in%c('data.frame','tibble')){
    write.table(to_log_data,file=log_file,append = T,row.names = F,sep='\t')
  }else{
    write(to_log_data,file=log_file,append = T)
  }
}
