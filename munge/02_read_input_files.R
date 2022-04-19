read_input_bam_files<-function(){
  if (analysis_or_report=='report'){return()}
  if (is.na(bam_folder)){
    bam_folder<-'./data/bam_analysis/input'
  }
  group_dirs<-list.dirs(glue('{bam_folder}'),recursive = F)
  input_bams_df<-NULL
  for (group_path in group_dirs){
    group_name<-basename(group_path)
    sample_path<-list.files(group_path,pattern = '*.bam$', recursive = T,full.names = T)
    if (length(sample_path)==0){
      message(glue('no bam files in {group_path}, skipping'))
      next
    }
    sample_name<-sample_path%>%basename()%>%stringr::str_replace('.bam','')
    input_bams_df<-input_bams_df%>%bind_rows(data.frame(group_name=group_name,
                                                        group_path=group_path,
                                                        sample_name,
                                                        sample_path))
  }
  return(input_bams_df)
}

input_bams_df<-read_input_bam_files()

read_input_bed_files<-function(){
  group_dirs<-list.dirs('./data/bed_intersect_analysis/input',recursive = F)
  input_beds_df<-NULL
  for (group_path in group_dirs){
    group_name<-basename(group_path)
    sample_path<-list.files(group_path,pattern = '*.bed$', recursive = T,full.names = T)
    sample_name<-sample_path%>%basename()%>%stringr::str_replace('.bed','')
    input_beds_df<-input_beds_df%>%bind_rows(data.frame(group_name=group_name,
                                                        group_path=group_path,
                                                        sample_name,
                                                        sample_path))
  }
  return(input_beds_df)
}
input_beds_df<-read_input_bed_files()
