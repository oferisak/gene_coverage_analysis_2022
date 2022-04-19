picard_path<-'/media/SSD/Bioinformatics/Tools/Picard/picard.jar'

run_BedToIntervalList<-function(input_bed,output_folder,reference_dict){
  output_intervals_file_name<-stringr::str_replace(basename(input_bed),'\\.bed','\\.intervals')
  output_intervals<-glue('{output_folder}/{output_intervals_file_name}')
  BedToIntervalList_command<-glue('java -jar {picard_path} BedToIntervalList -I "{input_bed}" -O {output_intervals} -SD {reference_dict}')
  system(BedToIntervalList_command)
  return(output_intervals)
}


run_CollectHsMetrics<-function(input_bam,output_metrics_file,reference_fasta,bait_intervals_bed,target_intervals_bed,per_target_file=NA){
  collectHsMetrics_command<-glue('java -jar {picard_path} CollectHsMetrics -I {input_bam} -O {output_metrics_file} -R {reference_fasta} -BAIT_INTERVALS {bait_intervals_bed} -TARGET_INTERVALS {target_intervals_bed}')
  if (!is.na(per_target_file)){
    collectHsMetrics_command<-glue('{collectHsMetrics_command} -PER_TARGET_COVERAGE {per_target_file}')
  }
  message(glue('running: {collectHsMetrics_command}..'))
  system(collectHsMetrics_command)
}


parse_picard_output<-function(picard_output){
  tmp <- readLines(picard_output)
  # start of each 'file'
  sof <- grep("#", tmp)
  
  # the actual start is one past that
  real_start <- sof + 1
  
  # figure out the end of each unique df
  real_end <- c(sof[-1] - 1, length(tmp))
  
  # the number of rows to read in
  to_read <- real_end - real_start
  to_read[to_read<0]<-1
  
  # a list to store your data.frames
  my_dfs <- vector("list", length = length(real_start))
  for(i in 1:length(my_dfs)){
    # suppressing warnings, as there are a lot 
    #  that come up when a column header is not
    #  in a specific data.frame
    my_dfs[[i]] <- data.table::fread(picard_output,
                                     sep = '\t',
                        skip = sof[i],
                        nrows = to_read[i],
                        fill = FALSE,
                        check.names = FALSE,
                        data.table = FALSE
    )
  }
  output<-list()
  output[['comand']]<-my_dfs[[1]]
  output[['picard_coverage_metrics']]<-my_dfs[[5]]
  output[['picard_coverage']]<-my_dfs[[6]]
  return(output)
}



# z<-all_pt_df%>%group_by(group,sample,transcript_name)%>%summarize(total_cds=n(),
#                                                cds_with_0x=sum(pct_0x>0),pct_cds_with_0x=cds_with_0x/n(),
#                                                cds_with_10pct_0x=sum(pct_0x>0.1),pct_cds_10pct_0x=cds_with_10pct_0x/n())%>%
#   ungroup()%>%group_by(group,sample)%>%summarise(n_transcripts=n(),
#                                                  n_transcripts_with_cds_not_fully_covered=sum(cds_with_0x>0),
#                                                  n_transcripts_with_most_cds_not_fully_covered=sum(pct_cds_with_0x>0.5))

parse_picard_output_folder<-function(picard_analysis_output_folder){
  coverage_analysis_results_file<-grep('analysis_results',list.files(picard_analysis_output_folder,full.names = T),value=T)
  if (length(coverage_analysis_results_file)>0){
    coverage_df<-readr::read_delim(coverage_analysis_results_file)
  }else{
    coverage_df<-NULL
    group_folders<-list.dirs(picard_analysis_output_folder,recursive = F)
    for (group_folder in group_folders){
      group_name<-basename(group_folder)
      sample_folders<-list.dirs(group_folder,recursive = F)
      for (sample_folder in sample_folders){
        sample_name<-basename(sample_folder)
        picard_cov_output_files<-grep('picard_cov_metrics',list.files(sample_folder,full.names = T),value=T)
        for (target_cov_file in picard_cov_output_files){
          target_cov<-parse_picard_output(target_cov_file)[['picard_coverage_metrics']]
          coverage_df<-coverage_df%>%bind_rows(target_cov%>%mutate(group_name=group_name,sample_name=sample_name,.before=1))
        }
      }
    }
  }
  coverage_df<-coverage_df%>%mutate(across(where(~class(.)=='integer64'),as.numeric))
  return(coverage_df)
}


parse_picard_per_target_output_folder<-function(picard_analysis_output_folder){
  coverage_analysis_results_file<-grep('analysis_results',list.files(picard_analysis_output_folder,full.names = T),value=T)
  
  coverage_df<-NULL
  group_folders<-list.dirs(picard_analysis_output_folder,recursive = F)
  for (group_folder in group_folders){
    group_name<-basename(group_folder)
    sample_folders<-list.dirs(group_folder,recursive = F)
    for (sample_folder in sample_folders){
      sample_name<-basename(sample_folder)
      picard_cov_output_files<-grep('picard_per_target_cov_metrics',list.files(sample_folder,full.names = T),value=T)
      for (target_cov_file in picard_cov_output_files){
        verify_ucsc_bed_based_coverage_analysis(target_cov_file)
        per_target_cov<-readr::read_delim(target_cov_file)
        coverage_df<-coverage_df%>%bind_rows(per_target_cov%>%mutate(group_name=group_name,sample_name=sample_name,.before=1))
      }
    }
  }
  coverage_df<-coverage_df%>%mutate(across(where(~class(.)=='integer64'),as.numeric))
  return(coverage_df)
}

# verify that the per-target coverage analysis was performed against a ucsc generated bed file 
verify_ucsc_bed_based_coverage_analysis<-function(per_target_output_file){
  per_target_output<-readr::read_delim(per_target_output_file)
  sample_name<-per_target_output%>%slice_sample(n=10)%>%pull(name)
  if (sum(grepl('^NM_',sample_name))!=length(sample_name)){stop(glue('The names in the target file {per_target_output_file} do not match the expected ucsc bed file format'))}
  cds_or_exon<-stringr::str_split(sample_name,'_',simplify = T)[,3]
  if (sum(grepl('cds|exon',cds_or_exon))!=length(sample_name)){stop(glue('The names in the target file {per_target_output_file} do not match the expected ucsc bed file format'))}
  
}

# ucsc bed per-transcript metrics
calculate_per_target_coverage<-function(coverage_df){
  # first split the target name into the different transcripts (picard joins same coordinates transcripts together)
  cov_data<-coverage_df%>%separate_rows(name,sep='\\|')
  # add the gene name to each transcript
  cov_data<-cov_data%>%mutate(transcript_name=stringr::str_extract(name,'(.+)_(exon|cds)')%>%stringr::str_replace('_(exon|cds)',''),
                              exon_num=as.numeric(stringr::str_extract(name,'(exon|cds)_\\d+')%>%stringr::str_extract('\\d+')))
  # join with gene name
  cov_data<-cov_data%>%left_join(refseq2gene)
  # get gene summary
  gene_cov_summary<-cov_data%>%
    left_join(cov_data%>%group_by(transcript_name,sample_name)%>%
                summarize(transcript_len=sum(length))%>%
                select(transcript_name,transcript_len)%>%
                distinct())%>%
    group_by(group_name,gene_name,sample_name,transcript_name)%>%
    summarize(across(c(mean_coverage,
                       pct_0x),
                     ~sum(.x*length/transcript_len)))%>%ungroup()%>%
    group_by(group_name,gene_name,transcript_name)%>%
    summarize(across(where(is.numeric),mean))%>%#now add min and max coverage
    left_join(cov_data%>%group_by(group_name,transcript_name)%>%
                summarize(min_coverage=min(min_coverage),
                          max_coverage=max(max_coverage)))
  
  return(gene_cov_summary)
}
