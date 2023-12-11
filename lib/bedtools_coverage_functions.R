sort_input_bed_file<-function(input_bed_file,output_folder,genome_file){
  message(glue('Sorting {input_bed_file}'))
  bed_file_name<-basename(input_bed_file)%>%stringr::str_replace('.bed','')
  input_bed<-readr::read_delim(input_bed_file,col_names = F)
  genome<-readr::read_delim(genome_file,col_names = F)
  # remove regions not found in the genome file
  input_bed_fixed<-input_bed%>%filter(X1%in%genome$X1)
  # sort bed file by genome file
  input_bed_fixed<-input_bed_fixed%>%group_by(X1)%>%arrange(X2)
  input_bed_fixed<-input_bed_fixed[order(match(input_bed_fixed$X1,genome$X1)),]
  input_bed_fixed_file<-glue('{output_folder}/{bed_file_name}.sorted.bed')
  options(scipen=10)
  write.table(input_bed_fixed,file=input_bed_fixed_file,row.names = F,col.names = F,sep='\t',quote = F)
  # sort_command<- ifelse(grepl('.gz',input_bed),
  #                       sprintf('gunzip -c %s',input_bed),
  #                       sprintf('cat %s',input_bed))
  #sort_command<-glue('bedtools sort -g {genome_file} -i {input_bed} > {output_folder}/{bed_file_name}.sorted.bed')
  
  #sort_command<-glue('{sort_command} | sort -k1,1 -k2,2n - > {output_folder}/{bed_file_name}.sorted.bed')
  #system(sort_command)
  return(input_bed_fixed_file)
}

prepare_genome_file<-function(input_bam,output_folder){
  samtools_command<-glue('samtools view -H {input_bam}')
  sample_name<-basename(output_folder)
  bam_headers<-system(samtools_command,intern=T)
  region_headers<-as.character(bam_headers%>%stringr::str_extract_all('SN:.+LN:\\d+',simplify = T))
  region_headers<-region_headers[region_headers!='']%>%stringr::str_replace_all('SN:|LN:','')
  
  genome_file<-glue('{output_folder}/{sample_name}.genome')
  #write.table(sorted_region_headers,file=genome_file,row.names = F,col.names = F,sep='\t',quote = F)
  writeLines(region_headers,genome_file)
  return(genome_file)
}

run_bedtools_coverage_by_bam<-function(input_bam,genome_file,target_region,output_folder,remove_duplicates=TRUE){
  bam_file_name<-basename(input_bam)%>%stringr::str_replace('.bam','')
  coverage_output<-glue('{output_folder}/{bam_file_name}.coverage')
  if (file.exists(coverage_output)){
    message(glue('{bam_file_name} already has a coverage file, will skip re-calculation.'))
    return(coverage_output)
  }
  #bed_coverage_command<-glue('bedtools bamtobed -i {input_bam} | bedtools coverage -g {genome_file} -a {target_region} -b - -hist -sorted > {coverage_output}')
  if (remove_duplicates){
    bed_coverage_command<-glue('samtools view -uF 0x400 {input_bam} | bedtools coverage -g {genome_file} -hist -b - -a {target_region} -sorted > {coverage_output}')
  }else{
    bed_coverage_command<-glue('bedtools bamtobed -i {input_bam} | bedtools coverage -g {genome_file} -a {target_region} -b - -hist -sorted > {coverage_output}')
  }
  message(sprintf('Running: %s',bed_coverage_command))
  info(logger,bed_coverage_command)
  system(bed_coverage_command)
  return(coverage_output)
}

# bed coverage by bed - check how one bed file (bed_b) covers a second bed (bed a)
run_bedtools_coverage_by_bed<-function(bed_a,bed_b,bed_name,output_folder){
  test_bed_name<-stringr::str_replace(basename(bed_b),'.bed','')
  output_dir<-glue('./{output_folder}/{bed_name}/')
  if (!dir.exists(output_dir)){dir.create(output_dir)}
  fix_input_bed_awk<-'awk -F\'\\t\' \'{print $1"\\t"$2"\\t"$3}\''
  # after grabing only the first three columns, merge the regions so there will be a maximum of one overlap
  bed_coverage_command<-glue('{fix_input_bed_awk} "{bed_b}"| sort -k 1,1 -k2,2n - | bedtools merge -i - | bedtools coverage -mean -a "{bed_a}" -b - > {output_dir}/{bed_name}.coverage')
  message(glue('Running {bed_coverage_command}'))
  system(bed_coverage_command)
  return(glue('./{output_dir}/{bed_name}.coverage'))
}

read_bed_coverage_file<-function(bedtools_coverage_output,remove_alt=T){
  cov_data <-
    readr::read_tsv(
      bedtools_coverage_output,
      show_col_types = FALSE,
      col_names = c(
        'chr',
        'start',
        'end',
        'interval_name',
        'zero',
        'strand',
        'coverage'
      )
    )
  
  # remove coverage stats for all the regions together
  cov_data<-cov_data%>%filter(chr!='all')
  
  # if remove_alt remove all chr names with alt in them
  cov_data<-cov_data%>%filter(!grepl('alt|fix|hap|chrUn|random',chr))
  
  # add the gene name to each transcript
  cov_data<-cov_data%>%mutate(transcript_name=stringr::str_extract(interval_name,'(.+)_(exon|cds)')%>%stringr::str_replace('_(exon|cds)',''),
                      exon_num=as.numeric(stringr::str_extract(interval_name,'(exon|cds)_\\d+')%>%stringr::str_extract('\\d+')))
  # join with gene name
  cov_data<-cov_data%>%left_join(refseq2gene)
  
  return(cov_data)
}

read_bed_coverage_hist_file<-function(bedtools_coverage_hist_output,remove_alt=T){
  message(glue('Reading coverage file {bedtools_coverage_hist_output}'))
  cov_data <-
    readr::read_tsv(
      bedtools_coverage_hist_output,
      show_col_types = FALSE,
      col_names = c(
        'chr',
        'start',
        'end',
        'interval_name',
        'zero',
        'strand',
        'coverage',
        'bases_at_coverage',
        'interval_length',
        'perc_at_coverage'
      )
    )
  
  # remove coverage stats for all the regions together
  cov_data<-cov_data%>%filter(chr!='all')
  
  # if remove_alt remove all chr names with alt in them
  cov_data<-cov_data%>%filter(!grepl('alt|fix|hap|chrUn|random',chr))
  
  # add the gene name to each transcript
  cov_data<-cov_data%>%mutate(transcript_name=stringr::str_extract(interval_name,'(.+)_(exon|cds)')%>%stringr::str_replace('_(exon|cds)',''),
                              exon_num=as.numeric(stringr::str_extract(interval_name,'(exon|cds)_\\d+')%>%stringr::str_extract('\\d+')))
  # join with gene name
  #cov_data<-cov_data%>%left_join(refseq2gene)
  
  return(cov_data)
}



# Given a bedtools coverage output (using the command with -mean), calculate the coverage per exon
summarize_coverage_per_exon<-function(bed_coverage_output,cov_threshs=c(0,10,20,50,100,200),gene_list_file=NA){
  message('Summarizing coverage per exon..')
  # Read the coverage bed file
  message(glue('Reading coverage file'))
  cov_data<-read_bed_coverage_file(bed_coverage_output)
  
  return(cov_data)
}

# Given a (pre-read) bedtools coverage output (using the command with -hist), calculate the coverage per exon
summarize_coverage_hist_per_exon<-function(cov_data){
  message('Summarizing coverage hist per exon..')
  # Read the coverage bed file
  
  #cov_data<-read_bed_coverage_hist_file(bed_coverage_hist_output)
  
  #test<-cov_data%>%slice(1:10000000)
  cov_per_exon<-cov_data%>%
    mutate(coverage_cat=cut(coverage,breaks=c(0,1,21,51,101,999999),right=F,labels=c('0x','1-20x','20-50x','50-100x','>100x')))%>%
    group_by(transcript_name,exon_num,chr,coverage_cat)%>%
    summarize(perc_at_cov_cat=sum(perc_at_coverage))
  cov_per_exon<-cov_per_exon%>%left_join(refseq2gene)
  return(cov_per_exon)
}

# Given a (pre-read) bedtools coverage output (using the command with -hist), calculate the coverage per exon
summarize_coverage_hist_per_transcript<-function(cov_data){
  message('Summarizing coverage hist per transcript..')
  # Read the coverage bed file
  
  #cov_data<-read_bed_coverage_hist_file(bed_coverage_hist_output)
  transcript_lengths<-cov_data%>%select(transcript_name,exon_num,chr,interval_length)%>%distinct()%>%
    group_by(transcript_name,chr)%>%summarize(transcript_length=sum(interval_length))
  #test<-cov_data%>%slice(1:10000000)
  
  cov_per_transcript<-cov_data%>%
    mutate(coverage_cat=cut(coverage,breaks=c(0,1,21,51,101,999999),right=F,labels=c('0x','1-20x','20-50x','50-100x','>100x')))%>%
    group_by(transcript_name,chr,coverage_cat)%>%
    summarize(bases_at_coverage=sum(bases_at_coverage))%>%
    left_join(transcript_lengths)%>%
    mutate(perc_at_cov_cat=bases_at_coverage/transcript_length)
    
  cov_per_transcript<-cov_per_transcript%>%left_join(refseq2gene)
  return(cov_per_transcript)
}

summarize_coverage_per_transcript<-function(bedtools_coverage_output,cov_threshs=c(0,10,20,50,100,200),gene_list_file = NA,output_path=NA){
  # check if transcript coverage analysis already performed
  if (!is.na(output_path) & file.exists(output_path)){
    message(glue('{output_path} already exists. Reading it..'))
    cov_per_transcript<-readr::read_delim(output_path)
    return(cov_per_transcript)
  }
  cov_data<-read_bed_coverage_file(bed_coverage_output)
  if (!is.na(gene_list_file)){
    gene_list<-parse_gene_list(gene_list_file)
    # now remove all transcripts not in the gene list
    cov_data<-cov_data%>%filter(transcript_name %in% gene_list$transcript_name)
    message(glue('Gene list provided, will analyze {length(unique(cov_data$transcript_name))} transcripts'))
  }

  message('Summarizing coverage per transcript')
  cov_per_transcript<-cov_data%>%mutate(interval_len=end-start+1)%>%
    group_by(gene_name,transcript_name)%>%
    summarise(transcript_size=sum(interval_len),
              transcript_cov=sum(coverage*interval_len/transcript_size),
              exons_without_full_cov=sum(coverage<1))
  
  if (!is.na(output_path)){write.table(cov_per_transcript,sep='\t',file=glue('{output_path}'),row.names = F)}
  return(cov_per_transcript)
}

collect_per_transcript_coverage_summaries<-function(bam_analysis_table,gene_list_file=NA){
  bam_analysis_table<-readr::read_delim(glue('{main_output_folder}/bam_analysis_table.csv'),show_col_types = FALSE)
  per_transcript_data<-NULL
  for (i in 1:nrow(bam_analysis_table)){
    bam_line<-bam_analysis_table%>%slice(i)
    group_name<-bam_line%>%pull(group_name)
    sample_name<-bam_line%>%pull(sample_name)
    message(glue('Collecting coverage metrics for {group_name}: {sample_name}'))
    cov_per_transcript<-summarize_coverage_per_transcript(bedtools_coverage_output = bam_line%>%pull(coverage_output),
                                                          gene_list_file=gene_list_file)
    per_transcript_data<-per_transcript_data%>%bind_rows(data.frame(group_name=group_name,
                                                                    sample_name,
                                                                    cov_per_transcript))
  }
  return(per_transcript_data)
}

plot_coverage<-function(bed_coverage_ouput){
  transcript_cov<-summarize_coverage_per_transcript(bed_coverage_ouput)
  # remove coverage stats for all the regions together
  g<-transcript_cov%>%ggplot(aes(x=transcript_name,y=perc_above_20,fill=perc_above_20))+
    geom_col(alpha=0.5)+facet_grid(gene_name~.,scales='free_y',space='free',switch="y")+
    geom_hline(yintercept = 0.9,linetype=2)+
    scale_y_continuous(labels = scales::percent)+
    scale_fill_gradient2(low='darkred',high='darkgreen',mid='orange',limits=c(0.5,1),midpoint=0.75,breaks=c(0,0.5,0.75,1))+
    theme_minimal()+
    labs(x=NULL,y='Percent above x20')+guides(fill='none')+
    theme(strip.placement = "outside",strip.text.y.left = element_text(face='bold',angle = 0))+
    coord_flip()
  g
  return(g)
}

get_target_bed_metrics<-function(cov_per_transcript){
  metrics<-data.frame(
    # count the number of transcripts without full coverage - since this analysis is based on an input bed file, 
    # there is no meaning to the mean coverage or perc above >Xx - only for the percent above 0x as a measure of 
    # how one bed file overlaps with the other
    num_o_transcripts_with_non_full_cov=sum(cov_per_transcript$transcript_cov<1),
    # per transcript full coverage metrics
    mean_coverage=cov_per_transcript%>%pull(transcript_cov)%>%mean()
  )
  return(metrics)
}

