# fix ucsc table clinvar tsv - some of the variants have no type.. will add . instead
fix_clinvar_ucsc_file<-function(clinvar_file,main_output_folder,only_p_or_lp=T){
  out_clinvar_file_name<-stringr::str_replace(basename(clinvar_file),'.(csv|tsv|txt).*',ifelse(only_p_or_lp,'_plp.bed','.bed'))
  out_clinvar_file<-glue('{main_output_folder}/{out_clinvar_file_name}')
  clinvar_data<-readr::read_delim(clinvar_file)
  clinvar_bed<-clinvar_data%>%select(`#chrom`,chromStart,chromEnd,clinSign,molConseq)%>%
    mutate(molConseq=ifelse(is.na(molConseq),'unspecified',molConseq))%>%
    mutate(`#chrom`=stringr::str_replace(`#chrom`,'chrMT','chrM'))
  if (only_p_or_lp){clinvar_bed<-clinvar_bed%>%filter(clinSign%in%c('Pathogenic','Likely pathogenic'))}
  write.table(clinvar_bed,sep='\t',row.names = F,col.names = F,file=out_clinvar_file,quote = F)
  return(out_clinvar_file)
}

# this function takes in a clinvar tsv file (downloaded from ucsc table browser *not bed*) and a bam file and calculates the coverage
# only_p_or_lp - analyzes only pathogenic/likely pathogenic variants
calculate_clinvar_coverage<-function(clinvar_bed,bam_file,genome_file,clinvar_cov_file){
  clinvar_coverage_command<-glue('bedtools sort -g {genome_file} -i {clinvar_bed} | bedtools coverage -a - -b {bam_file} -sorted -g {genome_file} > {clinvar_cov_file}')
  message(glue('Running {clinvar_coverage_command}'))
  system(clinvar_coverage_command)
  return(clinvar_cov_file)
}

intersect_target_with_clinvar<-function(clinvar_file,target_file,clinvar_vs_target_file,preprocess_input_clinvar=FALSE){
  if (preprocess_input_clinvar){
    awk_command<-'awk -F\'\\t\' \'{print $1"\\t"$2"\\t"$3"\\t"$14"\\t"$18}\''
    clinvar_coverage_command<-glue('gunzip -c {clinvar_file} | {awk_command} | grep "Pathogenic\\|Likely pathogenic" | bedtools intersect -c -a - -b {target_file} > {clinvar_vs_target_file}')
  }else{
    if (grepl('gz',clinvar_file)){
      clinvar_coverage_command<-glue('gunzip -c {clinvar_file} | bedtools intersect -c -a - -b {target_file} > {clinvar_vs_target_file}')
    }else{
      clinvar_coverage_command<-glue('bedtools intersect -c -a {clinvar_file} -b {target_file} > {clinvar_vs_target_file}')
    }
  }
  message(glue('Running {clinvar_coverage_command}'))
  system(clinvar_coverage_command)
  return(clinvar_vs_target_file)
}

read_clinvar_coverage_file<-function(clinvar_cov_file){
  clinvar_cov_raw<-readr::read_delim(clinvar_cov_file,
                                     col_names = c('chr','start','end','clin_sig','var_type','depth','number_of_bases_with_non_zero_cov','var_len','var_fract_covered'))
  return(clinvar_cov_raw)
}

read_clinvar_vs_target_file<-function(clinvar_vs_target_file){
  clinvar_cov_raw<-readr::read_delim(clinvar_vs_target_file,
                                     col_names = c('chr','start','end','clin_sig','var_type','num_o_overlap'))
  return(clinvar_cov_raw)
}
