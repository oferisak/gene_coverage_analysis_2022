bedtools_intersect <- function(exons_file, annot_file,annotation_colnames,transcript_id_parsed=T) {
  # -wo - Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported.
  intersect_command<-sprintf('bedtools intersect -wo -a %s -b %s',exons_file,annot_file)
  intersect_output<-read.table(text=system(intersect_command, intern = TRUE),sep='\t')
  if(transcript_id_parsed){
    colnames(intersect_output)<-c('chr','start','end','gene_symbol','transcript_id','exon_num',
                                  annotation_colnames,
                                  'overlap')
    # remove non-curated transcripts
    intersect_output<-intersect_output%>%filter(grepl('(NM|NR)_',transcript_id))
  }else{
    colnames(intersect_output)<-c('chr','start','end','interval_name','zero','strand',
                                  annotation_colnames,
                                  'overlap')
  }

  return(intersect_output)
}



# a function that takes in a bed file and a vcf file and returns all the variants that are not found in the 
bedtools_intersect_with_clinvar_vcf<-function(vcf_file,bed_file,get_missing=T,pad=50,output_dir){
  bed_file_name<-stringr::str_replace(basename(bed_file),'.[^.]+$','')
  vcf_file_name<-stringr::str_replace_all(basename(vcf_file),'.vcf|.gz','')
  intersect_file_name<-glue('{vcf_file_name}_vs_{bed_file_name}.bed.gz')
  if (get_missing){
    intersect_command<-glue('bedtools intersect -v -a {vcf_file} -b {bed_file} -bed | gzip - > {output_dir}/{intersect_file_name}')
  }else{
    intersect_command<-glue('bedtools intersect -a {vcf_file} -b {bed_file} -bed | gzip - > {output_dir}/{intersect_file_name}')
  }
  message(glue('Running {intersect_command}'))
  system(intersect_command)
  
  
}




