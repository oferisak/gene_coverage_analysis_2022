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
