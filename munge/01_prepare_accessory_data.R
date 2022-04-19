parse_refseq_transcript_to_gene_name<-function(){
  refseq2gene<-readr::read_tsv('./data/refseq_transcript_to_gene.txt.gz',comment = '')
  colnames(refseq2gene)<-c('transcript_name','gene_name')
  return(refseq2gene%>%distinct())
}

refseq2gene<-parse_refseq_transcript_to_gene_name()