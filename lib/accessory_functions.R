parse_gene_list<-function(gene_list_file,target_region){
  gene_list<-readr::read_delim(gene_list_file,col_names = 'gene_name',col_select = 1)
  # check if its hgnc or refseq_transcript
  if (sum(grepl('(NM|NR)_',gene_list$gene_name))==0){
    gene_list<-gene_list%>%left_join(refseq2gene)
  }
}

get_bam_sequences_from_header<-function(input_bam){
  headers_command<-glue('samtools view -H {input_bam}')
  bam_headers<-system(headers_command,intern = T)
  bam_headers_df<-read.table(textConnection(bam_headers),fill = NA,sep = '\t')
  bam_headers_seqs<-bam_headers_df%>%filter(V1=='@SQ')%>%pull(V2)%>%stringr::str_replace('SN:','')
  return(bam_headers_seqs)
}

# get the bam headers and build the reference fasta according to the order of the chr in the bam files
build_ref_according_to_bam<-function(input_bam,reference_chrs_path,reference_name){
  bam_headers_seqs<-get_bam_sequences_from_header(input_bam)
  setwd(reference_chrs_path)
  # verify chrs in bam headers are found in reference folder
  seqs_in_folder<-list.files(reference_chrs_path,pattern = '*.fa')
  for (bam_seq in bam_headers_seqs){
    if (!(glue('{bam_seq}.fa') %in% seqs_in_folder)){message(glue('could not find {bam_seq}.fa in folder'))}
  }
  join_command<-glue("cat {paste0(paste0(bam_headers_seqs,'.fa'),collapse=' ')} > {reference_name}")
  message(glue('Joining sequences with command: {join_command}'))
  system(join_command)
}

numbers_to_ranges<-function(numbers){
  f <- function(x){
    if (length(x) == 1) x else paste(x[1], x[length(x)], sep = "-")
  }  
  ranges<-tapply(numbers[order(numbers)], cumsum(c(TRUE, diff(numbers[order(numbers)]) != 1)), f)
  return(paste0(ranges,collapse=','))
}

