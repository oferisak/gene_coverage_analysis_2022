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

read_input_bam_files<-function(bam_folder=NA){
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

parse_refseq_transcript_to_gene_name<-function(){
  refseq2gene<-readr::read_tsv('./data/refseq_transcript_to_gene.txt.gz',comment = '')
  colnames(refseq2gene)<-c('transcript_name','gene_name')
  return(refseq2gene%>%distinct())
}