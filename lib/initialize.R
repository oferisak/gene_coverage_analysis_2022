input_files<-readr::read_delim('./config/input_files.txt',delim='\t')
# analysis_setup
analysis_setup<-readr::read_delim(analysis_setup_file,delim='\t')
# analysis type
analysis_type<-analysis_setup%>%filter(param=='analysis_type')%>%pull(value)
# bam folder
bam_folder<-analysis_setup%>%filter(param=='bam_folder')%>%pull(value)
input_bams_df<-read_input_bam_files(bam_folder)
# clinvar files
clinvar_file_name<-analysis_setup%>%filter(param=='clinvar_file')%>%pull(value)
clinvar_file<-input_files%>%filter(type=='clinvar_snv',name==clinvar_file_name)%>%pull(file)
# target regions
target_regions_names<-analysis_setup%>%filter(param=='target_regions')%>%pull(value)%>%strsplit(',')%>%unlist()
target_regions_table<-input_files%>%inner_join(data.frame(name=target_regions_names))
target_regions<-target_regions_table%>%pull(file)
target_reference_regions<-target_regions_table%>%pull(ref)
per_target_analysis_req<-analysis_setup%>%filter(param=='per_target_analysis')%>%pull(value)%>%strsplit(',')%>%unlist()
# reference regions
reference_fasta<-input_files%>%filter(type=='reference_fasta')%>%inner_join(data.frame(ref=target_reference_regions))%>%pull(file)
reference_dict<-input_files%>%filter(type=='reference_dict')%>%inner_join(data.frame(ref=target_reference_regions))%>%pull(file)
# remove duplicates
remove_duplicates<-ifelse(analysis_setup%>%filter(param=='remove_duplicates')%>%pull(value)==1,TRUE,FALSE)
# refseq transcript to gene name
refseq2gene<-parse_refseq_transcript_to_gene_name()