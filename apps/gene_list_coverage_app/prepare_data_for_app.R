project_dir<-'/media/SSD/Bioinformatics/Projects/gene_coverage_analysis_2022'
app_folder<-'/media/SSD/Bioinformatics/Projects/gene_coverage_analysis_2022/apps/gene_list_coverage_app'
library(ProjectTemplate)
setwd(project_dir)
load.project()

# Prepare gene name table
refseq_curated_to_gene<-readr::read_delim('/media/SSD/Bioinformatics/Databases/refseq/refseq_curated_to_gene.txt.gz',
                                          col_names = c('transcript_name','chrom','gene_name'),
                                          comment = '#')
refseq_select_to_gene<-readr::read_delim('/media/SSD/Bioinformatics/Databases/refseq/refseq_select_to_gene.txt.gz',
                                          col_names = c('transcript_name','chrom','gene_name'),
                                          comment = '#')

refseq_to_gene<-refseq_curated_to_gene%>%
  mutate(cannonical=ifelse(transcript_name%in%refseq_select_to_gene$transcript_name,T,F))%>%
  rename(gene_symbol='gene_name')
write.table(refseq_to_gene,
            file=glue('{app_folder}/accessory_data/refseq_to_gene.txt'),
            row.names = F,
            quote = F,
            sep='\t')

# Generate segmental duplications vs exons file
# merge segmental duplications and use the maximal dup level as the annotation
# the seg dups file should be in the format: chr star end dup_perc
exons_file<-'/media/SSD/Bioinformatics/Databases/refseq/refseq_hg19_curated_cds_20220403.bed.gz'
seg_dups_file<-'/media/SSD/Bioinformatics/Databases/hg19/genomicSuperDups_hg19_20220414.bed.gz'
seg_dups_merged<-'/media/SSD/Bioinformatics/Databases/hg19/genomicSuperDups_hg19_20220414_merged.bed'
system(glue('bedtools sort -i {seg_dups_file} | bedtools merge -c 4 -o max > {seg_dups_merged}'))


#annot_file<-'/media/SSD/Bioinformatics/Databases/hg38/segmental_duplications/genomicSuperDups_hg38_20211115_merged.bed'
annot_file<-'/media/SSD/Bioinformatics/Databases/hg19/genomicSuperDups_hg19_20220414_merged.bed'
annotation_colnames<-c('ann_chr','ann_start','ann_end','dup_perc')

intersect_output<-bedtools_intersect(exons_file,annot_file,annotation_colnames,transcript_id_parsed = F)
intersect_output<-intersect_output%>%
  mutate(transcript_name=stringr::str_extract(interval_name,'(.+)_(exon|cds)')%>%stringr::str_replace('_(exon|cds)',''),
         exon_num=as.numeric(stringr::str_extract(interval_name,'(exon|cds)_\\d+')%>%stringr::str_extract('\\d+')))

intersect_output<-intersect_output%>%mutate(dup_perc=round(dup_perc,3))
# for each exon get the maximal level of duplication
intersect_output<-intersect_output%>%group_by(transcript_name,exon_num)%>%top_n(n=1,wt=dup_perc)%>%
  distinct(transcript_name,chr,start,end,exon_num,.keep_all = T)

exons_vs_segdup<-intersect_output%>%group_by(transcript_name,start,end,exon_num,overlap)%>%
  top_n(1,dup_perc)%>%
  distinct(gene_symbol,transcript_id,start,end,exon_num,.keep_all = T)

exons_vs_segdup<-intersect_output%>%ungroup()%>%
  group_by(transcript_name,chr,start,end,exon_num,dup_perc)%>%
  summarize(overall_overlap=sum(overlap))%>%
  mutate(exon_len=end-start,overlap_perc=overall_overlap/exon_len)%>%
  ungroup()%>%
  group_by(transcript_name)%>%arrange(transcript_name,exon_num)
summary(exons_vs_segdup)

write.table(exons_vs_segdup,file=glue('{app_folder}/accessory_data/refseq_hg19_curated_cds_vs_segdup.csv'),row.names = F,sep ='\t')

# Get all the variants not covered in the given target file ####
clinvar_vcf<-'/media/SSD/Bioinformatics/Databases/clinvar/clinvar_20221224.plp.hg19.vcf'

# IDT
target_file<-'/media/SSD/Bioinformatics/Databases/idt/xgen_hg19_exome_with_mt_targets.pad50.bed'
bedtools_intersect_with_clinvar_vcf(clinvar_vcf,bed_file = target_file,get_missing = T,
                                    output_dir = glue('{app_folder}/accessory_data'))
# Twist
target_file<-'/media/SSD/Bioinformatics/Databases/twist/twist_hg19_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_0.pad50.chrM.bed'
bedtools_intersect_with_clinvar_vcf(clinvar_vcf,bed_file = target_file,get_missing = T,
                                    output_dir = glue('{app_folder}/accessory_data'))

# coverage data ####

# IDT
library("openxlsx")
idt_coverage_file<-'/media/SSD/Bioinformatics/Projects/gene_coverage_analysis_2022/apps/gene_list_coverage_app/data_for_prepare/coverage_data/idt_exnos_coverage_genoox_202302.xlsx'
openxlsx::getSheetNames(idt_coverage_file)
# in the genoox file, the exon numbers are 1 based and in the refseq bed files they are 0 based
idt_coverage_file<-openxlsx::read.xlsx(idt_coverage_file,sheet = 'Exons (canonical)')
# Checking this issue reducing one from the exon number in the genoox file results in a higher number of retained exons
# when joining with the refseq file. so thats the better way
idt_coverage_file<-idt_coverage_file%>%mutate(exon=exon-1)
# add the location information for each exon
refseq_curated<-readr::read_delim('./accessory_data/refseq_hg19_curated_exons_20220410.bed.gz',
                                  col_names = c('chr','start','end','interval_name','zero','strand'))%>%
  filter(!grepl('_',chr))%>%
  mutate(transcript_name=stringr::str_extract(interval_name,'(.+)_(exon|cds)')%>%stringr::str_replace('_(exon|cds)',''),
         exon_num=as.numeric(stringr::str_extract(interval_name,'(exon|cds)_\\d+')%>%stringr::str_extract('\\d+')))

idt_coverage_file<-idt_coverage_file%>%rename(transcript_name=transcript,exon_num=exon)

# fix missing transcript names (convert them)
missing_transcripts<-idt_coverage_file%>%
  filter(!transcript_name%in%refseq_curated$transcript_name)%>%
  select(transcript_name)%>%
  distinct()%>%
  mutate(no_iso_transcript_name=stringr::str_replace(transcript_name,'\\..+',''))%>%
  left_join(refseq_curated%>%
              mutate(no_iso_transcript_name=stringr::str_replace(transcript_name,'\\..+',''))%>%
              select(no_iso_transcript_name,refseq_transcript=transcript_name)%>%distinct())

idt_coverage_file<-idt_coverage_file%>%
  left_join(missing_transcripts%>%select(transcript_name,refseq_transcript))%>%
  mutate(transcript_name=ifelse(is.na(refseq_transcript),transcript_name,refseq_transcript))%>%
  select(-refseq_transcript)

# now join with the refseq file (and remove the missing transcripts)
# if there are exons without covereage for the transcript add them
gene_to_transcript<-idt_coverage_file%>%select(gene,transcript_name)%>%distinct()
idt_coverage_with_refseq<-refseq_curated%>%
  inner_join(gene_to_transcript)%>%
  filter(transcript_name %in% idt_coverage_file$transcript_name)%>%
  left_join(idt_coverage_file,by=c('gene','transcript_name','exon_num'))%>%
  mutate(exonLength=ifelse(is.na(exonLength),end-start,))

idt_coverage_file<-idt_coverage_file%>%
  left_join(refseq_curated,by=c('transcript_name','exon_num'))
write.table(idt_coverage_file,file=glue('./accessory_data/idt_genoox_coverage.{Sys.Date()}.csv'),sep='\t',row.names = F,quote = F)

