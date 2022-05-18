app_folder<-'/media/SSD/Bioinformatics/Projects/gene_coverage_analysis_2022/apps/gene_list_coverage_app'

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
