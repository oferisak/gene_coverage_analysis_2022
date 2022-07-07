library(glue)
app_folder<-'/media/SSD/Bioinformatics/Projects/gene_coverage_analysis_2022/apps/gene_list_coverage_app'
# set the padding around target region to use
pad<-50

# in order to generate the founder variants file i used the supp table from Bella davidovs paper
# i then copied the transcript:coding sequence column to the variantvalidator batch tool and used the output

target_region<-'/media/SSD/Bioinformatics/Databases/idt/xgen-exome-research-panel-v2-targets-hg19.bed'
# first create a padded version of the target region (50bp)
original_tr<-readr::read_delim(target_region,col_names = c('chr','start','end'))
padded_tr<-original_tr%>%mutate(start=as.character(start-pad),
                                end=as.character(end+pad))


padded_tr_name<-stringr::str_replace(basename(target_region),'.bed','_pad50.bed')
padded_tr_file<-glue('{app_folder}/accessory_data/founder_variants/{padded_tr_name}')

write.table(file=padded_tr_file,
            padded_tr,
            col.names = F,
            row.names = F,
            quote = F,
            sep='\t')

# founder varaints genomic coordinates obtained from https://variantvalidator.org/
original_founder_file<-glue('{app_folder}/accessory_data/founder_variants/davidov_et_al_clinical_genetics_s8.csv')
vv_founder_file<-glue('{app_folder}/accessory_data/founder_variants/davidov_variant_validator_output.csv')

original_founder_variants<-readr::read_delim(original_founder_file,col_names = c('gene_symbol','omim','Input'),skip = 1)
vv_founder_variants<-readr::read_delim(vv_founder_file,comment = '#')

# join the original file with the disease name and the vv output
founder_variants<-original_founder_variants%>%left_join(vv_founder_variants)


# now generate a bed file from the genomic positions -
# bed files are 0 based so i need to reduce 1 from the coordinates (the nucleotide with the coordinate 1 in a genome will have a value of 0 in column 2 and a value of 1 in column 3.)
chr_order<-paste0('chr',c(as.character(1:22),'X'))
# add "chr" to the chromsome names (should be ok)
founder_bed_file<-founder_variants%>%select(chr=GRCh37_CHR,end=GRCh37_POS,gene_symbol=Gene_Symbol,hgvs=Input,omim)%>%
  filter(!is.na(chr))%>%
  mutate(start=end-1,
         chr=fct_relevel(paste0('chr',chr),chr_order))

bed_file_name<-glue('{app_folder}/accessory_data/founder_variants/founder_variants_davidov.bed')
write.table(file=bed_file_name,
            founder_bed_file%>%select(chr,start,end,gene_symbol,hgvs,omim)%>%
              arrange(chr,desc(start))%>%
              mutate(start=as.character(start),
                     end=as.character(end)),
            col.names = F,
            row.names = F,
            quote = F,
            sep='\t')

# now intersect the founder bed with the target file
bedtools_command<-glue('bedtools intersect -a {bed_file_name} -b {padded_tr_file} -c')
bedtools_intersect_output<-system(bedtools_command,intern = T)
fv_vs_target<-read.table(text=bedtools_intersect_output,sep='\t')
# write the missing variants 
write.table(fv_vs_target%>%filter(V7==0),
            file=glue('{app_folder}/accessory_data/founder_variants/idt_missing_variants.csv'),
            row.names = F,
            col.names = F,
            quote = F,
            sep = '\t')
