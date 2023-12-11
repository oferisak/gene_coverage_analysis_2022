
############# functions #############

prepare_genome_file<-function(input_bam,output_folder){
  samtools_command<-glue('samtools view -H {input_bam}')
  sample_name<-basename(input_bam)
  bam_headers<-system(samtools_command,intern=T)
  region_headers<-as.character(bam_headers%>%stringr::str_extract_all('SN:.+LN:\\d+',simplify = T))
  region_headers<-region_headers[region_headers!='']%>%stringr::str_replace_all('SN:|LN:','')
  
  genome_file<-glue('{output_folder}/{sample_name}.genome')
  writeLines(region_headers,genome_file)
  return(genome_file)
}

get_coverage_on_file = function(input_bam,output_folder, target_region, genome_file){
  genome_file<-prepare_genome_file(input_bam,output_folder = output_folder)
  
  bam_file_name = basename(bam)
  coverage_output <-glue('{output_folder}/{bam_file_name}.coverage')
  bed_coverage_command<-glue('samtools view -uF 0x400 {input_bam} | bedtools bamtobed -i - | bedtools coverage -g {genome_file} -a {target_region} -b - -hist -sorted > {coverage_output}')
  
  message(sprintf('Running: %s',bed_coverage_command))
  system(bed_coverage_command)
  return(coverage_output)
}


process_coverage = function(coverage_output){
  
  cov = read.delim(coverage_output, header = F)
  coverage_one_sample = cov[cov$V1 != "all",c(1:5,11)] 
  coverage_one_sample_col = aggregate(coverage_one_sample$V11, by = list(coverage_one_sample$V1,
                                                                         coverage_one_sample$V2,
                                                                         coverage_one_sample$V3,
                                                                         coverage_one_sample$V4,
                                                                         coverage_one_sample$V5), min)
  sample_name = gsub(".bam.coverage", "", basename(coverage_output))
  
  common_col = c("chr","pos","dbsnp","ref","alt")
  names(coverage_one_sample_col) = c(common_col ,sample_name)
  return(list(coverage_one_sample_col, common_col))
}


################ analysis ##############

project_dir = "/data/gene_coverage_analysis_2022/"


# Gets in_dir and target file (vcf)
# Go over all bam files within in_dir and get coverage files 
# Write results to in_dir

in_dir = "/data/PGx_variants//bam_files_IDT/"
target_region = "/data/PGx_variants/pharmcat_positions.vcf"
out_file = "/data/PGx_variants/results/pharmcat_positions_coverage_IDT.csv"

library(ProjectTemplate)
setwd(project_dir)
load.project()

output_folder = in_dir

bam_files = list.files(in_dir, pattern = ".bam$")



for (bam in bam_files){
  input_bam = glue("{in_dir}/{bam}")
  
  genome_file = prepare_genome_file(input_bam,output_folder)
  
  coverage_output = get_coverage_on_file(input_bam,output_folder, target_region, genome_file)
  
  
  res = process_coverage(coverage_output)
  coverage_one_sample_col = res[[1]]
  common_col = res[[2]]
  
  if (!exists("all_cov")){
    all_cov = coverage_one_sample_col
  }else{
    all_cov = merge(all_cov, coverage_one_sample_col, by = common_col)
  }
}


number_col = ncol(all_cov)

all_cov$min_cov = apply(all_cov[,6:number_col], 1, min)
all_cov$avg_cov = apply(all_cov[,6:number_col], 1, mean)
all_cov$median_cov = apply(all_cov[,6:number_col], 1, median)

write.csv(all_cov, out_file)


stop()

########### merge with non covered ####


res_dir = "/data/PGx_variants/results/"

setwd(res_dir)

common_col = c("chr","pos","dbsnp","ref","alt")

not_cov = read.table("not_covered_idt.bed")
not_cov = not_cov[,1:5]
names(not_cov) = common_col
not_cov$IDT_target_covered = "IDT not covered"

idt  = read.csv("pharmcat_positions_coverage_IDT.csv", row.names = 1)
names(idt)[6:8] = gsub("X", "", names(idt)[6:8])
names(idt)[6:11] = paste0(names(idt)[6:11], "_IDT")

idt = merge(idt,not_cov, all = T)


twist = read.csv("pharmcat_positions_coverage_twist.csv", row.names = 1)
names(twist) = gsub("X", "", names(twist))
names(twist)[12:14] = paste0(names(twist)[12:14], "_TW")

idt_wt = merge(twist, idt,  all =T)

idt_wt = idt_wt[order(idt_wt$IDT_target_covered),]

write.csv(idt_wt,"pharmcat_positions_coverage_sum_results.csv", 
          row.names = F)
