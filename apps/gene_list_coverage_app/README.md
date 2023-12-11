## Data preparation description

the prepare_data_for_app.R script includes a sequence of functions to prepare different inputs for the app

### clinvar p/lp vs target bed

includes only P/LP clinvar variants that are not covered in the given bed file.

**input**:

a clinvar P/LP file: created by downloading the vcf from <https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/> and then running:

``` bash
gunzip -c clinvar_20221224.vcf.gz | grep -i "CLNSIG=Pathogenic\|CLNSIG=Likely_Pathogen\|#" > clinvar_20221224.plp.vcf
# convert to hg19
gunzip -c clinvar_20221224.plp.vcf.gz | awk '$1 ~ /^#/ { print; next } { print "chr" $0 }' - > clinvar_20221224.plp.hg19.vcf
```

a target bed file - in the case of IDT, ive also added some padding to the file

**output**:

the clinvar vs target bed file - this file is created by running bed intersect (in the prep_data_for_app.R script), with the option -v that returns only the vcf values that have *no* overlap with the target bed regions. The required input includes a clinvar file (make sure the reference (e.g hg19) matches your needs) and a target bed file.
