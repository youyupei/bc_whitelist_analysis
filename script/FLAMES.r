library(devtools)
library(FLAMES)

############
#Running FLAMES for quatifiaction 
##############

########
#Parametrs 
# edit distance of 2 
# annotation gencodev31
# genome hg38
#######

####Inputs that change depending on sampLe
###whitelists
Whitelist = "path_to_whitelsit"

#####fastq_files 
fastq = "path_to_fastq"

###out_dirs
output = "apth_to_output_dir"

####Inputs that don't change
####referecne files 
GTF = "gencode.v31.annotation.gtf"
genome = "hg38.fa"

minimap2_dir="path_to_minimap"

####Comand for Q20 -> BLAZE
sce <- sc_long_pipeline(fastq=fastq, outdir=output, reference_csv=Whitelist, annot=GTF, genome_fa=genome, match_barcode=TRUE, MAX_DIST=2, has_UMI=TRUE, minimap2_dir=minimap2_dir)
