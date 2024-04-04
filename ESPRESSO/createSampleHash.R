setwd('/rsrch5/home/epi/bhattacharya_lab/download/01-aws-download/bam')
bam_files = list.files()
bam_files = bam_files[!grepl('.bai',bam_files)]

df = data.frame(bam = file.path('/rsrch5/home/epi/bhattacharya_lab/download/01-aws-download/bam',
                                bam_files),
                sample = sapply(strsplit(bam_files,'.sort'),
                                function(x) x[1]))

data.table::fwrite(df,
                   '/rsrch5/home/epi/bhattacharya_lab/projects/placenta_mapqtl/LR_RNAseq/ESPRESSO/sampleHash.tsv',
                   sep='\t',
                   col.names=T,
                   row.names=F,
                   quote=F)