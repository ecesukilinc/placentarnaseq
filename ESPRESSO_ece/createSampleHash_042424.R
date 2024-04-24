.libPaths('/rsrch5/home/neuro_rsrch/ekilinc/R/libs')
library('data.table')

setwd('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/bam')
bam_files = list.files()
bam_files = bam_files[!grepl('.bai',bam_files)]

df = data.frame(bam = file.path('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/bam',
                                bam_files),
                sample = sapply(strsplit(bam_files,'.sort'),
                                function(x) x[1]))

data.table::fwrite(df,
                   '/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/Espresso/sampleHash.tsv',
                   sep='\t',
                   col.names=T,
                   row.names=F,
                   quote=F)