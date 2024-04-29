.libPaths('/rsrch5/home/neuro_rsrch/ekilinc/R/libs')
library("optparse")
library('data.table')

option_list <- list(
  make_option(c("-s", "--sample"), action="store_true", default=TRUE,
              help="barcode",type='numeric')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

sampleID = opt$sample
setwd('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/bam')
hash = read.table('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/Espresso/sampleHashTrial.tsv')
fff = hash$V1

bam_file = fff[sampleID]
prefix = strsplit(fff[sampleID],'.sorted')[[1]][1]
prefix = strsplit(prefix,'/')[[1]][9]

### BAM TO BED12
bed12_file = file.path(getwd(),
                       paste0(prefix,'.bed12'))
system(paste('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/software/bedtools bamtobed -bed12 -i',
             bam_file,'>',bed12_file))

### FLAIR CORRECT
ref_folder='/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/ref'
gtf_file = file.path(ref_folder,
                     'gencode.v44.annotation.gtf')
assembly_file = file.path(ref_folder,
                         'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna')
out_flair = '/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/Flair/flair_out'

system(paste('flair correct -q',bed12_file,
             '-f',gtf_file,
             '-g',assembly_file,
             '--output',file.path(out_flair,
                                  prefix)))
