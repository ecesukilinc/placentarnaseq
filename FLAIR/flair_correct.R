library("optparse")
option_list <- list(
  make_option(c("-s", "--sample"), action="store_true", default=TRUE,
              help="barcode",type='numeric')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

sampleID = opt$sample
setwd('/rsrch5/home/epi/bhattacharya_lab/download/01-aws-download/bam')
hash = read.table('/rsrch5/home/epi/bhattacharya_lab/projects/placenta_mapqtl/LR_RNAseq/ESPRESSO/sampleHash.tsv')
fff = hash$V1

bam_file = fff[sampleID]
prefix = strsplit(fff[sampleID],'.sorted')[[1]][1]
prefix = strsplit(prefix,'/')[[1]][9]

### BAM TO BED12
bed12_file = file.path(getwd(),
                       paste0(prefix,'.bed12'))
system(paste('/rsrch5/home/epi/bhattacharya_lab/software/bedtools bamtobed -bed12 -i',
             bam_file,'>',bed12_file))

### FLAIR CORRECT
ref_folder='/rsrch5/home/epi/bhattacharya_lab/projects/placenta_mapqtl/reference'
gtf_file = file.path(ref_folder,
                     'gencode.v44.annotation.gtf')
assembly_file = file.path(ref_folder,
                         'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna')
out_flair = '/rsrch5/home/epi/bhattacharya_lab/download/01-aws-download/flair_out'

system(paste('flair correct -q',bed12_file,
             '-f',gtf_file,
             '-g',assembly_file,
             '--output',file.path(out_flair,
                                  prefix)))
