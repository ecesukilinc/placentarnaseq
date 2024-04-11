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

setwd('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/')
hash = data.table::fread('lrReadBarcode_SampleIDHash_092623.csv')
hash$Barcode = paste0('barcode',as.numeric(sapply(strsplit(hash$`Barcode (PCB111.24)`,
                                                           'ode'),
                                                  function(x) x[2])))
all_samples = unique(hash$`Sample ID`)
all_samples = all_samples[!is.na(all_samples)]

hash = subset(hash,`Sample ID` == all_samples[sampleID])
muxID = unique(hash$MUX)
bc = as.numeric(strsplit(hash$Barcode[1],'ode')[[1]][2])

file_list = 
  data.table::fread('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/updatedFileList.tsv')
  setnames(file_list, c('V1','V2','V3','V4'), c('base_folder', 'level1', 'barcode', 'file'))

NanoFilt = 'NanoFilt'

file_list = subset(file_list,level1 != '')

barcode_all = bc

file_list = subset(file_list,
                   barcode == bc &
                     level1 %in% 
                     paste0('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/output/',
                            muxID))
all_files = file.path(file_list$level1,
                      file_list$file)
all_files = all_files[grepl('fastq_pass/',
                            all_files)]

if (!dir.exists('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/pass/')) {
  dir.create('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/pass/')
}

setwd('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/pass/')
out = paste0('SampleID_',all_samples[sampleID][1],'.fastq')

if (!file.exists(paste0(out, '.gz'))) {
  for (fq in all_files) {
    system(paste('zcat',fq,'>> ',out))
}
system(paste('gzip', out))
} else {
  message(paste('Merged FASTQ file', paste0(out, '.gz'), 'already exists'))
}

if (!dir.exists('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/trim/')) {
  dir.create('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/trim/')
}

out_trim = paste0('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/trim/',
                  paste0('SampleID_',all_samples[sampleID][1],'.fastq.gz'))

if (!file.exists(out_trim)) {
  system(paste('zcat',paste0(out,'.gz |'),
              NanoFilt,'-l 200  -q 7 >> ',
              out_trim))
  system(paste('gzip',paste0(out,'.gz |'), out_trim))
} else {
  message(paste('Trimmed FASTQ file', out_trim, 'already exists'))
}

# system(paste('zcat',paste0(out,'.gz |'),
#               NanoFilt,'-l 200  -q 7 >> ',
#               paste0('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/trim/',
#                      out)))

# system(paste('gzip',paste0(out,'.gz |'),
#               paste0('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/trim/',
#                      out)))

# out_trim = paste0('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/trim/',
#                      out,'.gz')

if (!dir.exists('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/bam/')){
  dir.create('/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/bam/')
} 
bam_folder='/rsrch5/scratch/neuro_rsrch/ekilinc/lnrnatrial/placenta_trial/bam/'

ref='/rsrch5/home/epi/bhattacharya_lab/projects/placenta_mapqtl/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'

aln_out=file.path(bam_folder,paste0('SampleID_',all_samples[sampleID][1],'.sam'))
bam_out=file.path(bam_folder,paste0('SampleID_',all_samples[sampleID][1],'.bam'))
bam_sort=file.path(bam_folder,paste0('SampleID_',all_samples[sampleID][1],'.sorted.bam'))
system(paste('minimap2 -ax map-ont --sam-hit-only',
             ref,out_trim,' > ',aln_out))

# module load samtools
system(paste('samtools view -S -b',
             aln_out,' > ',bam_out))
file.remove(aln_out)
system(paste('samtools sort',
             bam_out,' -o ',bam_sort))
file.remove(bam_out)
system(paste('samtools index',bam_sort))

