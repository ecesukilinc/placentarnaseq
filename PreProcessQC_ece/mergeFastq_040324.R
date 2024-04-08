library("optparse")
option_list <- list(
  make_option(c("-s", "--sample"), action="store_true", default=TRUE,
              help="barcode",type='numeric')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

sampleID = opt$sample

setwd('/rsrch5/home/epi/bhattacharya_lab/download/01-aws-download')
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
  data.table::fread('/rsrch5/home/epi/bhattacharya_lab/download/fileList.tsv')

NanoFilt = 'NanoFilt'

file_list = subset(file_list,level1 != '')

barcode_all = bc

file_list = subset(file_list,
                   barcode == bc &
                     level1 %in% 
                     paste0('/rsrch5/home/epi/bhattacharya_lab/download/01-aws-download/output/',
                            muxID))
all_files = file.path(file_list$level1,
                      file_list$file)
all_files = all_files[grepl('fastq_pass/',
                            all_files)]

dir.create('/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/pass')
setwd('/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/pass')
out = paste0('SampleID_',all_samples[sampleID][1],'.fastq')
for (fq in all_files){
  
  system(paste('zcat',fq,'>> ',out))
  
}

system(paste('gzip ',out))

system(paste('zcat',paste0(out,'.gz |'),
              NanoFilt,'-l 200  -q 7 >> ',
              paste0('/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/trim/',
                     out)))

system(paste('gzip',paste0(out,'.gz |'),
              paste0('/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/trim/',
                     out)))

out_trim = paste0('/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/trim/',
                     out,'.gz')

dir.create('/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/bam/')
bam_folder = '/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/bam/'

ref='/rsrch5/home/epi/bhattacharya_lab/projects/placenta_mapqtl/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'

aln_out=file.path(bam_folder,paste0('SampleID_',all_samples[sampleID][1],'.sam'))
bam_out=file.path(bam_folder,paste0('SampleID_',all_samples[sampleID][1],'.bam'))
bam_sort=file.path(bam_folder,paste0('SampleID_',all_samples[sampleID][1],'.sorted.bam'))
system(paste('minimap2 -ax map-ont --sam-hit-only',
              ref,out_trim,' > ',aln_out))
 
system(paste('samtools view -S -b',
           aln_out,' > ',bam_out))
file.remove(aln_out)
system(paste('samtools sort',
              bam_out,' -o ',bam_sort))
file.remove(bam_out)
system(paste('samtools index',bam_sort))


