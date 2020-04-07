library(sequenza)
library(optparse)

## For now, test if commands are in original, trailing format, or new opt-parse format
option_list = list(
  make_option(c("-s", "--snp_file"), type="character", default=NULL, 
              help="varscan snp file name", metavar="character"),
  make_option(c("-c", "--cnv_file"), type="character", default=NULL,
              help="varscan copy number file name [default= %default]", metavar="character"),
  make_option(c("-y", "--remove_Y_vars"), type="logical", default=TRUE,
            help="Remove Y chromosome variants? [default= %default]", metavar="logical"),
  make_option(c("-p", "--prefix"), type="character", default=NULL,
            help="Output prefix for processed data file [default= %default]", metavar="character")
 ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

snp.file <- opt$snp_file
cnv.file <- opt$cnv_file
remove_Y <- opt$remove_Y_vars
prefix   <- opt$prefix

print("Getting SeqZ file for Sequenza:")

########
# MAIN #
########

## read the snp and cnv files
snp <- read.table(snp.file,sep="\t",header = TRUE)
cnv <- read.table(cnv.file,sep="\t",header = TRUE)

## remove chromosome M because it seems to cause problems
snp <-subset(snp,chrom!="chrM")
cnv <-subset(cnv,chrom!="chrM")

## remove chromosome Y because it seems to cause problems
if(remove_Y){
snp<-subset(snp,chrom!="chrY")
cnv<-subset(cnv,chrom!="chrY")
}

##### Adjust for filtered data, if needed
cnv<-cnv[,1:8]
colnames(cnv)[7]="log2_ratio"

##### prepare sequenza data file
seqz.data <- VarScan2seqz(varscan.somatic = snp, 
                          varscan.copynumber = cnv)

file <- paste0(prefix, ".seqz")
write.table(seqz.data, file, col.names = TRUE, row.names = FALSE, sep = "\t")
