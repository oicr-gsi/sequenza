library(sequenza)
library(optparse)

## For now, test if commands are in original, trailing format, or new opt-parse format
option_list = list(
  make_option(c("-s", "--seqz_file"), type="character", default=NULL,
              help="varscan snp file name", metavar="character"),
  make_option(c("-l", "--ploidy_file"), type="character", default=NULL,
              help="ploidy file name", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default=NULL,
            help="Output prefix for result files [default= %default]", metavar="character"),
  make_option(c("-g", "--gamma"), type="integer", default=80,
            help="Gamma parameter for data extraction [default= %default]", metavar="integer"),
  make_option(c("-f", "--female"), type="logical", default=TRUE,
            help="Female Sex flag [default= %default]", metavar="logical"),
  make_option(c("-w", "--window"), type="integer", default=50,
            help="Window for detection of CNVs [default= %default]", metavar="integer"),
  make_option(c("-t", "--type"), type="character", default="PCGP",
            help="Window for detection of CNVs [default= %default]", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

FILE   <- opt$seqz_file
SAMPLE <- opt$prefix
PLOIDY <- opt$ploidy_file
GAMMA  <- opt$gamma
FEMALE <- opt$female
WINDOW <- opt$window
TYPE   <- opt$type

# ======================= PREPROCESSING ===================================
EXTR  <- sequenza.extract(FILE, gz = FALSE,
                          breaks.method = "het",
                          window = WINDOW,       # expose
                          gamma  = GAMMA,        # expose
                          min.reads.normal = 10, 
                          min.reads.baf = 1)     

print("Extract Ok")
ratio_priority = FALSE

load(file=PLOIDY) # goes to module
priors <- subset(ploidy_table,cancer_type==toupper(TYPE))[c("CN","value")]
CP.example <- sequenza.fit(EXTR,
                           priors.table = priors,
                           method = "mufreq",          
                           ratio.priority = ratio_priority,
                           chromosome.list = 1:22)

print("Fit Ok")
sequenza.results(EXTR, 
                 cp.table = CP.example, 
                 SAMPLE, 
                 out.dir = getwd(),
                 female = FEMALE,     # expose
                 CNt.max = 20,
                 ratio.priority = ratio_priority, 
                 XY = c(X = "chrX", Y = "chrY"),
                 chromosome.list = 1:22)
print("Results Ok")
