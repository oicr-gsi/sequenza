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
  make_option(c("-t", "--type"), type="character", default="PCPG",
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

# ======================= OVERRIDE THIS SEQUENZA FUNCTION =================

sequenza.extract <- function(file, gz = TRUE, window = 1e6, overlap = 1, gamma = 80, kmin = 10,
                             gamma.pcf = 140, kmin.pcf = 40, mufreq.treshold = 0.10, min.reads = 40,
                             min.reads.normal = 10, min.reads.baf = 1, max.mut.types = 1,
                             min.type.freq = 0.9, min.fw.freq = 0, verbose = TRUE, chromosome.list = NULL,
                             breaks = NULL, breaks.method = "het", assembly = "hg19", weighted.mean = TRUE,
                             normalization.method = "mean", gc.stats = NULL) {
  
  if (is.null(gc.stats)) {
    gc.stats <- gc.sample.stats(file, gz = gz)
  }
  chr.vect <- as.character(gc.stats$file.metrics$chr)
  if (normalization.method != "mean") {
    gc.vect  <- setNames(gc.stats$raw.median, gc.stats$gc.values)
  } else {
    gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
  }
  windows.baf   <- list()
  windows.ratio <- list()
  mutation.list <- list()
  segments.list <- list()
  coverage.list <- list()
  if (is.null(dim(breaks))) {
    breaks.all <- NULL
  } else {
    breaks.all <- breaks
  }
  if (is.null(chromosome.list)) {
    chromosome.list <- chr.vect
  } else {
    chromosome.list <- chromosome.list[chromosome.list %in% chr.vect]
  }
  for (chr in chromosome.list){
    if (verbose){
      message("Processing ", chr, ": ", appendLF = FALSE) 
    }
    file.lines <- gc.stats$file.metrics[which(chr.vect == chr), ]
    seqz.data   <- read.seqz(file, gz = gz, n.lines = c(file.lines$start, file.lines$end))
    seqz.data$adjusted.ratio <- round(seqz.data$depth.ratio / gc.vect[as.character(seqz.data$GC.percent)], 3)
    seqz.hom <- seqz.data$zygosity.normal == 'hom'
    seqz.het <- seqz.data[!seqz.hom, ]
    het.filt <- seqz.het$good.reads >= min.reads.baf
    seqz.r.win <- windowValues(x = seqz.data$adjusted.ratio,
                               positions = seqz.data$position,
                               chromosomes = seqz.data$chromosome,
                               window = window, overlap = overlap,
                               weight = seqz.data$depth.normal)
    if (nrow(seqz.het) > 0) {
      breaks.method.i <- breaks.method
      
      seqz.b.win <- windowValues(x = seqz.het$Bf[het.filt],
                                 positions = seqz.het$position[het.filt],
                                 chromosomes = seqz.het$chromosome[het.filt],
                                 window = window, overlap = overlap,
                                 weight = seqz.het$good.reads[het.filt])
      if (is.null(breaks.all)){
        if (breaks.method.i == "full") {
          breaks <- find.breaks(seqz.data, gamma = gamma.pcf, assembly = assembly, 
                                kmin = kmin.pcf, seg.algo = "pcf")
          breaks.het <- try(find.breaks(seqz.het, gamma = gamma, assembly = assembly,
                                        kmin = kmin, baf.thres = c(0, 0.5)),
                            silent = FALSE)
          if (!is.null(breaks.het)) {
            merge.breaks <- function (breaks, breaks.het) {
              merged.breaks <- unique(sort(c(breaks$start.pos, breaks$end.pos, breaks.het$start.pos, breaks.het$end.pos)))
              merged.breaks <- merged.breaks[diff(merged.breaks) > 1]
              merged.start <- merged.breaks
              merged.start[-1] <- merged.start[-1]+1
              breaks <- data.frame(chrom = unique(breaks$chrom),
                                   start.pos = merged.start[-(length(merged.start))],
                                   end.pos = merged.breaks[-1])
            }
            chr.p <- merge.breaks(breaks[breaks$arm == "p",], breaks.het[breaks.het$arm == "p",])
            chr.q <- merge.breaks(breaks[breaks$arm == "q",], breaks.het[breaks.het$arm == "q",])
            breaks <- rbind(chr.p, chr.q)
          }
        } else if (breaks.method.i == "het"){
          breaks <- try(find.breaks(seqz.het, gamma = gamma, assembly = assembly,
                                    kmin = kmin, baf.thres = c(0, 0.5)),
                        silent = FALSE)               
        } else if (breaks.method.i == "fast"){
          BAF <- data.frame(chrom = chr, pos = c(seqz.b.win[[1]]$start, tail(seqz.b.win[[1]]$end, n = 1)),
                            s1 = c(seqz.b.win[[1]]$mean, tail(seqz.b.win[[1]]$mean, n = 1)))
          BAF$s1[is.na(BAF$s1)] <- 0
          logR <- data.frame(chrom = chr, pos = c(seqz.r.win[[1]]$start, tail(seqz.r.win[[1]]$end, n = 1)),
                             s1 = c(log2(seqz.r.win[[1]]$mean), log2(tail(seqz.r.win[[1]]$mean, n = 1))))
          not.cover <- is.na(logR$s1)
          BAF  <- BAF[!not.cover, ]
          logR <- logR[!not.cover, ]
          logR.wins <- copynumber::winsorize(logR, verbose = FALSE)
          allele.seg <- copynumber::aspcf(logR = logR.wins, BAF = BAF, baf.thres = c(0, 0.5),
                                          verbose = FALSE, gamma = gamma, kmin = kmin)
          if (length(grep("chr", chr)) > 0) {
            allele.seg$chrom <- paste("chr", allele.seg$chrom, sep = "")
          }
          breaks   <- allele.seg[, c("chrom", "start.pos", "end.pos")]
          not.uniq <- which(breaks$end.pos == c(breaks$start.pos[-1],0))
          breaks$end.pos[not.uniq] <- breaks$end.pos[not.uniq] - 1
        } else {
          stop("The implemented segmentation methods are \'full\', \'het\' and \'fast\'.")
        }
      } else {
        breaks <- breaks.all[breaks.all$chrom == chr, ]
      }
      if (!is.null(breaks) && class(breaks) == "data.frame" && nrow(breaks) > 0){
        seg.s1    <- segment.breaks(seqz.tab = seqz.data, breaks = breaks,
                                    min.reads.baf = min.reads.baf, weighted.mean = weighted.mean)            
      } else {
        seg.s1 <- segment.breaks(seqz.data,
                                 breaks = data.frame(chrom = chr,
                                                     start.pos = min(seqz.data$position, na.rm = TRUE),
                                                     end.pos = max(seqz.data$position, na.rm = TRUE)),
                                 weighted.mean = weighted.mean)
      }
      
    } else {
      seqz.b.win <- list()
      seqz.b.win[[1]] <- data.frame(start = min(seqz.data$position, na.rm = TRUE),
                                    end = max(seqz.data$position, na.rm = TRUE), mean = 0.5,
                                    q0 = 0.5,  q1 = 0.5, N = 1)
      if (breaks.method == "full") {
        breaks <- find.breaks(seqz.data, gamma = gamma.pcf, assembly = assembly,
                              kmin = kmin.pcf, seg.algo = "pcf")
      } else {
        breaks = data.frame(chrom = chr,
                            start.pos = min(seqz.data$position, na.rm = TRUE),
                            end.pos = max(seqz.data$position, na.rm = TRUE))
      }
      seg.s1 <- segment.breaks(seqz.data,
                               breaks = breaks,
                               weighted.mean = weighted.mean)
    }
    mut.tab   <- mutation.table(seqz.data, mufreq.treshold = mufreq.treshold,
                                min.reads = min.reads, min.reads.normal = min.reads.normal,
                                max.mut.types = max.mut.types, min.type.freq = min.type.freq,
                                min.fw.freq = min.fw.freq, segments = seg.s1)
    windows.baf[[which(chromosome.list == chr)]]   <- seqz.b.win[[1]]
    windows.ratio[[which(chromosome.list == chr)]] <- seqz.r.win[[1]]
    mutation.list[[which(chromosome.list == chr)]] <- mut.tab
    segments.list[[which(chromosome.list == chr)]] <- seg.s1
    coverage.list[[which(chromosome.list == chr)]] <- data.frame(sum = sum(as.numeric(seqz.data$depth.tumor),
                                                                           na.rm = TRUE),
                                                                 N  = length(seqz.data$depth.tumor) )
    if (verbose){
      
      message(nrow(mut.tab), ' variant calls; ',
              nrow(seqz.het), ' heterozygous positions; ',
              sum(seqz.hom), ' homozygous positions.') 
    }
  }
  names(windows.baf)   <- chromosome.list
  names(windows.ratio) <- chromosome.list
  names(mutation.list) <- chromosome.list
  names(segments.list) <- chromosome.list
  coverage.list <- do.call(rbind, coverage.list)
  coverage <- sum(coverage.list$sum) / sum(coverage.list$N)
  return(list(BAF = windows.baf, ratio = windows.ratio, mutations = mutation.list,
              segments = segments.list, chromosomes = chromosome.list, gc = gc.stats,
              avg.depth = round(coverage,0)))
}


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
