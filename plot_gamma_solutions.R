library(ggplot2)
library(optparse)


## For now, test if commands are in original, trailing format, or new opt-parse format
option_list = list(
  make_option(c("-f", "--solutionsFile"), type="character", default=NULL,
              help="json file with gamma solutions from sequenza", metavar="character"),
  make_option(c("-o", "--outputFile"), type="character", default=NULL,
              help="output png file", metavar="character"),
  make_option(c("-w", "--width"), type="integer", default="1440",
              help="width of png file", metavar="integer"),
  make_option(c("-h", "--height"), type="integer", default="440",
              help="height of png file", metavar="integer"))

opt_parser = OptionParser(option_list = option_list, add_help_option = FALSE)
opt = parse_args(opt_parser);

solutionsFile = opt$solutionsFile
outFile = opt$outputFile
gamma.solutions<-read.csv(solutionsFile)

df.cellularity<-data.frame(gamma = gamma.solutions$gamma, value = gamma.solutions$cellularity, type="cellularity")
df.ploidy<-data.frame(gamma = gamma.solutions$gamma, value = gamma.solutions$ploidy, type="ploidy")
df.segmentation<-data.frame(gamma = gamma.solutions$gamma, value=gamma.solutions$no_segments, type="no. segments")
df<-rbind(df.segmentation, df.cellularity, df.ploidy)
p<-ggplot(df, aes(x = gamma, y = value))+
  geom_line(linetype = "dashed", aes(colour = type, group = type), size = 1)+
  geom_point(aes(colour = type, group = type),size = 4)+
  scale_x_continuous(breaks = as.numeric(gamma.solutions$gamma)) +
  facet_grid(rows = vars(type), scales="free")+
  theme_bw(base_size = 18)+
  theme(legend.position = "none")

png(filename = file.path(outFile), type="cairo", width = opt$width, height = opt$height)
print(p)   
dev.off() 
