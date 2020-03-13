### Generate sequenza ploidy priors table based on PANCAN Ascat results
### Script from Jeffrey Bruce, modified by Peter Ruzanov

library(openxlsx)
cmd_args=commandArgs(trailingOnly = TRUE)
FILE    = cmd_args[1]
purity_ploidy<-read.xlsx(
  FILE,
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
  detectDates = FALSE,
  skipEmptyRows = TRUE,
  skipEmptyCols = TRUE,
  rows = NULL,
  cols = NULL,
  check.names = FALSE,
  namedRegion = NULL,
  na.strings = "NA",
  fillMergedCells = FALSE
)

cancer_types <- unique(purity_ploidy$tumor_type)

# LAML 	Acute Myeloid Leukemia
# ACC 	Adrenocortical carcinoma
# BLCA 	Bladder Urothelial Carcinoma
# LGG 	Brain Lower Grade Glioma
# BRCA 	Breast invasive carcinoma
# CESC 	Cervical squamous cell carcinoma and endocervical adenocarcinoma
# CHOL 	Cholangiocarcinoma
# LCML 	Chronic Myelogenous Leukemia
# COAD 	Colon adenocarcinoma
# CNTL 	Controls
# ESCA 	Esophageal carcinoma
# FPPP 	FFPE Pilot Phase II
# GBM 	Glioblastoma multiforme
# HNSC 	Head and Neck squamous cell carcinoma
# KICH 	Kidney Chromophobe
# KIRC 	Kidney renal clear cell carcinoma
# KIRP 	Kidney renal papillary cell carcinoma
# LIHC 	Liver hepatocellular carcinoma
# LUAD 	Lung adenocarcinoma
# LUSC 	Lung squamous cell carcinoma
# DLBC 	Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
# MESO 	Mesothelioma
# MISC 	Miscellaneous
# OV 	Ovarian serous cystadenocarcinoma
# PAAD 	Pancreatic adenocarcinoma
# PCPG 	Pheochromocytoma and Paraganglioma
# PRAD 	Prostate adenocarcinoma
# READ 	Rectum adenocarcinoma
# SARC 	Sarcoma
# SKCM 	Skin Cutaneous Melanoma
# STAD 	Stomach adenocarcinoma
# TGCT 	Testicular Germ Cell Tumors
# THYM 	Thymoma
# THCA 	Thyroid carcinoma
# UCS 	Uterine Carcinosarcoma
# UCEC 	Uterine Corpus Endometrial Carcinoma
# UVM 	Uveal Melanoma

ploidy_table<- data.frame(CN = "", 
                          value = "",
                          cancer_type="")[0,]

for (cancer in cancer_types[which(cancer_types != "NA")]) {
  tab <- purity_ploidy[purity_ploidy$tumor_type==cancer,]
  num <- nrow(tab)
  
  tmp_frame <- data.frame(CN = 0:7, 
                          value = rep(0,8),
                          cancer_type=rep(cancer,8))
  
  for(i in 0:7){
    tmp_frame[i+1,"value"] <- sum(tab$ploidy>(i-0.5) & 
                                  tab$ploidy <= (i+0.5),na.rm=TRUE)/num
  }
  
  ploidy_table <- rbind(ploidy_table,tmp_frame)
}

## add the percentages for all cancers
tmp_frame <- data.frame(CN = 0:7, 
                        value = rep(0,8),
                        cancer_type=rep("all",8))

num <- nrow(purity_ploidy)
for(i in 0:7){
  tmp_frame[i+1,"value"] <- sum(purity_ploidy$ploidy>(i-0.5) & 
                                purity_ploidy$ploidy <= (i+0.5),na.rm=TRUE)/num
}

ploidy_table <- rbind(ploidy_table,tmp_frame)


save(ploidy_table,
     file="PANCAN_ASCAT_ploidy_prob.Rdata")

