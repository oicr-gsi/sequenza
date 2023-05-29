version 1.0


struct GenomeResources {
    Int genomeSize
    String ploidyFile
}

workflow sequenza {
input {
    # Normally we need only tumor bam, normal bam may be used when availablea
    File snpFile
    File cnvFile
    Array[String] gammaRange = ["50","100","200","300","400","500","600","700","800","900","1000","1250","1500","2000"]
    String outputFileNamePrefix = ""
    String reference
}

Map[String, GenomeResources] resources = {
  "hg19": {
    "genomeSize": 23,
    "ploidyFile": "$SEQUENZA_RES_ROOT/PANCAN_ASCAT_ploidy_prob.Rdata"
  },
  "hg38": {
    "genomeSize": 23,
    "ploidyFile": "$SEQUENZA_RES_ROOT/PANCAN_ASCAT_ploidy_prob.Rdata"
  },
  "mm9": {
    "genomeSize": 20,
    "ploidyFile": ""
  },
  "mm10": {
    "genomeSize": 20,
    "ploidyFile": ""
  }
}

String sampleID = if outputFileNamePrefix=="" then basename(snpFile, ".snp") else outputFileNamePrefix

# Preprocess VarScan data
call preprocessInputs { input: snpFile = snpFile, cnvFile = cnvFile, prefix = sampleID }
# Configure and run Sequenza
scatter (g in gammaRange) {
  call runSequenza { input: seqzFile = preprocessInputs.seqzFile, prefix = sampleID, gamma = g, reference = reference, ploidyFile = resources[reference].ploidyFile, genomeSize = resources[reference].genomeSize }
}
# Format combined json and re-zip results in a single zip
call formatJson { input: txtPaths = select_all(runSequenza.altSolutions), zips = runSequenza.outZip,  gammaValues = runSequenza.gammaOut, prefix = sampleID }

parameter_meta {
  snpFile: "File (data file with CNV calls from Varscan)."
  cnvFile: " File (data file with SNV calls from Varscan)."
  gammaRange: "List of gamma parameters for tuning Sequenza seqmentation step, used by copynumber package."
  outputFileNamePrefix: "Output prefix to prefix output file names with."
  reference: "Version of genome reference"
}

meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "Sequenza workflow, Given a pair of cellularity and ploidy parameters, the function returns the most likely allele-specific copy numbers with the corresponding log-posterior probability of the fit, for given values of B-allele frequency and depth ratio.\nSequenza workflow, Given a pair of cellularity and ploidy parameters, the function returns the most likely allele-specific copy numbers with the corresponding log-posterior probability of the fit, for given values of B-allele frequency and depth ratio.\n## Overview\n![sequenza outputs](docs/Screenshot_Sequenza_PDFs.png)"
  dependencies: [
      {
        name: "sequenza/2.1.2m",
        url: "https://sequenzatools.bitbucket.io"
      },
      {
        name: "sequenza-scripts/2.1.5m",
        url: "https://github.com/oicr-gsi/sequenza"
      },
      {
        name: "sequenza-res/2.1.2",
        url: "http://api.gdc.cancer.gov/data/dea893cd-9189-4091-9611-e761a1d31ebe"
      }
  ]
  output_meta: {
    resultZip: "All results from sequenza runs using gamma sweep.",
    resultJson: "Combined json file with ploidy and contamination data.",
    gammaSummaryPlot: "png for summary plot showing the effect of different gamma values",
    gammaMarkdownPdf: "rmarkdown pdf with all gamma-specific panels along with gamma effect summary plot"
  }
}

output {
 File resultZip = formatJson.combinedZip
 File? resultJson = formatJson.outputJson
 File gammaSummaryPlot = formatJson.gammaSummaryPlot
 File gammaMarkdownPdf = formatJson.gammaMarkdownPdf
}

}

# ==========================================
#  Prepare Sequenza data
# ==========================================
task preprocessInputs {
input {
  File snpFile
  File cnvFile
  String prefix = "SEQUENZA"
  String rScript = "$RSTATS_ROOT/bin/Rscript"
  String preprocessScript = "$SEQUENZA_SCRIPTS_ROOT/bin/SequenzaPreProcess_v2.2.R"
  String modules = "sequenza/2.1.2m sequenza-scripts/2.1.5m"
  Int  timeout = 20
  Int jobMemory = 38
}

parameter_meta {
  snpFile: ".snp file from VarScan"
  cnvFile: ".copynumber file from Varscan"
  prefix: "prefix for the output file name"
  rScript: "path to Rscript"
  preprocessScript: "Path to the preprocessing .R script"
  timeout: "timeout for this step in Hr, default is 20"
  modules: "modules needed to run preprocessing step"
  jobMemory: "Memory allocated for this job"
}

command <<<
  ~{rScript} ~{preprocessScript} -s ~{snpFile} -c ~{cnvFile} -y TRUE -p ~{prefix}
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File seqzFile = "~{prefix}.seqz"
}
}


# ==========================================
#  configure and run Sequenza
# ==========================================
task runSequenza {
input {
  File seqzFile
  # Parameters
  String gamma = "80"
  String rScript = "$RSTATS_ROOT/bin/Rscript"
  String prefix = "SEQUENZA"
  String reference 
  String sequenzaScript = "$SEQUENZA_SCRIPTS_ROOT/bin/SequenzaProcess_v2.2.R"
  String? ploidyFile
  String modules = "sequenza/2.1.2m sequenza-scripts/2.1.5m sequenza-res/2.1.2"
  String? female
  String? cancerType
  Float? minReadsNormal
  Int? minReadsBaf
  Int windowSize = 100000
  Int genomeSize = 23
  Int timeout = 20
  Int jobMemory = 24
}

parameter_meta {
 seqzFile:  ".seqz file from preprocessing step"
 gamma: "parameter for tuning Sequenza seqmentation step, used by copynumber package"
 windowSize: "parameter to define window size for segmentation"
 female: "logical, TRUE or FALSE. default is TRUE"
 cancerType: "acronym for cancer type (from ploidy table)"
 minReadsNormal: "threshold of minimum number of observation of depth ratio in a segment"
 minReadsBaf: "threshold of minimum number of observation of B-allele frequency in a segment"
 rScript: "Path to Rscript"
 sequenzaScript: "Sequenza wrapper script, instructions for running the pipeline"
 ploidyFile: "Resource used by sequenza to infer ploidy value"
 prefix: "prefix for the output file name"
 reference: "genome assembly, hg38 etc. the default is hg19"
 genomeSize: "number of chromosomes in haploid genome. Default is 23"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
 jobMemory: "Memory allocated for this job"
}

command <<<
 set -euo pipefail
 ~{rScript} ~{sequenzaScript} -s ~{seqzFile} -r ~{reference} -z ~{genomeSize} -w ~{windowSize} -g ~{gamma} -p ~{prefix} \
            ~{"-l " + ploidyFile} ~{"-f " + female} ~{"-t " + cancerType} ~{"-n " + minReadsNormal} ~{"-a " + minReadsBaf}
 zip -qr ~{prefix}_results.zip sol* ~{prefix}*
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File outZip = "~{prefix}_results.zip"
  String gammaOut = "~{gamma}"
  File? altSolutions = "~{prefix}_alternative_solutions.txt"
}
}

# ============================================================
#   .json - generating code, produce a combined report here
# ============================================================
task formatJson {
input {
  String prefix = "SEQUENZA"
  Array[File] txtPaths
  Array[File] zips
  Array[String] gammaValues
  Int jobMemory = 8
  Int width = 1200
  Int height = 400
  String modules = "sequenza-scripts/2.1.5m rmarkdown/0.1m"
  String summaryPlotScript = "$SEQUENZA_SCRIPTS_ROOT/bin/plot_gamma_solutions.R"
  String sequenzaRmd = "$SEQUENZA_SCRIPTS_ROOT/bin/SequenzaSummary.Rmd"
  String rScript = "$RSTATS_ROOT/bin/Rscript"
}

parameter_meta {
 prefix: "a String for prefixing all result files"
 txtPaths: "List of files which need to be processed"
 zips: "List of zip files from runSequenza"
 gammaValues: "List of gamma values for the used range"
 jobMemory: "Memory allocated for this job"
 width: "width of the summary plot, default is 1200"
 height: "height of the summary plot, default is 400"
 summaryPlotScript: "service script for plotting data from gamma solutions file, summary plot"
 sequenzaRmd: "Path to rmarkdown file for producing a .pdf report"
 rScript: "Path to Rscript"
 modules: "Names and versions of modules"
}

command <<<
 set -euo pipefail
 python3<<CODE
 import json
 import os
 import pandas as pd
 from os.path import isfile
 pts = "~{sep=' ' txtPaths}"
 paths = pts.split()
 gms = "~{sep=' ' gammaValues}"
 gammas = gms.split()
 zps = "~{sep=' ' zips}"
 zips = zps.split()
 json_name = "~{prefix}_alternative_solutions.sequenza.json"
 jsonDict = {}

 for g in range(0, len(gammas)):
   gamma = gammas[g]
   os.system("mkdir -p gammas/" + gamma)
   os.system("unzip " + zips[g] + " -d gammas/" + gamma + "/")
   chunk = {
     'cellularity': [],
     'ploidy': [],
     'SLPP': []
   }
   if not isfile(paths[g]):
     continue
   with open(paths[g]) as f:
     for line in f:
       if line.find("cellularity") > 0:
         continue
       line = line.rstrip()
       tmp = line.split("\t")
       chunk['cellularity'].append(tmp[0])
       chunk['ploidy'].append(tmp[1])
       chunk['SLPP'].append(tmp[2])
   f.close()
   jsonDict[gamma] = chunk
 
 if len(jsonDict.keys()) > 0: 
   with open(json_name, 'w') as json_file:
     json.dump(jsonDict, json_file)

 cellularity = []
 ploidy = []
 no_segments = []

 for g in gammas:
   print(g)
   solutions = pd.read_table(os.path.join("gammas", g, "~{prefix}" + '_alternative_solutions.txt'))
   row = solutions.loc[solutions['SLPP'].idxmax()]
   cellularity.append(float(row['cellularity']))
   ploidy.append(float(row['ploidy']))
   path_seg = os.path.join("gammas", g, "~{prefix}" + '_Total_CN.seg')
   no_segments.append(len(open(path_seg).readlines()) - 1)
 
 gamma_solutions = pd.DataFrame({"gamma": gammas,
                                 "cellularity": cellularity,
                                 "ploidy": ploidy,
                                 "no_segments": no_segments})
 gamma_solutions.to_csv('gamma_solutions.csv', index=False)
 CODE

 ~{rScript} ~{summaryPlotScript} -f gamma_solutions.csv -o ~{prefix}_summary.sequenza.png -w ~{width} -h ~{height}
 cp ~{sequenzaRmd} .
 ~{rScript} -e "rmarkdown::render('SequenzaSummary.Rmd', params = list(sample = '~{prefix}',summaryImage = '~{prefix}_summary.sequenza.png'))"
 mv SequenzaSummary.pdf ~{prefix}_summary.sequenza.pdf
 zip -qr ~{prefix}_results.sequenza.zip gammas/*
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File? outputJson = "~{prefix}_alternative_solutions.sequenza.json"
  File gammaSummaryPlot = "~{prefix}_summary.sequenza.png"
  File gammaMarkdownPdf = "~{prefix}_summary.sequenza.pdf"
  File combinedZip = "~{prefix}_results.sequenza.zip"
}
}
