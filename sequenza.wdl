version 1.0

workflow sequenza {
input {
    # Normally we need only tumor bam, normal bam may be used when availablea
    File snpFile
    File cnvFile
    Array[String] gammaRange = ["50","100","200","300","400","500","600","700","800","900","1000","1250","1500","2000"]
    String outputFileNamePrefix = ""
}

String sampleID = if outputFileNamePrefix=="" then basename(snpFile, ".snp") else outputFileNamePrefix

# Preprocess VarScan data
call preprocessInputs { input: snpFile = snpFile, cnvFile = cnvFile, prefix = sampleID }
# Configure and run Sequenza
scatter (g in gammaRange) {
  call runSequenza { input: seqzFile = preprocessInputs.seqzFile, gamma = g }
}
# Format combined json and re-zip results in a single zip
call formatJson { input: txtPaths = runSequenza.altSolutions, zips = runSequenza.outZip,  gammaValues = runSequenza.gammaOut }

parameter_meta {
  snpFile: "File (data file with CNV calls from Varscan)."
  cnvFile: " File (data file with SNV calls from Varscan)."
  gammaRange: "List of gamma parameters for tuning Sequenza seqmentation step, used by copynumber package."
  outputFileNamePrefix: "Output prefix to prefix output file names with."
}

meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "Sequenza workflow, Given a pair of cellularity and ploidy parameters, the function returns the most likely allele-specific copy numbers with the corresponding log-posterior probability of the fit, for given values of B-allele frequency and depth ratio."
  dependencies: [
      {
        name: "sequenza/2.1.2",
        url: "https://sequenzatools.bitbucket.io"
      },
      {
        name: "sequenza-scripts/2.1.2",
        url: "https://github.com/oicr-gsi/sequenza"
      },
      {
        name: "sequenza-res/2.1.2",
        url: "http://api.gdc.cancer.gov/data/dea893cd-9189-4091-9611-e761a1d31ebe"
      }
  ]
  output_meta: {
    resultZip: "All results from sequenza runs using gamma sweep.",
    resultJson: "Combined json file with ploidy and contamination data."
  }
}

output {
 File resultZip = formatJson.combinedZip
 File resultJson = formatJson.outputJson
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
  String modules = "sequenza/2.1.2 sequenza-scripts/2.1.2"
  Int  timeout = 20
  Int jobMemory = 10
}

parameter_meta {
  snpFile: ".snp file from VarScan"
  cnvFile: ".copynumber file from Varscan"
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
  String sequenzaScript = "$SEQUENZA_SCRIPTS_ROOT/bin/SequenzaProcess_v2.2.R"
  String ploidyFile = "$SEQUENZA_RES_ROOT/PANCAN_ASCAT_ploidy_prob.Rdata"
  String modules = "sequenza/2.1.2 sequenza-scripts/2.1.2 sequenza-res/2.1.2"
  Int  timeout = 20
  Int jobMemory = 10
}

parameter_meta {
 seqzFile:  ".seqz file from preprocessing step"
 gamma: "parameter for tuning Sequenza seqmentation step, used by copynumber package"
 rScript: "Path to Rscript"
 sequenzaScript: "Sequenza wrapper script, instructions for running the pipeline"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
 jobMemory: "Memory allocated for this job"
}

command <<<
 ~{rScript} ~{sequenzaScript} -s ~{seqzFile} -l ~{ploidyFile} -p ~{prefix}
 zip -qj "~{prefix}_results.zip" ~{prefix}*
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File outZip = "~{prefix}_results.zip"
  String gammaOut = "~{gamma}"
  File altSolutions = "~{prefix}_alternative_solutions.txt"
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
}

parameter_meta {
 txtPaths: "List of files which need to be processed"
 zips: "List of zip files from runSequenza"
 gammaValues: "List of gamma values for the used range"
 jobMemory: "Memory allocated for this job"
}

command <<<
 python <<CODE
 import json
 import os
 pts = "~{sep=' ' txtPaths}"
 paths = pts.split()
 gms = "~{sep=' ' gammaValues}"
 gammas = gms.split()
 zps = "~{sep=' ' zips}"
 zips = zps.split()
 json_name = "~{prefix}_alternative_solutions.json"
 jsonDict = {}

 for p in range(0, len(paths)):
     gamma = gammas[p]
     os.popen("mkdir -p gammas/" + gamma)
     os.popen("unzip " + zips[p] + " -d gammas/" + gamma + "/")
     chunk = {
         'cellularity': [],
         'ploidy': [],
         'SLPP': []
     }
     with open(paths[p]) as f:
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

 with open(json_name, 'w') as json_file:
     json.dump(jsonDict, json_file)
 CODE
 zip -qr "~{prefix}_results.zip" gammas/*
>>>

runtime {
  memory:  "~{jobMemory} GB"
}

output {
  File outputJson = "~{prefix}_alternative_solutions.json"
  File combinedZip = "~{prefix}_results.zip"
}
}
