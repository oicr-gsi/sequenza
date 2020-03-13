version 1.0

workflow sequenza {
input {
    # Normally we need only tumor bam, normal bam may be used when availablea
    File snpFile
    File cnvFile
    String? outputFileNamePrefix = ""
}

String? sampleID = if outputFileNamePrefix=="" then basename(snpFile, ".snp") else outputFileNamePrefix

# Preprocess VarScan data
call preprocessInputs { input: snpFile = snpFile, cnvFile = cnvFile, prefix = sampleID }
# Configure and run Sequenza
call runSequenza { input: seqzFile = preprocessInputs.seqzFile }

meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "Sequenza 2.0"
}

output {
 File resultZip = runSequenza.outZip
}

}

# ==========================================
#  Prepare Sequenza data
# ==========================================
task preprocessInputs {
input {
  File snpFile
  File cnvFile
  String? prefix = "SEQUENZA"
  String? rScript = "$RSTATS_ROOT/bin/Rscript"
  String preprocessScript = "$SEQUENZA_SCRIPTS/bin/SequenzaPreProcess_v2.2.R"
  String? modules = "sequenza/2.1.2 sequenza-scripts/2.1.2"
  Int  timeout = 20
  Int? jobMemory = 10
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
  module load ~{modules} 
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
  Int? gamma = 80
  String? rScript = "$RSTATS_ROOT/bin/Rscript"
  String? prefix = "SEQUENZA"
  String sequenzaScript = "$SEQUENZA_SCRIPTS/bin/SequenzaProcess_v2.2.R"
  String ploidyFile = "$SEQUENZA_RES_ROOT/PANCAN_ASCAT_ploidy_prob.Rdata"
  String? modules = "sequenza/2.1.2 sequenza-scripts/2.1.2 sequenza-res/2.1.2"
  Int  timeout = 20
  Int? jobMemory = 10
}

parameter_meta {
 seqzFile:  ".seqz file from preprocessing step"
 gamma: "parameter for tuning Sequenza seqmentation step, used by copynumber package"
 rScript: "Path to Rscript"
 sequenzaScript: "Sequenza wrapper script, instructions for running the pipeline"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
 module load ~{modules} 
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
}
}

# ================================================
#    a placeholder for .json - generating code ()
# ================================================


