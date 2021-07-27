# sequenza

Sequenza workflow, Given a pair of cellularity and ploidy parameters, the function returns the most likely allele-specific copy numbers with the corresponding log-posterior probability of the fit, for given values of B-allele frequency and depth ratio.

## Overview

## Dependencies

* [sequenza 2.1.2m](https://sequenzatools.bitbucket.io)
* [sequenza-scripts 2.1.5m](https://github.com/oicr-gsi/sequenza)
* [sequenza-res 2.1.2](http://api.gdc.cancer.gov/data/dea893cd-9189-4091-9611-e761a1d31ebe)


## Usage

### Cromwell
```
java -jar cromwell.jar run sequenza.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`snpFile`|File|File (data file with CNV calls from Varscan).
`cnvFile`|File| File (data file with SNV calls from Varscan).


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`gammaRange`|Array[String]|["50", "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", "1250", "1500", "2000"]|List of gamma parameters for tuning Sequenza seqmentation step, used by copynumber package.
`outputFileNamePrefix`|String|""|Output prefix to prefix output file names with.


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`preprocessInputs.rScript`|String|"$RSTATS_CAIRO_ROOT/bin/Rscript"|path to Rscript
`preprocessInputs.preprocessScript`|String|"$SEQUENZA_SCRIPTS_ROOT/bin/SequenzaPreProcess_v2.2.R"|Path to the preprocessing .R script
`preprocessInputs.modules`|String|"sequenza/2.1.2m sequenza-scripts/2.1.5m"|modules needed to run preprocessing step
`preprocessInputs.timeout`|Int|20|timeout for this step in Hr, default is 20
`preprocessInputs.jobMemory`|Int|38|Memory allocated for this job
`runSequenza.rScript`|String|"$RSTATS_CAIRO_ROOT/bin/Rscript"|Path to Rscript
`runSequenza.reference`|String|"hg19"|genome assembly, hg38 etc. the default is hg19
`runSequenza.sequenzaScript`|String|"$SEQUENZA_SCRIPTS_ROOT/bin/SequenzaProcess_v2.2.R"|Sequenza wrapper script, instructions for running the pipeline
`runSequenza.ploidyFile`|String?|None|Resource used by sequenza to infer ploidy value
`runSequenza.modules`|String|"sequenza/2.1.2m sequenza-scripts/2.1.5m sequenza-res/2.1.2"|Names and versions of modules
`runSequenza.female`|String?|None|logical, TRUE or FALSE. default is TRUE
`runSequenza.cancerType`|String?|None|acronym for cancer type (from ploidy table)
`runSequenza.minReadsNormal`|Float?|None|threshold of minimum number of observation of depth ratio in a segment
`runSequenza.minReadsBaf`|Int?|None|threshold of minimum number of observation of B-allele frequency in a segment
`runSequenza.windowSize`|Int|100000|parameter to define window size for segmentation
`runSequenza.genomeSize`|Int|23|number of chromosomes in haploid genome. Default is 23
`runSequenza.timeout`|Int|20|Timeout in hours, needed to override imposed limits
`runSequenza.jobMemory`|Int|24|Memory allocated for this job
`formatJson.jobMemory`|Int|8|Memory allocated for this job
`formatJson.width`|Int|1200|width of the summary plot, default is 1200
`formatJson.height`|Int|400|height of the summary plot, default is 400
`formatJson.modules`|String|"sequenza-scripts/2.1.5m rmarkdown/0.1m"|Names and versions of modules
`formatJson.summaryPlotScript`|String|"$SEQUENZA_SCRIPTS_ROOT/bin/plot_gamma_solutions.R"|service script for plotting data from gamma solutions file, summary plot
`formatJson.sequenzaRmd`|String|"$SEQUENZA_SCRIPTS_ROOT/bin/SequenzaSummary.Rmd"|Path to rmarkdown file for producing a .pdf report
`formatJson.rScript`|String|"$RSTATS_CAIRO_ROOT/bin/Rscript"|Path to Rscript


### Outputs

Output | Type | Description
---|---|---
`resultZip`|File|All results from sequenza runs using gamma sweep.
`resultJson`|File?|Combined json file with ploidy and contamination data.
`gammaSummaryPlot`|File|png for summary plot showing the effect of different gamma values
`gammaMarkdownPdf`|File|rmarkdown pdf with all gamma-specific panels along with gamma effect summary plot


## Commands

This section lists command(s) run by sequenza workflow

* Running sequenza

Sequenza produces the most likely allele-specific copy numbers for given values of B-allele frequency and depth ratio

Preprocessing:

```
  Rscript PREPROCESS_SCRIPT -s VARSCAN_SNP_FILE -c VARSCAN_CNV_FILE -y TRUE -p PREFIX

```
Prepearing data file using Varscan results:

```
 set -euo pipefail
 Rscript SEQUENZA_SCRIPT -s SEQZ_FILE -r REFERENCE -z GENOME_SIZE -w WINDOW_SIZE -g GAMMA -p PREFIX \
            -l PLOIDY_FILE (Optional) -f FEMALE_FLAG (Optional) -t CANCER_TYPE (Optional) -n MIN_READS_NORMAL (Optional) -a MIN_READS_BAF (OPtional)
 zip -qr PREFIX_results.zip sol* PREFIX*

```

Running analysis:


```
 ...

 In this section Sequenza runs for a range of gamma values (fragment shown):

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

 ...

```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
