# sequenza
Sequenza workflow, Given a pair of cellularity and ploidy parameters, the function returns the most likely allele-specific
copy numbers with the corresponding log-posterior probability of the fit, for given values of B-allele
frequency and depth ratio.

## Usage

## Cromwell

``` 
 java -jar cromwell.jar run sequenza.wdl --inputs inputs.json 

```

## Running Pipeline

Pipeline is run with two wrapper scripts (one for pre-processing data and another - for running the analysis)

```
 Rscript SequenzaPreProcess_v2.2.R -s [.snp file] -c [.copynumber file] -y [female flag] -p [prefix]

 Rscript SequenzaProcess_v2.2.R -s [.seqz file] -l [ploidy file] -p [prefix]

```

