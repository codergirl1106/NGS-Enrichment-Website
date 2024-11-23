## ----setup, echo=FALSE, results="hide", include = FALSE-------------------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=FALSE,
                      #dev="png",
                      message=FALSE, error=FALSE, warning=TRUE)
options(width=100)

## ----getPackage, eval=FALSE-----------------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("Rfastp")

## ----eval = FALSE---------------------------------------------------------------------------------
# BiocManager::install("remotes")
# BiocManager::install("RockefellerUniversity/Rfastp")

## -------------------------------------------------------------------------------------------------
library(Rfastp)

## -------------------------------------------------------------------------------------------------
se_read1 <- system.file("extdata","Fox3_Std_small.fq.gz",package="Rfastp")
pe_read1 <- system.file("extdata","reads1.fastq.gz",package="Rfastp")
pe_read2 <- system.file("extdata","reads2.fastq.gz",package="Rfastp")
outputPrefix <- tempfile(tmpdir = tempdir())

## -------------------------------------------------------------------------------------------------
se_json_report <- rfastp(read1 = se_read1, 
    outputFastq = paste0(outputPrefix, "_se"), thread = 4)

## -------------------------------------------------------------------------------------------------
pe_json_report <- rfastp(read1 = pe_read1, read2 = pe_read2,
    outputFastq = paste0(outputPrefix, "_pe"))

## -------------------------------------------------------------------------------------------------
pe_merge_json_report <- rfastp(read1 = pe_read1, read2 = pe_read2, merge = TRUE,
    outputFastq = paste0(outputPrefix, '_unpaired'),
    mergeOut = paste0(outputPrefix, "_merged.fastq.gz"))

## -------------------------------------------------------------------------------------------------
umi_json_report <- rfastp(read1 = pe_read1, read2 = pe_read2, 
    outputFastq = paste0(outputPrefix, '_umi1'), umi = TRUE, umiLoc = "read1",
    umiLength = 16)

## ----umi------------------------------------------------------------------------------------------
umi_json_report <- rfastp(read1 = pe_read1, read2 = pe_read2, 
    outputFastq = paste0(outputPrefix, '_umi2'), umi = TRUE, umiLoc = "read1",
    umiLength = 16, umiPrefix = "#", umiNoConnection = TRUE, 
    umiIgnoreSeqNameSpace = TRUE)

## -------------------------------------------------------------------------------------------------
clipr_json_report <- rfastp(read1 = se_read1, 
    outputFastq = paste0(outputPrefix, '_clipr'),
    disableTrimPolyG = TRUE,
    cutLowQualFront = TRUE,
    cutFrontWindowSize = 29,
    cutFrontMeanQual = 20,
    cutLowQualTail = TRUE,
    cutTailWindowSize = 1,
    cutTailMeanQual = 5,
    minReadLength = 29,
    adapterSequenceRead1 = 'GTGTCAGTCACTTCCAGCGG'
)

## -------------------------------------------------------------------------------------------------
pe001_read1 <- system.file("extdata","splited_001_R1.fastq.gz",
    package="Rfastp")
pe002_read1 <- system.file("extdata","splited_002_R1.fastq.gz",
    package="Rfastp")
pe003_read1 <- system.file("extdata","splited_003_R1.fastq.gz",
    package="Rfastp")
pe004_read1 <- system.file("extdata","splited_004_R1.fastq.gz",
    package="Rfastp")
inputfiles <- c(pe001_read1, pe002_read1, pe003_read1, pe004_read1)
cat_rjson_report <- rfastp(read1 = inputfiles, 
    outputFastq = paste0(outputPrefix, "_merged1"))

## -------------------------------------------------------------------------------------------------
pe001_read2 <- system.file("extdata","splited_001_R2.fastq.gz",
    package="Rfastp")
pe002_read2 <- system.file("extdata","splited_002_R2.fastq.gz",
    package="Rfastp")
pe003_read2 <- system.file("extdata","splited_003_R2.fastq.gz",
    package="Rfastp")
pe004_read2 <- system.file("extdata","splited_004_R2.fastq.gz",
    package="Rfastp")
inputR2files <- c(pe001_read2, pe002_read2, pe003_read2, pe004_read2)
catfastq(output = paste0(outputPrefix,"_merged2_R2.fastq.gz"), 
    inputFiles = inputR2files)

## -------------------------------------------------------------------------------------------------
dfsummary <- qcSummary(pe_json_report)

## -------------------------------------------------------------------------------------------------
p1 <- curvePlot(se_json_report)
p1

## -------------------------------------------------------------------------------------------------
p2 <- curvePlot(se_json_report, curve="content_curves")
p2

## -------------------------------------------------------------------------------------------------
dfTrim <- trimSummary(pe_json_report)

## -------------------------------------------------------------------------------------------------
?rfastp

## -------------------------------------------------------------------------------------------------
?catfastq

## -------------------------------------------------------------------------------------------------
?qcSummary

## -------------------------------------------------------------------------------------------------
?trimSummary

## -------------------------------------------------------------------------------------------------
?curvePlot

## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

