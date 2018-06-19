if (!any(rownames(installed.packages()) == "rjson")) install.packages("rjson")
require(rjson)

samples <- readLines("samples.txt")
hts_dir <- "01-HTS_Preproc"

hts_json <- lapply(samples, function(s) {
  fromJSON(paste(readLines(file.path(hts_dir,s,paste0(s,"_htsStats.log"))),collapse=""))
})

if (length(hts_json) != length(samples)) stop("not all samples have json log files")
if (!all(apply(sapply(hts_json,names),1, function(x) length(unique(sub("_[0-9]+$","",x)))) == 1)) stop("not all json log files have the same htsteam sub apps")

# sub("_[0-9]+$","",x)
apps <- apply(sapply(hts_json,names),1, function(x) unique(sub("_[0-9]+$","",x)))
napps <- length(apps)
apps <- make.names(apps, unique = T)

check_apps <- c("hts_Stats", "hts_SeqScreener", "hts_SuperDeduper", "hts_SeqScreener.1", "hts_AdapterTrimmer", "hts_QWindowTrim", "hts_NTrimmer", "hts_CutTrim", "hts_Stats.1")
if (!(length(check_apps) == length(apps) && all(check_apps == apps))) stop("expected sequence of apps not found")

## Table of total reads output after each stage
#totalReadsOut <- as.data.frame(sapply(hts_json,function(x) sapply(x, "[[","totalFragmentsOutput")))
#rownames(totalReadsOut) <- apps
#colnames(totalReadsOut) <- samples

### Super Deduper Saturation Plot
# sd <- which(apps == "hts_SuperDeduper")
# if(length(sd) > 1) stop("More than one hts super deduper found, which makes little sense")
#
# if (length(sd) == 1){
#   ds <- lapply(hts_json,function(x) x[[sd]]$duplicate_saturation)
#   plot(cumsum(diff(sapply(ds[[2]],"[[",1L))),diff(sapply(ds[[2]],"[[",2L))/diff(sapply(ds[[2]],"[[",1L)),type='l')
# }

## Pre-Stats
statsTotalReadsIn <- sapply(hts_json,function(x) sapply(x[1], "[[","totalFragmentsInput"))
statsTotalReadsCG <- sapply(hts_json,function(x) sum(unlist(sapply(x[1],"[[","Base_composition")[c("C","G"),])))
statsTotalReadsN <- sapply(hts_json,function(x) unlist(sapply(x[1],"[[","Base_composition")[c("N"),]))
statsTotalReadsR1_BpLen <- sapply(hts_json,function(x) unlist(sapply(x[1],"[[","Paired_end")[c("R1_bpLen"),]))
statsTotalReadsR1_BQ30 <- sapply(hts_json,function(x) unlist(sapply(x[1],"[[","Paired_end")[c("R1_bQ30"),]))
statsTotalReadsR2_BpLen <- sapply(hts_json,function(x) unlist(sapply(x[1],"[[","Paired_end")[c("R2_bpLen"),]))
statsTotalReadsR2_BQ30 <- sapply(hts_json,function(x) unlist(sapply(x[1],"[[","Paired_end")[c("R2_bQ30"),]))

## PhiX screen
seqScrTotalReadsIn <- sapply(hts_json,function(x) sapply(x[2], "[[","totalFragmentsInput"))
seqScrTotalReadsHits <- sapply(hts_json,function(x) unlist(sapply(x[2],"[[","Paired_end")[c("PE_hits"),]))

## Super Deduper PCR duplicates
SdTotalReadsIn <- sapply(hts_json,function(x) sapply(x[3], "[[","totalFragmentsInput"))
SdTotalReadsIgnored <- sapply(hts_json,function(x) sapply(x[3],"[[","ignored"))
SdTotalReadsDuplicate <- sapply(hts_json,function(x) sapply(x[3],"[[","duplicate"))

## Seq Screen rRNA
seqScr2TotalReadsIn <- sapply(hts_json,function(x) sapply(x[4], "[[","totalFragmentsInput"))
seqScr2TotalReadsHits <- sapply(hts_json,function(x) unlist(sapply(x[4],"[[","Paired_end")[c("PE_hits"),]))

## Adapter Trimmer
AdaptTotalReadsIn <- sapply(hts_json,function(x) sapply(x[5], "[[","totalFragmentsInput"))
AdaptTotalAdaptTrim <- sapply(hts_json,function(x) unlist(sapply(x[5],"[[","Paired_end")[c("PE_adapterTrim"),]))
AdaptTotalBpTrim <- sapply(hts_json,function(x) unlist(sapply(x[5],"[[","Paired_end")[c("PE_adapterBpTrim"),]))

## Q windowtrim
QwinTotalReadsIn <- sapply(hts_json,function(x) sapply(x[6], "[[","totalFragmentsInput"))
QwinTotalR1LT <- sapply(hts_json,function(x) unlist(sapply(x[6],"[[","Paired_end")[c("R1_leftTrim"),]))
QwinTotalR1RT <- sapply(hts_json,function(x) unlist(sapply(x[6],"[[","Paired_end")[c("R1_rightTrim"),]))
QwinTotalR2LT <- sapply(hts_json,function(x) unlist(sapply(x[6],"[[","Paired_end")[c("R2_leftTrim"),]))
QwinTotalR2RT <- sapply(hts_json,function(x) unlist(sapply(x[6],"[[","Paired_end")[c("R2_rightTrim"),]))
QwinTotalDiscard <- sapply(hts_json,function(x) unlist(sapply(x[6],"[[","Paired_end")[c("PE_discarded"),]))

## N trim
NwinTotalReadsIn <- sapply(hts_json,function(x) sapply(x[7], "[[","totalFragmentsInput"))
NwinTotalR1LT <- sapply(hts_json,function(x) unlist(sapply(x[7],"[[","Paired_end")[c("R1_leftTrim"),]))
NwinTotalR1RT <- sapply(hts_json,function(x) unlist(sapply(x[7],"[[","Paired_end")[c("R1_rightTrim"),]))
NwinTotalR2LT <- sapply(hts_json,function(x) unlist(sapply(x[7],"[[","Paired_end")[c("R2_leftTrim"),]))
NwinTotalR2RT <- sapply(hts_json,function(x) unlist(sapply(x[7],"[[","Paired_end")[c("R2_rightTrim"),]))
NwinTotalDiscard <- sapply(hts_json,function(x) unlist(sapply(x[7],"[[","Paired_end")[c("PE_discarded"),]))

## CutTrim
CutTrimTotalReadsIn <- sapply(hts_json,function(x) sapply(x[8], "[[","totalFragmentsInput"))
CutTrimTotalDiscard <- sapply(hts_json,function(x) unlist(sapply(x[8],"[[","Paired_end")[c("PE_discarded"),]))


## Post-Stats
stats2TotalReadsIn <- sapply(hts_json,function(x) sapply(x[9], "[[","totalFragmentsInput"))
stats2TotalReadsCG <- sapply(hts_json,function(x) sum(unlist(sapply(x[9],"[[","Base_composition")[c("C","G"),])))
stats2TotalReadsN <- sapply(hts_json,function(x) unlist(sapply(x[9],"[[","Base_composition")[c("N"),]))
stats2TotalReadsR1_BpLen <- sapply(hts_json,function(x) unlist(sapply(x[9],"[[","Paired_end")[c("R1_bpLen"),]))
stats2TotalReadsR1_BQ30 <- sapply(hts_json,function(x) unlist(sapply(x[9],"[[","Paired_end")[c("R1_bQ30"),]))
stats2TotalReadsR2_BpLen <- sapply(hts_json,function(x) unlist(sapply(x[9],"[[","Paired_end")[c("R2_bpLen"),]))
stats2TotalReadsR2_BQ30 <- sapply(hts_json,function(x) unlist(sapply(x[9],"[[","Paired_end")[c("R2_bQ30"),]))



LongTable <- data.frame(
    Raw_Reads=statsTotalReadsIn,
    Raw_Bp=(statsTotalReadsR1_BpLen+statsTotalReadsR2_BpLen),
    Raw_R1_PercentQ30=statsTotalReadsR1_BQ30/statsTotalReadsR1_BpLen*100,
    Raw_R2_PercentQ30=statsTotalReadsR2_BQ30/statsTotalReadsR2_BpLen*100,
    Raw_Percent_CG=statsTotalReadsCG/(statsTotalReadsR1_BpLen+statsTotalReadsR2_BpLen)*100,
    Raw_Chars_N = statsTotalReadsN,
    PhiX_IN = seqScrTotalReadsIn,
    PhiX_Discard = seqScrTotalReadsHits,
    PhiX_Percent_Discard = seqScrTotalReadsHits/seqScrTotalReadsIn*100,
    SuperDeduper_IN = SdTotalReadsIn,
    SuperDeduper_Ignored = SdTotalReadsIgnored,
    SuperDeduper_Duplicate = SdTotalReadsDuplicate,
    SuperDeduper_Percent_Duplicate = SdTotalReadsDuplicate/SdTotalReadsIn*100,
    rRNA_In = seqScr2TotalReadsIn,
    rRNA_Identified = seqScr2TotalReadsHits,
    rRNA_Percent_Identified = seqScr2TotalReadsHits/seqScr2TotalReadsIn*100,
    AdapterTrimmed_IN = AdaptTotalReadsIn,
    AdapterTrimmed_Reads = AdaptTotalAdaptTrim,
    AdapterTrimmed_Percent_Reads = AdaptTotalAdaptTrim/AdaptTotalReadsIn*100,
    AdapterTrimmed_BP = AdaptTotalBpTrim,
    QwindowTrimmed_IN = QwinTotalReadsIn,
    QwindowTrimmed_R1_LeftBpTrim = QwinTotalR1LT,
    QwindowTrimmed_R1_RightBpTrim = QwinTotalR1RT,
    QwindowTrimmed_R2_LeftBpTrim = QwinTotalR2LT,
    QwindowTrimmed_R2_RightBpTrim = QwinTotalR2RT,
    QwindowTrimmed_Discard = QwinTotalDiscard,
    NcharTrimmed_IN = NwinTotalReadsIn,
    NcharTrimmed_R1_LeftBpTrim = NwinTotalR1LT,
    NcharTrimmed_R1_RightBpTrim = NwinTotalR1RT,
    NcharTrimmed_R2_LeftBpTrim = NwinTotalR2LT,
    NcharTrimmed_R2_RightBpTrim = NwinTotalR2RT,
    NcharTrimmed_Discard = NwinTotalDiscard,
    MinLen_IN = CutTrimTotalReadsIn,
    MinLen_Discard = CutTrimTotalDiscard,
    Proc_Reads=stats2TotalReadsIn,
    Proc_Bp=(stats2TotalReadsR1_BpLen+stats2TotalReadsR2_BpLen),
    Proc_R1_PercentQ30=stats2TotalReadsR1_BQ30/stats2TotalReadsR1_BpLen*100,
    Proc_R2_PercentQ30=stats2TotalReadsR2_BQ30/stats2TotalReadsR2_BpLen*100,
    Proc_Percent_CG=stats2TotalReadsCG/(stats2TotalReadsR1_BpLen+stats2TotalReadsR2_BpLen)*100,
    Proc_Chars_N = stats2TotalReadsN,
    Final_Percent_Read=stats2TotalReadsIn/statsTotalReadsIn*100,
    Final_Percent_Bp=(stats2TotalReadsR1_BpLen+stats2TotalReadsR2_BpLen)/(statsTotalReadsR1_BpLen+statsTotalReadsR2_BpLen)*100
)

rownames(LongTable) <- samples
write.table(round(LongTable,2),"summary_hts.txt",sep="\t",row.names=T,col.names=T,quote=F)
