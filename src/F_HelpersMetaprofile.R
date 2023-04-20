### INFO: Helper function to prepare genomic coordinates for metagene profiles
### DATE: 20.04.2023
### AUTHOR: Internet, small modifications

# Custom function to generate plotting table
get_metaprofile_tab <- function (txGTF = NULL, txGFF = NULL, txGenomeVer = NULL, txTxdb = NULL, 
                                 txGuitarTxdb = NULL, txGuitarTxdbSaveFile = NA, stBedFiles = NULL, 
                                 stGRangeLists = NULL, stGroupName = NULL, stAmblguity = 5, 
                                 stSampleNum = 10, stSampleModle = "Equidistance", txfiveutrMinLength = 100, 
                                 txcdsMinLength = 100, txthreeutrMinLength = 100, txlongNcrnaMinLength = 100, 
                                 txlncrnaOverlapmrna = FALSE, txpromoterLength = 1000, txtailLength = 1000, 
                                 txAmblguity = 5, txPrimaryOnly = FALSE, txTxComponentProp = NULL, 
                                 txMrnaComponentProp = NULL, txLncrnaComponentProp = NULL, 
                                 mapFilterTranscript = TRUE, headOrtail = TRUE, enableCI = TRUE, 
                                 pltTxType = c("tx", "mrna", "ncrna"), overlapIndex = 1, 
                                 siteLengthIndex = 1, adjust = 1, CI_ResamplingTime = 1000, 
                                 CI_interval = c(0.025, 0.975), miscOutFilePrefix = NA) 
{
  genomeVersion2Txdb <- list(hg18 = "TxDb.Hsapiens.UCSC.hg18.knownGene", 
                             hg19 = "TxDb.Hsapiens.UCSC.hg19.knownGene", hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene", 
                             mm9 = "TxDb.Mmusculus.UCSC.mm9.knownGene", mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene")
  if (headOrtail) {
    txpromoterLength <- txpromoterLength
    txtailLength <- txtailLength
  }
  else {
    txpromoterLength <- 0
    txtailLength <- 0
  }
  print(format(Sys.time(), "%Y%m%d%H%M%S"))
  guitarTxdb <- Guitar:::.getGuitarTxdb(txGTF = txGTF, txGFF = txGFF, 
                                        txGenomeVer = txGenomeVer, txTxdb = txTxdb, txGuitarTxdb = txGuitarTxdb, 
                                        txfiveutrMinLength = txfiveutrMinLength, txcdsMinLength = txcdsMinLength, 
                                        txthreeutrMinLength = txthreeutrMinLength, txlongNcrnaMinLength = txlongNcrnaMinLength, 
                                        txlncrnaOverlapmrna = txlncrnaOverlapmrna, txpromoterLength = txpromoterLength, 
                                        txtailLength = txtailLength, txAmblguity = txAmblguity, 
                                        txTxComponentProp = txTxComponentProp, txMrnaComponentProp = txMrnaComponentProp, 
                                        txLncrnaComponentProp = txLncrnaComponentProp, txPrimaryOnly = txPrimaryOnly, 
                                        pltTxType = pltTxType, genomeVersion2Txdb)
  print(format(Sys.time(), "%Y%m%d%H%M%S"))
  if (!(is.na(txGuitarTxdbSaveFile))) {
    txGuitarTxdbSaveFile <- paste("GuitarTxdb", txGuitarTxdbSaveFile, 
                                  format(Sys.time(), "%Y%m%d"), sep = "-")
    save(guitarTxdb, file = txGuitarTxdbSaveFile)
  }
  sitesGroup <- Guitar:::.getStGroup(stBedFiles = stBedFiles, stGRangeLists = stGRangeLists, 
                                     stGroupName = stGroupName)
  GroupNames <- names(sitesGroup)
  sitesGroupNum <- length(sitesGroup)
  sitesPointsNormlize <- list()
  sitesPointsRelative <- list()
  pointWeight <- list()
  for (i in seq_len(sitesGroupNum)) {
    GroupName = GroupNames[[i]]
    print(paste("sample", stSampleNum, "points for", GroupName, 
                sep = " "))
    sitesPoints <- samplePoints(sitesGroup[i], stSampleNum = stSampleNum, 
                                stAmblguity = stAmblguity, pltTxType = pltTxType, 
                                stSampleModle = stSampleModle, mapFilterTranscript = mapFilterTranscript, 
                                guitarTxdb)
    for (txType in pltTxType) {
      sitesPointsNormlize[[txType]][[GroupName]] <- normalize(sitesPoints, 
                                                              guitarTxdb, txType, overlapIndex, siteLengthIndex)
      sitesPointsRelative[[txType]][[GroupName]] <- sitesPointsNormlize[[txType]][[GroupName]][[1]]
      pointWeight[[txType]][[GroupName]] <- sitesPointsNormlize[[txType]][[GroupName]][[2]]
    }
  }
  for (txType in pltTxType) {
    if (!(txType %in% guitarTxdb$txTypes)) {
      print(paste("Warning: Cannot plot distribution for", 
                  txType))
      next
    }
    print(paste("start figure plotting for", txType, "..."))
    if (txType == "mrna") {
      txType_name = "mRNA"
    }
    else if (txType == "ncrna") {
      txType_name = "ncRNA"
    }
    else {
      txType_name = "Transcript"
    }
    title <- paste("Distribution on", txType_name)
    #
    densityDataframe_CI <- Guitar:::.generateDensity_CI(sitesPointsRelative[[txType]], 
                                                        pointWeight[[txType]], 
                                                        CI_ResamplingTime, 
                                                        adjust = adjust, 
                                                        enableCI = enableCI)
    #
    comp_width <- guitarTxdb[[txType]]$componentWidthAverage_pct
    #
    if (enableCI) {
      peak <- max(densityDataframe_CI$confidenceUp)
    } 
    else {
      peak <- max(densityDataframe_CI$density)
    }
    #
    out <- list(densityDataframe_CI, 
                comp_width, 
                peak) %>% 
      magrittr::set_names(c("data_tab", "component_width", "peak_height"))
    return(out)
  }
}