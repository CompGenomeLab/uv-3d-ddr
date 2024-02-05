suppressMessages(library(LOLA))
suppressMessages(library(GenomicRanges))

options <- commandArgs(trailingOnly = TRUE)
firstinputbed = options[1]
secondinputbed = options[2]
outputfolder = options[3]

regionDB = readRDS("~/Desktop/projects/lola/hg38_robust_UniBind_LOLA.RDS")

firstregions = readBed(firstinputbed)
secondregions = readBed(secondinputbed)

thesets = GRangesList(firstregions, secondregions)
universe = disjoin(unlist(thesets))
locresults = runLOLA(firstregions, universe, regionDB, cores=8)

writeCombinedEnrichment(locresults, outFolder= outputfolder,
                        includeSplits=TRUE)

# make a new folder under the output folder, called "extracted_regions"
dir.create(file.path(outputfolder, "extracted_regions"))
options(scipen=999)
for (i in 1:nrow(locresults)) {
    pvaluelog = locresults[i,]$pValueLog
    pvalue = 10^-pvaluelog
    if (pvalue < 0.05) {
        gr = extractEnrichmentOverlaps(locresults[i,], firstregions, regionDB)
        df <- data.frame(
            seqnames=seqnames(gr),
            starts=start(gr)-1,
            ends=end(gr),
            names=c(rep(".", length(gr))),
            scores=c(rep(".", length(gr))),
            strands=strand(gr))
        write.table(df, file=file.path(outputfolder, "extracted_regions", paste0(i, ".bed")), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
}
