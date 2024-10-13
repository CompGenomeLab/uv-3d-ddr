comp <- c("00", "12")
hypo = "lessAbs"
alpha = 0.05
lfc = 0.5

library(tximport)
load("/home/carlos/oldies/projects/rna-seq/tx2gene.RDA")
t_comp1 <- paste0("t", comp[1])
t_comp2 <- paste0("t", comp[2])
samples <-data.frame (
  sample  = c(paste0(c("1","2","3"), comp[1]), paste0(c("1","2","3"), comp[2])),
  condition = c(rep(t_comp1, 3), rep(t_comp2, 3)),
  batch = c(rep(c(1,2,3), 2))
)
samples$condition <- relevel(as.factor(samples$condition), t_comp1)
samples$batch <- as.factor(samples$batch)
files <- file.path("quant", paste0("SU_" ,samples$sample), "quant.sf")
names(files) <- paste0(samples$sample)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

library(DESeq2)
dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~condition)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", t_comp2, t_comp1), alpha = alpha, lfcThreshold = lfc, altHypothesis = hypo)
res <- na.omit(res)
res$ensembl_gene_id <- gsub("\\.[0-9]*$", "", rownames(res))
degs <- res[res$padj <= alpha,]

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand"),
  filter="ensembl_gene_id",
  values=degs$ensembl_gene_id,
  uniqueRows=TRUE)

degs_df <- merge(as.data.frame(degs), annotLookup, by = "ensembl_gene_id")

saveTo <- paste0("/home/carlos/oldies/projects/rna-seq/",  "t0-t", comp[2], ".", hypo, ".tsv")
colnames(degs_df)[c(9:11)] <- c("chrom", "start", "end")
write.table(degs_df, file = saveTo, row.names = FALSE, sep = "\t", quote = FALSE)