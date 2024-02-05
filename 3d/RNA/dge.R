library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/home/carlos/Desktop/projects/rna-seq/GRCh38.gtf")
k <- keys(txdb, keytype = "GENEID")
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- tx2gene[c(2,1)]

comp <- c("00", "12")
t_comp1 <- paste0("t", comp[1])
t_comp2 <- paste0("t", comp[2])
degs_df_save_as <- paste0(comp[2], "_vs_", comp[1], ".tsv")

library(tximport)

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

# vsd <- vst(dds, blind = FALSE)
# rld <- rlog(dds, blind = FALSE)
# 
# library("dplyr")
# library("ggplot2")
# 
# dds <- estimateSizeFactors(dds)
# 
# df <- bind_rows(
#   as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
#     mutate(transformation = "log2(x + 1)"),
#   as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
#   as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
# 
# colnames(df)[1:2] <- c("x", "y")  
# 
# lvls <- c("log2(x + 1)", "vst", "rlog")
# df$transformation <- factor(df$transformation, levels=lvls)
# 
# ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
#   coord_fixed() + facet_grid( . ~ transformation)  
# 
# sampleDists <- dist(t(assay(vsd)))
# 
# library("pheatmap")
# library("RColorBrewer")
# 
# sampleDistMatrix <- as.matrix( sampleDists )
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows = sampleDists,
#          clustering_distance_cols = sampleDists,
#          col = colors)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", t_comp2, t_comp1), alpha = .05, lfcThreshold=0.5)
#res <- lfcShrink(dds = dds, res = res, coef=paste("condition", t_comp2, "vs", t_comp1, sep="_"), type="apeglm")
res <- na.omit(res)
res$ensembl_gene_id <- gsub("\\.[0-9]*$", "", rownames(res))
degs <- res[res$padj <= .05,]

#degs to Datafrane
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

saveTo <- paste0("/home/carlos/Desktop/projects/rna-seq/",  "t0-t12.degs.tsv")
colnames(degs_df)[c(9:11)] <- c("chrom", "start", "end")
write.table(degs_df, file = saveTo, row.names = FALSE, sep = "\t", quote = FALSE)

## Entrez
library("AnnotationDbi")
library("org.Hs.eg.db")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=res$ensembl_gene_id,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

entrez_df <- na.omit(res)


degs_entrez <- entrez_df[entrez_df$padj <= .05,]

geneList <- as.vector(degs_entrez$log2FoldChange)
names(geneList) <- degs_entrez$entrez

universe <- as.vector(entrez_df$log2FoldChange)
names(universe) <- entrez_df$entrez

#####
library("clusterProfiler")

#ego <- enrichGO(gene         = degs$ENSEMBL, universe = res$ENSEMBL, OrgDb         = org.Hs.eg.db, keyType       = 'ENSEMBL', ont           = "BP")
#cnetplot(ego, foldChange=geneList)

library(msigdbr)

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

em <- enricher(gene=degs_entrez$entrez, universe=entrez_df$entrez, TERM2GENE=m_t2g)
em_r <- setReadable(em, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
saveTo <- paste0("/home/carlos/Desktop/projects/rna-seq/time_course_Degs/",   paste0(comp[2], "_vs_", comp[1], ".shrunk.pdf"))
pdf(saveTo, 30, 30)
cnetplot(em_r , foldChange=universe)
dev.off()