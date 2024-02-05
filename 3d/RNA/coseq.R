outPath = "/home/carlos/Desktop/manuscripts/notebooks/RNA/" #with / at end


library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

OrgDb <- org.Hs.eg.db

library(tximport)
library(GenomicFeatures)

# txdb <- makeTxDbFromGFF("/home/carlos/Desktop/projects/rna-seq/GRCh38.gtf",)
# k <- keys(txdb, keytype = "GENEID")
# tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
# tx2gene <- tx2gene[c(2,1)]
# save(tx2gene, file='/home/carlos/Desktop/projects/rna-seq/tx2gene.RDA')
load("/home/carlos/Desktop/projects/rna-seq/tx2gene.RDA")

sampleNames=c("SU_100", "SU_200", "SU_300", "SU_112", "SU_212", "SU_312", "SU_130", "SU_230", "SU_330", "SU_160", "SU_260", "SU_360")
files <- file.path("/home/carlos/Desktop/projects/rna-seq/quant", sampleNames, "quant.sf")
names(files) <- paste0(sampleNames)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
conds <- c(rep("Control",3),
           rep("12_mins",3),
           rep("30_mins",3),
           rep("60_mins",3)
           )

conds <- factor(conds, levels = c('Control', 
                                  '12_mins', 
                                  '30_mins', 
                                  '60_mins'))

library(DESeq2)
dds <- DESeqDataSetFromTximport(txi.salmon, DataFrame(time= conds),design = ~ time)
dds <- dds[ rowSums(counts(dds)) > 3, ]
dds2 <- DESeq(dds, test=c("LRT"), reduced =~1)
res <- results(dds2, independentFiltering=FALSE, alpha = 0.05)

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotate_genes <- function(mart, df, ensLookup){
  annotLookup <- getBM(
    mart=mart,
    attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name", "entrezgene_id", "description", "chromosome_name", "start_position","end_position", "strand"),
    values=ensLookup, uniqueRows = T, filters = 'ensembl_gene_id')
  
  colnames(annotLookup) <- c('gene_id', colnames(annotLookup)[2:9])
  
  df <- as.data.frame(df)
  df$gene_id = ensLookup
  df = df[match(df$gene_id, annotLookup$gene_id),]
  total <- merge(df,annotLookup,by="gene_id")
  return(total)
}


res_df <- as.data.frame(subset(res, padj <= 0.05))
all_genes <- sub("\\.\\d+", "", rownames(res_df))
ensLookup <- all_genes
total <- annotate_genes(mart, res_df, ensLookup = ensLookup)
write.table(total, file = paste0(outPath, "all_deseq_lrt.tsv"), row.names=FALSE, quote = F, col.names=TRUE, sep="\t")


library(coseq)
coseqCounts = counts(dds2, norm=TRUE)[which(res$padj < 0.05),]
all_genes <- sub("\\.\\d+", "", rownames(coseqCounts))
# counts_df <- annotate_genes(mart, coseqCounts, all_genes)
# write.table(counts_df, file = paste0(outPath, "counts.all_deseq_lrt.tsv"), row.names=FALSE, quote = F, col.names=TRUE, sep="\t")

### annot
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens")
m_t2g_h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, ensembl_gene)
###
seed = 43
set.seed(seed)
runArcsin <- coseq(coseqCounts,
                   normFactors = "none",
                   K=5:15,
                   model="Normal",
                   transformation="arcsin",
                   parallel = TRUE)

filterTrsh=0.90
res = runArcsin
nameofRes = paste0("runArcsin_seed",seed)

cls = clusters(res)
genesAssay = assay(res)
#universe = sub("\\.\\d+", "", rownames(genesAssay))
universe = sub("\\.\\d+", "", rownames(coseqCounts))

pdf(paste0(outPath,"clusters_", nameofRes, ".pdf"),
    width=14, height=14)
print(plot(res, 
           graphs="boxplots", 
           conds=conds, 
           collapse_reps = "average"))

all_clusters = list()
for (i in sort(unique(cls))){
  cluster = i
  genes = rownames(genesAssay)[which(genesAssay[,i] >= filterTrsh)]
  genes = sub("\\.\\d+", "", genes)
  
  cname <- paste0('c',i)
  all_clusters[[cname]] <- genes
  
  print(i)
  outPathNow = paste0(outPath, cluster)
  df=data.frame(genes=genes, cluster=rep(i ,length(genes)))
  
  ensLookup <- df$genes
  annotLookup <- getBM(
    mart=mart,
    attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name", "entrezgene_id", "description", "chromosome_name", "start_position","end_position", "strand"),
    values=ensLookup, uniqueRows = T, filters = 'ensembl_gene_id')
  
  annotLookup <- data.frame(
    df$genes[match(annotLookup$ensembl_gene_id, genes)],
    annotLookup)
  
  colnames(annotLookup) <- c(
    "original_id",
    c("gene_id", "gene_biotype", "external_gene_name", "entrezgene_id", "description", "chromosome_name", "start_position","end_position", "strand"))
  
  df$gene_id = ensLookup
  total <- merge(df,annotLookup,by="gene_id")
  
  write.table(total, file = paste0(outPathNow,"_", nameofRes, ".geneIDs.", filterTrsh,".tsv"), row.names=FALSE, quote = F, col.names=TRUE, sep="\t")
  
  
  ego = enrichGO(gene          = genes,
                 OrgDb         = OrgDb,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 universe=universe)
  
  
  if (nrow(as.data.frame(ego)) > 0){
    ego = setReadable(ego, OrgDb, keyType = "ENSEMBL")
    print(paste0(cluster, ' - ego'))
    write.table(as.data.frame(ego), file = paste0(outPathNow,"_", nameofRes, ".ego.tsv"),row.names=FALSE, quote = F, col.names=TRUE, sep="\t")
    
    p1 = dotplot(ego, showCategory=20, orderBy="GeneRatio", font.size = 8, title = paste0("EGO - Cluster: ", cluster))
    p2 = enrichplot::cnetplot(ego, showCategory=20, colorEdge=TRUE, cex_label_category = 0.6, node_label="gene" ) + theme(text=element_text(size=6))
    print(p1)
    print(p2)
  }
  
  msig_h = enricher(genes, 
                  TERM2GENE = m_t2g_h,
                  universe = universe
                  )
  if (nrow(as.data.frame(msig_h)) > 0){
    msig_h = setReadable(msig_h, OrgDb, keyType = "ENSEMBL")
    print(paste0(cluster, ' - h - msigdb'))
    write.table(as.data.frame(msig_h), file = paste0(outPathNow,"_", nameofRes, ".msig_h.tsv"),row.names=FALSE, quote = F, col.names=TRUE, sep="\t")
    
    p1 = dotplot(msig_h, showCategory=20, orderBy="GeneRatio", font.size = 8, title = paste0("MSigDB - H - Cluster: ", cluster))
    p2 = enrichplot::cnetplot(msig_h, showCategory=20, colorEdge=TRUE, cex_label_category = 0.6, node_label="gene" ) + theme(text=element_text(size=6))
    print(p1)
    print(p2)
  }
}

ck <- compareCluster(geneCluster = all_clusters, 
                     fun = enrichGO, 
                     OrgDb         = OrgDb,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     universe = universe)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENSEMBL")
p = dotplot(ck, showCategory=20, font.size = 8, title = "GO - BP - All Cluster")
print(p)

ck <- compareCluster(geneCluster = all_clusters, fun = enricher, TERM2GENE = m_t2g_h,
                     universe = universe)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENSEMBL")
p = dotplot(ck, showCategory=20, font.size = 8, title = "MSigDB - H - All Cluster")
print(p)

write.table(table(cls), file = paste0(outPath, nameofRes, ".cls.tsv"),row.names=FALSE, quote = F, col.names=TRUE, sep="\t")
tcounts_df <- as.data.frame(tcounts(res))
tcounts_df$gene_id = rownames(tcounts_df)

write.table(tcounts_df, file = paste0(outPath, nameofRes, ".tcounts.tsv"),row.names=FALSE, quote = F, col.names=TRUE, sep="\t")

dev.off()
