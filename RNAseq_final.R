library(edgeR)
library(limma)
library(dplyr)
library(fgsea)
library(corrplot)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(gplots)
library(enrichplot)
library(clusterProfiler)
library(ggnewscale)
library(gplots)
library(readxl)
library(ggExtra)
source("Program/my_netgraph_functions.R")


###########################################################################################################
################################# Processing RNA-seq data #################################################
###########################################################################################################

##### read raw count and gene annotation
count <- read.csv("Data/myfullTable.csv")
gene.anno <- read.delim2("Data/gene_annotation_info.txt")

##### remove version number in Ensembl_Gene_ID
count$Ensembl_Gene_ID <- gsub("\\.[0-9]*", "", count$Ensembl_Gene_ID)
gene.anno$Ensembl_Gene_ID <- gsub("\\.[0-9]*", "", gene.anno$Ensembl_Gene_ID)

##### creat a DGEList object
group <- factor(c(rep('Adult',9),rep('CB',11)))
dglist <- DGEList(counts = count[,3:22], genes = gene.anno, group = group)

##### filtering
myrpkm <- rpkm(count[,3:22], gene.length = count$Length)
cri<- 0.3
k1 <- gene.anno$Chromosome!="chrM" & gene.anno$Gene_type %in% c("protein_coding","pseudogene")
k2 <- rowSums(myrpkm[,1:20] > cri) > 2
dglist <- dglist[k1&k2,,keep.lib.sizes=FALSE]

##### normalization
dglist <- calcNormFactors(dglist,method = "TMM")
# tmm <- cpm(dglist)

##### DE test
design = model.matrix(~group)
dglist = estimateDisp(dglist, design)

fit = glmFit(dglist,design) # likelihood ratio test
lrt = glmLRT(fit,coef = 2)
topTags(lrt)
summary(decideTests(lrt))
DE_res = topTags(lrt,nrow(lrt)) %>% as.data.frame()
DE_res$DEG_type <- ifelse(DE_res$FDR >= 0.05, 
                          "Not-DEG", 
                          ifelse(DE_res$logFC > 2,
                                 "Up-DEG",
                                 ifelse(DE_res$logFC < -2,
                                        "Down-DEG",
                                        "Not-DEG")))
write.csv(DE_res,"Results/RNAseq_edgeR_DE_result.csv", row.names = FALSE)
saveRDS(dglist, "Data/RNAseq_DGEList.RDS")


#### GSEA (MSigDB Hallmark)
rank_metric <- DE_res$logFC
rna_ranks <- sort(rank_metric, decreasing = TRUE)
names(rna_ranks) <- DE_res$Gene_name[order(rank_metric, decreasing = TRUE)]
# hallmark <- gmtPathways("~/lzy/RA/Single Cell/hspc/h.all.v7.0.symbols.gmt")
# rna_fgseaRes <- fgsea(hallmark, rna_ranks, nperm=100000)
# rna_fgseaRes_df <- apply(rna_fgseaRes,2,as.character)
# write.csv(rna_fgseaRes, "Results/RNAseq_GSEA_Hallmark_result.csv")
hallmark <- read.gmt("~/lzy/RA/Single Cell/hspc/h.all.v7.0.symbols.gmt")
hallmark$ont <- gsub("HALLMARK_","", hallmark$ont)
gseaResult <- GSEA(geneList = rna_ranks,
                   nPerm = 100000,
                   minGSSize = 1,
                   maxGSSize = Inf,
                   TERM2GENE = hallmark,
                   by = 'fgsea')
saveRDS(gseaResult, "Data/RNAseq_fgsea_result.RDS")


#### GSEA (KEGG)
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

ids<-bitr(gene.anno$Ensembl_Gene_ID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

df2 = DE_res[DE_res$Ensembl_Gene_ID %in% dedup_ids$ENSEMBL,]
df2$Y = dedup_ids$ENTREZID[match(df2$Ensembl_Gene_ID, dedup_ids$ENSEMBL)]

kegg_gene_list <- df2$logFC
names(kegg_gene_list) <- df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               nPerm        = 100000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none")
kk2_df <- kk2@result
write.csv(kk2_df, "Results/GSEA_KEGG_result.csv")





##### save rpkm data ready for plots
# myrpkm_filtered <- rpkm(dglist, gene.length = count$Length[k1&k2]) %>% as.data.frame()
# myrpkm_filtered$RNA_A_avglog2 <- rowMeans(log2(myrpkm_filtered[,1:9]+1))
# myrpkm_filtered$RNA_CB_avglog2 <- rowMeans(log2(myrpkm_filtered[,10:20]+1))
myrpkm_filtered <- myrpkm[k1&k2,] %>% as.data.frame()
myrpkm_filtered$RNA_A_avglog2 <- log2(rowMeans(myrpkm_filtered[,1:9]))
myrpkm_filtered$RNA_CB_avglog2 <- log2(rowMeans(myrpkm_filtered[,10:20]))
myrpkm_filtered$RNA_A_avglog2_rank <- rank(-myrpkm_filtered$RNA_A_avglog2)
myrpkm_filtered$RNA_CB_avglog2_rank <- rank(-myrpkm_filtered$RNA_CB_avglog2)
myrpkm_filtered$Gene <- dglist$genes$Gene_name
myrpkm_filtered$Ensembl_ID <- dglist$genes$Ensembl_Gene_ID
myrpkm_filtered$RNA_A_avglog2_adj <- log2(rowMeans(myrpkm_filtered[,1:9])+1)
myrpkm_filtered$RNA_CB_avglog2_adj <- log2(rowMeans(myrpkm_filtered[,10:20])+1)
write.xlsx(myrpkm_filtered, "Data/RPKM_data.xlsx", row.names = FALSE, col.names = TRUE)

##### get MT genes data
mt.keep <- gene.anno$Chromosome =="chrM"
mt <- myrpkm[mt.keep,] %>% as.data.frame()
mt$Gene_name <- gene.anno$Gene_name[mt.keep]
write.csv(mt,"Data/MTgenes_rpkm.csv", row.names = FALSE)

##### get ribosomal genes
ribosomal.keep <- grepl("^RP[SL]", gene.anno$Gene_name, ignore.case = FALSE)
ribosomal <- myrpkm[ribosomal.keep&k1&k2,] %>% as.data.frame()
ribosomal$Gene_name <- gene.anno$Gene_name[ribosomal.keep&k1&k2]
write.csv(ribosomal,"Data/Ribosomalgenes_rpkm.csv", row.names = FALSE)

###########################################################################################################
#################################### Generating plots #####################################################
###########################################################################################################

# read intermediate data
dglist <- readRDS("Data/RNAseq_DGEList.RDS")
myrpkm_filtered <- read_xlsx("Data/RPKM_data.xlsx", 1)
DE_res <- read.csv("Results/RNAseq_edgeR_DE_result.csv")
gseaResult <- readRDS("Data/RNAseq_fgsea_result.RDS")
mt <- read.csv("Data/MTgenes_rpkm.csv")
ribosomal <- read.csv("Data/Ribosomalgenes_rpkm.csv")
hallmark <- read.gmt("~/lzy/RA/Single Cell/hspc/h.all.v7.0.symbols.gmt")
hallmark$ont <- gsub("HALLMARK_","", hallmark$ont)

##### density curve
jpeg(filename = "Results/RNAseq_density_curve.jpeg", res=500, width = 8, height = 6, units = "in")
xl = expression(log[2]("RPKM"))
plot(density(log2(myrpkm_filtered[,1])),col="black",
     # main = "Density of RPKM distribution",
     main = "",
     xlab = xl,
     ylim = c(0,0.2),
     xlim = c(-20,20))
for(i in 2:9){
  lines(density(log2(myrpkm_filtered[,i])),col="black")
}
for(j in 10:20){
  lines(density(log2(myrpkm_filtered[,j])),col="red")
}
legend("topright",legend=c("Adult","CB"),lty=c(1,1),col=c("black","red"),bty="n")
dev.off()


##### abundance rank plot
ggplot()+
  geom_point(data = subset(myrpkm_filtered, RNA_A_avglog2 > log2(0.3)), 
             aes(x = RNA_A_avglog2_rank, y = RNA_A_avglog2, color = 'a'))+
  geom_point(data = subset(myrpkm_filtered, RNA_A_avglog2 <= log2(0.3)), 
             aes(x = RNA_A_avglog2_rank, y = RNA_A_avglog2), color='grey70')+
  geom_point(data = subset(myrpkm_filtered, RNA_CB_avglog2 > log2(0.3)), 
             aes(x = RNA_CB_avglog2_rank, y = RNA_CB_avglog2, color = 'b'))+
  geom_point(data = subset(myrpkm_filtered, RNA_CB_avglog2 <= log2(0.3)), 
             aes(x = RNA_CB_avglog2_rank, y = RNA_CB_avglog2), color='pink')+
  theme_classic()+
  geom_hline(yintercept=log2(0.3), linetype="dashed", color = "blue")+
  # geom_text(aes(x=100, y=0.6, label="RPKM=0.3"), color='blue')+
  scale_color_manual(
    values = c("black", 'red'),
    labels = c("Adult", 'CB')
  )+
  labs(x="Rank of gene", y="log2-avergae-RPKM", colour='group')
ggsave("Results/RNAseq_rank_unlabeled.jpeg", dpi = 700, height = 6,width = 8)




##### Supplementary Table3
# Adult-restricted genes
N_restricted_df <- myrpkm_filtered[which(myrpkm_filtered$RNA_A_avglog2 > log2(0.3) & 
                                         myrpkm_filtered$RNA_CB_avglog2 <= log2(0.3)),]
write.xlsx(N_restricted_df, "Data/Adult_restricted_genes_data.xlsx", row.names = FALSE, col.names = TRUE)

# N_restricted_df_melt <- data.frame(type = rep(c("Adult", "CB"),each=88),
#                                    value = c(N_restricted_df$RNA_A_avglog2, N_restricted_df$RNA_CB_avglog2))
# ggplot(N_restricted_df_melt, aes(x=value, fill=type)) +
#   geom_histogram(binwidth=0.5, alpha=.6, position="identity")+
#   theme_classic()+
#   xlab("CB RNAseq (log2 RPKM)")+
#   ylab("frequency")+
#   scale_x_continuous(breaks = seq(-6,4,1), limits = c(-6.5,4.5))
# ggsave("Results/Adult_restricted_genes_dist.jpeg", dpi = 700, height = 6,width = 8)

N_restricted_df <- merge(N_restricted_df, DE_res, by.x='Gene', by.y='Gene_name')
table(N_restricted_df$DEG_type)

# ggplot(N_restricted_df, aes(x=logFC)) +
#   geom_histogram(binwidth=0.5, fill="grey", alpha=0.4, color='black')+
#   theme_classic()+
#   xlab("log2 fold-change (CB:Adult)")+
#   ylab("frequency")+
#   scale_x_continuous(breaks = seq(-9,-1,1), limits = c(-9,-1))+
#   scale_y_continuous(breaks = seq(0,20,5), limits = c(0,25), expand = c(0, 0))
# ggsave("Results/Adult_restricted_genes_logfc_hist.jpeg", dpi = 700, height = 6,width = 8)

# myrpkm_filtered$N_restricted <- ifelse(myrpkm_filtered$Gene %in% N_restricted_df$Gene, 
#                                        "Adult-restricted genes (N=88)", 
#                                        "other genes")
# ggplot(myrpkm_filtered, aes(x = RNA_A_avglog2_adj, y = RNA_CB_avglog2_adj)) +
#   geom_point(aes(color = factor(N_restricted))) + 
#   scale_color_manual(values = c("red","black")) +
#   scale_size_manual(values = c(1,3))+
#   theme_classic() +
#   geom_text(x=9, y=14, label="R = 0.92, p < 2.2e-16")+
#   geom_text(x=9, y=13.5, label="#genes = 11741")+
#   xlab("Adult RNASeq (log2 RPKM)")+
#   ylab("CB RNASeq (log2 RPKM)")+
#   labs(color='gene type') +
#   scale_x_continuous(breaks = seq(0, 16, by = 2))+
#   scale_y_continuous(breaks = seq(0, 16, by = 2))+
#   theme(
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank())
# ggsave('Results/RNAseq_corrplot_with_Adult_restricted.jpeg', dpi = 700, height = 8, width = 10)

p1 <- ggplot(N_restricted_df, aes(x = RNA_CB_avglog2_adj, y = RNA_A_avglog2_adj, color=DEG_type, label=Gene))+
  geom_point()+
  scale_color_manual(values = c("Not-DEG"="grey", "Down-DEG"="blue"))+
  theme_classic()+
  theme(legend.position="bottom")+
  xlab("CB RNASeq (adjusted log2 RPKM)")+
  ylab("Adult RNASeq (adjusted log2 RPKM)")+
  # scale_x_continuous(breaks = seq(-2, 4, by = 1), limits = c(-2, 4.5))+
  # scale_y_continuous(breaks = seq(-6, -2, by = 1), limits = c(-6.5, -1.5))+
  scale_x_continuous(breaks = seq(0, 0.4, by = 0.1), limits = c(0, 0.4), expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0, 5), expand = c(0, 0))+
  geom_text_repel(data = subset(N_restricted_df, Gene %in% c("ZNF385D", "L3MBTL4")),
                  nudge_y = 5 - subset(N_restricted_df, Gene %in% c("ZNF385D", "L3MBTL4"))$RNA_A_avglog2_adj,
                  hjust=0,
                  vjust=0,
                  size=3)
p1 <- ggMarginal(p1, type = "box")
ggsave(plot=p1,'Results/Adult_restricted_genes_corrplot.jpeg', dpi = 700, height = 5, width = 5)


# CB-restricted genes
CB_restricted_df <- myrpkm_filtered[which(myrpkm_filtered$RNA_A_avglog2_adj <= log2(1.3) & 
                                           myrpkm_filtered$RNA_CB_avglog2_adj > log2(1.3)),]
CB_restricted_df <- merge(CB_restricted_df, DE_res, by.x='Gene', by.y='Gene_name')

p2 <- ggplot(CB_restricted_df, aes(x = RNA_A_avglog2_adj, y = RNA_CB_avglog2_adj, color=DEG_type, label=Gene))+
  geom_point()+
  scale_color_manual(values = c("Not-DEG"="grey", "Up-DEG"="red"))+
  theme_classic()+
  theme(legend.position="bottom")+
  xlab("Adult RNASeq (adjusted log2 RPKM)")+
  ylab("CB RNASeq (adjusted log2 RPKM)")+
  # scale_x_continuous(breaks = seq(-12, -1, by = 1), limits = c(-12, -1), expand = c(0, 0))+
  # scale_y_continuous(breaks = seq(-2, 5, by = 1), limits = c(-2, 5), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, 0.4, by = 0.1), limits = c(0, 0.4), expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0, 5), expand = c(0, 0))+
  geom_text_repel(data = subset(CB_restricted_df, Gene %in% c('PAICS',"IGF2BP1", "COL4A5")),
                  nudge_y = 5- subset(CB_restricted_df, Gene %in% c('PAICS',"IGF2BP1", "COL4A5"))$RNA_CB_avglog2_adj,
                  hjust=0.2,
                  vjust=1,
                  size=3)
p2 <- ggMarginal(p2, type = "box")
ggsave(plot=p2,'Results/CB_restricted_genes_corrplot.jpeg', dpi = 700, height = 5, width = 5)


ggplot(data = CB_restricted_df, aes(x=logCPM, y=logFC, color=DEG_type, label=Gene)) + 
  geom_point()+
  scale_color_manual(values = c("Not-DEG"="grey", "Up-DEG"="red"))+
  theme_classic()+
  xlab("Average logCPM")+
  ylab("log2 fold-change")+
  scale_x_continuous(breaks = seq(-10, 20, by = 5))+
  scale_y_continuous(breaks = seq(-10, 15, by = 5))+
  labs(colour = "gene type")





##### Venn diagram
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
venn.diagram(x=list(myrpkm_filtered$Gene[which(myrpkm_filtered$RNA_A_avglog2 > log2(0.3))],
                    myrpkm_filtered$Gene[which(myrpkm_filtered$RNA_CB_avglog2 > log2(0.3))]), 
             filename = 'Results/RNAseq_rank_Venn.jpeg', 
             category.names = c('Adult', 'CB'), 
             fill= c('black','red'),
             cex = 1.5,
             cat.cex = 2,
             units = 'in',
             height = 6,
             width = 6)


##### MDS plot
jpeg(filename = "Results/RNAseq_MDS_plot.jpeg", res=300, width = 6, height = 6, units = "in")
plotMDS(dglist, col=c(rep("black",9), rep("red",11)), pch=19) # plotMDS uses log2(cpm) to calculate distance
legend("right", c("Adult", "CB"),pch = 19, col = c("black", "red"), cex = 0.8)
dev.off()


##### dendrogram
d <- myrpkm_filtered[,1:20] %>% as.matrix() %>% t() %>% scale() %>% dist(., method = "euclidean")
hc <- hclust(d, method = "ward.D2")
dhc <- as.dendrogram(hc)
cols <- rep(c('black','red'),c(9,11))
cols <- cols[order.dendrogram(dhc)]

jpeg(filename = "Results/RNAseq_dendrogram.jpeg", res=300, width = 8, height = 6, units = "in")
dhc %>% set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 1.5) %>%  # node point size
  set("leaves_col", cols) %>% #node point color
  plot(main = "Dendrogram of RNA-seq RPKM", leaflab = "none")
legend("topright", c("Adult", "CB"),pch = 19, col = c("black", "red"), cex = 0.8)
dev.off()


##### MA plot
ggplot(data = DE_res, aes(x=logCPM, y=logFC, color=DEG_type, label=Gene_name)) + 
  geom_point()+
  scale_color_manual(values = c("Not-DEG"="grey", "Up-DEG"="red", 'Down-DEG'="blue"))+
  theme_classic()+
  xlab("Average logCPM")+
  ylab("log2 fold-change")+
  scale_x_continuous(breaks = seq(-10, 20, by = 5))+
  scale_y_continuous(breaks = seq(-10, 15, by = 5))+
  labs(colour = "gene type")
ggsave("Results/RNAseq_MAplot_unlabeled.jpeg", dpi = 700, height = 6,width = 8)

ggplot(data = DE_res, aes(x=logCPM, y=logFC, color=DEG_type, label=Gene_name)) + 
  geom_point()+
  scale_color_manual(values = c("Not-DEG"="grey", "Up-DEG"="red", 'Down-DEG'="blue"))+
  theme_classic()+
  xlab("Average logCPM")+
  ylab("log2 fold-change")+
  scale_x_continuous(breaks = seq(-10, 20, by = 5))+
  scale_y_continuous(breaks = seq(-10, 15, by = 5))+
  geom_text_repel(data = subset(DE_res, logFC > 6),
                  nudge_y = 12 - subset(DE_res, logFC > 6)$logFC,
                  hjust=0,
                  vjust=0,
                  size=2.5,
                  direction = 'both')+
  geom_text_repel(data = subset(DE_res, logFC < -5 | logCPM > 13),
                  nudge_y = -10 - subset(DE_res, logFC < -5 | logCPM > 13)$logFC,
                  hjust=0,
                  vjust=0,
                  size=2.5,
                  direction = 'both')+
  labs(colour = "gene type")
ggsave("Results/RNAseq_MAplot_labeled.jpeg", dpi = 700, height = 6,width = 8)


ma_df <- data.frame(Ensembl_Gene_ID = myrpkm_filtered$Ensembl_ID, 
                    logRPKM = log2(rowMeans(myrpkm_filtered[,1:20]))) %>% 
  merge(., DE_res[,c('Ensembl_Gene_ID', 'logFC','DEG_type','Gene_name')])
ggplot(data = ma_df, aes(x=logRPKM, y=logFC, color=DEG_type, label=Gene_name)) + 
  geom_point()+
  scale_color_manual(values = c("Not-DEG"="grey", "Up-DEG"="red", 'Down-DEG'="blue"))+
  theme_classic()+
  xlab("log2-average-RPKM")+
  ylab("log2-fold-change")+
  xlim(-5,15)+
  ylim(-10,15)+
  labs(colour = "gene type")+
  geom_text_repel(data = subset(ma_df, logFC > 6),
                  nudge_y = 12 - subset(ma_df, logFC > 6)$logFC,
                  hjust=0,
                  vjust=0,
                  size=2.5,
                  direction = 'both')+
  geom_text_repel(data = subset(ma_df, logFC < -5 | logRPKM > 13),
                  nudge_y = -10 - subset(ma_df, logFC < -5 | logRPKM > 13)$logFC,
                  hjust=0,
                  vjust=0,
                  size=2.5,
                  direction = 'both')
ggsave("Results/RNAseq_MAplot_RPKM_labeled.jpeg", dpi = 700, height = 6,width = 8)


##### volcano plot
DE_res$log10padj <- -log10(DE_res$FDR)
ggplot(data = DE_res, aes(x=logFC, y=log10padj, color=DEG_type)) + 
  geom_point()+
  scale_color_manual(values = c("Not-DEG"="grey", "Up-DEG"="red", 'Down-DEG'="blue"))+
  theme_classic()+
  xlab("log2 fold-change")+
  ylab("-log10(padj)")+
  scale_x_continuous(breaks = seq(-12, 14, by = 2))+
  scale_y_continuous(breaks = seq(0, 130, by = 20))+
  geom_text_repel(aes(label= ifelse(logFC < -4 | logFC > 6 | log10padj > 50 | (logFC > 5 & log10padj > 40) | (logFC < -3 & log10padj > 30),
                                    as.character(Gene_name),'')),hjust=0.5,vjust=0,size=3)+
  labs(colour = "gene type")
ggsave("Results/RNAseq_volcano_labeled.jpeg", dpi = 700, height = 8,width = 10)

ggplot(data = DE_res, aes(x=logFC, y=log10padj, color=DEG_type)) + 
  geom_point()+
  scale_color_manual(values = c("Not-DEG"="grey", "Up-DEG"="red", 'Down-DEG'="blue"))+
  theme_classic()+
  xlab("log2 fold-change")+
  ylab("-log10(padj)")+
  scale_x_continuous(breaks = seq(-12, 14, by = 2))+
  scale_y_continuous(breaks = seq(0, 130, by = 20))+
  labs(colour = "gene type")
ggsave("Results/RNAseq_volcano_unlabeled.jpeg", dpi = 700, height = 6,width = 8)


#### correlation scatter plot
# within group
# rna_n_cor <- cor(myrpkm_filtered[,grep("^N", colnames(myrpkm_filtered))])
# jpeg(file = 'Results/Corr_matrix_RNA_rpkm_Adult.jpeg', units = "in", width = 5, height = 5, res = 500)
# corrplot(rna_n_cor,method = "color", type = "lower",
#          addCoef.col = "white", tl.col = "black", tl.srt = 0,
#          diag = FALSE, number.cex = 0.8, number.digits = 3)
# dev.off()
# 
# rna_cb_cor <- cor(myrpkm_filtered[,grep("^CB", colnames(myrpkm_filtered))])
# jpeg(file = 'Results/Corr_matrix_RNA_rpkm_CB.jpeg', units = "in", width = 5, height = 5, res = 500)
# corrplot(rna_cb_cor,method = "color", type = "lower",
#          addCoef.col = "white", tl.col = "black", tl.srt = 0,
#          diag = FALSE, number.cex = 0.8, number.digits = 3)
# dev.off()


# rna-seq: n vs cb (Figure4)
cor.test(myrpkm_filtered$RNA_A_avglog2_adj, myrpkm_filtered$RNA_CB_avglog2_adj)

ggplot(myrpkm_filtered, aes(x = RNA_A_avglog2_adj, y = RNA_CB_avglog2_adj, label=Gene)) +
  geom_point() + 
  theme_classic() +
  geom_text(x=9, y=14, label="R = 0.92, p < 2.2e-16")+
  geom_text(x=9, y=13.5, label="#genes = 11741")+
  xlab("RNA-seq_A_log2-average-RPKM")+
  ylab("RNA-seq_CB_log2-average-RPKM")+
  scale_x_continuous(breaks = seq(0, 16, by = 2))+
  scale_y_continuous(breaks = seq(0, 16, by = 2))+
  geom_text_repel(data          = subset(myrpkm_filtered, (0.90*RNA_A_avglog2_adj-RNA_CB_avglog2_adj+0.82)/sqrt(0.90^2+1) > 2),
                  nudge_y       = -1 - subset(myrpkm_filtered, (0.90*RNA_A_avglog2_adj-RNA_CB_avglog2_adj+0.82)/sqrt(0.90^2+1) > 2)$RNA_CB_avglog2_adj,
                  size          = 2.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x")+
  geom_text_repel(data          = subset(myrpkm_filtered, (0.90*RNA_A_avglog2_adj-RNA_CB_avglog2_adj+0.82)/sqrt(0.90^2+1) < -2),
                  nudge_y       = 13 - subset(myrpkm_filtered, (0.90*RNA_A_avglog2_adj-RNA_CB_avglog2_adj+0.82)/sqrt(0.90^2+1) < -2)$RNA_CB_avglog2_adj,
                  size          = 2.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "both")
  # geom_text_repel(aes(label = ifelse(Gene %in% c('HBG1', 'HBG2','DEFA3','ZNF385D','CPT1A'), 
  #                                    as.character(Gene), '')))
ggsave('Results/RNAseq_corrplot_labeled.jpeg', dpi = 700, height = 8, width = 8)


ggplot(myrpkm_filtered, aes(x = RNA_A_avglog2_adj, y = RNA_CB_avglog2_adj, label=Gene)) +
  geom_point() + 
  theme_classic() +
  geom_text(x=9, y=14, label="R = 0.92, p < 2.2e-16")+
  geom_text(x=9, y=13.5, label="#genes = 11741")+
  xlab("RNA-seq_A_avg-log2-RPKM")+
  ylab("RNA-seq_CB_avg-log2-RPKM")+
  scale_x_continuous(breaks = seq(0, 16, by = 2))+
  scale_y_continuous(breaks = seq(0, 16, by = 2))
ggsave('Results/RNAseq_corrplot_unlabeled.jpeg', dpi = 700, height = 8, width = 8)




#### DEG Heatmap
rna_de_up <- DE_res$Gene_name[DE_res$DEG_type == 'Up-DEG']
rna_de_down <- DE_res$Gene_name[DE_res$DEG_type == 'Down-DEG']

rna_heatmap_df <- rbind(myrpkm_filtered[which(myrpkm_filtered$Gene %in% rna_de_up),c(1:20,25)], 
                        myrpkm_filtered[which(myrpkm_filtered$Gene %in% rna_de_down),c(1:20,25)])
row.names(rna_heatmap_df) <- rna_heatmap_df$Gene
jpeg(filename = "Results/RNAseq_DEG_heatmap.jpeg", units = "in", width = 6, height = 7, res = 300)
rna_heatmap_df[,1:20] %>% as.matrix() %>% t() %>% scale() %>% t() %>%
  heatmap.2(Rowv = FALSE, Colv = TRUE, labRow = FALSE, col = "bluered",
            scale = "none", dendrogram='column', cexRow = 0.8,
            trace = "none", density.info = "none", key = FALSE,
            RowSideColors = c(rep("indianred1", length(rna_de_up)), rep("skyblue1", length(rna_de_down))),
            hclustfun = function(x) hclust(x,method = 'ward.D'))
dev.off()



##### GSEA network graph
my_emapplot(gseaResult, min_edge=0.01, cex_label_category = 0.8)
ggsave('Results/RNAseq_GSEA_hallmark_emapplot.jpeg', dpi = 700, width = 8, height = 6)

rank_metric <- DE_res$logFC
rna_ranks <- sort(rank_metric, decreasing = TRUE)
names(rna_ranks) <- DE_res$Gene_name[order(rank_metric, decreasing = TRUE)]
DE_res$log10padj <- -log10(DE_res$FDR)
volcano_genes <- subset(DE_res, logFC < -4 | logFC > 6 | log10padj > 50 | (logFC > 5 & log10padj > 40) | (logFC < -3 & log10padj > 30))
my_cnetplot(gseaResult, cex_label_category = 1, cex_gene = 0.8,cex_label_gene = 0.8,
            genelist = c(as.character(volcano_genes$Gene_name),"TOM1L1","VAMP3","AP1G1","DNM1L","CLCN3", "SNAP23"),
            foldChange = rna_ranks,
            layout = 'nicely')
ggsave('Results/RNAseq_GSEA_hallmark_cnetplot.jpeg', dpi = 700, width = 11, height = 11)


# # over-represented analysis for up-regulted DEGs against GO
# up_deg <- DE_res$Ensembl_Gene_ID[which(DE_res$DEG_type=='Up-DEG')] %>% as.vector()
# ego1 <- enrichGO(gene = up_deg, 
#                  ont ="BP", 
#                  keyType = "ENSEMBL", 
#                  minGSSize = 3, 
#                  maxGSSize = 800, 
#                  pvalueCutoff = 0.05, 
#                  OrgDb = organism, 
#                  pAdjustMethod = "BH",
#                  readable = TRUE)
# emapplot(ego1, showCategory = 50)
# ggsave('Results/RNAseq_upDEG_GO_network.jpeg', dpi = 700, height = 10, width = 10)
# # cnetplot(ego1, node_label='category', categorySize="pvalue", showCategory = 10,foldChange=genelist)
# 
# # over-represented analysis for down-regulted DEGs against GO
# down_deg <- DE_res$Ensembl_Gene_ID[which(DE_res$DEG_type=='Down-DEG')] %>% as.vector()
# ego2 <- enrichGO(gene = down_deg, 
#                  ont ="all", 
#                  keyType = "ENSEMBL", 
#                  minGSSize = 1, 
#                  maxGSSize = 300, 
#                  pvalueCutoff = 0.05, 
#                  OrgDb = organism, 
#                  pAdjustMethod = "BH",
#                  readable = TRUE)
# emapplot(ego2, showCategory = 50)




##### Heatmap of "protein secretion" pathway genes
ps_all_genes <- hallmark$gene[which(hallmark$ont == "PROTEIN_SECRETION")] # All genes in the pathway
ps_in_genes <- intersect(ps_all_genes, myrpkm_filtered$Gene) # pathway genes that in our data
rna_fgseaRes_df <- gseaResult@result
ps_core_genes <- rna_fgseaRes_df["PROTEIN_SECRETION",]$core_enrichment %>% # core enrichment genes
  strsplit(.,"/") %>% 
  unlist()

ps_df <- cpm(dglist) %>% as.data.frame()
ps_df <- ps_df[match(ps_core_genes,myrpkm_filtered$Gene),]
row.names(ps_df) <- ps_core_genes

lmat <- rbind(c(0, 4,5), c(0,1,0), c(3,2,0))
lhei <- c(1.5, 0.2, 4)
lwid <- c(0.5,4,1)

rowcol <- c("black", "red")[(row.names(ps_df) %in% c("CLCN3","SNAP23","DNM1L","AP1G1","VAMP3","TOM1L1")) +1]

jpeg(filename = "Results/RNAseq_Protein_Secretion_core_genes_heatmap.jpeg", 
     units = "in", width = 6, height = 5, res = 500)
par(mar=c(0.1,0.1,0.1,0.1))
ps_df %>% as.matrix() %>% t() %>% scale() %>% t() %>%
  heatmap.2(Rowv = FALSE, Colv = TRUE,
            scale = "none", dendrogram='column', cexCol = 1,
            trace = "none", density.info = "none", cexRow = 0.6,
            ColSideColors = c(rep("black", 9), rep("red", 11)),
            hclustfun = function(x) hclust(x,method = 'ward.D'),
            col="bluered", labCol=FALSE, key.ylab=NA, key.xlab="scaled value",
            lmat = lmat, lwid = lwid, lhei = lhei, keysize = 0.1,
            rowsep = c(0,40:44,46:50), 
            colsep = c(0,11,20), sepwidth = c(0.03,0.03), sepcolor = 'black',
            colRow = c(rowcol), margins = c(1,3))
# coords <- locator(1)  # coords is a list recording the place to put legend (coords$x = 0.8733997, coords$y = 0.7199757)
coords <- list(x=0.8733997,y = 0.7199757)
legend(coords,      
       legend = c("Adult", "CB"),
       col = c("black","red"), 
       lty= 1,             
       lwd = 5,           
       cex=.6
)
dev.off()






##### Heatmap of MT genes
row.names(mt) = mt$Gene_name
mt <- mt[,1:20]

jpeg(filename = "Results/Riboseq_MTgenes_heatmap.jpeg", units = "in", width = 6, height = 6, res = 300)
mt %>% as.matrix() %>% t() %>% scale() %>% t() %>%
  heatmap.2(Rowv = FALSE, Colv = TRUE,
            scale = "none", dendrogram='column', cexCol = 1,
            trace = "none", density.info = "none", cexRow = 0.5,
            # ColSideColors = c(rep("red", 9), rep("black", 11)),
            hclustfun = function(x) hclust(x,method = 'ward.D'))
dev.off()


##### ranked scatter plot of MT genes
mt_ttest_df <- data.frame(ratio  = rowMeans(mt[,1:9])/rowMeans(mt[,10:20]),
                          pvalue = apply(mt, 1, function(x) t.test(x[1:9],x[10:20], alternative = 'greater')$p.value))
mt_ttest_df$adjp <- p.adjust(mt_ttest_df$pvalue)
mt_ttest_df$rank <- rank(-mt_ttest_df$ratio)
mt_ttest_df$type <- ifelse(mt_ttest_df$adjp<0.05, "adjp < 0.05", "adjp >= 0.05")
ggplot(data = mt_ttest_df, mapping = aes(x=rank, y=ratio, color=type))+
  geom_point()+
  scale_color_manual(values=c("adjp < 0.05"="red", "adjp >= 0.05"='black'), name = "t.test result")+
  theme_classic() +
  xlab("Ordinal of mitochondrial genes")+
  ylab("Ratio of average RPKM (Adult:CB)")+
  scale_x_continuous(breaks = seq(1, 37, by = 1))+
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), limits=c(0,3))+
  geom_hline(yintercept=1, linetype="dashed", color = "blue")+
  theme(axis.text.x = element_text(size=5))
  
ggsave("Results/Riboseq_MTgenes_rankplot.jpeg", dpi = 700, width = 6,height = 4)



##### Heatmap of Ribosomal genes
row.names(ribosomal) <- ribosomal$Gene_name

jpeg(filename = "Results/Riboseq_RIBOSOMALgenes_heatmap.jpeg", units = "in", width = 6, height = 6, res = 300)
ribosomal[,1:20] %>% as.matrix() %>% t() %>% scale() %>% t() %>%
  heatmap.2(Rowv = FALSE, Colv = TRUE, labRow = FALSE,
            scale = "none", dendrogram='column', cexCol = 1,
            trace = "none", density.info = "none", cexRow = 0.5,
            # ColSideColors = c(rep("red", 9), rep("black", 11)),
            hclustfun = function(x) hclust(x,method = 'ward.D'))
dev.off()
