library(dplyr)
library(ggplot2)
library(ggrepel)


###########################################################################################################
####################################### Merging 2 datasets ################################################
###########################################################################################################

##### read RNA-seq data
rna_rpkm <- read.csv("/Users/zhaoyan/lzy/RA/Platelet/Data/RPKM_data.csv")
rna_de_res <- read.csv("/Users/zhaoyan/lzy/RA/Platelet/Results/RNAseq_edgeR_DE_result.csv")[,c(8,5,10:15)]
colnames(rna_de_res)[1:2] <- c('Gene','Ensembl_ID')
rna_de_res <- merge(rna_de_res, rna_rpkm[,c("Ensembl_ID", "RNA_A_avglog2_adj", "RNA_CB_avglog2_adj")],
                    by = 'Ensembl_ID')

##### read Ribo-seq data
ribo_rpkm <- read.csv("2/Data/RPKM_data.csv")
ribo_de_res <- read.csv("2/Results/Riboseq_edgeR_sva_DE_result.csv")[,-3]
ribo_de_res <- merge(ribo_de_res, ribo_rpkm[,c("Ensembl_ID", "RIBO_A_avglog2_adj", "RIBO_CB_avglog2_adj")],
                    by = 'Ensembl_ID')

##### merge two data
merge_df <- merge(rna_de_res, ribo_de_res, by = 'Ensembl_ID', suffixes = c(".RNA", ".RIBO"))
write.csv(merge_df, "2/Data/RNAseq_Riboseq_merged_data.csv",row.names = FALSE)



###########################################################################################################
#################################### Generating plots #####################################################
###########################################################################################################

#### read merged data
merge_df <- read.csv("2/Data/RNAseq_Riboseq_merged_data.csv")
ribo_rpkm <- read.csv("2/Data/RPKM_data.csv")
rna_rpkm <- read.csv("/Users/zhaoyan/lzy/RA/Platelet/Data/RPKM_data.csv")


#### Venn
rna_id <- rna_rpkm$Ensembl_ID
ribo_id <- ribo_rpkm$Ensembl_ID
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
venn.diagram(x=list(rna_id, ribo_id), 
             imagetype = 'png',
             filename = '2/Results/Integrated_gene_Venn.png', 
             category.names = c('RNA-seq', 'Ribo-seq'),
             fill=c('steelblue1','lightcoral'), 
             na='remove',
             cex = 2,
             cat.cex = 1.3)

venn.diagram(x=list(merge_df$Gene.RNA[which(merge_df$RNA_A_avglog2_adj > log2(0.3+1))],
                    merge_df$Gene.RNA[which(merge_df$RIBO_A_avglog2_adj > log2(0.3+1))]), 
             filename = '2/Results/Integrated_Adult_expressed_Venn.jpeg', 
             category.names = c('RNA-seq Adult', 'Ribo-seq Adult'), 
             fill= c('blue','yellow'),
             cex = 1.4,
             cat.cex = 1,
             units = 'in',
             height = 6,
             width = 6)
venn.diagram(x=list(merge_df$Gene.RNA[which(merge_df$RNA_CB_avglog2_adj > log2(0.3+1))],
                    merge_df$Gene.RNA[which(merge_df$RIBO_CB_avglog2_adj > log2(0.3+1))]), 
             filename = '2/Results/Integrated_CB_expressed_Venn.jpeg', 
             category.names = c('RNA-seq CB', 'Ribo-seq CB'), 
             fill= c('blue','yellow'),
             cex = 1.4,
             cat.cex = 1,
             units = 'in',
             height = 6,
             width = 6)


#### correlation
# for Adult group
cor.test(merge_df$RNA_A_avglog2_adj, merge_df$RIBO_A_avglog2_adj)

ggplot(merge_df, aes(x = RNA_A_avglog2_adj, y = RIBO_A_avglog2_adj, label=Gene.RNA)) +
  geom_point() + 
  theme_classic() +
  geom_text(x=9, y=16, label="R = 0.87, p < 2.2e-16")+
  geom_text(x=9, y=15.5, label="#genes = 7558")+
  xlab("RNA-seq_A_log2-average-RPKM")+
  ylab("Ribo-seq_A_log2-average-RPKM")+
  scale_x_continuous(breaks = seq(0, 16, by = 2))+
  scale_y_continuous(breaks = seq(0, 16, by = 2))+
  geom_text_repel(data          = subset(merge_df, (1.08*RNA_A_avglog2_adj-RIBO_A_avglog2_adj+0.81)/sqrt(1.08^2+1) > 3),
                  nudge_y       = -1 - subset(merge_df, (1.08*RNA_A_avglog2_adj-RIBO_A_avglog2_adj+0.81)/sqrt(1.08^2+1) > 3)$RIBO_A_avglog2_adj,
                  size          = 2.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x")+
  geom_text_repel(data          = subset(merge_df, (1.08*RNA_A_avglog2_adj-RIBO_A_avglog2_adj+0.81)/sqrt(1.08^2+1) < -3),
                  nudge_y       = 14 - subset(merge_df, (1.08*RNA_A_avglog2_adj-RIBO_A_avglog2_adj+0.81)/sqrt(1.08^2+1) < -3)$RIBO_A_avglog2_adj,
                  size          = 2.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "both")
ggsave('2/Results/Integrated_Adult_corrplot_labeled.jpeg', dpi = 700, height = 10, width = 10)

ggplot(merge_df, aes(x = RNA_A_avglog2_adj, y = RIBO_A_avglog2_adj)) +
  geom_point() + 
  theme_classic() +
  geom_text(x=9, y=16, label="R = 0.87, p < 2.2e-16")+
  geom_text(x=9, y=15.5, label="#genes = 7558")+
  xlab("RNA-seq_A_log2-average-RPKM")+
  ylab("Ribo-seq_A_log2-average-RPKM")+
  scale_x_continuous(breaks = seq(0, 16, by = 2))+
  scale_y_continuous(breaks = seq(0, 16, by = 2))
ggsave('2/Results/Integrated_Adult_corrplot_unlabeled.jpeg', dpi = 700, height = 8, width = 8)




# for CB
cor.test(merge_df$RNA_CB_avglog2_adj, merge_df$RIBO_CB_avglog2_adj)

ggplot(merge_df, aes(x = RNA_CB_avglog2_adj, y = RIBO_CB_avglog2_adj, label=Gene.RNA)) +
  geom_point() + 
  theme_classic() +
  geom_text(x=7, y=15.5, label="R = 0.85, p < 2.2e-16")+
  geom_text(x=7, y=15, label="#genes = 7558")+
  xlab("RNA-seq_CB_avg-log2-RPKM")+
  ylab("Ribo-seq_CB_avg-log2-RPKM")+
  scale_x_continuous(breaks = seq(0, 16, by = 2))+
  scale_y_continuous(breaks = seq(0, 16, by = 2))+
  geom_text_repel(data          = subset(merge_df, (1.08*RNA_CB_avglog2_adj-RIBO_CB_avglog2_adj+0.75)/sqrt(1.08^2+1) > 3),
                  nudge_y       = -1 - subset(merge_df, (1.08*RNA_CB_avglog2_adj-RIBO_CB_avglog2_adj+0.75)/sqrt(1.08^2+1) > 3)$RIBO_CB_avglog2_adj,
                  size          = 2.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "x")+
  geom_text_repel(data          = subset(merge_df, (1.08*RNA_CB_avglog2_adj-RIBO_CB_avglog2_adj+0.75)/sqrt(1.08^2+1) < -3),
                  nudge_y       = 14 - subset(merge_df, (1.08*RNA_CB_avglog2_adj-RIBO_CB_avglog2_adj+0.75)/sqrt(1.08^2+1) < -3)$RIBO_CB_avglog2_adj,
                  size          = 2.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "both")
ggsave('2/Results/Integrated_CB_corrplot_labeled.jpeg', dpi = 700, height = 10, width = 10)

ggplot(merge_df, aes(x = RNA_CB_avglog2_adj, y = RIBO_CB_avglog2_adj, label=Gene.RNA)) +
  geom_point() + 
  theme_classic() +
  geom_text(x=7, y=15.5, label="R = 0.85, p < 2.2e-16")+
  geom_text(x=7, y=15, label="#genes = 7558")+
  xlab("RNA-seq_CB_avg-log2-RPKM")+
  ylab("Ribo-seq_CB_avg-log2-RPKM")+
  scale_x_continuous(breaks = seq(0, 16, by = 2))+
  scale_y_continuous(breaks = seq(0, 16, by = 2))
ggsave('2/Results/Integrated_CB_corrplot_unlabeled.jpeg', dpi = 700, height = 8, width = 8)



#### log2FC histogram
# lfc_df <- data.frame(Log2FC = c(merge_df$RNA_LFC, merge_df$LFC), 
#                      tech = rep(c('RNA-seq', 'Ribo-seq'), c(nrow(merge_df), nrow(merge_df))))
# mean(merge_df$RNA_LFC) #1.06
# sd(merge_df$RNA_LFC) #1.14
# mean(merge_df$LFC) #0.73
# sd(merge_df$LFC) #2.46
# ggplot(lfc_df, aes(x=Log2FC, fill=tech, color=tech)) +
#   geom_histogram(alpha = 0.3, binwidth = 0.5) +
#   theme_classic()
# ggsave('Log2FC_hist_01.jpeg')



#### DEG comparison
# DE Venn diagram
rna_de_gene <- alltags2$Ensembl_ID[which(alltags2$RNA_FDR < 0.05 & abs(alltags2$RNA_logFC) > 2)]
rna_de_up <- alltags2$Ensembl_ID[which(alltags2$RNA_FDR < 0.05 & alltags2$RNA_logFC > 2)]
rna_de_down <- alltags2$Ensembl_ID[which(alltags2$RNA_FDR < 0.05 & alltags2$RNA_logFC < -2)]

ribo_de_gene <- ribo_rpkm_2$Ensembl_ID[which(ribo_rpkm_2$sva_RIBO_FDR < 0.05 & abs(ribo_rpkm_2$sva_RIBO_logFC) > 1.5)]
ribo_de_up <- ribo_rpkm_2$Ensembl_ID[which(ribo_rpkm_2$sva_RIBO_FDR < 0.05 & ribo_rpkm_2$sva_RIBO_logFC > 1.5)]
ribo_de_down <- ribo_rpkm_2$Ensembl_ID[which(ribo_rpkm_2$sva_RIBO_FDR < 0.05 & ribo_rpkm_2$sva_RIBO_logFC < -1.5)]

# venn.diagram(x=list(rna_de_gene, ribo_de_gene), filename = '2/Venn_all_DE_genes_ver2.png', 
#              category.names = c('RNA-seq', 'Ribo-seq'),
#              fill=c('steelblue1','lightcoral'),
#              cex = 1.5,
#              cat.cex = 2,
#              margin = 0.1,
#              cat.pos = c(-90,90),
#              cat.dist = c(0.1,0.1))
# venn.diagram(x=list(rna_de_up, ribo_de_up), filename = '2/Venn_up_DE_genes_ver2.png', 
#              category.names = c('RNA-seq', 'Ribo-seq'),
#              fill=c('steelblue1','lightcoral'),
#              cex = 1.5,
#              cat.cex = 2,
#              margin = 0.1,
#              cat.pos = c(-90,90),
#              cat.dist = c(0.1,0.1))
# venn.diagram(x=list(rna_de_down, ribo_de_down), filename = '2/Venn_down_DE_genes_ver2.png', 
#              category.names = c('RNA-seq', 'Ribo-seq'),
#              fill=c('steelblue1','lightcoral'),
#              cex = 1.5,
#              cat.cex = 2,
#              margin = 0.1,
#              cat.pos = c(-90,90),
#              cat.dist = c(0.1,0.1))
venn.diagram(x=list(rna_rpkm$Ensembl_ID,ribo_rpkm_2$Ensembl_ID, rna_de_gene, ribo_de_gene), 
             filename = '2/Venn_all.png', 
             category.names = c('RNA-seq', 'Ribo-seq', "RNA-seq_DE","Ribo-seq_DE"), 
             fill= c('steelblue1','lightcoral','#B3E2CD','#FDCDAC'))



# DE log2FC scatter plot
# merge_df$DE_type <- factor('not DE/DT',
#                            levels = c("DE only","DT only","homo-directional DE&DT","opposite-directional DE&DT","not DE/DT"))
# for(i in 1:nrow(merge_df)){
#   if(merge_df$DEG_type.RNA[i] == merge_df$DEG_type.RIBO[i] & merge_df$DEG_type.RNA[i] != 'Not-DEG'){
#     merge_df$DE_type[i] <- 'homo-directional DE&DT'
#   }
#   else if(merge_df$DEG_type.RNA[i] != 'Not-DEG' & merge_df$DEG_type.RIBO[i] != 'Not-DEG' & merge_df$DEG_type.RNA[i] != merge_df$DEG_type.RIBO[i]){
#     merge_df$DE_type[i] <- 'opposite-directional DE&DT'
#   }
#   else{
#     if(merge_df$DEG_type.RNA[i] != 'Not-DEG'){
#       merge_df$DE_type[i] <- 'DE only'
#     }
#     else if(merge_df$DEG_type.RIBO[i] != 'Not-DEG'){
#       merge_df$DE_type[i] <- 'DT only'
#     }
#   }
# }
merge_df$DE_type <- factor('other genes',
                           levels = c("up-regulated DE&DT", "down-regulated DE&DT", "DT only", "other genes"))
for(i in 1:nrow(merge_df)){
  if(merge_df$DEG_type.RNA[i] == merge_df$DEG_type.RIBO[i] & merge_df$DEG_type.RNA[i] == 'Up-DEG'){
    merge_df$DE_type[i] <- 'up-regulated DE&DT'
  }
  if(merge_df$DEG_type.RNA[i] == merge_df$DEG_type.RIBO[i] & merge_df$DEG_type.RNA[i] == 'Down-DEG'){
    merge_df$DE_type[i] <- 'down-regulated DE&DT'
  }
  if(merge_df$DEG_type.RNA[i] == 'Not-DEG' & merge_df$DEG_type.RIBO[i] != 'Not-DEG'){
    merge_df$DE_type[i] <- 'DT only'
  }
}
cols <- c("up-regulated DE&DT" = "red", 
          "down-regulated DE&DT" = "green", 
          "DT only" = "#56B4E9", 
          "other genes" = "#E1E0DF")
ggplot(mapping = aes(x=logFC.RNA, y=logFC.RIBO, label=Gene.RNA, color=DE_type)) + 
  geom_point(data=subset(merge_df, DE_type=='other genes'),size=2) +
  geom_point(data=subset(merge_df, DE_type=='DT only'), size=2) +
  geom_point(data=subset(merge_df, DE_type=='up-regulated DE&DT'), size=3) +
  geom_point(data=subset(merge_df, DE_type=='down-regulated DE&DT'), size=3) +
  # scale_colour_manual(values=c("#D55E00", "#56B4E9", "#009E73", "#F0E442", "#E1E0DF")) +
  # scale_colour_manual(values=c("#56B4E9", "#009E73","#D55E00", "#E1E0DF")) +
  # scale_colour_discrete(limits = c("up-regulated DE&DT", "down-regulated DE&DT", "DT only", "other genes"))+
  scale_colour_manual(values = cols)+
  theme_classic() +
  # geom_hline(yintercept=2.46, linetype="dashed", color = "black")+
  # geom_hline(yintercept=-2.46, linetype="dashed", color = "black")+
  # geom_vline(xintercept=1.14, linetype="dashed", color = "black")+
  # geom_vline(xintercept=-1.14, linetype="dashed", color = "black")
  geom_hline(yintercept=2, linetype="dashed", color = "black")+
  geom_hline(yintercept=-2, linetype="dashed", color = "black")+
  geom_vline(xintercept=2, linetype="dashed", color = "black")+
  geom_vline(xintercept=-2, linetype="dashed", color = "black")+
  # geom_abline(slope = 1, intercept = 0, color='black', linetype="dashed")+
  xlim(-11,11)+
  ylim(-11,11)+
  xlab("RNA-seq_logFC")+
  ylab("Ribo-seq_logFC")+
  geom_text(x=-2, y=-10, label="x=-2", color='black')+
  geom_text(x=2, y=-10, label="x=2", color='black')+
  geom_text(x=-10, y=2, label="y=2", color='black')+
  geom_text(x=-10, y=-2, label="y=-2", color='black')
  # geom_text(aes(label= ifelse(Gene.RNA %in% c('HBG1', 'HBG2','DEFA3','ZNF385D','CPT1A'),
  #                             as.character(Gene.RNA),'')),hjust=0,vjust=0, size=3)
  # geom_text_repel(data          = subset(merge_df, logFC.RIBO>0 & DE_type == 'up-regulated DE&DT'),
  #                 nudge_y       = 11 - subset(merge_df, logFC.RIBO>0 & DE_type == 'up-regulated DE&DT')$logFC.RIBO,
  #                 size          = 2.5,
  #                 segment.size  = 0.2,
  #                 segment.color = "grey50",
  #                 direction     = "both")+
  # geom_text_repel(data          = subset(merge_df, logFC.RIBO<0 & DE_type == 'down-regulated DE&DT'),
  #                 nudge_y       = -8 - subset(merge_df, logFC.RIBO<0 & DE_type == 'down-regulated DE&DT')$logFC.RIBO,
  #                 size          = 2.5,
  #                 segment.size  = 0.2,
  #                 segment.color = "grey50",
  #                 direction     = "x")
ggsave('2/Results/Integrated_logFC_unlabeled.jpeg', height = 4, width = 6, dpi = 700)



##### GSEA comparison
# es_df <- merge(rna_fgseaRes[,c('pathway','padj','NES')], ribo_fgseaRes[,c('pathway','padj','NES')],
#                suffixes = c('_RNA_seq', '_Ribo_seq'),
#                by = 'pathway')
# es_df[,2:5] <- apply(es_df[,2:5], MARGIN = 2, as.numeric)
# es_df$type <- factor('non enriched',
#                      levels = c("RNA-seq enriched only","Ribo-seq enriched only","both enriched","non enriched"))
# for(i in 1:nrow(es_df)){
#   if(es_df$padj_RNA_seq[i] < 0.05 & es_df$padj_Ribo_seq[i] >= 0.05){es_df$type[i] <- "RNA-seq enriched only"}
#   if(es_df$padj_RNA_seq[i] >= 0.05 & es_df$padj_Ribo_seq[i] < 0.05){es_df$type[i] <- "Ribo-seq enriched only"}
#   if(es_df$padj_RNA_seq[i] < 0.05 & es_df$padj_Ribo_seq[i] < 0.05){es_df$type[i] <- "both enriched"}
# }
# 
# ggplot(es_df, aes(x=NES_RNA_seq, y=NES_Ribo_seq, color=type))+
#   geom_point()+
#   theme_classic()+
#   geom_hline(yintercept=0, linetype="dashed", color = "black")+
#   geom_vline(xintercept=0, linetype="dashed", color = "black")+
#   xlab("RNA-seq_NES")+
#   ylab("Ribo-seq_NES")
# ggsave('2/GSEA_scatter_plot_ver2.jpeg', height = 6, width = 8)
# 
# 
# 
#### Translation efficiency
te_df <- merge(data.frame(Ensembl_ID = rna_rpkm$Ensembl_ID, Gene = rna_rpkm$Gene, RNA_A_mean = rowMeans(rna_rpkm[,1:9]), RNA_CB_mean = rowMeans(rna_rpkm[,10:20])),
               data.frame(Ensembl_ID = ribo_rpkm$Ensembl_ID, RIBO_A_mean = rowMeans(ribo_rpkm[,1:3]), RIBO_CB_mean = rowMeans(ribo_rpkm[,4:6]), RIBO_all_mean = rowMeans(ribo_rpkm[,1:6])),
               by = "Ensembl_ID")
te_df$A_TE <- te_df$RIBO_A_mean / te_df$RNA_A_mean
te_df$CB_TE <- te_df$RIBO_CB_mean / te_df$RNA_CB_mean

te_df$logFC <- log2((te_df$CB_TE+1e-9)/(te_df$A_TE+1e-9))
te_df$RIBO_log2_all_mean <- log2(te_df$RIBO_all_mean)

ggplot(data = te_df, aes(x=RIBO_log2_all_mean, y=logFC,label=Gene)) + 
  geom_point()+
  # scale_color_manual(values = c(""="grey", "express only in Ribo-seq"="red", 'express only in RNA-seq'="blue"))+
  theme_classic()+
  xlab("Ribo-seq log2-average-RPKM")+
  ylab("TE log2FC (CB:Adult)")+
  xlim(-5,16)+
  ylim(-40, 40)+
  # labs(colour = "gene type")
  geom_text_repel(data          = subset(te_df, logFC>33),
                  nudge_y       = 40 - subset(te_df, logFC>33)$logFC,
                  size          = 2.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "both")+
  geom_text_repel(data          = subset(te_df, logFC < -30),
                  nudge_y       = -40 - subset(te_df, logFC < -30)$logFC,
                  size          = 2.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "both")+
  geom_text_repel(aes(label = ifelse(RIBO_log2_all_mean>15, as.character(Gene),'')), hjust=0, vjust=0.5, cex=2.5)
ggsave("2/Results/Integrated_TE_MAplot_labeled.jpeg", dpi = 700, height = 10,width = 10)


# cor.test(te_df$A_TE, te_df$CB_TE)
# ggplot(te_df, aes(x=A_TE, y=CB_TE, label = Gene))+
#   geom_point()+
#   theme_classic() +
#   xlab("TE_N")+
#   ylab("TE_CB")+
#   # ylim(-1,8)+
#   geom_abline(slope = 1, intercept = 0, color='blue', linetype="dashed")+
#   # geom_text(x=7, y=7, label="y=x", color='blue')+
#   geom_text(x=4, y=7, label="R = 0.55, p < 2.2e-16", color='black')
#   # geom_text(aes(label=ifelse(abs(CB_TE_log-N_TE_log)/sqrt(2) > 3,
#   #                            as.character(Gene.x),'')),hjust=0, vjust=0)
#   geom_text_repel(data          = subset(merge_df, Gene.x %in% g4),
#                   nudge_y       = -1 - subset(merge_df, Gene.x %in% g4)$CB_TE_log,
#                   size          = 2.5,
#                   segment.size  = 0.2,
#                   segment.color = "grey50",
#                   direction     = "x")+
#   geom_text_repel(data          = subset(merge_df, Gene.x %in% g5),
#                   nudge_y       = 8 - subset(merge_df, Gene.x %in% g5)$CB_TE_log,
#                   size          = 2.5,
#                   segment.size  = 0.2,
#                   segment.color = "grey50",
#                   direction     = "x")+
#   geom_point(data=subset(merge_df, Gene.x == "DEFA3"), colour="red") +
#   geom_text(aes(label= ifelse(Gene.x == 'DEFA3',as.character(Gene.x),'')),hjust=0,vjust=0, size=3, colour = "red")
# ggsave('2/TE_scatter_plot_ver2.jpeg')
