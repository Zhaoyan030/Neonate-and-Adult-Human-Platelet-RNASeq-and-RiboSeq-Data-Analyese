library(edgeR)
library(limma)
library(dplyr)
library(ggplot2)
library(fgsea)
library(readxl)
library(ggrepel)
library(sva)
library(corrplot)
library(gplots)
library(VennDiagram)


###########################################################################################################
################################### Process Ribo-seq data ###############################################
###########################################################################################################

##### read Ribo-seq data
ribo <- read_xlsx('2/raw/All_gene_expressions_table.xlsx')
colnames(ribo) <- gsub('-', '_', colnames(ribo))
ribo_rpkm <- ribo[,c(17:22,1,3)]
ribo_count <- ribo[,c(11:16,1,3,5)]

##### filtering
k1 <- rowSums(ribo_count[,1:6] >= 5) > 0
k2 <- !is.na(ribo_count$Ensembl_ID)
ribo_count_filtered <- ribo_count[k1&k2,] %>% subset(., !duplicated(Ensembl_ID))
ribo_rpkm_filtered <- ribo_rpkm[k1&k2,] %>% subset(., !duplicated(Ensembl_ID))

##### creat a DGEList object
group <- factor(c('A','A','A','CB','CB','CB'), levels = c('A','CB'))
y <- DGEList(ribo_count_filtered[,1:6], group=group, genes = ribo_count_filtered[,7:9])

##### normalization
y <- calcNormFactors(y, method = 'TMM')

##### DE test
y_logcpm <- cpm(y, log = TRUE)
mod <- model.matrix(~ group, data=group)
mod0 <- model.matrix(~ 1, data=group)
svobj <- svaseq(y_logcpm, mod, mod0)
design <- cbind(mod, svobj$sv)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit, coef = 2)
ribo_de_res <- topTags(lrt,nrow(lrt)) %>% as.data.frame()
ribo_de_res$DEG_type <- ifelse(ribo_de_res$FDR >= 0.05, 
                          "Not-DEG", 
                          ifelse(ribo_de_res$logFC > 1.5,
                                 "Up-DEG",
                                 ifelse(ribo_de_res$logFC < -1.5,
                                        "Down-DEG",
                                        "Not-DEG")))
write.csv(ribo_de_res, "2/Results/Riboseq_edgeR_sva_DE_result.csv", row.names = FALSE)
saveRDS(y, "2/Data/Riboseq_DGEList.RDS")


#### GSEA
rank_metric <- ribo_de_res$logFC
ribo_ranks <- sort(rank_metric, decreasing = FALSE)
names(ribo_ranks) <- ribo_de_res$Gene[order(rank_metric, decreasing = FALSE)]

hallmark <- gmtPathways("~/lzy/RA/Single Cell/hspc/h.all.v7.0.symbols.gmt")
ribo_fgseaRes <- fgsea(hallmark, ribo_ranks, nperm=100000)
ribo_fgseaRes <- apply(ribo_fgseaRes,2,as.character)
write.csv(ribo_fgseaRes, "2/Results/RIBOseq_GSEA_Hallmark_result.csv")


##### save rpkm data ready for plots
ribo_rpkm_filtered$RIBO_A_avglog2 <- log2(rowMeans(ribo_rpkm_filtered[,1:3]))
ribo_rpkm_filtered$RIBO_CB_avglog2 <- log2(rowMeans(ribo_rpkm_filtered[,4:6]))
ribo_rpkm_filtered$RIBO_A_avglog2_adj <- log2(rowMeans(ribo_rpkm_filtered[,1:3])+1)
ribo_rpkm_filtered$RIBO_CB_avglog2_adj <- log2(rowMeans(ribo_rpkm_filtered[,4:6])+1)
write.csv(ribo_rpkm_filtered, "2/Data/RPKM_data.csv", row.names = FALSE)



###########################################################################################################
#################################### Generating plots #####################################################
###########################################################################################################

##### read intermediate data
ribo_rpkm_filtered <- read.csv("2/Data/RPKM_data.csv")
ribo_de_res <- read.csv("2/Results/Riboseq_edgeR_sva_DE_result.csv")
y <- readRDS("2/Data/Riboseq_DGEList.RDS")


##### MDS plot
jpeg(filename = "2/Results/Riboseq_MDS_plot.jpeg", res=300, width = 6, height = 6, units = "in")
plotMDS(y, 
        col=c(rep("black",3), rep("red",3)), 
        # labels= c(paste("Adult",c('A','B','C'),sep = "_"), paste("Cord",c('A','B','C'),sep = "_")),
        xlim = c(-2,3), 
        ylim = c(-2,2),
        pch=19)
legend("right", c("Adult", "CB"),pch = 19, col = c("black", "red"), cex = 0.8)
dev.off()


##### Venn diagram
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
venn.diagram(x=list(ribo_rpkm_filtered$Gene[which(ribo_rpkm_filtered$RIBO_A_avglog2 > log2(0.3+1))],
                    ribo_rpkm_filtered$Gene[which(ribo_rpkm_filtered$RIBO_CB_avglog2 > log2(0.3+1))]), 
             filename = '2/Results/Riboseq_expressed_Venn.jpeg', 
             category.names = c('Adult', 'CB'), 
             fill= c('black','red'),
             cex = 1.4,
             cat.cex = 2,
             units = 'in',
             height = 6,
             width = 6)



#### correlation plots
# within group
# ribo_n_cor <- cor(ribo_rpkm_2[,grep("^Adult", colnames(ribo_rpkm_2))])
# row.names(ribo_n_cor) <- paste("Adult",c('A','B','C'), sep ='_')
# colnames(ribo_n_cor) <- paste("Adult",c('A','B','C'), sep ='_')
# jpeg(file = '2/Corr_matrix_RIBO_rpkm_N_ver2.jpeg', units = "in", width = 5, height = 5, res = 500)
# corrplot(ribo_n_cor,method = "color", type = "lower",
#          addCoef.col = "white", tl.col = "black", tl.srt = 0, tl.cex = 1.5,
#          diag = FALSE, number.cex = 1.5, number.digits = 3)
# dev.off()
# 
# ribo_cb_cor <- cor(ribo_rpkm_2[,grep("^Cord", colnames(ribo_rpkm_2))])
# row.names(ribo_cb_cor) <- paste("Cord",c('A','B','C'), sep ='_')
# colnames(ribo_cb_cor) <- paste("Cord",c('A','B','C'), sep ='_')
# jpeg(file = '2/Corr_matrix_RIBO_rpkm_CB_ver2.jpeg', units = "in", width = 5, height = 5, res = 500)
# corrplot(ribo_cb_cor,method = "color", type = "lower",
#          addCoef.col = "white", tl.col = "black", tl.srt = 0, tl.cex = 1.5,
#          diag = FALSE, number.cex = 1.5, number.digits = 3)
# dev.off()

# ribo-seq: n vs cb (Figure10)
cor.test(ribo_rpkm_filtered$RIBO_A_avglog2, ribo_rpkm_filtered$RIBO_CB_avglog2)

ggplot(ribo_rpkm_filtered, aes(x = RIBO_A_avglog2, y = RIBO_CB_avglog2, label=Gene)) +
  geom_point() + 
  theme_classic() +
  geom_text(x=9, y=16, label="R = 0.90, p < 2.2e-16")+
  geom_text(x=9, y=15.5, label="#genes = 8230")+
  xlab("Ribo-seq_A_log2-average-RPKM")+
  ylab("Ribo-seq_CB_log2-average-RPKM")+
  scale_x_continuous(breaks = seq(0, 16, by = 2))+
  scale_y_continuous(breaks = seq(0, 16, by = 2))+
  geom_text_repel(data          = subset(ribo_rpkm_filtered, (0.87*RIBO_A_avglog2-RIBO_CB_avglog2+1.15)/sqrt(0.87^2+1) > 2 & RIBO_A_avglog2 > 2),
                  nudge_x       = 11 - subset(ribo_rpkm_filtered, (0.87*RIBO_A_avglog2-RIBO_CB_avglog2+1.15)/sqrt(0.87^2+1) > 2 & RIBO_A_avglog2 > 2)$RIBO_CB_avglog2,
                  size          = 2.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "both")+
  geom_text_repel(data          = subset(ribo_rpkm_filtered, (0.87*RIBO_A_avglog2-RIBO_CB_avglog2+1.15)/sqrt(0.87^2+1) < -2),
                  nudge_y       = 13 - subset(ribo_rpkm_filtered, (0.87*RIBO_A_avglog2-RIBO_CB_avglog2+1.15)/sqrt(0.87^2+1) < -2)$RIBO_CB_avglog2,
                  size          = 2.5,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "both")
ggsave('2/Results/Riboseq_corrplot_labeled.jpeg', dpi = 700, height = 10, width = 10)


ggplot(ribo_rpkm_filtered, aes(x = RIBO_A_avglog2, y = RIBO_CB_avglog2, label=Gene)) +
  geom_point() + 
  theme_classic() +
  geom_text(x=9, y=16, label="R = 0.90, p < 2.2e-16")+
  geom_text(x=9, y=15.5, label="#genes = 8230")+
  xlab("Ribo-seq_A_log2-average-RPKM")+
  ylab("Ribo-seq_CB_log2-average-RPKM")+
  scale_x_continuous(breaks = seq(0, 16, by = 2))+
  scale_y_continuous(breaks = seq(0, 16, by = 2))
ggsave('2/Results/Riboseq_corrplot_unlabeled.jpeg', dpi = 700, height = 8, width = 8)



###### volcano plot
ribo_de_res$log10padj <- -log10(ribo_de_res$FDR)
ggplot(data = ribo_de_res, aes(x=logFC, y=log10padj, color=DEG_type)) + 
  geom_point()+
  scale_color_manual(values = c("Not-DEG"="grey", "Up-DEG"="red", 'Down-DEG'="blue"))+
  theme_classic()+
  xlab("log2 fold-change")+
  ylab("-log10(padj)")+
  scale_x_continuous(breaks = seq(-12, 14, by = 2))+
  scale_y_continuous(breaks = seq(0, 130, by = 20))+
  geom_text_repel(aes(label= ifelse(DEG_type != 'Not-DEG' & (logFC < -5 | logFC > 7 | log10padj > 8 ),
                                    as.character(Gene),'')),hjust=0.5,vjust=0,size=3)+
  labs(colour = "gene type")
ggsave("2/Results/Riboseq_volcano_labeled.jpeg", dpi = 700, height = 8,width = 10)

ggplot(data = ribo_de_res, aes(x=logFC, y=log10padj, color=DEG_type)) + 
  geom_point()+
  scale_color_manual(values = c("Not-DEG"="grey", "Up-DEG"="red", 'Down-DEG'="blue"))+
  theme_classic()+
  xlab("log2 fold-change")+
  ylab("-log10(padj)")+
  scale_x_continuous(breaks = seq(-12, 14, by = 2))+
  scale_y_continuous(breaks = seq(0, 130, by = 20))+
  labs(colour = "gene type")
ggsave("2/Results/Riboseq_volcano_unlabeled.jpeg", dpi = 700, height = 6,width = 8)


##### MA plot
ggplot(data = ribo_de_res, aes(x=logCPM, y=logFC, color=DEG_type, label=Gene)) + 
  geom_point()+
  scale_color_manual(values = c("Not-DEG"="grey", "Up-DEG"="red", 'Down-DEG'="blue"))+
  theme_classic()+
  xlab("Average logCPM")+
  ylab("log2 fold-change")+
  scale_x_continuous(breaks = seq(-10, 20, by = 5))+
  scale_y_continuous(breaks = seq(-10, 15, by = 5))+
  labs(colour = "gene type")
ggsave("2/Results/Riboseq_MAplot_unlabeled.jpeg", dpi = 700, height = 8,width = 10)

ggplot(data = ribo_de_res, aes(x=logCPM, y=logFC, color=DEG_type, label=Gene)) + 
  geom_point()+
  scale_color_manual(values = c("Not-DEG"="grey", "Up-DEG"="red", 'Down-DEG'="blue"))+
  theme_classic()+
  xlab("Average logCPM")+
  ylab("log2 fold-change")+
  scale_x_continuous(breaks = seq(-10, 20, by = 5))+
  scale_y_continuous(breaks = seq(-10, 15, by = 5))+
  geom_text_repel(data = subset(ribo_de_res, logFC > 7),
                  nudge_y = 12 - subset(ribo_de_res, logFC > 7)$logFC,
                  hjust=0,
                  vjust=0,
                  size=2.5,
                  direction = 'both')+
  geom_text_repel(data = subset(ribo_de_res, logFC < -5 | logCPM > 15),
                  nudge_y = -10 - subset(ribo_de_res, logFC < -5 | logCPM > 15)$logFC,
                  hjust=0,
                  vjust=0,
                  size=2.5,
                  direction = 'both')+
  labs(colour = "gene type")
ggsave("2/Results/Riboseq_MAplot_labeled.jpeg", dpi = 700, height = 8,width = 10)


#### DE result: Heatmap
ribo_de_up <- ribo_de_res$Gene[which(ribo_de_res$DEG_type == "Up-DEG")]
ribo_de_down <- ribo_de_res$Gene[which(ribo_de_res$DEG_type == "Down-DEG")]

ribo_heatmap_df <- rbind(ribo_rpkm_filtered[which(ribo_rpkm_filtered$Gene %in% ribo_de_up),c(1:7)], 
                         ribo_rpkm_filtered[which(ribo_rpkm_filtered$Gene %in% ribo_de_down),c(1:7)])
row.names(ribo_heatmap_df) <- ribo_heatmap_df$Gene
jpeg(filename = "2/Results/Riboseq_DEG_heatmap.jpeg", units = "in", width = 6, height = 6, res = 300)
ribo_heatmap_df[,1:6] %>% as.matrix() %>% t() %>% scale() %>% t() %>%
  heatmap.2(Rowv = FALSE, Colv = TRUE, labRow = FALSE, col = "bluered",
            scale = "none", dendrogram='column', cexCol = 1.2,
            trace = "none", density.info = "none", key = FALSE,
            RowSideColors = c(rep("indianred1", length(ribo_de_up)), rep("skyblue1", length(ribo_de_down))),
            hclustfun = function(x) hclust(x,method = 'ward.D'))
dev.off()



