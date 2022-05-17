library(data.table)
library(magrittr)
library(dplyr)
library(reshape)
library(ggplot2)
library(ggpubr)
library(car)
library(corrplot)
library(drc)
library(gridExtra)
library(grid)
library(cowplot)
source("Program/my_plot_functions.R")

#######################################################################################################
#################################### Processing data ##################################################
#######################################################################################################

df <- fread("Data/my_functional_data.csv",nrow=20)
pselection_df <- df[,-grep("AnnexinV", names(df), value = TRUE), with=FALSE]
annexinv_df <- df[,- grep("Pselection", names(df), value = TRUE), with=FALSE]

pselection_df %<>%
  melt(., measure.vars = grep("Pselection", names(df), value = TRUE),value.name = "protein_level") 
cols <- colsplit(pselection_df$variable, split = "_", names = c('Agonist', 'Dose', 'Protein'))
pselection_df <- cbind(pselection_df, cols)
pselection_df$Dose_cat <- as.factor(pselection_df$Dose)

annexinv_df %<>% 
  melt(., measure.vars = grep("AnnexinV", names(df), value = TRUE),value.name = "protein_level") 
cols <- colsplit(annexinv_df$variable, split = "_", names = c('Agonist', 'Dose', 'Protein'))
annexinv_df <- cbind(annexinv_df, cols)
annexinv_df$Dose_cat <- as.factor(annexinv_df$Dose)


#### remove outliers
# pselection_outliers <- get_outliers(pselection_df, select_dose)
# annexinv_outliers <- get_outliers(annexinv_df, select_dose)

# outID <- c("NO95 6/15/17 F", "5/4/16 F", "12/3/15 #1 M")
outID <- c("NO95 6/15/17 F")
pselection_df <- pselection_df[! sampleID %in% outID,]
annexinv_df <- annexinv_df[! sampleID %in% outID,]


agonist <- pselection_df$Agonist %>% unique() %>% as.character()

#######################################################################################################
#################################### LOESS curve fitting ##############################################
#######################################################################################################

pselection_list <- get_dr_curves(pselection_df)
annexinv_list <- get_dr_curves(annexinv_df)
# curve_plot_list <- list()
# for(i in 1:length(pselection_list)){
#   curve_plot_list[[2*i-1]] <- pselection_list[[i]]
#   curve_plot_list[[2*i]] <- annexinv_list[[i]]
# }
dr1 <- ggarrange(plotlist=pselection_list, ncol = 1, nrow = 4,align = 'hv',common.legend = FALSE)
dr1 <- annotate_figure(dr1,top = text_grob("P-Selectin", face = "bold"))

dr2 <- ggarrange(plotlist=annexinv_list, ncol = 1, nrow = 4,align = 'hv',common.legend = FALSE)
dr2 <- annotate_figure(dr2,top = text_grob("Annexin V", face = "bold"))

shared_legend <- get_shared_legend(pselection_list)
p1 <- my_arrange(list(dr1,dr2), shared_legend, legend_pos = 'right', nrow = 1, ncol = 2)
# ggsave("Results/Dose_response_curves.jpeg", dpi = 700, height = 5.5, width = 4)


#######################################################################################################
################################## Correlation scatter plot ###########################################
#######################################################################################################

select_dose <- c("Thrombin"=0.25, "TRAP"=3, "ADP"=1.6, "PAR4"=45)
scatter_plot_list <- get_pa_scatters(pselection_df, annexinv_df, select_dose)
scatter_plot <- ggarrange(plotlist=scatter_plot_list, ncol = 2, nrow = 4,align = 'hv',
                        common.legend = TRUE, legend="right")
scatter_plot <- scatter_plot+theme(plot.margin = margin(0, 20, 0, 20))
p2 <- annotate_figure(scatter_plot,
                left = text_grob("Annexin V", rot = 90, face = 'bold'),
                bottom = text_grob("P-Selectin", face = 'bold'),
                fig.lab = "B", fig.lab.face = "bold", fig.lab.size=20)
# p2 <- p2+theme(plot.margin = margin(0, 20, 0, 20))
# ggsave("Results/Dose_response_scatters_AnnV&Psel.jpeg", dpi = 700, height = 5.5, width = 4)


scatter_plot_list2 <- get_pp_scatters(pselection_df, select_dose)
scatter_plot2 <- ggarrange(plotlist=scatter_plot_list2, ncol = 6, nrow = 2,align = 'hv',
                          common.legend = TRUE, legend="right")
scatter_plot2 <- scatter_plot2+theme(plot.margin = margin(20, 20, 20, 0))
p3 <- annotate_figure(scatter_plot2,
                left = text_grob("P-Selectin", rot = 90, face = 'bold'),
                bottom = text_grob("P-Selectin", face = 'bold'),
                fig.lab = "C", fig.lab.face = "bold", fig.lab.size=20)
# p3 <- p3+theme(plot.margin = margin(20, 20, 20, 0))
# ggsave("Results/Dose_response_scatters_Psel&Psel.jpeg", dpi = 700, height = 3, width = 10)


scatter_plot_list3 <- get_pp_scatters(annexinv_df, select_dose)
scatter_plot3 <- ggarrange(plotlist=scatter_plot_list3, ncol = 6, nrow = 2,align = 'hv',
                           common.legend = TRUE, legend="right")
scatter_plot3 <- scatter_plot3+theme(plot.margin = margin(20, 20, 20, 0))
p4 <- annotate_figure(scatter_plot3,
                left = text_grob("Annexin V", rot = 90, face = 'bold'),
                bottom = text_grob("Annexin V", face = 'bold'),
                fig.lab = "D", fig.lab.face = "bold", fig.lab.size=20)
# p4 <- p4+theme(plot.margin = margin(20, 20, 20, 0))
# ggsave("Results/Dose_response_scatters_AnnV&AnnV.jpeg", dpi = 700, height = 3, width = 10)


#######################################################################################################
############################### Arrange 4 panels of Figure1 ###########################################
#######################################################################################################
# lay <- rbind(c(1,1,2,2),
#              c(1,1,2,2),
#              c(3,3,3,3),
#              c(4,4,4,4))
# g <- grid.arrange(p1,p2,p3,p4, layout_matrix = lay)
top_row <- plot_grid(p1, p2, ncol = 2)
g <- plot_grid(
  top_row, p3, p4,
  ncol = 1,
  align = 'hv',
  rel_heights = c(2,1,1)
)
# g <- plot_grid(g,
#                labels = "Z. Liu et. al., Fig. 1",
#                label_fontfamily = "Arial",
#                label_size = 10,
#                label_x = 0.83,
#                label_y = 0.02)
ggsave("Results/test.jpeg",g, dpi = 700, height = 14, width = 12)


#######################################################################################################
#################################### hierarchical clustering ##########################################
#######################################################################################################


row.names(df) <- df$sampleID
dist <- df[,4:43] %>% scale() %>% dist()
hclust <- hclust(dist, method='ward.D')
plot(hclust, hang=-1, labels = df$sampleID)


#######################################################################################################
########################################### PCA #######################################################
#######################################################################################################


library(factoextra)
pca <- prcomp(df[,4:43], center = TRUE,scale. = TRUE)
fviz_pca_ind(pca)



#######################################################################################################
#################################### Sigmoid curve fitting ############################################
#######################################################################################################

jpeg("Results/drc_curves2.jpeg", res = 300, height = 10, width = 7,units = 'in')
par(mfcol=c(4,2))
ymax <- c(100,100,60,100)
# ticks <- list(c(0,0.1,0.25,0.5,1),
#               c(0,0.75,1,1.5,3,6),
#               c(0,0.4,0.8,1,1.6,3.6),
#               c(0,10,25,45,65,85,100))
ec50_df <- data.frame(Protein = rep(c("P-selection","Annexin V"),4), 
                      Agonist = rep(agonist, each = 2),
                      Adult = NA,
                      CB = NA,
                      pvalue = NA)

for(j in 1:2){
  data <- list(pselection_df,annexinv_df)[[j]]
  for(i in 1:length(agonist)){
    ag <- agonist[i]
    d <- subset(data, Agonist==ag)
    # out <- d[,boxplot.stats(protein_level)$out, by=.(Dose,group)]$V1
    # out_id <- d[protein_level %in% out, sampleID]
    
    DR <- drm(protein_level ~ Dose, 
              data= d,
              group,
              robust = 'mean',
              fct = LL.4(names = c("Hill slope", "Min", "Max", "EC50")))
    ec50_df$Adult[2*(i-1)+j] <- DR$coefficients['EC50:Adult']
    ec50_df$CB[2*(i-1)+j] <- DR$coefficients['EC50:CB']
    ec50_df$pvalue[2*(i-1)+j] <- compParm(DR, "EC50", "-", display = FALSE)[4]
    
    DR.cb <- drm(protein_level ~ Dose, 
                 # data= subset(d, group=='CB' & (! sampleID %in% out_id)),
                 data= subset(d, group=='CB'),
                 robust = 'mean',
                 fct = LL.4(names = c("Hill slope", "Min", "Max", "EC50")))
    DR.n <- drm(protein_level ~ Dose, 
                # data= subset(d, group=='Adult' & (! sampleID %in% out_id)),
                data= subset(d, group=='Adult'),
                robust = 'mean',
                fct = LL.4(names = c("Hill slope", "Min", "Max", "EC50")))
    plot(DR.n,
         col = 'black',
         # xlim = c(0, 1000),
         ylim = c(0, ymax[i]),
         type='average',
         legend = FALSE,
         xlab= "",
         ylab = "",
         pch=16,
         lwd=1)
    plot(DR.n,
         col = 'black',
         # xlim = c(0, 1000),
         ylim = c(0, ymax[i]),
         add=TRUE,
         type='confidence',
         confidence.level = 0.90)
    segments(x0 = DR$coefficients['EC50:Adult'], y0=0, 
             x1 = DR$coefficients['EC50:Adult'], 
             y1 = predict(DR.n,newdata = data.frame(Dose=DR$coefficients['EC50:Adult'])),
             col = 'black', lty = 'dashed')
    par(new=TRUE)
    plot(DR.cb,
         col = 'red',
         # xlim = c(0, upper),
         ylim = c(0, ymax[i]),
         type='average',
         pch=16,
         lwd=1,
         xlab= "",
         ylab = "")
    plot(DR.cb,
         col = 'red',
         # xlim = c(0, upper),
         ylim = c(0, ymax[i]),
         add=TRUE,
         type='confidence',
         confidence.level = 0.90)
    segments(x0 = DR$coefficients['EC50:CB'], y0=0, 
             x1 = DR$coefficients['EC50:CB'], 
             y1 = predict(DR.cb,newdata = data.frame(Dose=DR$coefficients['EC50:CB'])),
             col = 'red', lty = 'dashed')
  }
}
dev.off()



#######################################################################################################
#################################### Correlation matrix ###############################################
#######################################################################################################

a <- df[group=='Adult',.(Thrombin_0.25_Pselection,TRAP_3_Pselection,ADP_1.6_Pselection,PAR4_45_Pselection, 
                         Thrombin_0.25_AnnexinV,  TRAP_3_AnnexinV,ADP_1.6_AnnexinV, PAR4_45_AnnexinV)] %>% as.matrix()
colnames(a) <- lapply(colnames(a), 
                     FUN=function(x) paste(unlist(strsplit(x, "_"))[c(3,1)],collapse = '_')) %>% unlist()
b <- df[group=='CB',.(Thrombin_0.25_Pselection,TRAP_3_Pselection,ADP_1.6_Pselection,PAR4_45_Pselection, 
                         Thrombin_0.25_AnnexinV,  TRAP_3_AnnexinV,ADP_1.6_AnnexinV, PAR4_45_AnnexinV)] %>% as.matrix()
colnames(b) <- lapply(colnames(b), 
                      FUN=function(x) paste(unlist(strsplit(x, "_"))[c(3,1)],collapse = '_')) %>% unlist()
M.n <-cor(a)
P.n <- cor.mtest(a, conf.level = 0.95)[[1]]
M.cb <-cor(b)
P.cb <- cor.mtest(b, conf.level = 0.95)[[1]]
for(i in 1:nrow(M.n)){
  j <- which(! grepl(paste(unlist(strsplit(row.names(M.n)[i], "_")),collapse = '|'),colnames(M.n)))
  M.n[i,j] <- 0
  P.n[i,j] <- 0
  
  M.cb[i,j] <- 0
  P.cb[i,j] <- 0
}
corrplot(M.n, type="lower", order="original", diag = FALSE, p.mat = P.n, sig.level = 0.05,
         tl.cex = 0.8, tl.srt = 45)
corrplot(M.cb, type="lower", order="original", diag = FALSE, p.mat = P.cb, sig.level = 0.05,
         tl.cex = 0.8, tl.srt = 45)

jpeg("Results/corrplot_Adult_pselection.jpeg", res = 300, height = 4.5, width = 4.5,units = 'in')
corrplot(M.n[1:4,1:4], type="lower", order="original", diag = FALSE, p.mat = P.n[1:4,1:4], sig.level = 0.05,
         tl.cex = 0.8, tl.srt = 45, tl.pos = 'n', pch.cex = 10, cl.pos = 'n')
dev.off()

jpeg("Results/corrplot_Adult_annexinv.jpeg", res = 300, height = 4.5, width = 4.5,units = 'in')
corrplot(M.n[5:8,5:8], type="lower", order="original", diag = FALSE, p.mat = P.n[5:8,5:8], sig.level = 0.05,
         tl.cex = 0.8, tl.srt = 45, tl.pos = 'n', pch.cex = 10, cl.pos = 'n')
dev.off()

jpeg("Results/corrplot_Adult_p&a.jpeg", res = 300, height = 6, width = 6,units = 'in')
corrplot(M.n[1:4,5:8], type="full", order="original", diag = TRUE, p.mat = P.n[1:4,5:8], sig.level = 0.05,
         tl.cex = 0.8, tl.srt = 45, tl.pos = 'n', pch.cex = 10, cl.pos = 'b')
dev.off()

jpeg("Results/corrplot_CB_pselection.jpeg", res = 300, height = 4.5, width = 4.5,units = 'in')
corrplot(M.cb[1:4,1:4], type="lower", order="original", diag = FALSE, p.mat = P.cb[1:4,1:4], sig.level = 0.05,
         tl.cex = 0.8, tl.srt = 45, tl.pos = 'n', pch.cex = 10, cl.pos = 'n')
dev.off()

jpeg("Results/corrplot_CB_annexinv.jpeg", res = 300, height = 4.5, width = 4.5,units = 'in')
corrplot(M.cb[5:8,5:8], type="lower", order="original", diag = FALSE, p.mat = P.cb[5:8,5:8], sig.level = 0.05,
         tl.cex = 0.8, tl.srt = 45, tl.pos = 'n', pch.cex = 10, cl.pos = 'n')
dev.off()

jpeg("Results/corrplot_CB_p&a.jpeg", res = 300, height = 6, width = 6,units = 'in')
corrplot(M.cb[1:4,5:8], type="full", order="original", diag = TRUE, p.mat = P.cb[1:4,5:8], sig.level = 0.05,
         tl.cex = 0.8, tl.srt = 45, tl.pos = 'n', pch.cex = 10, cl.pos = 'b')
dev.off()








#######################################################################################################
######################################## Two-way ANOVA ################################################
#######################################################################################################
anova_res_df <- data.frame(PSelectin = rep(NA,4), AnnexinV = rep(NA,4))
row.names(anova_res_df) <- agonist
for(ag in agonist){
  d <- pselection_df[Agonist == ag, .(group, Dose, protein_level)]
  d$Dose <- factor(d$Dose, 
                   levels = sort(unique(d$Dose)),
                   labels = paste("D",sort(unique(d$Dose)), sep = ''))
  aov.res <- aov(protein_level~group+Dose, data=d)
  anova_res_df[ag, 1] <- summary(aov.res)[[1]][["Pr(>F)"]][1]
}
for(ag in agonist){
  d <- annexinv_df[Agonist == ag, .(group, Dose, protein_level)]
  d$Dose <- factor(d$Dose, 
                   levels = sort(unique(d$Dose)),
                   labels = paste("D",sort(unique(d$Dose)), sep = ''))
  aov.res <- aov(protein_level~group+Dose, data=d)
  anova_res_df[ag, 2] <- summary(aov.res)[[1]][["Pr(>F)"]][1]
}

write.csv(anova_res_df, "Results/ANOVA_result.csv")


