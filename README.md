# Neonate-and-Adult-Human-Platelet-RNASeq-and-RiboSeq-Data-Analyese
The R code used to analyze neonate and adult human platelet RNA-seq, Rib-seq and dose-response data. 
- Workflow:

![worflow](https://github.com/Zhaoyan030/Neonate-and-Adult-Human-Platelet-RNASeq-and-RiboSeq-Data-Analyese/blob/main/Workflow.png)

- [RNAseq_final.R](https://github.com/Zhaoyan030/Neonate-and-Adult-Human-Platelet-RNASeq-and-RiboSeq-Data-Analyese/blob/main/RNAseq_final.R): the R script to preprocess RNA-seq data, perform DE analysis, GSEA and other statistical analyses, and visualize results
- [Ribo_seq.R](https://github.com/Zhaoyan030/Neonate-and-Adult-Human-Platelet-RNASeq-and-RiboSeq-Data-Analyese/blob/main/Ribo_seq.R): the R script to preprocess Ribo-seq data, perform DE analysis, GSEA and other statistical analyses, and visualize results
- [RNAseq_Riboseq_integrated.R](https://github.com/Zhaoyan030/Neonate-and-Adult-Human-Platelet-RNASeq-and-RiboSeq-Data-Analyese/blob/main/RNAseq_Riboseq_integrated.R): the R script to merge RNA-seq and Ribo-seq data, and perform comparison study on the integrated data
- [functional_data_plots.R](https://github.com/Zhaoyan030/Neonate-and-Adult-Human-Platelet-RNASeq-and-RiboSeq-Data-Analyese/blob/main/functional_data_plots.R): the R script to analyze the platelet dose-response data, and visualize the results
- [my_netgraph_functions.R](https://github.com/Zhaoyan030/Neonate-and-Adult-Human-Platelet-RNASeq-and-RiboSeq-Data-Analyese/blob/main/my_netgraph_functions.R): the R script including the modifed functions (from [R package clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)) to generate network graphs for RNA-seq GSEA results
- [my_plot_functions.R](https://github.com/Zhaoyan030/Neonate-and-Adult-Human-Platelet-RNASeq-and-RiboSeq-Data-Analyese/blob/main/my_plot_functions.R): the R script including functions to generate dose-response curves and correlation scatter plots for dose-response data
