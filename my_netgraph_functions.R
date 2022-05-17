gseaScores <- function(geneList, geneSet, exponent=1, fortify=FALSE) {
  geneSet <- intersect(geneSet, names(geneList))
  
  N <- length(geneList)
  Nh <- length(geneSet)
  
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet ## logical
  
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  
  Pmiss[!hits] <-  1/(N-Nh)
  Pmiss <- cumsum(Pmiss)
  
  runningES <- Phit - Pmiss
  
  ## ES is the maximum deviation from zero of Phit-Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if( abs(max.ES) > abs(min.ES) ) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  
  df <- data.frame(x=seq_along(runningES),
                   runningScore=runningES,
                   position=as.integer(hits)
  )
  
  if(fortify==TRUE) {
    return(df)
  }
  
  df$gene = names(geneList)
  res <- list(ES=ES, runningES = df)
  return(res)
}

get_gseaResult <- function(fgseaResult, collection, genelist){
  res <- data.frame(
    ID = as.character(fgseaResult$pathway),
    setSize = fgseaResult$size,
    enrichmentScore = fgseaResult$ES,
    NES = fgseaResult$NES,
    pvalue = fgseaResult$pval,
    p.adjust = fgseaResult$padj,
    stringsAsFactors = FALSE
  )
  res <- res[!is.na(res$pvalue),]
  res <- res[ res$p.adjust <= 0.05, ]
  idx <- order(res$pvalue, decreasing = FALSE)
  res <- res[idx, ]
  
  row.names(res) <- res$ID
  observed_info <- lapply(collection[res$ID], function(gs)
    gseaScores(geneSet=gs,
               geneList=genelist,
               exponent=1)
  )
  ledge <- leading_edge(observed_info)
  
  res$rank <- ledge$rank
  res$leading_edge <- ledge$leading_edge
  res$core_enrichment <- sapply(ledge$core_enrichment, paste0, collapse='/')
  
  params <- list(pvalueCutoff = 0.05,
                 nPerm = 100000,
                 pAdjustMethod = 0.05,
                 exponent = 1,
                 minGSSize = 1,
                 maxGSSize = Inf
  )
  
  gseaResult <- new("gseaResult",
                    result     = res,
                    geneSets   = hallmark,
                    geneList   = rna_ranks,
                    params     = params,
                    readable   = FALSE,
                    keytype    = "SYMBOL",
                    setType    = "DO")
  
  return(gseaResult)
}

leading_edge <- function(observed_info) {
  core_enrichment <- lapply(observed_info, function(x) {
    runningES <- x$runningES
    runningES <- runningES[runningES$position == 1,]
    ES <- x$ES
    if (ES >= 0) {
      i <- which.max(runningES$runningScore)
      leading_gene <- runningES$gene[1:i]
    } else {
      i <- which.min(runningES$runningScore)
      leading_gene <- runningES$gene[-c(1:(i-1))]
    }
    return(leading_gene)
  })
  
  rank <- sapply(observed_info, function(x) {
    runningES <- x$runningES
    ES <- x$ES
    if (ES >= 0) {
      rr <- which.max(runningES$runningScore)
    } else {
      i <- which.min(runningES$runningScore)
      rr <- nrow(runningES) - i + 1
    }
    return(rr)
  })
  
  tags <- sapply(observed_info, function(x) {
    runningES <- x$runningES
    runningES <- runningES[runningES$position == 1,]
    ES <- x$ES
    if (ES >= 0) {
      i <- which.max(runningES$runningScore)
      res <- i/nrow(runningES)
    } else {
      i <- which.min(runningES$runningScore)
      res <- (nrow(runningES) - i + 1)/nrow(runningES)
    }
    return(res)
  })
  
  ll <- sapply(observed_info, function(x) {
    runningES <- x$runningES
    ES <- x$ES
    if (ES >= 0) {
      i <- which.max(runningES$runningScore)
      res <- i/nrow(runningES)
    } else {
      i <- which.min(runningES$runningScore)
      res <- (nrow(runningES) - i + 1)/nrow(runningES)
    }
    return(res)
  })
  
  N <- nrow(observed_info[[1]]$runningES)
  setSize <- sapply(observed_info, function(x) sum(x$runningES$position))
  signal <- tags * (1-ll) * (N / (N - setSize))
  
  tags <- paste0(round(tags * 100), "%")
  ll <- paste0(round(ll * 100), "%")
  signal <- paste0(round(signal * 100), "%")
  leading_edge <- paste0('tags=', tags, ", list=", ll, ", signal=", signal)
  
  res <- list(rank = rank,
              tags = tags,
              list = ll,
              signal = signal,
              leading_edge = leading_edge,
              core_enrichment = core_enrichment)
  return(res)
}

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}

get_w <- function(y, geneSets){
  id <- y[, "ID"]
  geneSets <- geneSets[id]
  n <- nrow(y)
  w <- matrix(NA, nrow=n, ncol=n)
  colnames(w) <- rownames(w) <- y$ID
  for (i in seq_len(n-1)) {
    for (j in (i+1):n) {
      w[i,j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]])
    }
  }
  return(w)
}

emap_graph_build <- function(y, geneSets, color, cex_line, min_edge){
  w <- get_w(y, geneSets)
  wd <- reshape::melt(w)
  wd <- wd[wd[,1] != wd[,2],]
  wd <- wd[!is.na(wd[,3]),]
  
  g <- igraph::graph.data.frame(wd[, -3], directed=FALSE)
  igraph::E(g)$width <- sqrt(wd[, 3] * 5) * cex_line
  igraph::E(g)$weight <- wd[, 3]
  g <- igraph::delete.edges(g, igraph::E(g)[wd[, 3] < min_edge])
  idx <- unlist(sapply(igraph::V(g)$name, function(x) which(x == y$ID)))
  cnt <- sapply(geneSets[idx], length)
  igraph::V(g)$size <- cnt
  colVar <- y[idx, color]
  igraph::V(g)$color <- colVar
  return(g)
}

my_emapplot <- function(x, showCategory = 30, color="p.adjust",
                       layout = "nicely", min_edge=0.2,
                       cex_label_category  = 1, cex_category = 1,
                       cex_line = 1){
  y <- as.data.frame(x)
  geneSets <- geneInCategory(x)
  g <- emap_graph_build(y, geneSets, color=color, cex_line=cex_line, min_edge=min_edge)
  
  p <- ggraph::ggraph(g, layout=layout)
  if(length(igraph::E(g)$width) > 0) {
    p <- p + ggraph::geom_edge_link(alpha=.8, aes_(width=~I(width)),
                            colour='darkgrey')
  }
  p <- p + ggraph::geom_node_point(aes_(color=~color, size=~size))
  label_category <- 5
  p <- p + ggraph::geom_node_text(aes_(label=~name), repel=TRUE,
                          size = label_category * cex_label_category)

  p + theme_void() +
    scale_color_continuous(low="red", high="blue", name = color,
                           guide=guide_colorbar(reverse=TRUE)) +
    scale_size(range=c(3, 8) * cex_category)
}

list2graph <- function(inputList) {
  x <- list2df(inputList)
  g <- igraph::graph.data.frame(x, directed=FALSE)
  return(g)
}
list2df <- function(inputList) {
  ldf <- lapply(seq_len(length(inputList)), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })
  
  do.call('rbind', ldf)
}

my_cnetplot <- function(x,
                        showCategory = 5,
                        layout = "kk",
                        foldChange = NULL,
                        colorEdge = FALSE,
                        color = "p.adjust",
                        circular = FALSE,
                        genelist = NULL,
                        cex_category = 1,
                        cex_gene = 1,
                        cex_label_category = 1,
                        cex_label_gene = 1,
                        gene_label_repel = TRUE) {

  
  label_category <- 5
  label_gene <- 5
  if (circular) {
    layout <- "linear"
    geom_edge <- ggraph::geom_edge_arc
  } else {
    geom_edge <- ggraph::geom_edge_link
  }
  
  geneSets <- geneInCategory(x)
  
  g <- list2graph(geneSets)
  
  size <- sapply(geneSets, length)
  igraph::V(g)$size <- min(size)/2
  n <- length(geneSets)
  igraph::V(g)$size[1:n] <- size
  node_scales <- c(rep(cex_category, n), rep(cex_gene, (length(igraph::V(g)) - n)))
  if (colorEdge) {
    igraph::E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
    edge_layer <- geom_edge(aes_(color = ~category), alpha=.8)
  } else {
    edge_layer <- geom_edge(alpha=.8, colour='darkgrey')
  }
  
  y <- as.data.frame(x)
  idx <- unlist(sapply(igraph::V(g)$name, function(x) which(x == y$ID)))
  colVar <- y[idx, color]
  igraph::V(g)$color <- NA
  igraph::V(g)$color[1:n] <- colVar
  
  if (!is.null(foldChange)) {
    fc <- foldChange[igraph::V(g)$name[(n+1):length(igraph::V(g))]]
    igraph::V(g)$color[(n+1):length(igraph::V(g))] <- fc
    show_legend <- c(TRUE, FALSE)
    names(show_legend) <- c("color", "size")
    p <- ggraph::ggraph(g, layout=layout, circular = circular)
    p <- p + edge_layer +
      ggraph::geom_node_point(aes_(color=~as.numeric(color), size=~size),
                      data = p$data[1:n, ]) +
      scale_color_continuous(low="red", high="blue", name = color,
                             guide=guide_colorbar(reverse=TRUE))+
      scale_size(range=c(3, 8) * cex_category) +
      ggnewscale::new_scale("size") +
      ggnewscale::new_scale_color() +
      ggraph::geom_node_point(aes_(color=~as.numeric(color), size=~size),
                      data = p$data[(n+1):length(igraph::V(g)), ], show.legend = show_legend) +
      scale_size(range=c(3, 3) * cex_gene) +
      scale_colour_gradient2(name = "log2FC", low = "blue",
                             mid = "white", high = "red")+
      labs(colour = 'p.adjust')
  } else {
    p <- ggraph::ggraph(g, layout=layout, circular=circular)
    p <- p + edge_layer +
      ggraph::geom_node_point(aes_(color=~as.numeric(color), size=~size), data = p$data[1:n, ]) +
      scale_color_continuous(low="red", high="blue", name = color,
                             guide=guide_colorbar(reverse=TRUE))+
      scale_size(range=c(3, 8) * cex_category) +
      ggnewscale::new_scale("size") +
      ggraph::geom_node_point(aes_(size=~size),
                      data = p$data[(n+1):length(igraph::V(g)), ], show.legend = FALSE,
                      color="#B3B3B3") +
      scale_size(range=c(3, 3) * cex_gene) +
      labs(colour = 'p.adjust')
  }
  
  p <- p + theme_void()
  
  p <- p + ggraph::geom_node_text(aes_(label=~name), data = p$data[1:n,],repel=TRUE,
                          size = label_category * cex_label_category,
                          fontface = "bold")
  
  if (!is.null(genelist)){
    gene_idx <- which(igraph::V(g)$name %in% genelist)
    p <- p + ggraph::geom_node_text(aes_(label=~name), data = p$data[gene_idx,],
                            repel=gene_label_repel, size = label_gene * cex_label_gene)
  }
  else{
    p <- p + ggraph::geom_node_text(aes_(label=~name), data = p$data[-c(1:n),],
                            repel=gene_label_repel, size = label_gene * cex_label_gene)
  }
  
  return(p)
}

