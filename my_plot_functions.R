get_outliers <- function(data, select_dose){
  agonist <- data$Agonist %>% unique() %>% as.character()
  outlier_list <- list()
  for(i in 1:length(agonist)){
    ag <- agonist[i]
    d <- subset(data, Agonist==ag & Dose==select_dose[ag])
    d[,protein_level:=log(protein_level/(100-protein_level))]
    out <- d[,boxplot.stats(protein_level)$out, by=.(group)]$V1
    out_id <- d[protein_level %in% out, sampleID]
    outlier_list[[i]] <- out_id
  }
  names(outlier_list) <- agonist
  return(outlier_list)
}



ggdrc <- function(data){
  agonist <- data$Agonist %>% unique() %>% as.character()
  plot_list <- list()
  for(i in 1:length(agonist)){
    ag <- agonist[i]
    d <- subset(data, Agonist==ag)
    out <- d[,boxplot.stats(protein_level)$out, by=.(Dose,group)]$V1
    out_id <- d[protein_level %in% out, sampleID]
    ticks <- d[Agonist==ag, unique(Dose)]
    
    p <- ggplot()+
      geom_point(data=d, mapping = aes(x=Dose, y=protein_level, color=group))+
      theme_classic()+
      theme(legend.position = "none", axis.title.y=element_blank())+
      scale_y_continuous(breaks = c(0,25,50,75,100), limits = c(-5,110))+
      scale_color_manual(values = c("Adult"='black', "CB"="red"))
    
    for(k in c("Adult","CB")){
      DR.mi <- drm(protein_level ~ Dose, 
                   data= subset(d, group==k & (! sampleID %in% out_id)),
                   robust = 'mean',
                   fct = LL.4(names = c("Hill slope", "Min", "Max", "EC50")))
      upper <- DR.mi$coefficients[3] %>% as.numeric() %>% ceiling()
      fits <- expand.grid(conc=seq(0, upper, length=100))
      pm <- predict(DR.mi, newdata=fits, interval="confidence") 
      fits$p <- pm[,1]
      fits$pmin <- pm[,2]
      fits$pmax <- pm[,3]
      p <- p+
        geom_ribbon(data=fits, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2, fill='red') +
        geom_line(data=fits, aes(x=conc, y=p), color='red')+
        scale_x_continuous(breaks = ticks, limits = c(0,upper))
    }
    
    plot_list[[i]] <- p
  }
}

plotdrc <- function(data){
  agonist <- data$Agonist %>% unique() %>% as.character()
  for(i in 1:length(agonist)){
    ag <- agonist[i]
    d <- subset(data, Agonist==ag)
    out <- d[,boxplot.stats(protein_level)$out, by=.(Dose,group)]$V1
    out_id <- d[protein_level %in% out, sampleID]
    
    DR.cb <- drm(protein_level ~ Dose, 
                 data= subset(d, group=='CB' & (! sampleID %in% out_id)),
                 robust = 'mean',
                 fct = LL.4(names = c("Hill slope", "Min", "Max", "EC50")))
    DR.n <- drm(protein_level ~ Dose, 
                data= subset(d, group=='Adult' & (! sampleID %in% out_id)),
                robust = 'mean',
                fct = LL.4(names = c("Hill slope", "Min", "Max", "EC50")))
    plot(DR.n,
         col = 'black',
         # xlim = c(0, 1000),
         ylim = c(0, 100),
         type='average',
         legend = FALSE,
         xlab= ag,
         ylab = "",
         pch=16,
         lwd=1)
    plot(DR.n,
         col = 'black',
         # xlim = c(0, 1000),
         ylim = c(0, 100),
         add=TRUE,
         type='confidence',
         confidence.level = 0.90)
    par(new=TRUE)
    plot(DR.cb,
         col = 'red',
         # xlim = c(0, upper),
         ylim = c(0, 100),
         type='average',
         pch=16,
         lwd=1,
         xlab= ag,
         ylab = "")
    plot(DR.cb,
         col = 'red',
         # xlim = c(0, upper),
         ylim = c(0, 100),
         add=TRUE,
         type='confidence',
         confidence.level = 0.90)
  }
}


findInt <- function(model, value) {
  function(x) {
    predict(model, data.frame(Dose=x)) - value
  }
}

get_dr_curves <- function(data){
  plot_list <- list()
  agonist <- data$Agonist %>% unique() %>% as.character()
  len <- length(agonist)
  
  for(i in 1:len){
    ag <- agonist[i]
    ticks <- data[Agonist==ag, unique(Dose)]
    
    
    # get EC50
    d_cb <- subset(data, Agonist==ag & group=="CB")
    model_cb <- loess(protein_level ~ Dose, data=d_cb, span = 0.95)
    ec50_cb <- d_cb[Dose==max(Dose), mean(protein_level)] / 2
    x_cb <- uniroot(findInt(model_cb, ec50_cb), range(d_cb$Dose))$root
    # assign(paste0('x_cb',i),uniroot(findInt(model_cb, ec50_cb), range(d_cb$Dose))$root)
    
    d_n <- subset(data, Agonist==ag & group=="Adult")
    model_n <- loess(protein_level ~ Dose, data=d_n, span = 0.95)
    ec50_n <- d_n[Dose==max(Dose), mean(protein_level)] / 2
    x_n <- uniroot(findInt(model_n, ec50_n), range(d_n$Dose))$root
    # assign(paste0('x_n',i),uniroot(findInt(model_n, ec50_n), range(d_n$Dose))$root)
    
    seg_df <- data.frame(x_cb = x_cb,
                         x_n = x_n,
                         y_cb = ec50_cb,
                         y_n = ec50_n)
    
    # get two-way ANOVA p-value
    d <- data[Agonist == ag, .(group, Dose, protein_level)]
    d$Dose <- factor(d$Dose, 
                     levels = sort(unique(d$Dose)),
                     labels = paste("D",sort(unique(d$Dose)), sep = ''))
    aov.res <- aov(protein_level~group+Dose, data=d)
    aov.p <- summary(aov.res)[[1]][["Pr(>F)"]][1]
    p.label <- ifelse(aov.p<0.01, 
                      paste0("=",formatC(aov.p, format = "e", digits = 2)), 
                      paste0("=",round(aov.p,2)))
    
    
    g <- ggplot()+
      geom_smooth(data=d_cb, aes(x=Dose, y=protein_level, color='CB', fill='CB'),
                  method = "loess", formula = y~x, span=0.95, level=0.9)+
      geom_segment(data=seg_df, aes(x=x_cb, xend=x_cb, y=0, yend=y_cb), linetype=2, color='red')+
      geom_smooth(data=d_n, aes(x=Dose, y=protein_level, color='Adult', fill='Adult'),
                  method = "loess", formula = y~x, span=0.95, level=0.9)+
      geom_segment(data=seg_df, aes(x= x_n, xend=x_n, y=0, yend=y_n), linetype=2)+
      theme_classic()+
      theme(legend.position = "none", axis.title.y=element_blank(),
            plot.margin = margin(0, 10, 0, 0),
            axis.title.x = element_text(size=8))+
      scale_x_continuous(breaks=ticks, expand = c(0, 0))+
      scale_y_continuous(breaks = c(0,25,50,75,100), limits = c(0,110), expand = c(0, 0))+
      scale_color_manual(values = c("Adult"='black', "CB"="red"), name='group')+
      scale_fill_manual(values = c("Adult"='grey50', "CB"="pink"), name='group')+
      annotate("text",label=bquote(EC50[Adult]~"="~.(round(x_n,2))~","~EC50[CB]~"="~.(round(x_cb,2))), 
               x=max(ticks)/2, y=100, color='black', size=3)+
      # annotate("text",label=bquote(EC50[CB]~"="~.(round(x_cb,2))), 
      #          x=max(ticks)/4, y=100-15, color='red',size=3)+
      annotate("text",label=bquote(p~.(p.label)), 
               x=max(ticks)/4, y=100-15, color='brown3',size=3)
    
    if(ag=="Thrombin"){
      g <- g+xlab("Thrombin (U/mL)")
    }else{
      g <- g+xlab(bquote(.(ag)~(mu*M)))
    }
    plot_list[[i]] <- g
    rm(g)
  }
  return(plot_list)
}


get_pa_scatters <- function(data1, data2, select_dose){
  plot_list <- list()
  agonist <- data1$Agonist %>% unique() %>% as.character()
  len <- length(agonist)
  
  for(i in 1:len){
    ag <- agonist[i]
    j <- 1
    for(k in c('Adult', 'CB')){
      d_n1 <- subset(data1, Agonist==ag & Dose==select_dose[ag] & group==k)
      d_n2 <- subset(data2, Agonist==ag & Dose==select_dose[ag] & group==k)
      d <- merge(d_n1[,.(sampleID, group, protein_level)], d_n2[,.(sampleID, protein_level)],
                 by = 'sampleID')
      d[,`:=`(x = log(protein_level.x/(100-protein_level.x)), 
              y = log(protein_level.y/(100-protein_level.y)))]
      
      cor <- cor.test(d$x, d$y)
      p.col <- ifelse(cor$p.value<0.1, "red", "black")
      p.label <- ifelse(cor$p.value<0.01, "<0.01", paste0("=",round(cor$p.value,2)))
      
      g <- ggplot(data = d, aes(x=x, y=y, color=group))+
        geom_point(size=1.5)+
        theme_classic()+
        theme(axis.title=element_text(size=8),
              legend.position = "none")+
        scale_x_continuous(limits = c(-4,3), breaks = c(-4:3))+
        scale_y_continuous(limits = c(-4,3), breaks = c(-4:3))+
        scale_color_manual(values = c("Adult"='black', "CB"="red"))+
        annotate("text", label=bquote(R~"="~.(round(cor$estimate,2))), x=-2, y=2.5, size=3)+
        annotate("text", label=bquote(p~.(p.label)), x=-2, y=1, color=p.col, size=3)
      
      if(ag=="Thrombin"){
        g <- g+
          xlab("Thrombin (0.25 U/mL)")+
          ylab("Thrombin (0.25 U/mL)")
      }else{
        g <- g+
          xlab(bquote(.(ag)~(.(select_dose[ag])~mu*M)))+
          ylab(bquote(.(ag)~(.(select_dose[ag])~mu*M)))
      }

      plot_list[[2*i-j]] <- g
      j <- j-1
    }
  }
  return(plot_list)
}


get_pp_scatters <- function(data1, select_dose){
  plot_list <- list()
  agonist <- data1$Agonist %>% unique() %>% as.character()
  len <- length(agonist)
  l <- 1
  for(k in c('Adult', 'CB')){
    for(i in 1:(len-1)){
      ag1 <- agonist[i]
      for(m in (i+1):len){
        ag2 <- agonist[m]
        
        d_n1 <- subset(data1, Agonist==ag1 & Dose==select_dose[ag1] & group==k)
        d_n2 <- subset(data1, Agonist==ag2 & Dose==select_dose[ag2] & group==k)
        d <- merge(d_n1[,.(sampleID, group, protein_level)], d_n2[,.(sampleID, protein_level)],
                   by = 'sampleID')
        d[,`:=`(x = log(protein_level.x/(100-protein_level.x)), 
                y = log(protein_level.y/(100-protein_level.y)))]
        
        cor <- cor.test(d$x, d$y)
        p.col <- ifelse(cor$p.value<0.1, "red", "black")
        p.label <- ifelse(cor$p.value<0.01, "<0.01", paste0("=",round(cor$p.value,2)))
        
        g <- ggplot(data = d, aes(x=x, y=y, color=group))+
          geom_point(size=1.5)+
          theme_classic()+
          theme(legend.position = "none", axis.title=element_text(size=8))+
          scale_x_continuous(limits = c(-4,3), breaks = c(-4:3))+
          scale_y_continuous(limits = c(-4,3), breaks = c(-4:3))+
          scale_color_manual(values = c("Adult"='black', "CB"="red"))+
          annotate("text", label=bquote(R~"="~.(round(cor$estimate,2))), x=-2, y=2.5, size=3)+
          annotate("text", label=bquote(p~.(p.label)), x=-2, y=1.5, color=p.col, size=3)+
          ylab(bquote(.(ag2)~(.(select_dose[ag2])~mu*M)))
        
        if(ag1=="Thrombin"){
          g <- g+xlab("Thrombin (0.25 U/mL)")
        }else{
          g <- g+xlab(bquote(.(ag1)~(.(select_dose[ag1])~mu*M)))
        }
        
        if(cor$p.value<0.1){
          g <- g+theme(panel.border = element_rect(color = "green",fill = NA,size = 2))
        }else{
          if(cor$p.value<0.15){
            g <- g+theme(panel.border = element_rect(color = "blue",fill = NA,size = 2))
          }
        }
        
        
        plot_list[[l]] <- g
        l <- l+1
      }
    }
  }
  return(plot_list)
}







get_shared_legend <-function(plots) {
  g <-ggplotGrob(plots[[1]] + theme(legend.position = 'right'))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  return(legend)
}

my_arrange <- function(plots,legend,legend_pos,nrow, ncol){
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  p <- ggarrange(plotlist=plots, ncol = ncol, nrow = nrow,align = 'hv',common.legend = FALSE)
  p <- annotate_figure(p, left = text_grob("Response (% Positive)", rot = 90, size = 10),
                       fig.lab = "A", fig.lab.face = "bold", fig.lab.size=20)
  gl <- c(list(p), ncol = 1, nrow = 1)
  
  combined <- switch(
    legend_pos,
    "bottom" = arrangeGrob(
      do.call(arrangeGrob, gl),
      legend,
      ncol = 1,
      heights = unit.c(unit(1, "npc") - lheight, lheight)
    ),
    "right" = arrangeGrob(
      do.call(arrangeGrob, gl),
      legend,
      ncol = 2,
      widths = unit.c(unit(1, "npc") - lwidth, lwidth)
    )
  )
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}