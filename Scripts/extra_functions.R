multi_track_GSEA<-function (all_pathways,
                            pathway_of_interest,
                            stats, 
                            gseaParam=1,
                            ticksSize = 0.2) 
{
  
  pathways<-all_pathways[pathway_of_interest]
  stats=stats[!is.na(stats)]
  
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  midpoint<-median(which(statsAdj==0))
  
  pathway_results<-list()
  for(pathway_index in 1:length(pathways)){         
    pathway<- pathways[[pathway_index]]
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    pathway_results[[pathway_index]] <- data.frame(x = c(0, xs, n + 1),
                                             "pathway"=names(pathways)[[pathway_index]])
  }
  toPlot<-do.call(rbind,pathway_results)%>%
            mutate(pathway=factor(pathway,levels=rev(pathway_of_interest)))
  
  g <- ggplot(toPlot, aes(x = as.numeric(x), y = pathway,color=x)) +#      
                geom_point_rast(pch = "|",cex=ticksSize)+
                geom_vline(xintercept = midpoint,lty=2,color="black")+
                scale_color_distiller(palette = "RdYlBu",direction = 1)+
                theme_classic(base_size=16)+
                xlab("Gene rank")+
                guides(color="none")
  return(g)
}


add_flag<-function(pheatmap,
         kept.labels,
         repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  return(heatmap)
}


multi_sample_enrichment_plot<-function (pathway_of_interest, sample_list, gseaParam = 0.5, ticksSize = 0.2) 
{
  results<- list()
  for(i in 1:length(sample_list)){
    stats <- sample_list[[i]]
    pathway <- pathway_of_interest
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                            returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms))/8
    x = y = NULL
      results[[i]] <- list("toPlot"=toPlot,
                           "pathway"=pathway,
                           "tops"=tops,
                           "bottoms"=bottoms,
                           "diff"=diff)
    }
    names(results)<-names(sample_list)
    return(results)
}