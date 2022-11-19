draw_h_v <- function(exprSet,need_DEG,n='DEseq2'){
  ## we only need two columns of DEG, which are log2FoldChange and pvalue
  ## heatmap
  #exprSet=sarc_choose
  #need_DEG=nrDEG
  library(pheatmap)
  exprSet=log(edgeR::cpm(exprSet)+1)
  need_DEG=need_DEG[order(need_DEG[,1]),]
  choose_gene=c(head(rownames(need_DEG),50),
                tail(rownames(need_DEG),50)) ## 50 maybe better
  choose_matrix=exprSet[choose_gene,]
  pheatmap(choose_matrix,
           filename = paste0(n,'_need_DEG_top100_raw.png'))
  dat=t(scale(t(choose_matrix)))
  dat[dat>2]=2 
  dat[dat< -2]= -2
  dat[1:4,1:4]  
  group_list=clinic$short.histo
  ac=data.frame(group_list=group_list)
  rownames(ac)=colnames(dat) 
  pheatmap(dat,show_colnames =F,show_rownames = F,
          annotation_col=ac,
           filename = paste0(n,'_need_DEG_top100_scale.png'))
  
  
  
  logFC_cutoff <- with(need_DEG,mean(abs( log2FoldChange)) + 2*sd(abs( log2FoldChange)) )
  # logFC_cutoff=1
  
  need_DEG$change = as.factor(ifelse(need_DEG$pvalue < 0.05 & abs(need_DEG$log2FoldChange) > logFC_cutoff,
                                     ifelse(need_DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                      '\nThe number of up gene is ',nrow(need_DEG[need_DEG$change =='UP',]) ,
                      '\nThe number of down gene is ',nrow(need_DEG[need_DEG$change =='DOWN',])
  )
  library(ggplot2)
  g = ggplot(data=need_DEG, 
             aes(x=log2FoldChange, y=-log10(pvalue), 
                 color=change)) +
    geom_point(alpha=0.4, size=1.75) +
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("log2 fold change") + ylab("-log10 p-value") +
    ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
  print(g)
  ggsave(g,filename = paste0(n,'_volcano.png'))
}
