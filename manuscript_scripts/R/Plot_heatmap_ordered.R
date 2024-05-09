Plot_heatmap_ordered <- function(x, ColSideColors, geneInfo, mainLab=schemaTmp) {
  if(ncol(as.matrix(ColSideColors))<=15) {
    ColSideColorsSize <- ncol(as.matrix(ColSideColors))
  } else {
    ColSideColorsSize <- 15
  }
  
  heatmap.3(main= mainLab
            ,x = as.matrix(x)
            ,Rowv = FALSE
            ,Colv = FALSE
            ,dendrogram = "none" ,trace="none" ,scale="row"
            ,labRow = geneInfo$geneSymbol
            ,labCol = colnames(x)
            ,col = colorRampPalette(c("darkblue","blue","white","red","darkred"))(99)
            ,ColSideColors = as.matrix(ColSideColors)
            ,ColSideColorsSize = ColSideColorsSize
            ,RowSideColors = t(as.matrix(geneInfo$Color))
            ,RowSideColorsSize=1
            ,lwid = c(1,5)
            ,lhei = c(1,3)
            ,margins =c(8,18)
            ,keysize=0.5
            ,cex.lab=0.4
            ,cexRow = 0.5
            ,cexCol = 0.6
  )
  legend(xy.coords(x=.8,y=.98),
         legend=c("DECODER_weight","High","Low",
                  "PurIST_score", "High","Low",
                  "PurIST_graded","Strong Classical","Likely Classical","Lean Classical","Lean Basal-like","Likely Basal-like","Strong Basal-like",
                  "PurIST","Basal-like","Classical",
                  "Moffitt","Basal-like","Classical",
                  "Collisson","Classical","Exocrine","QM",
                  "Bailey","ADEX","Immunogenic","Pancreatic Progenitor","Squamous" ,
                  "Chan-Seng-Yue","Basal-like A","Basal-like B","Hybrid","Classical A","Classical B",
                  "Puleo", "Pure Classical","Immune Classical","Desmoplastic","Stroma Activated","Pure Basal-like",
                  "MS","Absent","Activated","Normal",
                  "MSI","Activated","Immune","Normal",
                  "Elyada","apCAF","iCAF","myCAF",
                  "Maurer","ECM-rich","Immune-rich"),
         fill=c("white","black","snow2",
                "white","black","snow2",
                "white","#0000FF","#6666FF","#CCCCFF","#FFECCB","#FFC965","#FFA500",
                "white","orange","blue",
                "white","orange","blue",
                "white","steelblue1","forestgreen","magenta1" ,
                "white","darkorchid1","red","slategray2","gold1",
                "white","brown","orange","grey","blue3","skyblue",
                "white","dodgerblue4","green4","mediumpurple","darkorange2","red3",
                "white","black","brown","skyblue",
                "white","brown","red","skyblue",
                "white","darkslateblue","coral","darkgreen",     
                "white","purple3","forestgreen"),
         border=FALSE, bty="n",
         y.intersp = 0.85 , cex = 0.45)
}