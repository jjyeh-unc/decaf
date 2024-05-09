Plot_heatmap_CC <- function(tmp = cafSubtype[[tmp]], schemaTmp, ColSideColors, mainLab=schemaTmp) {
  i <- which(schemaList %in% schemaTmp)
  genesTmp <- subtypeGeneList[[i]]
  idxGene <- match(tmp$featInfo,genesTmp$geneSymbol)
  geneInfo <- genesTmp[idxGene,]
  
  if(ncol(as.matrix(ColSideColors))<=15) {
    ColSideColorsSize <- ncol(as.matrix(ColSideColors))
  } else {
    ColSideColorsSize <- 15
  }
  
  heatmap.3(main= mainLab
            ,x = as.matrix(tmp$ex)
            ,Rowv = tmp$Rowv
            ,Colv = tmp$Colv
            ,dendrogram = "both" ,trace="none" ,scale="row"
            ,labRow = tmp$featInfo
            ,labCol = tmp$sampID
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
                  "Maurer","ECM-rich","Immune-rich",
                  "MS","Absent","Activated","Normal",
                  "MSI","Activated","Immune","Normal",
                  "Elyada","iCAF","myCAF",
                  "SCISSORS","apCAF","iCAF","myCAF","panCAF","Perivascular"),
         fill=c("white","black","snow2",
                "white","black","snow2",
                "white","#0000FF","#6666FF","#CCCCFF","#FFECCB","#FFC965","#FFA500",
                "white","orange","blue",
                "white","orange","blue",
                "white","steelblue1","forestgreen","magenta1" ,
                "white","darkorchid1","red","slategray2","gold1",
                "white","brown","orange","grey","blue3","skyblue",
                "white","dodgerblue4","green4","mediumpurple","darkorange2","red3",
                "white","purple3","forestgreen",
                "white","black","brown","skyblue",
                "white","brown","red","skyblue",
                "white","coral","darkgreen",     
                "white","darkslateblue","coral","darkgreen","khaki","steelblue1"    ),
         border=FALSE, bty="n",
         y.intersp = 0.85 , cex = 0.45)
}
