datTmp <- data.frame(DeCAF = c("permCAF","permCAF","restCAF","restCAF"),
                     PurIST = c("Basal-like","Classical","Basal-like","Classical"),
                     Freq = c(167,529,92,536),
                     stringsAsFactors = FALSE)

bar_plots <- ggplot(datTmp, aes(x=factor(datTmp$PurIST, levels=c("Basal-like", "Classical")), y=Freq, fill=DeCAF)) +
  geom_bar(position="fill", stat = "identity")+
  scale_fill_manual(values = alpha(c("permCAF" = "#FE46A5", "restCAF" = "#008080"))) +
  geom_text(aes(label = paste0(Freq)), position = position_fill(vjust = 0.5),
            col="white", size=6)+
  #geom_text(aes(x = 1.5, y = 1.05,
  #              label = paste("Fisher's Exact, p =", signif(datStat$p.value, digits = 3))),
  #          hjust = 0, vjust = -0.5, color = "black", size = 5, family = "Helvetica")+
  theme_classic(base_size=15)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  )+
  labs(title = "Pooled public data",
       y = "Proportion of patients",
       x = "PurIST subtype"
  )


pdf("../figure_DeCAF/decaf_purist_bar.pdf")
bar_plots

dev.off()
