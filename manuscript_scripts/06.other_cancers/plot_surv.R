############################## Functions and libraries ##############################
# set dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
library(ConsensusClusterPlus) # R3.xx and R4.xx have different versions of ConsensusClusterPlus
library("RColorBrewer")
library(openxlsx)
library(stringr)
library(survival)
library(survminer)

pdf("../jdamrauer/Desktop/CAF_Expression/CAF_Yeh/")
# TCGA --------------------------------------------------------------------


tcga<-read.table("New.TCGA.Survival.txt", sep="\t", header=T, row.names=1)
TCGA<-tcga %>%
  subset(M_stage %in% c("M0"))


TCGA <- within(TCGA, {
  DeCAF.UNC <- factor(DeCAF.UNC, levels=c("restCAF_Luminal", "permCAF_Luminal", "restCAF_Basal", "permCAF_Basal"))
  DeCAF.UNC.rev <- factor(DeCAF.UNC, levels=c("restCAF_Basal","restCAF_Luminal", "permCAF_Basal", "permCAF_Luminal" ))
  DeCAF <- factor(DeCAF, levels=c("restCAF", "permCAF"))
  UNC <- factor(UNC, levels=c("Luminal", "Basal"))
  
})

TCGA.C<-subset(TCGA, DeCAF.C %in% c("restCAF_LumP", "permCAF_LumP", "restCAF_Ba/Sq", "permCAF_Ba/Sq"))

TCGA.C <- within(TCGA.C, {
  DeCAF.C <- factor(DeCAF.C, levels=c("restCAF_LumP", "permCAF_LumP", "restCAF_Ba/Sq", "permCAF_Ba/Sq"))
  DeCAF.C.rev <- factor(DeCAF.C, levels=c("restCAF_Ba/Sq","restCAF_LumP", "permCAF_Ba/Sq", "permCAF_LumP"))
  
})

pdf("survival.plots.stratified.v2.pdf", height = 7, width = 7.5)

fit<- survfit(Surv(OS.time.5, OS.5) ~DeCAF , data=TCGA)


ggsurvplot(font.legend = 10,  legend = "right",
           font.title = 16,        
           font.x = 16,  
           font.y = 16, 
           font.tickslab = 16,
           risk.table.fontsize = 7,
           risk.table.height = 0.35,
           surv.median.line = "hv", 
           size = 0.3, 
           censor.size=7, 
           pval.size = 8,
           fit, 
           conf.int = F, pval = T,
           legend.title="",break.time.by = 12,
           legend.labs=c("restCAF","permCAF"),
           palette = c("turquoise4", "violetred1"),
           xlab = "Time (months)", risk.table = T,
           title = "TCGA: Overall Survival")

pdf("TCGA.forest.OS.pdf", width = 7.5, height = 3)
ggforest(coxph(Surv(OS.time.5, OS.5) ~DeCAF, data=TCGA))
dev.off()

fit<- survfit(Surv(OS.time.5, OS.5) ~DeCAF.UNC, data=TCGA)

ggsurvplot(font.legend = 10,  legend = "right",
           font.title = 16,        
           font.x = 16,  
           font.y = 16, 
           font.tickslab = 16,
           risk.table.fontsize = 7,
           risk.table.height = 0.35,
           surv.median.line = "hv", 
           size = 0.3, 
           censor.size=7, 
           pval.size = 8,
           fit, 
           conf.int = F, pval = T,
           legend.labs=c("restCAF_Luminal", "permCAF_Luminal", "restCAF_Basal", "permCAF_Basal"),
           palette = c("turquoise3", "forestgreen", "violetred1","violetred4"),
           legend.title="",break.time.by = 12,
           linetype = c(2,1,2,1),
           xlab = "Time (months)", risk.table = T,
           title = "TCGA: Overall Survival")

ggforest(coxph(Surv(OS.time.5, OS.5) ~DeCAF.UNC, data=TCGA))
ggforest(coxph(Surv(OS.time.5, OS.5) ~DeCAF.UNC.rev, data=TCGA))
 

fit<- survfit(Surv(OS.time.5, OS.5) ~DeCAF.C, data=TCGA.C)

ggsurvplot(font.legend = 10,  legend = "right",
           font.title = 16,        
           font.x = 16,  
           font.y = 16, 
           font.tickslab = 16,
           risk.table.fontsize = 7,
           risk.table.height = 0.35,
           surv.median.line = "hv", 
           size = 0.3, 
           censor.size=7, 
           pval.size = 8,
           fit, 
           conf.int = F, pval = T,
           legend.labs=c("restCAF_LumP", "permCAF_LumP", "restCAF_Ba/Sq", "permCAF_Ba/Sq"),
           palette = c("turquoise3", "forestgreen", "violetred1","violetred4"),
           legend.title="",break.time.by = 12,
           linetype = c(2,1,2,1),
           xlab = "Time (months)", risk.table = T,
           title = "TCGA: Overall Survival")

ggforest(coxph(Surv(OS.time.5, OS.5) ~DeCAF.C, data=TCGA.C))
ggforest(coxph(Surv(OS.time.5, OS.5) ~DeCAF.C.rev, data=TCGA.C))

dev.off()


fit<- survfit(Surv(OS.time.5, OS.5) ~DeCAF.C, data=TCGA.C)

ggsurvplot(font.legend = 16, data=TCGA.C,
                 font.title = 16,        
                 font.x = 16,  
                 font.y = 16, 
                 font.tickslab = 16,
                 risk.table.fontsize = 7,
                 risk.table.height = 0.35,
                 surv.median.line = "hv", 
                 size = 0.3, 
                 censor.size=7, 
                 pval.size = 8,
                 fit, 
                 conf.int = F, pval = T,
                 legend.title="",break.time.by = 12,
                 legend.labs=c("Basal","Luminal"),
                 palette = c("violetred1","turquoise4"),
                 linetype = c("strata"),
                 xlab = "Time (months)", risk.table = T,
                 title = "TCGA: Overall Survival")
ggforest(coxph(Surv(OS.time.5, OS.5) ~DeCAF.C, data=TCGA.C))
ggforest(coxph(Surv(OS.time.5, OS.5) ~DeCAF.C.rev, data=TCGA.C))


dev.off()




fit<- survfit(Surv(OS.time.5, OS.5) ~UNC, data=TCGA)

ggsurvplot_facet(font.legend = 16, data=TCGA, facet.by = "DeCAF",
                 font.title = 16,        
                 font.x = 16,  
                 font.y = 16, 
                 font.tickslab = 16,
                 risk.table.fontsize = 7,
                 risk.table.height = 0.35,
                 surv.median.line = "hv", 
                 size = 0.3, 
                 censor.size=7, 
                 pval.size = 8,
                 fit, 
                 conf.int = F, pval = T,
                 legend.title="",break.time.by = 12,
                 legend.labs=c("permCAF","restCAF"),
                 palette = c("orange","blue"),
                 linetype = c("strata"),
                 xlab = "Time (months)", risk.table = T,
                 title = "TCGA: Overall Survival")


ggforest(coxph(Surv(OS.time.5, OS.5) ~UNC + DeCAF, data=TCGA))


fit<- survfit(Surv(OS.time.5, OS.5) ~DeCAF + UNC, data=TCGA)

ggsurvplot(font.legend = 16,
           font.title = 16,        
           font.x = 16,  
           font.y = 16, 
           font.tickslab = 16,
           risk.table.fontsize = 7,
           risk.table.height = 0.35,
           surv.median.line = "hv", 
           size = 0.3, 
           censor.size=7, 
           pval.size = 8,
           fit, 
           conf.int = F, pval = T,
           legend.title="",break.time.by = 12,
           legend.labs=c("permCAF:Basal","permCAF:Luminal","restCAF:Basal","restCAF:Luminal"),
           palette = c("violetred1", "pink3", "turquoise4", "forestgreen"),
           xlab = "Time (months)", risk.table = T,
           title = "TCGA: Overall Survival")


fit<- survfit(Surv(OS.time.5, OS.5) ~UNC + DeCAF, data=TCGA)

ggsurvplot(font.legend = 16,
           font.title = 16,        
           font.x = 16,  
           font.y = 16, 
           font.tickslab = 16,
           risk.table.fontsize = 7,
           risk.table.height = 0.35,
           surv.median.line = "hv", 
           size = 0.3, 
           censor.size=7, 
           pval.size = 8,
           fit, 
           conf.int = F, pval = T,
           legend.title="",break.time.by = 12,
           legend.labs=c("Basal:permCAF", "Basal:restCAF", "Luminal:permCAF","Luminal:restCAF"),
           palette = c("violetred1", "pink3", "turquoise4", "forestgreen"),
           xlab = "Time (months)", risk.table = T,
           title = "TCGA: Overall Survival")



fit<- survfit(Surv(OS.time.5, OS.5) ~UNC, data=TCGA)

ggsurvplot_facet(font.legend = 16, facet.by="DeCAF",
           font.title = 16,        
           font.x = 16,  
           font.y = 16, 
           font.tickslab = 16,
           risk.table.fontsize = 7,
           risk.table.height = 0.35,
           surv.median.line = "hv", 
           size = 0.3, 
           censor.size=7, 
           pval.size = 8,
           fit, 
           conf.int = F, pval = T,
           legend.title="",break.time.by = 12,
           legend.labs=c("Basal:permCAF", "Basal:restCAF", "Luminal:permCAF","Luminal:restCAF"),
           palette = c("violetred1", "pink3", "turquoise4", "forestgreen"),
           xlab = "Time (months)", risk.table = T,
           title = "TCGA: Overall Survival")


fit<- survfit(Surv(OS.time.5, OS.5) ~UNC, data=TCGA)

ggsurvplot_facet(font.legend = 16, facet.by = "DeCAF", data=TCGA,
           font.title = 16,        
           font.x = 16,  
           font.y = 16, 
           font.tickslab = 16,
           risk.table.fontsize = 7,
           risk.table.height = 0.35,
           surv.median.line = "hv", 
           size = 0.3, 
           censor.size=7, 
           pval.size = 8,
           fit, 
           conf.int = F, pval = T,
           legend.title="",break.time.by = 12,
           legend.labs=c("Basal", "Luminal"),
           palette = c("violetred1", "turquoise4"),
           xlab = "Time (months)", risk.table = T,
           title = "TCGA: Overall Survival")



# IMVIG -------------------------------------------------------------------

imvig<-read.table("IMvigor.CAF.20231209.txt", sep="\t", header=T, row.names=1)

imvigor <- within(imvig, {
  DeCAF <- factor(DeCAF, levels=c("restCAF", "permCAF"))
})

fit<- survfit(Surv(os, censOS) ~DeCAF , data=imvig)

pdf("IMvigor210.OS.plot.pdf", height = 7, width = 7.5)

ggsurvplot(font.legend = 16, 
           font.title = 16,        
           font.x = 16,  
           font.y = 16, 
           font.tickslab = 16,
           risk.table.fontsize = 7,
           risk.table.height = 0.35,
           surv.median.line = "hv", 
           size = 0.3, 
           censor.size=7, 
           pval.size = 8,
           fit, 
           conf.int = F, pval = T,
           legend.title="",break.time.by = 12,
           legend.labs=c("permCAF","restCAF"),
           palette = c("violetred1","turquoise4"),
           xlab = "Time (months)", risk.table = T,
           title="IMvigor210")
dev.off()


pdf("IMvigor210.forest.plot.pdf", width = 7.5, height = 3)
ggforest(coxph(Surv(os, censOS) ~DeCAF , data=imvigor))
dev.off()


summary(coxph(Surv(os, censOS) ~DeCAF , data=imvigor))

surv_median(fit)

library(patchwork)
pdf("IMvigor210.updated.Survival.pdf", width = 7, height = 7)
imvigor.plot
dev.off()

# UROMOL ----------------------------------------------------------------

UROMOL<-read.table("UROMOL.CAF.20231209.txt", sep="\t", header=T, row.names=1)

fit<- survfit(Surv(Prog.Time, Progression) ~DeCAF , data=UROMOL)
uromol.pfs<-ggsurvplot(font.legend = 10, 
           font.title = 16,        
           font.x = 16,  
           font.y = 16, 
           font.tickslab = 16,
           risk.table.fontsize = 7,
           risk.table.height = 0.35,
           surv.median.line = "hv", 
           size = 0.3, 
           censor.size=7, 
           pval.size = 8,
           #pval.coord = c(96, 0.75),
           #xlim = c(0,144),
           fit, 
           conf.int = F, pval = T,
           legend.title="",break.time.by = 12,
           legend.labs=c("permCAF","restCAF"),
           palette = c("violetred1","turquoise4"),
           xlab = "Time (months)", risk.table = T,
           title = "UROMOL: Progression Free Survival")
           
uromol.1 <- within(UROMOL, {
  DeCAF <- factor(DeCAF, levels=c("restCAF", "permCAF"))
})


pdf("UROMOL.forest.PFS.Survival.pdf", width = 7.5, height = 3)
ggforest(coxph(Surv(Prog.Time, Progression) ~DeCAF, data=uromol.1))
dev.off()


pdf("UROMOL.RFS.pdf", width = 7.5, height = 7)

fit<- survfit(Surv(Recur.Time, Recur) ~DeCAF , data=UROMOL)
ggsurvplot(font.legend = 10, 
           font.title = 16,        
           font.x = 16,  
           font.y = 16, 
           font.tickslab = 16,
           risk.table.fontsize = 7,
           risk.table.height = 0.35,
           surv.median.line = "hv", 
           size = 0.3, 
           censor.size=7, 
           pval.size = 8,
           #pval.coord = c(96, 0.75),
           #xlim = c(0,144),
           fit, 
           conf.int = F, pval = T,
           legend.title="",break.time.by = 12,
           legend.labs=c("permCAF","restCAF"),
           palette = c("violetred1","turquoise4"),
           xlab = "Time (months)", risk.table = T,
           title = "UROMOL: Recurrence Free Survival")
dev.off()


UROMOL <- within(UROMOL, {
  DeCAF <- factor(DeCAF, levels=c("restCAF", "permCAF"))
})

summary(coxph(Surv(Prog.Time, Progression) ~DeCAF, data=UROMOL))
surv_median(fit)


URO.Surv<-subset(UROMOL, Stage %in% c("T1", "Ta"))

fit<- survfit(Surv(Recur.Time, Recur) ~DeCAF, data=URO.Surv)
uro.1<-ggsurvplot_facet(font.legend = 10, facet.by = c("Staging"), data=URO.Surv,
           font.title = 16,        
           font.x = 16,  
           font.y = 16, 
           font.tickslab = 16,
           risk.table.fontsize = 7,
           risk.table.height = 0.35,
           surv.median.line = "hv", 
           size = 0.3, 
           censor.size=7, 
           pval.size = 8,
           pval.coord = c(40, 0.85),
           #xlim = c(0,144),
           fit, 
           conf.int = F, pval = T,
           legend.title="",break.time.by = 12,
          legend.labs=c("permCAF","restCAF"),
          palette = c("violetred1","turquoise4"),
           xlab = "Time (months)", risk.table = T,
          ncol=5,
           title = "UROMOL: Recurrence Free Survival")

fit<- survfit(Surv(Prog.Time, Progression) ~DeCAF, data=URO.Surv)
uro.3<-ggsurvplot_facet(font.legend = 10, facet.by = c("Staging"), data=URO.Surv,
                 font.title = 16,        
                 font.x = 16,  
                 font.y = 16, 
                 font.tickslab = 16,
                 risk.table.fontsize = 7,
                 risk.table.height = 0.35,
                 surv.median.line = "hv", 
                 size = 0.3, 
                 censor.size=7, 
                 pval.size = 8,
                 #pval.coord = c(96, 0.75),
                 #xlim = c(0,144),
                 fit, 
                 conf.int = F, pval = T,
                 legend.title="",break.time.by = 12,
                 legend.labs=c("permCAF","restCAF"),
                 palette = c("violetred1","turquoise4"),
                 xlab = "Time (months)", risk.table = T,
                 ncol=5,
                 title = "UROMOL: Progression Free Survival")


URO.Surv <- within(URO.Surv, {
  DeCAF <- factor(DeCAF, levels=c("restCAF", "permCAF"))
  
})


pdf("UROMOL.Prog.Recur.forest.pdf", width = 7, height = 3)

ggforest(coxph(Surv(Prog.Time, Progression) ~DeCAF, data=URO.Surv),
         main = "UROMOL: Progression Free Survival")

ggforest(coxph(Surv(Recur.Time, Recur) ~DeCAF, data=URO.Surv), 
         main = "UROMOL: Recurrence Free Survival")

dev.off()


library(patchwork)
pdf("UROMOL.Prog.Recur.Survival.pdf", width = 12, height = 7)
grid.arrange(uro.1, uro.3, ncol = 1)
dev.off()


roc.plot<-URO.Surv %>%
  roc(Progression, DeCAF.prob, smooth=TRUE, plot=TRUE, ci=TRUE, na.rm=TRUE, ci.se=TRUE, ci.sp=TRUE)

roc.plot.1<-URO.Surv %>%
  roc(Progression, Ba.Sq, smooth=TRUE, plot=TRUE, ci=TRUE, na.rm=TRUE, ci.se=TRUE, ci.sp=TRUE)

coords(roc.plot, "best", ret=c("threshold", "specificity", "1-npv"))
coords(roc.plot.1, "local maximas", ret=c("threshold", "sens", "spec", "ppv", "npv"))
power.roc.test(auc=0.8, power=0.9)
roc.test(roc.plot, roc.plot.1, reuse.auc=FALSE)
