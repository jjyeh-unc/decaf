#####################################################################################

##                       Survival Figure Generating Script                  ############

####################################################################################




######################################################################################
############## Step 1: Create Functions for Generating Figures ########################
######################################################################################




# load libraries
library(survival)
library(survminer)
library(gtable)
library(gridExtra)
library(grid)
library(gplots)
library(coin)
library(xlsx)
library(ggpattern)
library(patchwork)

## Colors for DeCAF
DeCAFList = c("permCAF","restCAF")
DeCAFCol = c("violetred1","turquoise4")


#### Load Clean Survival Data ####
Load_clean_survival <- function(rDataName) {
  
  dataSet =  read.csv(sprintf("clean_data/survival/%s_survival.csv", rDataName ))
  return(dataSet)
}



#### Plot Survival Curves

# plot survival curve
Plot_survival <- function(survDat, km, mainLabel, DeCAF_graph, save, individual, schemaList){
  pos = 96
  if(mainLabel == "CPTAC"){
    pos = 30
  }
  if(mainLabel == "Dijk"){
    pos = 84
  }
  if(mainLabel == "Grunwald"){
    pos = 42
  }
  
  if(mainLabel == "Linehan"){
    pos = 36
  }
  
  if(mainLabel == "Moffitt_GEO_array"){
    pos = 42
  }
  if(mainLabel == "PACA_AU_array"){
    pos = 36
  }
  if(mainLabel == "PACA_AU_seq"){
    pos = 36
  }
  if(mainLabel == "TCGA_PAAD"){
    pos = 52
  }
  
  splots <- list()
  smodels <- list()
  # Elyada
  index_elyada = which(schemaList == "Elyada_CAF")
  if(individual == F){
  p = try({ coxph(km ~ Elyada_CAF + strata(dataSet), data = survDat) })
  } else {
    p = try({ coxph(km ~ Elyada_CAF , data = survDat) })
  }
  smodels[[1]] = p
  if(is.character(p)){
    km_fit <- survfit(km ~ Elyada_CAF, data = survDat, type = "kaplan-meier")
    splots[[1]] <- ggsurvplot(font.legend = 16,                            
                              font.title = 16,                             
                              font.x = 16,                               
                              font.y = 16,                               
                              font.tickslab = 16,  
                              risk.table.fontsize = 5.5,
                              risk.table.height = 0.35,
                              surv.median.line = "hv",
                              size = 0.3,
                              censor.size = 3,
                              pval.size = 8,
                              pval.coord =c(pos,0.75),
                              km_fit, conf.int = F, pval = T,
                              legend.title="",break.time.by = 12,
                              legend.labs="myCAF",
                              palette = "darkgreen",
                              xlab = "Time (months)", risk.table = T,
                              title = paste(mainLabel," - Elyada ",sep=""))
    
    
  } else{
    bic = round(BIC(p),3)
    hr <- round(summary(p)$coefficients[2],3)
    km_fit <- survfit(km ~ Elyada_CAF, data = survDat, type = "kaplan-meier")
    splots[[1]] <- ggsurvplot(font.legend = 16,
                              font.title = 16,
                              font.x = 16, 
                              font.y = 16, 
                              font.tickslab = 16, 
                              risk.table.fontsize = 5.5,
                              risk.table.height = 0.35,
                              surv.median.line = "hv",
                              size = 0.3,
                              censor.size = 3,
                              pval.size = 8,
                              pval.coord =c(pos,0.75),
                              km_fit, conf.int = F, pval = T,
                              legend.title="",break.time.by = 12,
                              legend.labs=subtypeList[[index_elyada]],
                              palette = subtypeColList[[index_elyada]],
                              xlab = "Time (months)", risk.table = T,
                              title = paste(mainLabel," - Elyada \n myCAF vs iCAF HR=",hr," (BIC=",
                                            bic,")",
                                            sep=""))
    
  }
  
  
  if(!is.character(p)){
    hr_out = c(NA, round(summary(p)$coefficients[2],3))
    LCL_hr_0.95 = c(NA, round(summary(p)$conf.int[4],3))
    UCL_hr_0.95 = c(NA, round(summary(p)$conf.int[3],3))
    p_val_out = c(NA, round(summary(p)$coefficients[5],4))
    
  }
  
  if(is.character(p)){
    hr_out = NA
    LCL_hr_0.95 = NA
    UCL_hr_0.95 = NA
    p_val_out = NA
    
  }
  if(!all.equal(surv_pvalue(survfit(km ~ Elyada_CAF, data = survDat))$pval,
                1-pchisq(survdiff(km ~ Elyada_CAF , data = survDat)$chisq, df = 1))){
    stop("My two log rank p values dont match something is wrong!")
  }
  log_rank_p_out = round(surv_pvalue(survfit(km ~ Elyada_CAF, data = survDat))$pval,4)
  if(surv_pvalue(survfit(km ~ Elyada_CAF, data = survDat))$pval < 0.0001){
    log_rank_p_out = "<0.0001"
  }
  stable <- data.frame(N = km_fit$n, Event = summary(km_fit)$table[,"events"],
                       HR = round(hr_out,3),
                       LCL_hr_0.95 = round(LCL_hr_0.95,3),
                       UCL_hr_0.95 = round(UCL_hr_0.95,3),
                       Median = round(summary(km_fit)$table[,"median"],3),
                       LCL_0.95 = round(summary(km_fit)$table[,"0.95LCL"],3),
                       UCL_0.95 = round(summary(km_fit)$table[,"0.95UCL"],3),
                       Year_0.5 = round(summary(km_fit, times = 6, extend = T)$surv,3),
                       Year_1 = round(summary(km_fit, times = 12, extend = T)$surv,3),
                       Year_2 = round(summary(km_fit, times = 24, extend = T)$surv,3),
                       Year_5 = round(summary(km_fit, times = 5*12, extend = T)$surv,3),
                       p_val = p_val_out,
                       log_rank_p = c(NA, log_rank_p_out))
  stable = stable[rev(1:nrow(stable)),]
  if(!is.na(stable$p_val[1])){
    stable$p_val[1] = round(stable$p_val[1],4)
    if(
      stable$p_val[1] < 0.0001){
      stable$p_val[1] = "<0.0001"
    } else{
      stable$p_val[1] = round(stable$p_val[1],4)
    }
  }
  stable3 = stable
  
  

  
  
  # MoffittStroma
  index_moffit = which(schemaList == "MS_K2")
  if(individual == F){
    p = try({ coxph(km ~ MS_K2 + strata(dataSet), data = survDat) })
  } else {
    p = try({ coxph(km ~ MS_K2 , data = survDat) })
  }
  smodels[[2]] = p
  
  if(is.character(p)){
    km_fit <- survfit(km ~ MS_K2, data = survDat, type = "kaplan-meier")
    splots[[2]] <- ggsurvplot(font.legend = 16,                               
                              font.title = 16,                               
                              font.x = 16,                                
                              font.y = 16,                                
                              font.tickslab = 16,  
                              risk.table.fontsize = 5.5,
                              risk.table.height = 0.35,
                              surv.median.line = "hv",
                              size = 0.3,
                              censor.size = 3,
                              pval.size = 8,
                              pval.coord =c(pos,0.75),                              
                              km_fit, conf.int = F, pval = T,
                              legend.title="",break.time.by = 12,
                              legend.labs="Activated",
                              palette = "brown",
                              xlab = "Time (months)", risk.table = T,
                              title = paste(mainLabel," - Moffitt Stroma ",sep=""))
    
  } else{
    bic = round(BIC(p),3)
    hr <- round(1/summary(p)$coefficients[2],3)
    km_fit <- survfit(km ~ MS_K2, data = survDat, type = "kaplan-meier")
    splots[[2]] <- ggsurvplot(font.legend = 16,                               
                              font.title = 16,                               
                              font.x = 16,             
                              font.y = 16,                                
                              font.tickslab = 16,  
                              risk.table.fontsize = 5.5,
                              risk.table.height = 0.35,
                              surv.median.line = "hv",
                              size = 0.3,
                              censor.size = 3,
                              pval.size = 8,
                              pval.coord =c(pos,0.75),                              
                              km_fit, conf.int = F, pval = T,
                              legend.title="",break.time.by = 12,
                              legend.labs=subtypeList[[index_moffit]],
                              palette = subtypeColList[[index_moffit]],
                              xlab = "Time (months)", risk.table = T,
                              title = paste(mainLabel," - Moffitt Stroma \n Activated vs Normal HR=",hr," (BIC=",
                                            bic,")",
                                            sep=""))
    
  }
  
  
  if(!is.character(p)){
    hr_out = c(round(1/summary(p)$coefficients[2],3), NA)
    LCL_hr_0.95 = c(round(1/summary(p)$conf.int[4],3), NA)
    UCL_hr_0.95 = c(round(1/summary(p)$conf.int[3],3), NA)
    p_val_out = c(round(summary(p)$coefficients[5],4), NA)
    
  }
  
  if(is.character(p)){
    hr_out = NA
    LCL_hr_0.95 = NA
    UCL_hr_0.95 = NA
    p_val_out = NA
    
  }
  if(!all.equal(surv_pvalue(survfit(km ~ MS_K2, data = survDat))$pval,
                1-pchisq(survdiff(km ~ MS_K2 , data = survDat)$chisq, df = 1))){
    stop("My two log rank p values dont match something is wrong!")
  }
  log_rank_p_out = round(surv_pvalue(survfit(km ~ MS_K2, data = survDat))$pval,4)
  if(surv_pvalue(survfit(km ~ MS_K2, data = survDat))$pval < 0.0001){
    log_rank_p_out = "<0.0001"
  }
  stable <- data.frame(N = km_fit$n, Event = summary(km_fit)$table[,"events"],
                       HR = round(hr_out,3),
                       LCL_hr_0.95 = round(LCL_hr_0.95,3),
                       UCL_hr_0.95 = round(UCL_hr_0.95,3),
                       Median = round(summary(km_fit)$table[,"median"],3),
                       LCL_0.95 = round(summary(km_fit)$table[,"0.95LCL"],3),
                       UCL_0.95 = round(summary(km_fit)$table[,"0.95UCL"],3),
                       Year_0.5 = round(summary(km_fit, times = 6, extend = T)$surv,3),
                       Year_1 = round(summary(km_fit, times = 12, extend = T)$surv,3),
                       Year_2 = round(summary(km_fit, times = 24, extend = T)$surv,3),
                       Year_5 = round(summary(km_fit, times = 5*12, extend = T)$surv,3),
                       p_val = p_val_out,
                       log_rank_p = c(log_rank_p_out, NA))
  if(!is.na(stable$p_val[1])){
    stable$p_val[1] = round(stable$p_val[1],4)
    if(
      stable$p_val[1] < 0.0001){
      stable$p_val[1] = "<0.0001"
    } else{
      stable$p_val[1] = round(stable$p_val[1],4)
    }
  }
  stable2 = stable

  
  # SCISSORS_CAF K2 Model

  if(individual == F){
    p = try({ coxph(km ~ SCISSORS_CAF_K2 + strata(dataSet), data = survDat) })
  } else {
    p = try({ coxph(km ~ SCISSORS_CAF_K2 , data = survDat) })
  }
  smodels[[3]] = p
  
  bic = round(BIC(p),3)
  hr <- round(1/summary(p)$coefficients[2],3)
  km_fit <- survfit(km ~ SCISSORS_CAF_K2, data = survDat, type = "kaplan-meier")
  splots[[3]] <- ggsurvplot(font.legend = 16,                               
                            font.title = 16,                               
                            font.x = 16,                               
                            font.y = 16,                               
                            font.tickslab = 16,     
                            risk.table.fontsize = 5.5,
                            risk.table.height = 0.35,
                            surv.median.line = "hv",
                            size = 0.3,
                            censor.size = 3,
                            pval.size = 8,
                            pval.coord =c(pos,0.75),                        
                            km_fit, conf.int = F, pval = T,
                            legend.labs=DeCAFList,
                            legend.title="",break.time.by = 12,
                            palette = DeCAFCol,
                            xlab = "Time (months)", risk.table = T,
                            title = paste(mainLabel," - SCISSORS CAF \n permCAF vs restCAF HR=",hr,
                                          " (BIC=",
                                          bic,")",sep=""))
  
  
  if(!is.character(p)){
    hr_out = c(round(1/summary(p)$coefficients[2],3), NA)
    LCL_hr_0.95 = c(round(1/summary(p)$conf.int[4],3), NA)
    UCL_hr_0.95 = c(round(1/summary(p)$conf.int[3],3), NA)
    p_val_out = c(round(summary(p)$coefficients[5],4), NA)
    
  }
  
  if(is.character(p)){
    hr_out = NA
    LCL_hr_0.95 = NA
    UCL_hr_0.95 = NA
    p_val_out = NA
    
  }
  if(!all.equal(surv_pvalue(survfit(km ~ SCISSORS_CAF_K2, data = survDat))$pval,
                1-pchisq(survdiff(km ~ SCISSORS_CAF_K2 , data = survDat)$chisq, df = 1))){
    stop("My two log rank p values dont match something is wrong!")
  }
  log_rank_p_out = round(surv_pvalue(survfit(km ~ SCISSORS_CAF_K2, data = survDat))$pval,4)
  if(surv_pvalue(survfit(km ~ SCISSORS_CAF_K2, data = survDat))$pval < 0.0001){
    log_rank_p_out = "<0.0001"
  }
  stable <- data.frame(N = km_fit$n, Event = summary(km_fit)$table[,"events"],
                       HR = round(hr_out,3),
                       LCL_hr_0.95 = round(LCL_hr_0.95,3),
                       UCL_hr_0.95 = round(UCL_hr_0.95,3),
                       Median = round(summary(km_fit)$table[,"median"],3),
                       LCL_0.95 = round(summary(km_fit)$table[,"0.95LCL"],3),
                       UCL_0.95 = round(summary(km_fit)$table[,"0.95UCL"],3),
                       Year_0.5 = round(summary(km_fit, times = 6, extend = T)$surv,3),
                       Year_1 = round(summary(km_fit, times = 12, extend = T)$surv,3),
                       Year_2 = round(summary(km_fit, times = 24, extend = T)$surv,3),
                       Year_5 = round(summary(km_fit, times = 5*12, extend = T)$surv,3),
                       p_val = p_val_out,
                       log_rank_p = c(log_rank_p_out, NA))
  if(!is.na(stable$p_val[1])){
    stable$p_val[1] = round(stable$p_val[1],4)
    if(
      stable$p_val[1] < 0.0001){
      stable$p_val[1] = "<0.0001"
    } else{
      stable$p_val[1] = round(stable$p_val[1],4)
    }
  }
  stable1 = stable
  
  # DeCAF
  if(DeCAF_graph == TRUE){
    
    if(individual == F){
      p = try({ coxph(km ~ DeCAF + strata(dataSet), data = survDat) })
    } else {
      p = try({ coxph(km ~ DeCAF , data = survDat) })
    }
    smodels[[4]] = p
    
    if(is.character(p)){
      km_fit <- survfit(km ~ DeCAF, data = survDat, type = "kaplan-meier")
      
      
      
        splots[[4]] <- ggsurvplot(font.legend = 16,                               
                                font.title = 16,                              
                                font.x = 16,                                
                                font.y = 16,                               
                                font.tickslab = 16,    
                                risk.table.fontsize = 5.5,
                                risk.table.height = 0.35,
                                surv.median.line = "hv",
                                size = 0.3,
                                censor.size = 3,
                                pval.size = 8,
                                pval.coord =c(pos,0.75),                          
                                km_fit, conf.int = F, pval = T,
                                legend.title="",break.time.by = 12,
                                legend.labs="permCAF",
                                palette = "violetred1",
                                xlab = "Time (months)", risk.table = T,
                                title = paste(mainLabel," DeCAF",sep=""))
      
    } else{
      bic = round(BIC(p),3)
      hr <- round(1/summary(p)$coefficients[2],3)
      km_fit <- survfit(km ~ DeCAF, data = survDat, type = "kaplan-meier")
      splots[[4]] <- ggsurvplot(font.legend = 16,                
                                font.title = 16,               
                                font.x = 16,                   
                                font.y = 16,                    
                                font.tickslab = 16,  
                                risk.table.fontsize = 5.5,
                                risk.table.height = 0.35,
                                surv.median.line = "hv",
                                size = 0.3,
                                censor.size = 3,
                                pval.size = 8,
                                pval.coord =c(pos,0.75),            
                                km_fit, conf.int = F, pval = T,
                                legend.title="",break.time.by = 12,
                                legend.labs=DeCAFList,
                                palette = DeCAFCol,
                                xlab = "Time (months)", risk.table = T,
                                title = paste(mainLabel," - DeCAF \n permCAF vs restCAF HR=",hr,
                                              " (BIC=",
                                              bic,")",sep=""))
      
    }
    if(!is.character(p)){
      hr_out = c(round(1/summary(p)$coefficients[2],3), NA)
      LCL_hr_0.95 = c(round(1/summary(p)$conf.int[4],3), NA)
      UCL_hr_0.95 = c(round(1/summary(p)$conf.int[3],3), NA)
      p_val_out = c(round(summary(p)$coefficients[5],4), NA)
      
    }
    
    if(is.character(p)){
      hr_out = NA
      LCL_hr_0.95 = NA
      UCL_hr_0.95 = NA
      p_val_out = NA
      
    }
    if(!all.equal(surv_pvalue(survfit(km ~ DeCAF, data = survDat))$pval,
              1-pchisq(survdiff(km ~ DeCAF , data = survDat)$chisq, df = 1))){
      stop("My two log rank p values dont match something is wrong!")
    }
    log_rank_p_out = round(surv_pvalue(survfit(km ~ DeCAF, data = survDat))$pval,4)
    if(surv_pvalue(survfit(km ~ DeCAF, data = survDat))$pval < 0.0001){
      log_rank_p_out = "<0.0001"
    }
    stable <- data.frame(N = km_fit$n, Event = summary(km_fit)$table[,"events"],
                         HR = round(hr_out,3),
                         LCL_hr_0.95 = round(LCL_hr_0.95,3),
                         UCL_hr_0.95 = round(UCL_hr_0.95,3),
                         Median = round(summary(km_fit)$table[,"median"],3),
                         LCL_0.95 = round(summary(km_fit)$table[,"0.95LCL"],3),
                         UCL_0.95 = round(summary(km_fit)$table[,"0.95UCL"],3),
                         Year_0.5 = round(summary(km_fit, times = 6, extend = T)$surv,3),
                         Year_1 = round(summary(km_fit, times = 12, extend = T)$surv,3),
                         Year_2 = round(summary(km_fit, times = 24, extend = T)$surv,3),
                         Year_5 = round(summary(km_fit, times = 5*12, extend = T)$surv,3),
                         p_val = p_val_out,
                         log_rank_p = c(log_rank_p_out, NA))
    if(!is.na(stable$p_val[1])){
      stable$p_val[1] = round(stable$p_val[1],4)
      if(
        stable$p_val[1] < 0.0001){
        stable$p_val[1] = "<0.0001"
      } else{
        stable$p_val[1] = round(stable$p_val[1],4)
      }
    }
  }
  stable77 =  stable
  ## UNC SCISSORS_CAF ALL Model
  if(individual == F){
    p = try({ coxph(km ~ SCISSORS_CAF_ALL + strata(dataSet), data = survDat) })
  } else {
    p = try({ coxph(km ~ SCISSORS_CAF_ALL , data = survDat) })
  }
  smodels[[5]] = p
  
  bic = round(BIC(p),3)
  hr <- round(1/summary(p)$coefficients[2],3)
  survDat$SCISSORS_CAF_ALL = factor(survDat$SCISSORS_CAF_ALL, levels = c("Absent","restCAF", "Mixed.restCAF", "Mixed", "Mixed.permCAF", "permCAF"))
  km_fit <- survfit(km ~ SCISSORS_CAF_ALL, data = survDat, type = "kaplan-meier")
  if(nrow(km_fit) !=  8){
    splots[[5]] <- ggsurvplot(font.legend = 16, 
                              font.title = 16,   
                              font.x = 16,         
                              font.y = 16,              
                              font.tickslab = 16,  
                              risk.table.fontsize = 5.5,
                              risk.table.height = 0.35,
                              surv.median.line = "hv",
                              size = 0.3,
                              censor.size = 3,
                              pval.size = 8,
                              pval.coord =c(pos,0.75),     
                              km_fit, conf.int = F, pval = F,
                              legend.title="",break.time.by = 12,
                              xlab = "Time (months)", risk.table = T,
                              title = paste(mainLabel," - SCISSORS CAF \n (BIC=",
                                            bic,")",sep=""),
                              legend.labs = sub("SCISSORS_CAF_ALL=","",names(km_fit$strata)))
    
  } else{
    splots[[5]] <- ggsurvplot(font.legend = 16,      
                              font.title = 16,       
                              font.x = 16,             
                              font.y = 16,             
                              font.tickslab = 16,   
                              risk.table.fontsize = 5.5,
                              risk.table.height = 0.35,
                              surv.median.line = "hv",
                              size = 0.3,
                              censor.size = 3,
                              pval.size = 8,
                              pval.coord =c(pos,0.75),     
                              km_fit, conf.int = F, pval = F,
                              legend.title="",break.time.by = 12,
                              xlab = "Time (months)", risk.table = T,
                              title = paste(mainLabel," - SCISSORS CAF \n (BIC=",
                                            bic,")",sep=""),
                              legend.labs = sub("SCISSORS_CAF_ALL=","",names(km_fit$strata)))
    
  }
  
 
  
 
  ## SCISSORS_CAF K2 - Random Assignment Model
  set.seed(1)
  index_list = which(survDat$SCISSORS_CAF_ALL %in% c("Mixed","Absent"))
  index_subset = sample(index_list, floor(length(index_list)/2))
  
  survDat$SCISSORS_CAF_ALL = survDat$SCISSORS_CAF_ALL
  survDat$SCISSORS_CAF_ALL[survDat$SCISSORS_CAF_ALL == "Mixed.permCAF"] = "permCAF"
  survDat$SCISSORS_CAF_ALL[survDat$SCISSORS_CAF_ALL == "Mixed.restCAF"] = "restCAF"
  
  survDat$SCISSORS_CAF_ALL[index_subset] = "permCAF"
  survDat$SCISSORS_CAF_ALL[ survDat$SCISSORS_CAF_ALL %in% c("Mixed","Absent")] = "restCAF"
  
  survDat$SCISSORS_CAF_ALL = factor(
    survDat$SCISSORS_CAF_ALL, levels = c("restCAF", "permCAF"))
  
  
  if(individual == F){
    p = try({ coxph(km ~ SCISSORS_CAF_ALL + strata(dataSet), data = survDat) })
  } else {
    p = try({ coxph(km ~ SCISSORS_CAF_ALL , data = survDat) })
  }
  smodels[[6]] = p
  p_value_r = 1- pchisq(p$score, df = 1)
    if(p_value_r < 0.0001){
       p_value_r = "p < 0.0001"
     } else {
       p_value_r = round(p_value_r,5)
     }
 
  bic = round(BIC(p),3)
  hr <- round(1/summary(p)$coefficients[2],3)
  km_fit <- survfit(km ~ SCISSORS_CAF_ALL, data = survDat, type = "kaplan-meier")
  splots[[6]] <- ggsurvplot(font.legend = 16,           
                            font.title = 16,             
                            font.x = 16,                    
                            font.y = 16,                      
                            font.tickslab = 16,   
                            risk.table.fontsize = 5.5,
                            risk.table.height = 0.35,
                            surv.median.line = "hv",
                            size = 0.3,
                            censor.size = 3,
                            pval.size = 8,
                            pval.coord =c(pos,0.75),              
                            km_fit, conf.int = F, pval = p_value_r,
                            legend.title="",break.time.by = 12,
                            legend.labs=DeCAFList[2:1],
                            palette = DeCAFCol[2:1],
                            xlab = "Time (months)", risk.table = T,
                            
                            title = paste(mainLabel," - SCISSORS CAF \n permCAF vs restCAF HR=",hr,
                                          " (BIC=",
                                          bic,")",sep=""))
  

  index_maurer = which(schemaList == "Maurer")
  if(individual == F){
    p = try({ coxph(km ~ Maurer + strata(dataSet), data = survDat) })
  } else {
    p = try({ coxph(km ~ Maurer , data = survDat) })
  }
  smodels[[7]] = p
  if(is.character(p)){
    km_fit <- survfit(km ~ Maurer, data = survDat, type = "kaplan-meier")
    splots[[7]] <- ggsurvplot(font.legend = 16,                            
                              font.title = 16,                             
                              font.x = 16,                               
                              font.y = 16,                               
                              font.tickslab = 16,   
                              risk.table.fontsize = 5.5,
                              risk.table.height = 0.35,
                              surv.median.line = "hv",
                              size = 0.3,
                              censor.size = 3,
                              pval.size = 8,
                              pval.coord =c(pos,0.75),   
                              legend.labs="ECM-rich",
                              palette = "purple3",
                              km_fit, conf.int = F, pval = T,
                              legend.title="",break.time.by = 12,
                              xlab = "Time (months)", risk.table = T,
                              title = paste(mainLabel," - Maurer ",sep=""))
    
    
  } else{
    bic = round(BIC(p),3)
    hr <- round(1/summary(p)$coefficients[2],3)
    km_fit <- survfit(km ~ Maurer, data = survDat, type = "kaplan-meier")
    splots[[7]] <- ggsurvplot(font.legend = 16,
                              font.title = 16,
                              font.x = 16, 
                              font.y = 16, 
                              font.tickslab = 16, 
                              risk.table.fontsize = 5.5,
                              risk.table.height = 0.35,
                              surv.median.line = "hv",
                              size = 0.3,
                              censor.size = 3,
                              pval.size = 8,
                              pval.coord =c(pos,0.75),
                              km_fit, conf.int = F, pval = T,
                              legend.title="",break.time.by = 12,
                              legend.labs=subtypeList[[index_maurer]],
                              palette = subtypeColList[[index_maurer]],
                              xlab = "Time (months)", risk.table = T,
                              title = paste(mainLabel," - Maurer \n ECM-rich vs Immune-rich HR=",hr," (BIC=",
                                            bic,")",
                                            sep=""))
    
  }
  
  
  if(!is.character(p)){
    hr_out = c(round(1/summary(p)$coefficients[2],3), NA)
    LCL_hr_0.95 = c(round(1/summary(p)$conf.int[4],3), NA)
    UCL_hr_0.95 = c(round(1/summary(p)$conf.int[3],3), NA)
    p_val_out = c(round(summary(p)$coefficients[5],4), NA)
    
  }
  
  if(is.character(p)){
    hr_out = NA
    LCL_hr_0.95 = NA
    UCL_hr_0.95 = NA
    p_val_out = NA
    
  }
  if(!is.na(surv_pvalue(survfit(km ~ Maurer, data = survDat))$pval)){
    
    if(!all.equal(surv_pvalue(survfit(km ~ Maurer, data = survDat))$pval,
                  1-pchisq(survdiff(km ~ Maurer , data = survDat)$chisq, df = 1))){
      stop("My two log rank p values dont match something is wrong!")
    }
    log_rank_p_out = round(surv_pvalue(survfit(km ~ Maurer, data = survDat))$pval,4)
    if(surv_pvalue(survfit(km ~ Maurer, data = survDat))$pval < 0.0001){
      log_rank_p_out = "<0.0001"
    }
    
      }
  if(is.na(surv_pvalue(survfit(km ~ Maurer, data = survDat))$pval)){
    log_rank_p_out = NA
  }
  if("numeric" %in% class(summary(km_fit)$table) ){
    stable <- data.frame(N = km_fit$n, Event = summary(km_fit)$table["events"],
                         HR = round(hr_out,3)[1],
                         LCL_hr_0.95 = round(LCL_hr_0.95,3)[1],
                         UCL_hr_0.95 = round(UCL_hr_0.95,3)[1],
                         Median = round(summary(km_fit)$table["median"],3),
                         LCL_0.95 = round(summary(km_fit)$table["0.95LCL"],3),
                         UCL_0.95 = round(summary(km_fit)$table["0.95UCL"],3),
                         Year_0.5 = round(summary(km_fit, times = 6, extend = T)$surv,3)[1],
                         Year_1 = round(summary(km_fit, times = 12, extend = T)$surv,3)[1],
                         Year_2 = round(summary(km_fit, times = 24, extend = T)$surv,3)[1],
                         Year_5 = round(summary(km_fit, times = 5*12, extend = T)$surv,3)[1],
                         p_val = p_val_out[1],
                         log_rank_p = c(log_rank_p_out))
    
    label_m = survDat$Maurer[!is.na(survDat$Maurer)][1]
    rownames(stable) = paste("Maurer=",label_m,sep = "")
    }
  
  if("matrix" %in% class(summary(km_fit)$table) ){
    stable <- data.frame(N = km_fit$n, Event = summary(km_fit)$table[,"events"],
                         HR = round(hr_out,3),
                         LCL_hr_0.95 = round(LCL_hr_0.95,3),
                         UCL_hr_0.95 = round(UCL_hr_0.95,3),
                         Median = round(summary(km_fit)$table[,"median"],3),
                         LCL_0.95 = round(summary(km_fit)$table[,"0.95LCL"],3),
                         UCL_0.95 = round(summary(km_fit)$table[,"0.95UCL"],3),
                         Year_0.5 = round(summary(km_fit, times = 6, extend = T)$surv,3),
                         Year_1 = round(summary(km_fit, times = 12, extend = T)$surv,3),
                         Year_2 = round(summary(km_fit, times = 24, extend = T)$surv,3),
                         Year_5 = round(summary(km_fit, times = 5*12, extend = T)$surv,3),
                         p_val = p_val_out,
                         log_rank_p = c(log_rank_p_out, NA))
    
  }
  if(!is.na(stable$p_val[1])){
    if(
      stable$p_val[1] < 0.0001){
      stable$p_val[1] = "<0.0001"
    } else{
      stable$p_val[1] = round(stable$p_val[1],4)
    }
  }
  stable4 = stable
  
  
  if(save == TRUE){
    if(individual == TRUE){
      filename = paste("data/results/figures/Individual_Survival_Plots/", mainLabel,"_Survival.pdf", sep ="")
    } else {
      filename = "data/results/figures/Pooled_Survival_Plots/Pooled_Survival.pdf"
    }
    pdf(filename)
    arrange_ggsurvplots(splots, 
                        print = TRUE,
                        ncol = 1, nrow = 1)
    dev.off()
    
  } 
  s_output = list(splots, smodels,stable77,stable1,stable2,stable3,stable4)
  return(s_output)

  
}





# plot survival curve
Plot_survival_interaction <- function(survDat, km, mainLabel, DeCAF_graph, save, individual, schemaList){
  pos = 96
  splots <- list()
  smodels <- list()
  # Elyada
  if(individual == F){
    p = try({ coxph(km ~ Elyada_CAF + PurIST + strata(dataSet), data = survDat) })
  } else {
    p = try({ coxph(km ~ Elyada_CAF + PurIST , data = survDat) })
  }
  smodels[[1]] <- p
  
  
  
  if(!"try-error" %in% class(p)){
    bic = round(BIC(p),3)
    hr <- round(1/summary(p)$coefficients[2],3)
  }
    km_fit <- survfit(km ~ Elyada_CAF+PurIST, data = survDat, type = "kaplan-meier")
    if(!is.character(p)){
      hr_out = rbind(round(summary(p)$conf.int,3)[1,c(1,3,4)],
                     round(1/summary(p)$conf.int,3)[2,c(1,4,3)])
      hr_out = as.data.frame(hr_out)
      hr_out$p_val = (summary(p)$coefficients[,5])
      for(qq in 1:nrow(hr_out)){
        if(hr_out$p_val[qq]<0.0001){
          hr_out$p_val[qq]="<0.0001"
        } else {
          hr_out$p_val[qq] = round(summary(p)$coefficients[,5],4)[qq]
        }
      }
    }
    if(!all.equal(surv_pvalue(survfit(km ~ Elyada_CAF + PurIST, data = survDat), method = "survdiff")$pval,
                  1-pchisq(survdiff(km ~ Elyada_CAF + PurIST, data = survDat)$chisq, df = (length(survdiff(km ~ Elyada_CAF + PurIST, data = survDat)$n)-1)))){
      stop("My two log rank p values dont match something is wrong!")
    }
    log_rank_p_out = surv_pvalue(survfit(km ~ Elyada_CAF + PurIST, data = survDat), method = "survdiff")$pval
    if(log_rank_p_out < 0.0001){
      log_rank_p_out = "<0.0001"
    } else{
      log_rank_p_out = round(log_rank_p_out,4)
    }
    
    
    
    if(is.character(p)){
      hr_out = data.frame("exp(coef)" = NA, "upper .95" = NA, "lower .95" = NA, p_val = NA)
      
    }
    
    
    stable <- data.frame(N = round(km_fit$n,3), Event = round(summary(km_fit)$table[,"events"],3),
                         Median = round(summary(km_fit)$table[,"median"],3),
                         LCL_0.95 = round(summary(km_fit)$table[,"0.95LCL"],3),
                         UCL_0.95 = round(summary(km_fit)$table[,"0.95UCL"],3),
                         Year_0.5 = round(summary(km_fit, times = 6, extend = T)$surv,3),
                         Year_1 = round(summary(km_fit, times = 12, extend = T)$surv,3),
                         Year_2 = round(summary(km_fit, times = 24, extend = T)$surv,3),
                         Year_5 = round(summary(km_fit, times = 5*12, extend = T)$surv,3))
    
    stable = stable[c(grep("myCAF",row.names(stable)),    grep("iCAF",row.names(stable))),]
    stable$log_rank_p = NA
    stable$log_rank_p[1] = log_rank_p_out
    stables  = list()
    stables[[1]] = hr_out
    stables[[2]] = stable
    stables3 = stables
    splots[[1]] <- try(ggsurvplot(font.legend = 16,    
                              font.title = 16,    
                              font.x = 16,         
                              font.y = 16,           
                              font.tickslab = 16,   
                              risk.table.fontsize = 5.5,
                              risk.table.height = 0.35,
                              size = 0.3,
                              censor.size = 3,
                              pval.size = 8,
                              pval.coord =c(pos,0.75),  
                              km_fit, conf.int = F, pval = F,
                              legend.title="",break.time.by = 12,
                              xlab = "Time (months)", risk.table = T,
                               legend.labs = c("iCAF-Basal","iCAF-Classical", "myCAF-Basal","myCAF-Classical"),
                              title = paste(mainLabel," - Elyada PurIST Interaction", sep="")))
    if("try-error" %in% class(splots[[1]])){
      splots[[1]] <- try(ggsurvplot(font.legend = 16,    
                                  font.title = 16,    
                                  font.x = 16,         
                                  font.y = 16,           
                                  font.tickslab = 16,   
                                  risk.table.fontsize = 5.5,
                                  risk.table.height = 0.35,
                                  size = 0.3,
                                  censor.size = 3,
                                  pval.size = 8,
                                  pval.coord =c(pos,0.75),  
                                  km_fit, conf.int = F, pval = F,
                                  legend.title="",break.time.by = 12,
                                  xlab = "Time (months)", risk.table = T,
                                  
                                  title = paste(mainLabel," - Elyada PurIST Interaction", sep="")))
    }
  
  ## Moffit Stroma
    if(individual == F){
      p = try({ coxph(km ~ MS_K2 + PurIST + strata(dataSet), data = survDat) })
    } else {
      p = try({ coxph(km ~ MS_K2 + PurIST , data = survDat) })
    }
    smodels[[2]] <- p
    
    
    if(!"try-error" %in% class(p)){
      bic = round(BIC(p),3)
      hr <- round(1/summary(p)$coefficients[2],3)
    }
    km_fit <- survfit(km ~ MS_K2 + PurIST  , data = survDat, type = "kaplan-meier")
    if(!is.character(p)){
      hr_out = round(1/summary(p)$conf.int,3)[,c(1,4,3)]
      hr_out = as.data.frame(hr_out)
      hr_out$p_val = (summary(p)$coefficients[,5])
      for(qq in 1:nrow(hr_out)){
        if(hr_out$p_val[qq]<0.0001){
          hr_out$p_val[qq]="<0.0001"
        } else {
          hr_out$p_val[qq] = round(summary(p)$coefficients[,5],4)[qq]
        }
      }
    }
    if(!all.equal(surv_pvalue(survfit(km ~ MS_K2 + PurIST, data = survDat), method = "survdiff")$pval,
                  1-pchisq(survdiff(km ~ MS_K2 + PurIST, data = survDat)$chisq, df = (length(survdiff(km ~ MS_K2 + PurIST, data = survDat)$n)-1)))){
      stop("My two log rank p values dont match something is wrong!")
    }
    log_rank_p_out = surv_pvalue(survfit(km ~ MS_K2 + PurIST, data = survDat), method = "survdiff")$pval
    if(log_rank_p_out < 0.0001){
      log_rank_p_out = "<0.0001"
    } else{
      log_rank_p_out = round(log_rank_p_out,4)
    }
    
    
    
    if(is.character(p)){
      hr_out = data.frame("exp(coef)" = NA, "upper .95" = NA, "lower .95" = NA, p_val = NA)
      
    }
    
    
    stable <- data.frame(N = round(km_fit$n,3), Event = round(summary(km_fit)$table[,"events"],3),
                         Median = round(summary(km_fit)$table[,"median"],3),
                         LCL_0.95 = round(summary(km_fit)$table[,"0.95LCL"],3),
                         UCL_0.95 = round(summary(km_fit)$table[,"0.95UCL"],3),
                         Year_0.5 = round(summary(km_fit, times = 6, extend = T)$surv,3),
                         Year_1 = round(summary(km_fit, times = 12, extend = T)$surv,3),
                         Year_2 = round(summary(km_fit, times = 24, extend = T)$surv,3),
                         Year_5 = round(summary(km_fit, times = 5*12, extend = T)$surv,3))
    
    stable$log_rank_p = NA
    stable$log_rank_p[1] = log_rank_p_out
    stables  = list()
    stables[[1]] = hr_out
    stables[[2]] = stable
    stables2 = stables
    splots[[2]] <- try(ggsurvplot(font.legend = 16,             
                              font.title = 16,               
                              font.x = 16,                     
                              font.y = 16,                      
                              font.tickslab = 16,   risk.table.fontsize = 5.5,
                              risk.table.height = 0.35,
                              size = 0.3,
                              censor.size = 3,
                              pval.size = 8,
                              pval.coord =c(pos,0.75),  
                              km_fit, conf.int = F, pval = F,
                              legend.title="",break.time.by = 12,
                              xlab = "Time (months)", risk.table = T,
                          
                              legend.labs = c("Activated-Basal","Activated-Classical", "Normal-Basal","Normal-Classical"),
                          
                              title = paste(mainLabel," - Moffitt Stroma PurIST Interaction",sep="")))
    if("try-error" %in% class(splots[[2]])){
      splots[[2]] = try(ggsurvplot(font.legend = 16,             
                                   font.title = 16,               
                                   font.x = 16,                     
                                   font.y = 16,                      
                                   font.tickslab = 16,   risk.table.fontsize = 5.5,
                                   risk.table.height = 0.35,
                                   size = 0.3,
                                   censor.size = 3,
                                   pval.size = 8,
                                   pval.coord =c(pos,0.75),  
                                   km_fit, conf.int = F, pval = F,
                                   legend.title="",break.time.by = 12,
                                   xlab = "Time (months)", risk.table = T,
                                   
                                   
                                   
                                   title = paste(mainLabel," - Moffitt Stroma PurIST Interaction",sep="")))
    }
    
    ## UNC SCISSORS_CAF K2
    smodels[[3]] <- p
    if(individual == F){
      p = try({ coxph(km ~ SCISSORS_CAF_K2 + PurIST + strata(dataSet), data = survDat) })
    } else {
      p = try({ coxph(km ~ SCISSORS_CAF_K2 + PurIST , data = survDat) })
    }
    
    if(!"try-error" %in% class(p)){
      bic = round(BIC(p),3)
      hr <- round(1/summary(p)$coefficients[2],3)
    }
    km_fit <- survfit(km ~ SCISSORS_CAF_K2 + PurIST , data = survDat, type = "kaplan-meier")
    if(!is.character(p)){
      hr_out = round(1/summary(p)$conf.int,3)[,c(1,4,3)]
      hr_out = as.data.frame(hr_out)
      hr_out$p_val = (summary(p)$coefficients[,5])
      for(qq in 1:nrow(hr_out)){
        if(hr_out$p_val[qq]<0.0001){
          hr_out$p_val[qq]="<0.0001"
        } else {
          hr_out$p_val[qq] = round(summary(p)$coefficients[,5],4)[qq]
        }
      }
    }
    if(!all.equal(surv_pvalue(survfit(km ~ SCISSORS_CAF_K2 + PurIST, data = survDat), method = "survdiff")$pval,
                  1-pchisq(survdiff(km ~ SCISSORS_CAF_K2 + PurIST, data = survDat)$chisq, df = (length(survdiff(km ~ SCISSORS_CAF_K2 + PurIST, data = survDat)$n)-1)))){
      stop("My two log rank p values dont match something is wrong!")
    }
    log_rank_p_out = surv_pvalue(survfit(km ~ SCISSORS_CAF_K2 + PurIST, data = survDat), method = "survdiff")$pval
    if(log_rank_p_out < 0.0001){
      log_rank_p_out = "<0.0001"
    } else{
      log_rank_p_out = round(log_rank_p_out,4)
    }
    
    
    
    if(is.character(p)){
      hr_out = data.frame("exp(coef)" = NA, "upper .95" = NA, "lower .95" = NA, p_val = NA)
      
    }
    
    
    stable <- data.frame(N = round(km_fit$n,3), Event = round(summary(km_fit)$table[,"events"],3),
                         Median = round(summary(km_fit)$table[,"median"],3),
                         LCL_0.95 = round(summary(km_fit)$table[,"0.95LCL"],3),
                         UCL_0.95 = round(summary(km_fit)$table[,"0.95UCL"],3),
                         Year_0.5 = round(summary(km_fit, times = 6, extend = T)$surv,3),
                         Year_1 = round(summary(km_fit, times = 12, extend = T)$surv,3),
                         Year_2 = round(summary(km_fit, times = 24, extend = T)$surv,3),
                         Year_5 = round(summary(km_fit, times = 5*12, extend = T)$surv,3))
    
    stable$log_rank_p = NA
    stable$log_rank_p[1] = log_rank_p_out
    stables  = list()
    stables[[1]] = hr_out
    stables[[2]] = stable
    stables1 = stables
  splots[[3]] <- try(ggsurvplot(font.legend = 16,           
                            font.title = 16,           
                            font.x = 16,                
                            font.y = 16,                 
                            font.tickslab = 16,  risk.table.fontsize = 5.5,
                            risk.table.height = 0.35,
                            size = 0.3,
                            censor.size = 3,
                            pval.size = 8,
                            pval.coord =c(pos,0.75),            
                            km_fit, conf.int = F, pval = F,
                            legend.title="",break.time.by = 12,
                            legend.labs = c("permCAF-Basal","permCAF-Classical","restCAF-Basal","restCAF-Classical"),
                            xlab = "Time (months)", risk.table = T,
                            title = paste(mainLabel," - SCISSORS CAF PurIST Interaction",sep="")))
   if("try-error" %in% class(splots[[3]])){
     splots[[3]] <- try(ggsurvplot(font.legend = 16,           
                                   font.title = 16,           
                                   font.x = 16,                
                                   font.y = 16,                 
                                   font.tickslab = 16,  risk.table.fontsize = 5.5,
                                   risk.table.height = 0.35,
                                   size = 0.3,
                                   censor.size = 3,
                                   pval.size = 8,
                                   pval.coord =c(pos,0.75),            
                                   km_fit, conf.int = F, pval = F,
                                   legend.title="",break.time.by = 12,
                                   xlab = "Time (months)", risk.table = T,
                                   title = paste(mainLabel," - SCISSORS CAF PurIST Interaction",sep="")))
     
   }
  
  # DeCAF
  if(DeCAF_graph == TRUE){
    
    index_caf = which(schemaList == "SCISSORS_CAF_K2_top10")
    if(individual == F){
      p = try({ coxph(km ~ DeCAF + PurIST + strata(dataSet), data = survDat) })
    } else {
      p = try({ coxph(km ~ DeCAF + PurIST , data = survDat) })
    }  
    smodels[[4]] <- p
    if(!"try-error" %in% class(p)){
      bic = round(BIC(p),3)
      hr <- round(1/summary(p)$coefficients[2],3)
    }
      km_fit <- survfit(km ~ DeCAF + PurIST, data = survDat, type = "kaplan-meier")
      if(!is.character(p)){
        hr_out = round(1/summary(p)$conf.int,3)[,c(1,4,3)]
        hr_out = as.data.frame(hr_out)
        hr_out$p_val = (summary(p)$coefficients[,5])
        for(qq in 1:nrow(hr_out)){
          if(as.numeric(hr_out$p_val[qq])<0.0001){
            hr_out$p_val[qq]="<0.0001"
          } else {
            hr_out$p_val[qq] = round(summary(p)$coefficients[,5],4)[qq]
          }
        }
      }
      if(!all.equal(surv_pvalue(survfit(km ~ DeCAF + PurIST, data = survDat), method = "survdiff")$pval,
                    1-pchisq(survdiff(km ~ DeCAF + PurIST, data = survDat)$chisq, df = (length(survdiff(km ~ DeCAF + PurIST, data = survDat)$n)-1)))){
        stop("My two log rank p values dont match something is wrong!")
      }
      log_rank_p_out = surv_pvalue(survfit(km ~ DeCAF + PurIST, data = survDat), method = "survdiff")$pval
      if(log_rank_p_out < 0.0001){
        log_rank_p_out = "<0.0001"
      } else{
        log_rank_p_out = round(log_rank_p_out,4)
      }
      
      
      
      if(is.character(p)){
        hr_out = data.frame("exp(coef)" = NA, "upper .95" = NA, "lower .95" = NA, p_val = NA)
        
      }
     
      
      stable <- data.frame(N = round(km_fit$n,3), Event = round(summary(km_fit)$table[,"events"],3),
                           Median = round(summary(km_fit)$table[,"median"],3),
                           LCL_0.95 = round(summary(km_fit)$table[,"0.95LCL"],3),
                           UCL_0.95 = round(summary(km_fit)$table[,"0.95UCL"],3),
                           Year_0.5 = round(summary(km_fit, times = 6, extend = T)$surv,3),
                           Year_1 = round(summary(km_fit, times = 12, extend = T)$surv,3),
                           Year_2 = round(summary(km_fit, times = 24, extend = T)$surv,3),
                           Year_5 = round(summary(km_fit, times = 5*12, extend = T)$surv,3))
      
      stable$log_rank_p = NA
      stable$log_rank_p[1] = log_rank_p_out
      stables  = list()
      stables[[1]] = hr_out
      stables[[2]] = stable
      stables0 = stables
      
      
      splots[[4]] <- try(ggsurvplot(font.legend = 16,         
                                font.title = 16,              
                                font.x = 16,                     
                                font.y = 16,                        
                                font.tickslab = 16,  risk.table.fontsize = 5.5,
                                risk.table.height = 0.35,
                                size = 0.3,
                                censor.size = 3,
                                pval.size = 8,
                                pval.coord =c(pos,0.75),
                                km_fit, legend.labs = c("permCAF-Basal","permCAF-Classical","restCAF-Basal","restCAF-Classical"),
                                conf.int = F, pval = F  ,
                                legend.title="",break.time.by = 12,
                                 xlab = "Time (months)", risk.table = T,
                                title = paste(mainLabel," - DeCAF PurIST Interaction",sep="")))
      
      if("try-error" %in% class(splots[[4]])){
        splots[[4]] <- try(ggsurvplot(font.legend = 16,         
                                font.title = 16,              
                                font.x = 16,                     
                                font.y = 16,                        
                                font.tickslab = 16,  risk.table.fontsize = 5.5,
                                risk.table.height = 0.35,
                                size = 0.3,
                                censor.size = 3,
                                pval.size = 8,
                                pval.coord =c(pos,0.75),
                                km_fit,
                                conf.int = F, pval = F  ,
                                legend.title="",break.time.by = 12,
                                 xlab = "Time (months)", risk.table = T,
                                title = paste(mainLabel," - DeCAF PurIST Interaction",sep="")))
      
      }
      
      splots[[5]] <- try(ggsurvplot(font.legend = 16,         
                                font.title = 16,              
                                font.x = 16,                     
                                font.y = 16,                        
                                font.tickslab = 16,risk.table.fontsize = 5.5,
                                risk.table.height = 0.35,
                                size = 0.3,
                                censor.size = 3,
                                pval.size = 8,
                                pval.coord =c(pos,0.75), 
                                km_fit, legend.labs = c( "permCAF-Basal","permCAF-Classical","restCAF-Basal","restCAF-Classical"),
                                conf.int = F, pval = F  ,
                                legend.title="DeCAF",break.time.by = 12,
                                xlab = "Time (months)", risk.table = T,
                                title = paste(mainLabel," - DeCAF PurIST Interaction",sep=""),linetype = c("PurIST")))
      
      if("try-error" %in% class(splots[[5]])){
        splots[[5]] <- try(ggsurvplot(font.legend = 16,         
                                      font.title = 16,              
                                      font.x = 16,                     
                                      font.y = 16,                        
                                      font.tickslab = 16,risk.table.fontsize = 5.5,
                                      risk.table.height = 0.35,
                                      size = 0.3,
                                      censor.size = 3,
                                      pval.size = 8,
                                      pval.coord =c(pos,0.75), 
                                      km_fit, 
                                      conf.int = F, pval = F  ,
                                      legend.title="DeCAF",break.time.by = 12,
                                      xlab = "Time (months)", risk.table = T,
                                      title = paste(mainLabel," - DeCAF PurIST Interaction",sep=""),linetype = c("PurIST")))
        
      }
  }
      
    ## SCISSORS_CAF K4
    
      if(individual == F){
        p = try({ coxph(km ~ SCISSORS_CAF_K4 + PurIST + strata(dataSet), data = survDat) })
      } else {
        p = try({ coxph(km ~ SCISSORS_CAF_K4 + PurIST , data = survDat) })
      }
  smodels[[5]] <- p
  if(!"try-error" %in% class(p)){
    bic = round(BIC(p),3)
    hr <- round(1/summary(p)$coefficients[2],3)
  }
      splots[[6]] <- try(ggsurvplot(font.legend = 16,         
                                font.title = 16,          
                                font.x = 16,                 
                                font.y = 16,                       
                                font.tickslab = 16, risk.table.fontsize = 5.5,
                                risk.table.height = 0.35,
                                size = 0.3,
                                censor.size = 3,
                                pval.size = 8,
                                pval.coord =c(pos,0.75), 
                                km_fit, conf.int = F, pval = F,
                                legend.title="",break.time.by = 12,
                                xlab = "Time (months)", risk.table = T,
                                title = paste(mainLabel," - SCISSORS CAF PurIST Interaction",sep="")))
    
  
  

  ## SCISSORS_CAF ALL
  if(individual == F){
    p = try({ coxph(km ~ SCISSORS_CAF_ALL + PurIST + strata(dataSet), data = survDat) })
  } else {
    p = try({ coxph(km ~ SCISSORS_CAF_ALL + PurIST , data = survDat) })
  }
      smodels[[6]] <- p
      if(!"try-error" %in% class(p)){
        bic = round(BIC(p),3)
        hr <- round(1/summary(p)$coefficients[2],3)
      }
  splots[[7]] <- ggsurvplot(font.legend = 16,            
                            font.title = 16,             
                            font.x = 16,                     
                            font.y = 16,                         
                            font.tickslab = 16, risk.table.fontsize = 5.5,
                            risk.table.height = 0.35,
                            size = 0.3,
                            censor.size = 3,
                            pval.size = 8,
                            pval.coord =c(pos,0.75), 
                            km_fit, conf.int = F, pval = F,
                            legend.title="",break.time.by = 12,
                            xlab = "Time (months)", risk.table = T,
                            title = paste(mainLabel," - SCISSORS CAF PurIST Interaction",sep=""))
  




  # Elyada
  if(individual == F){
    p = try({ coxph(km ~ Maurer + PurIST + strata(dataSet), data = survDat) })
  } else {
    p = try({ coxph(km ~ Maurer + PurIST , data = survDat) })
  }
  smodels[[7]] <- p
  
  
  
  
  if(!"try-error" %in% class(p)){
    bic = round(BIC(p),3)
    hr <- round(1/summary(p)$coefficients[2],3)
  }
  km_fit <- survfit(km ~ Maurer+PurIST , data = survDat, type = "kaplan-meier")
  if(!is.character(p)){
    hr_out = round(1/summary(p)$conf.int,3)[,c(1,4,3)]
    hr_out = as.data.frame(hr_out)
    hr_out$p_val = (summary(p)$coefficients[,5])
    for(qq in 1:nrow(hr_out)){
      if(hr_out$p_val[qq]<0.0001){
        hr_out$p_val[qq]="<0.0001"
      } else {
        hr_out$p_val[qq] = round(summary(p)$coefficients[,5],4)[qq]
      }
    }
  }
  if(!all.equal(surv_pvalue(survfit(km ~ Maurer + PurIST, data = survDat), method = "survdiff")$pval,
                1-pchisq(survdiff(km ~ Maurer + PurIST, data = survDat)$chisq, df = (length(survdiff(km ~ Maurer + PurIST, data = survDat)$n)-1)))){
    stop("My two log rank p values dont match something is wrong!")
  }
  log_rank_p_out = surv_pvalue(survfit(km ~ Maurer + PurIST, data = survDat), method = "survdiff")$pval
  if(log_rank_p_out < 0.0001){
    log_rank_p_out = "<0.0001"
  } else{
    log_rank_p_out = round(log_rank_p_out,4)
  }
  
  
  
  if(is.character(p)){
    hr_out = data.frame("exp(coef)" = NA, "upper .95" = NA, "lower .95" = NA, p_val = NA)
    
  }
  
  
  stable <- data.frame(N = round(km_fit$n,3), Event = round(summary(km_fit)$table[,"events"],3),
                       Median = round(summary(km_fit)$table[,"median"],3),
                       LCL_0.95 = round(summary(km_fit)$table[,"0.95LCL"],3),
                       UCL_0.95 = round(summary(km_fit)$table[,"0.95UCL"],3),
                       Year_0.5 = round(summary(km_fit, times = 6, extend = T)$surv,3),
                       Year_1 = round(summary(km_fit, times = 12, extend = T)$surv,3),
                       Year_2 = round(summary(km_fit, times = 24, extend = T)$surv,3),
                       Year_5 = round(summary(km_fit, times = 5*12, extend = T)$surv,3))
  
  stable$log_rank_p = NA
  stable$log_rank_p[1] = log_rank_p_out
  stables  = list()
  stables[[1]] = hr_out
  stables[[2]] = stable
  stables4 = stables
  splots[[8]] <- try(ggsurvplot(font.legend = 16,    
                            font.title = 16,    
                            font.x = 16,         
                            font.y = 16,           
                            font.tickslab = 16,    risk.table.fontsize = 5.5,
                            risk.table.height = 0.35,
                            size = 0.3,
                            censor.size = 3,
                            pval.size = 8,
                            pval.coord =c(pos,0.75), 
                            km_fit, conf.int = F, pval = F,
                            legend.title="",break.time.by = 12,
                            xlab = "Time (months)", risk.table = T,
          
                            legend.labs = c("ECM-rich-Basal","ECM-rich-Classical", "Immune-rich-Basal","Immune-rich-Classical"),
                            title = paste(mainLabel," - Maurer PurIST Interaction", sep="")))
  if("try-error" %in% class(splots[[8]])){
    splots[[8]] <- ggsurvplot(font.legend = 16,    
                              font.title = 16,    
                              font.x = 16,         
                              font.y = 16,           
                              font.tickslab = 16,    risk.table.fontsize = 5.5,
                              risk.table.height = 0.35,
                              size = 0.3,
                              censor.size = 3,
                              pval.size = 8,
                              pval.coord =c(pos,0.75), 
                              km_fit, conf.int = F, pval = F,
                              legend.title="",break.time.by = 12,
                              xlab = "Time (months)", risk.table = T,
                              
                              
                              title = paste(mainLabel," - Maurer PurIST Interaction", sep=""))
  }
  if(save == TRUE){
    if(individual == TRUE){
      filename = paste("data/results/figures/Individual_Survival_Plots_Interaction/", mainLabel,"_Survival.pdf", sep ="")
    } else {
      filename = "data/results/figures/Pooled_Survival_Plots/Pooled_PurIST_Interaction_Survival.pdf"
    }
    pdf(filename)
    arrange_ggsurvplots(splots, 
                        print = TRUE,
                        ncol = 1, nrow = 1)
    dev.off()
    
  } else{
  
    
  }
  s_output = list(splots, smodels, stables0, stables1, stables2, stables3, stables4)
  return(s_output)
  
  
}



######################################################################################
######################### Step 2: Load Combined Data  ################################
######################################################################################

## Set working directory
## setwd(PATH TO DeCAF_Classifier_Clean folder)

## Load Survival
SurvDatCmb <- read.csv("data/clean_data/combined_data/combined_survival_DeCAF.csv")

## Load Classifier to obtain training groups
load("DeCAF_classifier/final_DeCAF_classifier.Rdata")
training = unique(final_DeCAF_classifier$anno$study)

## Load Schema List - Aka Color Codes for other classifiers
load("data/color_codes/cmbsubtypes.RData")

###################################################################################
######################### Step 3: Create Figures  #################################
####################################################################################

## First Save Individual Level Plots
for(i in 1:length(table(SurvDatCmb$dataSet))){
  
  study = names(table(SurvDatCmb$dataSet))[i]
  print(study)
  survDat = SurvDatCmb[SurvDatCmb$dataSet == study,]
  km = Surv(survDat$time, survDat$status)
  ## Plot Individual Survival Figure
  ind_plot = Plot_survival(survDat, km, mainLabel = study, DeCAF_graph = T, save = T,
                individual = T, schemaList = schemaList)
  int_plot = Plot_survival_interaction(survDat, km, mainLabel = study, DeCAF_graph = T, save = T,
                           individual = T, schemaList = schemaList)
  
}




## First Remove Duplicated Patients in "PACA_AU_seq" and "PACA_AU_array" - keep only those in PACA AU seq
idDup <- SurvDatCmb$ID[(duplicated(SurvDatCmb$ID))]
idDup_Paca_array <- which((SurvDatCmb$ID %in% idDup) & (SurvDatCmb$dataSet == "PACA_AU_array"))
SurvDatCmb <- SurvDatCmb[-idDup_Paca_array,]

## Next Generate Figures with Pooled Total Data
survDat = SurvDatCmb

## Compute Kaplan Meier and generate plots
km = Surv(survDat$time, survDat$status)
survival_figures <- Plot_survival(survDat, km, mainLabel = "Full Pooled Survival Data", DeCAF_graph = T, save = FALSE,
                                  individual = FALSE, schemaList = schemaList)


## Save Plots
pdf("data/results/figures/Pooled_Survival_Plots/AllData_MoffitStroma_Dec2023_noYehSeq_noAguirre.pdf", onefile = F)
survival_figures[[1]][2]
dev.off()

pdf("data/results/figures/Pooled_Survival_Plots/AllData_ElyadaCAF_Dec2023_noYehSeq_noAguirre.pdf", onefile = F)
survival_figures[[1]][1]
dev.off()



pdf("data/results/figures/Pooled_Survival_Plots/AllData_ScissorsCAFK2_Dec2023_noYehSeq_noAguirre.pdf", onefile = F)
survival_figures[[1]][3]
dev.off()

pdf("data/results/figures/Pooled_Survival_Plots/AllData_DeCAF_Dec2023_noYehSeq_noAguirre.pdf", onefile = F)
survival_figures[[1]][4]
dev.off()

pdf("data/results/figures/Pooled_Survival_Plots/AllData_ScissorsCAF_ALL_Dec2023_noYehSeq_noAguirre.pdf", onefile = F)
survival_figures[[1]][5]
dev.off()


pdf("data/results/figures/Pooled_Survival_Plots/AllData_Maurer_Dec2023_noYehSeq_noAguirre.pdf", onefile = F)
survival_figures[[1]][7]
dev.off()


##################################

## PurIST x DeCAF Interaction Models on All Data
survDat = SurvDatCmb
km = Surv(survDat$time, survDat$status)

## Generate Validation PurIST x DeCAF Plots
survival_interaction_all = Plot_survival_interaction(survDat, km, mainLabel = "Full Pooled Survival Data", DeCAF_graph = TRUE, save = F,
                                                     individual = FALSE, schemaList = schemaList)


pdf("data/results/figures/Pooled_Survival_Plots/AllData_DeCAF_PurIST_Dec2023_noYehSeq_noAguirre.pdf", onefile = F, width = 8, height = 8)
(survival_interaction_all[[1]][[4]]$plot + guides(colour = guide_legend(ncol = 2)) )/ (survival_interaction_all[[1]][[4]]$table )+
  plot_layout(ncol = 1, heights = c(3, 1))
dev.off()

pdf("data/results/figures/Pooled_Survival_Plots/AllData_DeCAF_PurIST_Dec2023_noYehSeq_noAguirre_COL.pdf", onefile = F, width = 8, height = 8)
(survival_interaction_all[[1]][[5]]$plot +
    scale_color_manual(values = c("violetred1", "violetred1","turquoise4", "turquoise4"), labels = c("permCAF", NA, "restCAF", NA)) +
    guides(colour = guide_legend(ncol = 2))  )/ (survival_interaction_all[[1]][[4]]$table +  scale_color_manual(values = c("turquoise4", "turquoise4", "violetred1", "violetred1")) ) + plot_layout(ncol = 1, heights = c(3, 1)) +
  theme(axis.text.y = ggtext::element_markdown(color = c("turquoise4", "turquoise4","violetred1", "violetred1")))
dev.off()




##################### Forest Plots ##########################################


SurvDatCmb = SurvDatCmb %>% mutate(`strata(dataSet)` = dataSet)
SurvDatCmb$PurIST <- factor(SurvDatCmb$PurIST, levels = c("Classical","Basal-like"))
SurvDatCmb$DeCAF <- factor(SurvDatCmb$DeCAF, levels = c("restCAF","permCAF"))
fit <- coxph(formula = Surv(SurvDatCmb$time, SurvDatCmb$status) ~ PurIST+DeCAF  + strata(dataSet), data = SurvDatCmb)

pdf("data/results/figures/Pooled_Survival_Plots/AllData_DeCAF_PurIST_Dec2023_noYehSeq_noAguirre_Forrest_Plot.pdf", onefile = F, width = 8, height = 8)
ggforest(fit, 
         data = SurvDatCmb,
         fontsize = 1.2,
         cpositions = c(0.0, 0.10, 0.32))
dev.off()


