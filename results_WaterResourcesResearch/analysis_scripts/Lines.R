####Average Sensitivity Deviation Line####
comp_sens_2d_avg <- function(cs, v, def_val=NULL,sepMetric=FALSE, sameFig=FALSE,
                             demRates=FALSE,badShocks=TRUE, badCutoff_PS=-0.5,
                             badCutoff_Q=-0.2)
{
  comp_df_g_avg <- sens_avg_df(cs, v, badShocks=badShocks, 
                               badCutoff_PS=badCutoff_PS, badCutoff_Q=badCutoff_Q)
  
  if(sepMetric)
  {
    if(length(cs)==1 && cs==3){
      plt_rel_min <- filter(comp_df_g_avg,Metric=="Rel_Min" | Metric=="Rel_Avg") %>% 
        ggplot(aes(x=Param2_Value,y=Dev,lty=Metric)) +
        geom_line(color="#7CAE00") +
        theme_classic() + geom_hline(yintercept=0) + xlab(p_labels$p[v]) +
        ylab("Dev. of Avg Sens") + theme(text=element_text(size=14))
    }else{
      plt_rel_min <- filter(comp_df_g_avg,Metric=="Rel_Min" | Metric=="Rel_Avg") %>% 
        ggplot(aes(x=Param2_Value,y=Dev,color=City,lty=Metric)) +
        geom_line() +
        theme_classic() + geom_hline(yintercept=0) + xlab(p_labels$p[v]) +
        ylab("Dev. of Avg Sens") + theme(text=element_text(size=14))
    }
    
    
    if(demRates){
      if(length(cs)==1 && cs==3){
        plt_demRates <- filter(comp_df_g_avg,
                               Metric=="Rates_End" | Metric=="d_End") %>% 
          ggplot(aes(x=Param2_Value,y=Dev,lty=Metric)) +
          geom_line(color="#7CAE00") + 
          theme_classic() +
          geom_hline(yintercept=0) + xlab(p_labels$p[v]) +
          ylab("Dev. of Avg Sens") + theme(text=element_text(size=14))
      }else{
        plt_demRates <- filter(comp_df_g_avg,
                               Metric=="Rates_End" | Metric=="d_End") %>% 
          ggplot(aes(x=Param2_Value,y=Dev,color=City,lty=Metric)) +
          geom_line() + 
          theme_classic() +
          geom_hline(yintercept=0) + xlab(p_labels$p[v]) +
          ylab("Dev. of Avg Sens") + theme(text=element_text(size=14))
      }
    }else{
      if(length(cs)==1 && cs==3){
        plt_fin <- filter(comp_df_g_avg,Metric=="Rates_End") %>% 
          ggplot(aes(x=Param2_Value,y=Dev)) +
          geom_line(color="#7CAE00") + 
          theme_classic() +
          geom_hline(yintercept=0) + xlab(p_labels$p[v]) +
          ylab("Dev. of Avg Sens") + theme(text=element_text(size=14))
        plt_d <- filter(comp_df_g_avg, Metric=="d_End") %>% 
          ggplot(aes(x=Param2_Value,y=Dev)) +
          geom_line(color="#7CAE00") + 
          theme_classic() +
          geom_hline(yintercept=0) + xlab(p_labels$p[v]) +
          ylab("Dev. of Avg Sens") + theme(text=element_text(size=14))
      }else{
        plt_fin <- filter(comp_df_g_avg,Metric=="Rates_End") %>% 
          ggplot(aes(x=Param2_Value,y=Dev,color=City)) +
          geom_line() + 
          theme_classic() +
          geom_hline(yintercept=0) + xlab(p_labels$p[v]) +
          ylab("Dev. of Avg Sens") + theme(text=element_text(size=14))
        plt_d <- filter(comp_df_g_avg, Metric=="d_End") %>% 
          ggplot(aes(x=Param2_Value,y=Dev,color=City)) +
          geom_line() + 
          theme_classic() +
          geom_hline(yintercept=0) + xlab(p_labels$p[v]) +
          ylab("Dev. of Avg Sens") + theme(text=element_text(size=14))
      }
    }
    
    if(v %in% c(7,8,9))
    {
      plt_rel_min <- plt_rel_min + scale_x_continuous(trans='log', breaks=c(4,22,110,220))
      if(demRates){
        plt_demRates <- plt_demRates + scale_x_continuous(trans='log', breaks=c(4,22,110,220))
      }else{
        plt_fin <- plt_fin + scale_x_continuous(trans='log', breaks=c(4,22,110,220))
        plt_d <- plt_d + scale_x_continuous(trans='log', breaks=c(4,22,110,220))
      }
    }
    
    if(!is.null(def_val))
    {
      plt_rel_min <- plt_rel_min + geom_vline(xintercept = def_val, lty="dashed")
      if(demRates){
        plt_demRates <- plt_demRates + geom_vline(xintercept = def_val, lty="dashed")
      }else{
        plt_fin <- plt_fin + geom_vline(xintercept = def_val, lty="dashed")
        plt_d <- plt_d + geom_vline(xintercept = def_val, lty="dashed")
      }
    }
    
    if(demRates){
      return(list(plt_rel_min,plt_demRates))
    }else{
      return(list(plt_rel_min,plt_fin,plt_d))
    }
  }else{
    if(sameFig)
    {
      if(length(cs)==1 && cs==3)
      {
        plt <- filter(comp_df_g_avg,Metric!="Rel_Avg") %>% 
          ggplot(aes(x=Param2_Value,y=Dev, lty=Metric)) +
          geom_line(color="#7CAE00") + 
          theme_classic() +
          geom_hline(yintercept=0,lty="dashed") + xlab(p_labels$p[v]) +
          ylab("Dev. of Avg Sens")
      }else{
        plt <- filter(comp_df_g_avg,Metric!="Rel_Avg") %>% 
          ggplot(aes(x=Param2_Value,y=Dev,color=City, lty=Metric)) +
          geom_line() + 
          theme_classic() +
          geom_hline(yintercept=0,lty="dashed") + xlab(p_labels$p[v]) +
          ylab("Dev. of Avg Sens")
      }
    }else{
      plt <- filter(comp_df_g_avg,Metric!="Rel_Avg") %>% #comp_df_g_avg %>% 
        ggplot(aes(x=Param2_Value,y=Dev,color=City)) +
        geom_line() +facet_grid(rows=vars(Metric),scales="free") + 
        theme_classic() +
        geom_hline(yintercept=0,lty="dashed") + xlab(p_labels$p[v]) +
        ylab("Dev. of Avg Sens")
    }
    
    
    
    if(v %in% c(7,8,9))
    {
      plt <- plt + scale_x_continuous(trans='log', breaks=c(4,22,110,220))
    }
    
    if(!is.null(def_val))
    {
      plt <- plt + geom_vline(xintercept = def_val, lty="dashed")
    }
  }
  
  return(plt)
}

####CAP Magnitude Robustness Line####
comp_sens_1d <- function(cs,v,metricSep=TRUE,noavg=TRUE)
{
  comp_df <- sens_Table(cs,v) 
  
  comp_df_g <- comp_df %>% gather(Metric, Value, 4:7)
  
  if(noavg){
    comp_df_filt <- filter(comp_df_g, Metric != "Rel_Avg")
  }else{
    comp_df_filt<-comp_df_g
  }
  
  if(!metricSep)
  {
    plt_comp <- filter(comp_df_filt) %>% 
      ggplot(aes(x=Param1_Value,y=Value,color=City)) +
      geom_line() + xlab(p_labels$p[v]) + ylab("Sensitivity") +
      theme_classic() + expand_limits(y=0) +
      facet_grid(rows = vars(Metric))
    
    if(v==1)
    {
      plt_comp <- plt_comp+scale_x_reverse(labels = scales::percent)
    }
    
    return(plt_comp)
  }
  
  plts_comp <- list()
  
  plts_comp[[1]] <- filter(comp_df_g,Metric=="d_End") %>% 
    ggplot(aes(x=Param1_Value,y=Value,color=City)) +
    geom_line() + xlab(p_labels$p[v]) + ylab("Demand Sensitivity") +
    theme_classic() + expand_limits(y=0)
  plts_comp[[2]] <- filter(comp_df_g,Metric=="Rates_End") %>% 
    ggplot(aes(x=Param1_Value,y=Value,color=City)) +
    geom_line() + xlab(p_labels$p[v]) + ylab("Financial Sensitivity") +
    theme_classic() + expand_limits(y=0) 
  plts_comp[[3]] <- filter(comp_df_g,Metric=="Rel_Min") %>% 
    ggplot(aes(x=Param1_Value,y=Value,color=City)) +
    geom_line() + xlab(p_labels$p[v]) + ylab("Min. Reliability Sensitivity") +
    theme_classic() + expand_limits(y=0) 
  
  if(!noavg){
    plts_comp[[4]] <- filter(comp_df_g,Metric=="Rel_Avg") %>% 
      ggplot(aes(x=Param1_Value,y=Value,color=City)) +
      geom_line() + xlab(p_labels$p[v]) + ylab("Avg. Reliability Sensitivity") +
      theme_classic() + expand_limits(y=0)
  }
  
  if(v==1)
    for(i in 1:length(plts_comp))
      plts_comp[[i]] <- plts_comp[[i]] + scale_x_reverse(labels = scales::percent)
  
  return(plts_comp)
}

plt_sens_CAP_default <- function(cs,metricSep=TRUE,noavg=TRUE)
{
  plts_CAP <- comp_sens_1d(cs,1,metricSep,noavg)
  
  ##Seperate Plot for QC Rate Sensitivity
  plt_QC_rates_data <- filter(plts_CAP[[2]]$data,City=="QC")
  plt_other_rates_data <- filter(plts_CAP[[2]]$data,City!="QC")
  plt_QC_rates <- plt_QC_rates_data %>% 
    ggplot(aes(x=Param1_Value,y=Value,color=City)) +
    geom_line() + xlab(p_labels$p[1]) + ylab("Financial Sensitivity") +
    scale_color_manual(values="#7CAE00") +
    theme_classic() + expand_limits(y=0)+ scale_x_reverse(labels = scales::percent) 
  plt_other_rates <- plt_other_rates_data %>% 
    ggplot(aes(x=Param1_Value,y=Value,color=City)) +
    geom_line() + xlab(p_labels$p[1]) + ylab("Financial Sensitivity") +
    theme_classic() + expand_limits(y=0) + scale_x_reverse(labels = scales::percent) +
    ylim(-0.1,1)
  
  ##Standardize Limits for Demand and Reliability Sensitivity Plots
  plts_CAP[[1]] <- plts_CAP[[1]] + ylim(-0.1,1)
  plts_CAP[[3]] <- plts_CAP[[3]] + ylim(-0.1,1)
  
  plts_CAP_new <- list(plts_CAP[[1]],plts_CAP[[3]],plt_other_rates,plt_QC_rates)
  
  ##Add Annotations for Representative Shortages
  for(i in 1:length(plts_CAP_new)){
    if(i==4){
      plts_CAP_new[[i]] <- plts_CAP_new[[i]] + geom_vline(xintercept=-0.172,lty="dashed") + 
        annotate(geom="text",angle=90,size=4, x=-0.152, y = 5.5, label="Tier 2a") +
        geom_vline(xintercept=-0.284,lty="dashed") + 
        annotate(geom="text",angle=90,size=4, x=-0.264, y = 5.5, label="Tier 3") +
        geom_vline(xintercept=-0.525,lty="dashed") + 
        annotate(geom="text",angle=90,size=4, x=-0.505, y = 5.5, label="Alt 2") +
        geom_vline(xintercept=-0.96,lty="dashed") + 
        annotate(geom="text",angle=90,size=4, x=-0.94, y = 5.5, label="Alt 1")
    }else{
      plts_CAP_new[[i]] <- plts_CAP_new[[i]] + geom_vline(xintercept=-0.172,lty="dashed") + 
        annotate(geom="text",angle=90,size=4, x=-0.152, y = 0.8, label="Tier 2a") +
        geom_vline(xintercept=-0.284,lty="dashed") + 
        annotate(geom="text",angle=90,size=4, x=-0.264, y = 0.8, label="Tier 3") +
        geom_vline(xintercept=-0.525,lty="dashed") + 
        annotate(geom="text",angle=90,size=4, x=-0.505, y = 0.8, label="Alt 2") +
        geom_vline(xintercept=-0.96,lty="dashed") + 
        annotate(geom="text",angle=90,size=4, x=-0.94, y = 0.8, label="Alt 1")
    }
  }
  
  plt_CAP <- do.call("grid.arrange",c(plts_CAP_new,ncol=2))
  
  return(plt_CAP)
}

