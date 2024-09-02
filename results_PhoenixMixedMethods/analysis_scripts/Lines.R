####Average Sensitivity Deviation Line for 2-D Results####
comp_sens_2d_avg <- function(cs, v, ms=c("Rates_End","Rel_Min"),def_val=NULL,sepMetric=FALSE, sameFig=FALSE,
                             demRates=FALSE,badShocks=TRUE, badCutoff_PS=-0.6,
                             badCutoff_Q=-0.2)
{
  comp_df_g_avg <- sens_avg_df(cs, v, badShocks=badShocks, 
                               badCutoff_PS=badCutoff_PS, badCutoff_Q=badCutoff_Q)
  
  if(sepMetric)
  {
    plt_rel_min <- filter(comp_df_g_avg,Metric=="Rel_Min" | Metric=="Rel_Avg") %>% 
        ggplot(aes(x=Param2_Value,y=Dev,color=Metric)) +
        geom_line(linewidth=0.75) +
        theme_classic(base_size=18,base_family="serif") + geom_hline(yintercept=0) + xlab(p_labels$p[v]) +
        ylab("Dev. of Avg Sens")
    
    
    if(demRates){
        plt_demRates <- filter(comp_df_g_avg,
                               Metric=="Rates_End" | Metric=="d_End") %>% 
          ggplot(aes(x=Param2_Value,y=Dev,color=City,lty=Metric)) +
          geom_line(linewidth=0.75) + 
          theme_classic(base_size=18,base_family="serif") +
          geom_hline(yintercept=0) + xlab(p_labels$p[v]) +
          ylab("Dev. of Avg Sens")
    }else{
        plt_fin <- filter(comp_df_g_avg,Metric=="Rates_End") %>% 
          ggplot(aes(x=Param2_Value,y=Dev,color=City)) +
          geom_line(linewidth=0.75) + 
          theme_classic(base_size=18,base_family = "serif") +
          geom_hline(yintercept=0) + xlab(p_labels$p[v]) +
          ylab("Dev. of Avg Sens") 
        plt_d <- filter(comp_df_g_avg, Metric=="d_End") %>% 
          ggplot(aes(x=Param2_Value,y=Dev,color=City)) +
          geom_line(linewidth=0.75) + 
          theme_classic(base_size=18,base_family = "serif") +
          geom_hline(yintercept=0) + xlab(p_labels$p[v]) +
          ylab("Dev. of Avg Sens") 
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
      plt <- filter(comp_df_g_avg,Metric %in% ms) %>% 
          ggplot(aes(x=Param2_Value,y=Dev,color=Metric)) +
          geom_line(linewidth=0.75) +
          theme_classic(base_size=18,base_family = "serif") +
          geom_hline(yintercept=0,lty="dashed") + xlab(p_labels$p[v]) +
          ylab("Deviation of Average Sensitivity")
    }else{
      plt <- filter(comp_df_g_avg,Metric!="Rel_Avg") %>% 
        ggplot(aes(x=Param2_Value,y=Dev,color=City)) +
        geom_line(linewidth=0.75) + facet_grid(rows=vars(Metric),scales="free") + 
        theme_classic(base_size=14,base_family = "serif") +
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

####Average Sensitivity Deviation Line for 3-D Results####
comp_sens_3d_avg <- function(cs, v1, v2, shock_type, avg_v=NULL, set_v_value=NULL, 
                             ms=c("Rates_End","Rel_Min"),
                             def_val=NULL,badShocks=TRUE, badCutoff_PS=-0.6,
                             avg_all=FALSE,ymin=NA, ymax=NA)
{
  comp_df_g_avg <- sens_avg_df_3d(cs,v1,v2,shock_type,avg_v=avg_v,
                               set_v_value=set_v_value,
                               badShocks=badShocks, badCutoff_PS=badCutoff_PS,
                               avg_all=avg_all)
  
  v<-ifelse(v1==avg_v,v2,v1)
  
  if("Shortage" %in% ms){
    ms_1 <- ms[-which(ms=="Shortage")]
    ms_1 <- append(ms_1,"Rel_Avg")
    df <- filter(comp_df_g_avg,Metric %in% ms_1)
    df <- df %>% mutate(Dev = ifelse(Metric=="Rel_Avg",Dev*50,Dev))
    df <- df %>% mutate(Metric=ifelse(Metric=="Rel_Avg","Shortage",Metric))
    
    df_s <- df[-c(5,6)] %>% spread(key=Metric,value=Dev)
    
    plt <- df_s %>% ggplot(aes(x=Param_Value)) +
      geom_line(aes(y=Rates_End),linewidth=0.75,color="darkgreen") +
      geom_line(aes(y=Shortage/10),linewidth=0.75,color="blue") +
      theme_classic(base_size=18,base_family = "serif") +
      geom_hline(yintercept=0,lty="dashed") + xlab(p_labels$p[v]) +
      scale_y_continuous(name = "Deviation of Average\nSensitivity (Rates)",
                         limits = c(ymin,ymax),
                         sec.axis = sec_axis(~.*10, 
                                             name="Dev. Avg. Sens. (Shortage)")) +
      theme(
        axis.title.y = element_text(color = "darkgreen"),
        axis.title.y.right = element_text(color = "blue")
      ) 
  }else{
    df <- filter(comp_df_g_avg,Metric %in% ms) 
    
    plt <- df %>% 
      ggplot(aes(x=Param_Value,y=Dev,color=Metric)) +
      geom_line(linewidth=0.75) +
      theme_classic(base_size=18,base_family = "serif") +
      geom_hline(yintercept=0,lty="dashed") + xlab(p_labels$p[v]) +
      ylab("Deviation of Average Sensitivity")
  }
  
  
    
    
  if(v %in% c(7,8,9))
  {
    plt <- plt + scale_x_continuous(trans='log', breaks=c(4,22,110,220))
  }
    
  if(!is.null(def_val))
  {
    plt <- plt + geom_vline(xintercept = def_val, lty="dashed")
  }
  
  return(plt)
}



####CAP Magnitude Robustness Line####
comp_sens_1d <- function(cs,v,metrics=c("Infrast_End","d_End","Rates_End","Rel_Min"),metricSep=TRUE)
{
  comp_df <- sens_Table(cs,v) 
  
  comp_df_g <- comp_df %>% gather(Metric, Value, 4:8)
  
  comp_df_filt <- filter(comp_df_g, Metric %in% metrics)
  
  if(!metricSep)
  {
    if(length(cs)==1){
      plt_comp <- filter(comp_df_filt) %>% 
        ggplot(aes(x=Param1_Value,y=Value,color=Metric)) +
        geom_line(linewidth=0.75) + xlab(p_labels$p[v]) + ylab("Sensitivity") +
        theme_classic(base_size=18,base_family = "serif") + expand_limits(y=0)
      
    }else{
      plt_comp <- filter(comp_df_filt) %>% 
        ggplot(aes(x=Param1_Value,y=Value,color=City)) +
        geom_line(linewidth=0.75) + xlab(p_labels$p[v]) + ylab("Sensitivity") +
        theme_classic(base_size=18,base_family = "serif") + expand_limits(y=0) +
        facet_grid(rows = vars(Metric))
    }
    
    if(v==1)
    {
      plt_comp <- plt_comp+scale_x_reverse(labels = scales::percent)
    }
    
    return(plt_comp)
  }
  
  plts_comp <- list()
  
  plts_comp[[1]] <- filter(comp_df_g,Metric=="d_End") %>% 
    ggplot(aes(x=Param1_Value,y=Value,color=City)) +
    geom_line(linewidth=0.75) + xlab(p_labels$p[v]) + ylab("Ending Demand Sensitivity") +
    theme_classic(base_size=18,base_family = "serif") + expand_limits(y=0)
  plts_comp[[2]] <- filter(comp_df_g,Metric=="Rates_End") %>% 
    ggplot(aes(x=Param1_Value,y=Value,color=City)) +
    geom_line(linewidth=0.75) + xlab(p_labels$p[v]) + ylab("Ending Rates Sensitivity") +
    theme_classic(base_size=18,base_family = "serif") + expand_limits(y=0) 
  plts_comp[[3]] <- filter(comp_df_g,Metric=="Rel_Min") %>% 
    ggplot(aes(x=Param1_Value,y=Value,color=City)) +
    geom_line(linewidth=0.75) + xlab(p_labels$p[v]) + ylab("Min. Reliability Sensitivity") +
    theme_classic(base_size=18,base_family = "serif") + expand_limits(y=0)
  plts_comp[[4]] <- filter(comp_df_g,Metric=="Infrast_End") %>%
    ggplot(aes(x=Param1_Value,y=Value,color=City)) +
    geom_line(linewidth=0.75) + xlab(p_labels$p[v]) + ylab("Ending Infrastructure Metric") +
    theme_classic(base_size=18,base_family = "serif") + expand_limits(y=0)
  
  if(!noavg){
    plts_comp[[5]] <- filter(comp_df_g,Metric=="Rel_Avg") %>% 
      ggplot(aes(x=Param1_Value,y=Value,color=City)) +
      geom_line(linewidth=0.75) + xlab(p_labels$p[v]) + ylab("Avg. Reliability Sensitivity") +
      theme_classic(base_size=18,base_family = "serif") + expand_limits(y=0)
  }
  
  if(v==1)
    for(i in 1:length(plts_comp))
      plts_comp[[i]] <- plts_comp[[i]] + scale_x_reverse(labels = scales::percent)
  
  return(plts_comp)
}

plt_sens_CAP_default <- function(cs,metrics=c("Infrast_End","d_End","Rates_End","Rel_Min"),metricSep=TRUE)
{
  plt_CAP <- comp_sens_1d(cs,1,metrics,metricSep) 
    
  if(length(metrics)==2){
    plt_CAP <- plt_CAP + ylim(-0.25,0.5) +
      #theme(legend.position = "top") +
      geom_vline(xintercept=-0.172,lty="dashed") +
      annotate(geom="text",family="serif",angle=90,size=6, 
               x=-0.152, y = 0.4, label="Tier 2a") +
      geom_vline(xintercept=-0.284,lty="dashed") + 
      annotate(geom="text",family="serif",angle=90,size=6, 
               x=-0.264, y = 0.4, label="Tier 3") +
      geom_vline(xintercept=-0.525,lty="dashed") + 
      annotate(geom="text",family="serif",angle=90,size=6, 
               x=-0.522, y = 0.4, label="D-SEIS \n Share") +
      geom_vline(xintercept=-0.96,lty="dashed") + 
      annotate(geom="text",family="serif",angle=90,size=6, 
               x=-0.957, y = 0.3, label="D-SEIS \n Priority")
    
  }else{
    plt_CAP <- plt_CAP + ylim(0,0.3) +
      #theme(legend.position = "top") +
      geom_vline(xintercept=-0.172,lty="dashed") +
      annotate(geom="text",family="serif",angle=90,size=6, 
               x=-0.152, y = 0.25, label="Tier 2a") +
      geom_vline(xintercept=-0.284,lty="dashed") + 
      annotate(geom="text",family="serif",angle=90,size=6, 
               x=-0.264, y = 0.25, label="Tier 3") +
      geom_vline(xintercept=-0.525,lty="dashed") + 
      annotate(geom="text",family="serif",angle=90,size=6, 
               x=-0.522, y = 0.25, label="D-SEIS \n Share") +
      geom_vline(xintercept=-0.96,lty="dashed") + 
      annotate(geom="text",family="serif",angle=90,size=6, 
               x=-0.957, y = 0.25, label="D-SEIS \n Priority")
  }
  return(plt_CAP)
}

