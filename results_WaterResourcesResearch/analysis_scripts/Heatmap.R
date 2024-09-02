####Sensitivity Heat Map (Raw)####
plot_Sens_heat <- function(cs,v1,v2,metric=NA,noavg=TRUE,v1_def=NULL,v2_def=NULL,
                           cmin=NA,cmax=NA)
{
  df_heat <- sens_Table(cs,v1,v2)
  
  df_heat_g <- df_heat %>% gather(Metric,Sens, 6:9)
  
  if(noavg)
  {
    df_heat_g <- filter(df_heat_g, Metric != "Rel_Avg")
  }
  
  if(is.na(metric)){
    plt <- df_heat_g %>% ggplot(aes(y=Param1_Value,x=Param2_Value,fill=Sens)) +
      geom_tile() + 
      theme_classic(base_size=14) + xlab(p_labels$p[v2]) + ylab(p_labels$p[v1]) +
      facet_grid(rows = vars(City), cols=vars(Metric))
  }else{
    df_heat_g <- filter(df_heat_g,Metric==metric)
    
    plt <- df_heat_g %>% 
      ggplot(aes(y=Param1_Value,x=Param2_Value,fill=Sens)) +
      geom_tile() + 
      theme_classic(base_size=14) + xlab(p_labels$p[v2]) + ylab(p_labels$p[v1]) +
      facet_grid(rows = vars(City))
  }
  if(is.na(cmin)){
    plt <- plt + scale_fill_distiller(type="seq",direction=-1,palette="RdBu",
                         limits=c(min(df_heat_g$Sens),max(df_heat_g$Sens)))
  }else{
    plt <- plt + scale_fill_distiller(type="seq",direction=-1,palette="RdBu",
                                      limits=c(cmin,cmax))
  }
  
  
  
  if(v1 %in% c(7,8,9))
  {
    plt <- plt + scale_y_continuous(trans='log', breaks=c(4,22,110,220)) 
  }
  else if(v2 %in% c(7,8,9))
  {
    plt <- plt + scale_x_continuous(trans='log', breaks=c(4,22,110,220))
  }
  if(v1 == 1)
  {
    plt <- plt + scale_y_reverse()
  }
  if(!is.null(v1_def)){
    plt <- plt+geom_vline(xintercept=v1_def,lty="dashed")
  }
  if(!is.null(v2_def)){
    plt <- plt+geom_hline(yintercept=v2_def,lty="dashed")
  }
  
  return(plt)
}