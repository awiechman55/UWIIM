####Background Plotting Functions##############
####Compare 1 variable across cities 
#v m and cs are indices 
comp_1d <- function(cs,v,ms)
{
  comp_df <- data.frame(City=character(),Var=numeric())
  colnames(comp_df)[2]=p_names$p[v]
  
  for(m in 1:length(ms))
  {
    comp_df[metrics$Metric[ms[m]]] <- numeric()
  }
  
  for(c in 1:length(cs))
  {
    rez_cv <- load_rez(cs[c], v)
    rez_cv <- fixColNames(rez_cv, 1)
    rez_cv$City = cities$City[c]
    colnames(rez_cv)[1] <- p_names$p[v]
    
    if(26 %in% ms)
    {
      rez_cv <- rez_cv %>% mutate(Avg_Water_Restrict = ifelse(Num_Water_Restrict_Events==0,0,Total_Water_Restrict/Num_Water_Restrict_Events))
    }
    
    rez_cv <- rez_cv[,c("City",p_names$p[v],metrics$Metric[ms])]
    
    comp_df <- rbind(comp_df,rez_cv)
  }
  
  comp_df_g <- comp_df %>% gather(Metric, Value, 3:(length(ms)+2))
  
  #return(comp_df_g)
  
  plts_comp <- list()
  
  for(m in 1:length(ms))
  {
    plt_comp_m <- filter(comp_df_g,Metric==metrics$Metric[ms[m]]) %>% 
      ggplot(aes(x=comp_df[,2],y=Value,color=City)) +
      geom_line() + xlab(p_names$p[v]) + ylab(metrics$Metric[ms[m]]) +
      theme_classic() + expand_limits(y=0)
    plts_comp[[m]] <- plt_comp_m
  }
  
  if(length(ms)%%3==0)
    n=length(ms)/3
  else
    n=length(ms)%/%3+1
  
  plt_comp <- do.call("grid.arrange",c(plts_comp,ncol=n))
  
  return(plt_comp)
}









