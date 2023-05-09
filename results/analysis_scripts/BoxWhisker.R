####Compare OAT Sensitivity Analysis with Box-Whiskers Plot of Avg Sensitivity Deviations####
comp_bw <- function(cs, vs, noavg=FALSE)
{
  comp_df_g_avg <- sens_avg_df(cs,vs)
  
  p_labels_temp <- p_labels
  p_labels_temp$p_name = p_names$p
  colnames(p_labels_temp) <- c("Param2_label","Param2_Name")
  comp_df_g_avg <- left_join(comp_df_g_avg, p_labels_temp, by="Param2_Name")
  
  
  if(length(cs)==1){
    plt_bw <- filter(comp_df_g_avg,Metric!="Rel_Avg") %>% ggplot(aes(x=Param2_label,y=Dev)) +
      geom_boxplot(fill="grey")
  }else{
    if(noavg){
      plt_bw <- comp_df_g_avg %>% 
        ggplot(aes(x=Param2_label,y=Dev,fill=City)) + geom_boxplot()
    }else{
      plt_bw <- filter(comp_df_g_avg,Metric!="Rel_Avg") %>% 
        ggplot(aes(x=Param2_label,y=Dev,fill=City)) + geom_boxplot()
    }
  }
  
  plt_bw <- plt_bw + 
    theme_classic() + xlab(NULL) + ylab("Deviation of Average Sensitivity from Default") +
    scale_x_discrete(guide = guide_axis(angle = 70)) +
    facet_grid(rows=vars(Metric)) + theme(text=element_text(size=16))
  
  return(plt_bw)
}
