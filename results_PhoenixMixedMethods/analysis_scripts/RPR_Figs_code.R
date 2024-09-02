rpr_fig2_line <- function(v1,v2,ms,set_v,set_v_values){
  cs <- c(1)
  
  rpr_fig2_maxflex <- comp_sens_3d_avg(cs,v1,v2,2,avg_v=set_v,set_v_value=set_v_values[1]) +
    ylim(-0.15,0.15) + ggtitle("Maximum Flexibility") + 
    theme(legend.position = "none",axis.title.y = element_blank())
  rpr_fig2_modflex <- comp_sens_3d_avg(cs,v1,v2,2,avg_v=set_v,set_v_value=set_v_values[2]) +
    ylim(-0.15,0.15) + ggtitle("Moderate Flexibility") + 
    theme(legend.position = "none",axis.title.y = element_blank())
  rpr_fig2_minflex <- comp_sens_3d_avg(cs,v1,v2,2,avg_v=set_v,set_v_value=set_v_values[3]) +
    ylim(-0.15,0.15) + ggtitle("Minimum Flexibility") + 
    theme(legend.position = "none") 
  
  rpr_fig2_legend <- get_legend(rpr_fig2_maxflex + theme(legend.position = "top"))
  
  rpr_fig2_plots <- plot_grid(rpr_fig2_minflex, rpr_fig2_modflex,rpr_fig2_maxflex,ncol=3)
  rpr_fig2 <- plot_grid(rpr_fig2_legend,rpr_fig2_plots,ncol=1,
                        rel_heights = c(0.15,0.85))
  
  return(rpr_fig2)
}

rpr_fig2 <- function(v1,v2,ms,set_v,set_v_values){
  cs <- c(1)
  n_i <- length(set_v_values)
  
  if(set_v==8)
    set_v_name = "IF_sens_l"
  else
    set_v_name = "IF_sens_r"
  
  figs_m=list()
  titles = list()
  
  for(m in 1:length(ms)){
    for(i in 1:n_i){
      #Create Metric Labels
      figs_m[[(m-1)*(n_i+1)+1]] <- ggdraw() + 
        draw_label(ms[m], fontfamily="serif",
                   x=0.5, size = 20)
      
      #Set Title Based on Flexibility Setting
      if(i==1){
        title <- paste("Max Flex")
      }else if(i==2){
        title <- paste("Mod Flex")
      }else{
        title <- paste("Min Flex")
      }
      #Create Plot
      plt_m_i <- plot_Sens_heat_1constant(cs, v1, v2, 1, set_v_name, set_v_values[i], 
                                           metric=ms[m], cmin=-0.5,cmax=0.5) +
        ggtitle(title) + theme(legend.position = "none")
      #Store Plot
      figs_m[[(m-1)*(n_i+1)+i+1]] <- plt_m_i
    }
  }
  
  #Get Legend
  legend <- get_legend(figs_m[[2]]+
                         theme(legend.position = "top",
                               legend.key.width = unit(2.5, "cm")))
  
  #Create Figure Title
  AS <- ggdraw() + draw_label(ifelse(v1==8, "Investment \n Action Situation", 
               "Rate-Making \n Action Situation"), 
               fontfamily="serif",x=0.5, size = 20)
  title <- plot_grid(AS,legend,rel_widths=c(0.33,0.66))
  
  #Arrange Plots
  fig2_plots <- plot_grid(plotlist=figs_m, ncol=n_i+1,
                          rel_widths = c(0.16,0.28,0.28,0.28))
  
  #Arrange Figure
  fig2 <- plot_grid(title,fig2_plots,ncol=1,
                    rel_heights=c(0.15,0.85))
  
  
  return(fig2)
}
