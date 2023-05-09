################## Set up Workspace
setwd("C:/Users/awiec/GitHub_local/UWIIM/results")

library(tidyverse)
library(gridExtra)

#This script loads the meta-data needed to analyze the sensitivity analysis
source("./analysis_scripts/GenGlobalVars.R") 

#This script contains the functions used to load and clean a results file
source("./analysis_scripts/LoadCleanRez.R") 

#This script contains the functions used to calculate the sensitivity metrics 
#and deviations of average sensitivity used to compare results
source("./analysis_scripts/SensMetric.R") 

#This script contains the function used to generate the box-whisker plot 
source("./analysis_scripts/BoxWhisker.R")

#This script contains the function used to generate the 2-d heatmaps
source("./analysis_scripts/Heatmap.R")

#This script contains the functions used to generate sensitivity curves and 
#average deviation curves
source("./analysis_scripts/Lines.R")

#Vector with PMA city indices (1=PHX, 2=Sc, 3=QC)
cs <- c(1,2,3)

####Sensitivity Analysis Fig 1: Default CAP Magnitude Sensitivity####
#Generate First Sensitivity Analysis Figure (Default s variation)
wrr_fig1 <- plt_sens_CAP_default(cs,noavg=TRUE)

####Sensitivity Analysis Fig 2: Average Deviation of IC Variation####
#Function to Generate the Second Sensitivity Analysis Figure (Institutional Costs)
fig2 <- function(badShocks=TRUE, badCutoff_PS=-0.5, badCutoff_Q=-0.2)
{
  ###Institutional Costs
  ##Investment AS
  cs<-c(1,2,3)
  #PHX&Sc
  plt_avg_ic_invest_wrr_PS <- comp_sens_2d_avg(c(1,2), 11, sepMetric=TRUE, 
                                             demRates=TRUE, badShocks=badShocks, 
                                             badCutoff_PS=badCutoff_PS,badCutoff_Q=badCutoff_Q)
  plt_avg_ic_invest_wrr_PS_rm <- plt_avg_ic_invest_wrr_PS[[1]] + ylim(-0.13,0.12) +
    theme(legend.position="top") + guides(lty="none") +
    ggtitle("Reliability Metrics")
  plt_avg_ic_invest_wrr_PS_dr <- plt_avg_ic_invest_wrr_PS[[2]] + ylim(-0.13,0.12) +
    theme(legend.position="top") + guides(lty="none") +
    ggtitle("Adjacent Metrics")

  #QC
  plt_avg_ic_invest_wrr_Q <- comp_sens_2d_avg(3, 11, sepMetric=TRUE, 
                                            demRates=TRUE, badShocks=badShocks, 
                                            badCutoff_PS=badCutoff_PS,badCutoff_Q=badCutoff_Q)
  plt_avg_ic_invest_wrr_Q_rm <- plt_avg_ic_invest_wrr_Q[[1]] + ylim(-2.6,2.4) +
    theme(legend.position="top") +
    ggtitle("") +
    annotate(geom="text",x=0.1,y=1.5,label="QC")
  plt_avg_ic_invest_wrr_Q_dr <- plt_avg_ic_invest_wrr_Q[[2]] + ylim(-2.6,2.4) +
    theme(legend.position="top") +
    ggtitle("")

  ##Rate-setting AS
  #PHX&Sc
  plt_avg_ic_rates_wrr_PS <- comp_sens_2d_avg(c(1,2), 12, sepMetric=TRUE, 
                                            demRates=TRUE, badShocks=TRUE, 
                                            badCutoff_PS=-0.5,badCutoff_Q=-0.2)
  plt_avg_ic_rates_wrr_PS_rm <- plt_avg_ic_rates_wrr_PS[[1]] + ylim(-0.21,0.2) +
    theme(legend.position="top") + guides(lty="none") +
    ggtitle("")
  plt_avg_ic_rates_wrr_PS_dr <- plt_avg_ic_rates_wrr_PS[[2]] + ylim(-0.21,0.2) +
    theme(legend.position="top") + guides(lty="none") +
    ggtitle("")

  #QC
  plt_avg_ic_rates_wrr_Q <- comp_sens_2d_avg(3, 12, sepMetric=TRUE, 
                                           demRates=TRUE, badShocks=TRUE, 
                                           badCutoff_PS=-0.5,badCutoff_Q=-0.2)
  plt_avg_ic_rates_wrr_Q_rm <- plt_avg_ic_rates_wrr_Q[[1]] + ylim(-21,20) +
    theme(legend.position="top")+
    ggtitle("")
  plt_avg_ic_rates_wrr_Q_dr <- plt_avg_ic_rates_wrr_Q[[2]] + ylim(-11,10) +
    theme(legend.position="top") +
    ggtitle("")

  ##All Together
  plt_avg_ic_invest_wrr <- grid.arrange(plt_avg_ic_invest_wrr_PS_rm,
                                      plt_avg_ic_invest_wrr_PS_dr,
                                      plt_avg_ic_invest_wrr_Q_rm,
                                      plt_avg_ic_invest_wrr_Q_dr,
                                      plt_avg_ic_rates_wrr_PS_rm,
                                      plt_avg_ic_rates_wrr_PS_dr,
                                      plt_avg_ic_rates_wrr_Q_rm,
                                      plt_avg_ic_rates_wrr_Q_dr,
                                      layout_matrix=rbind(c(1,3),c(5,7),c(2,4),c(6,8)))
  
  return(plt_avg_ic_invest_wrr)
}

#Generate the Second Sensitivity Analysis Figure (Institutional Costs)
wrr_fig2 <- fig2(badShocks=TRUE,badCutoff_PS=-0.6)

####Sensitivity Analysis Fig 3: Average Deviation of Ambig Variation####
#Function to Generate Third Sensitivity Analysis Figure (Institutional Sensitivity)
fig3 <- function(badShocks=TRUE, badCutoff_PS=-0.6, badCutoff_Q=-0.2,
                 noRel=FALSE, defVal=22, demRates=TRUE, sepQC=TRUE,
                 withHeat=TRUE)
{
  if(withHeat){
    QC_heat_invest_rel <- plot_Sens_heat(3,8,11,metric="Rel_Min",
                                         v1_def=0,v2_def=22,
                                         cmin=-1,cmax=1) + 
      ggtitle("Min Reliability") +
      annotate(geom="text",angle=0,size=4, x=0.1, y = 27, label="Default") +
      annotate(geom="text",angle=0,size=4, x=0.35, y = 8, label="High \n Ambiguity") +
      annotate(geom="text",angle=0,size=4, x=0.35, y = 90, label="Low \n Ambiguity") 
    QC_heat_invest_rates <- plot_Sens_heat(3,8,11,metric="Rates_End",
                                           v1_def=0,v2_def=22,
                                           cmin=-1,cmax=1) +
      ggtitle("Ending Rates") 
      #annotate(geom="text",angle=0,size=4, x=0.1, y = 17, label="Default")
    QC_heat_rates_rel <- plot_Sens_heat(3,9,12,metric="Rel_Min",
                                        v1_def=0,v2_def=22,
                                        cmin=-4,cmax=4) +
      ggtitle("")
      #annotate(geom="text",angle=0,size=4, x=0.18, y = 17, label="Default")
    QC_heat_rates_rates <- plot_Sens_heat(3,9,12,metric="Rates_End",
                                          v1_def=0,v2_def=22,
                                          cmin=-2.0,cmax=2.0) +
      ggtitle("")
      #annotate(geom="text",angle=0,size=4, x=0.18, y = 17, label="Default")
  }
  
  ###Institutional Ambiguity
  ##Investment AS
  cs<-c(1,2,3)
  
  if(sepQC){
    #PHX&Sc
    plt_avg_invest_wrr_PS <- comp_sens_2d_avg(c(1,2), 8, sepMetric=TRUE, def_val=defVal,
                                              demRates=demRates, badShocks=badShocks, 
                                              badCutoff_PS=badCutoff_PS,badCutoff_Q=badCutoff_Q)
    plt_avg_invest_wrr_PS_rm <- plt_avg_invest_wrr_PS[[1]] + ylim(-0.11,0.11) +
      theme(legend.position="top") + guides(lty="none") 
    plt_avg_invest_wrr_PS_dr <- plt_avg_invest_wrr_PS[[2]] + ylim(-0.11,0.11) +
      theme(legend.position="top") + guides(lty="none")
    
    #QC
    plt_avg_invest_wrr_Q <- comp_sens_2d_avg(3, 8, sepMetric=TRUE, def_val=defVal,
                                             demRates=demRates, badShocks=badShocks,
                                             badCutoff_PS=badCutoff_PS,badCutoff_Q=badCutoff_Q)
    plt_avg_invest_wrr_Q_rm <- plt_avg_invest_wrr_Q[[1]] + ylim(-0.11,0.11) +
      theme(legend.position="top") 
    plt_avg_invest_wrr_Q_dr <- plt_avg_invest_wrr_Q[[2]] + ylim(-0.11,0.11) +
      theme(legend.position="top")
    
    ##Rate-setting AS
    #PHX&Sc
    plt_avg_rates_wrr_PS <- comp_sens_2d_avg(c(1,2), 9, sepMetric=TRUE, def_val=defVal,
                                             demRates=demRates, badShocks=TRUE, 
                                             badCutoff_PS=-0.5,badCutoff_Q=-0.2)
    plt_avg_rates_wrr_PS_rm <- plt_avg_rates_wrr_PS[[1]] + ylim(-0.21,0.025) +
      theme(legend.position="top") + guides(lty="none") 
    plt_avg_rates_wrr_PS_dr <- plt_avg_rates_wrr_PS[[2]] + ylim(-0.21,0.025) +
      theme(legend.position="top") + guides(lty="none")
    
    #QC
    plt_avg_rates_wrr_Q <- comp_sens_2d_avg(3, 9, sepMetric=TRUE, def_val=defVal,
                                            demRates=demRates, badShocks=TRUE, 
                                            badCutoff_PS=-0.5,badCutoff_Q=-0.2)
    plt_avg_rates_wrr_Q_rm <- plt_avg_rates_wrr_Q[[1]] + ylim(-2.1,0.25) +
      theme(legend.position="top") 
    plt_avg_rates_wrr_Q_dr <- plt_avg_rates_wrr_Q[[2]] + ylim(-2.1,0.25) +
      theme(legend.position="top")
  }else{
    ##All PMA Cities
    #Invest AS
    plt_avg_invest_wrr <- comp_sens_2d_avg(cs, 8, sepMetric=TRUE, def_val=defVal,
                                              demRates=demRates, badShocks=badShocks, 
                                              badCutoff_PS=badCutoff_PS,badCutoff_Q=badCutoff_Q)
    plt_avg_invest_wrr_rm <- plt_avg_invest_wrr[[1]] + #ylim(-0.11,0.11) +
      theme(legend.position="top") 
    plt_avg_invest_wrr_dr <- plt_avg_invest_wrr[[2]] + #ylim(-0.11,0.11) +
      theme(legend.position="none")  +
      annotate(geom="text",angle=90,size=4, x=19, y = -0.058, label="Default") +
      ggtitle("Ending Rates (Avg Sens)")
    
    #Rate-making AS
    plt_avg_rates_wrr <- comp_sens_2d_avg(cs, 9, sepMetric=TRUE, def_val=defVal,
                                             demRates=demRates, badShocks=TRUE, 
                                             badCutoff_PS=-0.5,badCutoff_Q=-0.2)
    plt_avg_rates_wrr_rm <- plt_avg_rates_wrr[[1]] + #ylim(-0.21,0.025) +
      theme(legend.position="top") + guides(lty="none") 
    plt_avg_rates_wrr_dr <- plt_avg_rates_wrr[[2]] + #ylim(-0.21,0.025) +
      theme(legend.position="top") + guides(lty="none")
  }
  
  ##All Together
  if(noRel){
    if(sepQC){
      plt_avg_invest_wrr <- grid.arrange(#plt_avg_invest_wrr_PS_rm,
        plt_avg_invest_wrr_PS_dr,
        #plt_avg_invest_wrr_Q_rm,
        plt_avg_invest_wrr_Q_dr,
        #plt_avg_rates_wrr_PS_rm,
        plt_avg_rates_wrr_PS_dr,
        #plt_avg_rates_wrr_Q_rm,
        plt_avg_rates_wrr_Q_dr,
        layout_matrix=rbind(c(2,4),
                            c(6,8)))
    }else{
      plt_avg_invest_wrr <- grid.arrange(plt_avg_invest_wrr_dr,
                                         plt_avg_rates_wrr_dr,
                                         QC_heat_invest_rel,
                                         QC_heat_invest_rates,
                                         QC_heat_rates_rel,
                                         QC_heat_rates_rates,
                                         layout_matrix=rbind(c(1,2),
                                                             c(3,5),
                                                             c(4,6)))
    }
  }else{
    if(sepQC){
      plt_avg_invest_wrr <- grid.arrange(plt_avg_invest_wrr_PS_rm,
                                         plt_avg_invest_wrr_PS_dr,
                                         plt_avg_invest_wrr_Q_rm,
                                         plt_avg_invest_wrr_Q_dr,
                                         plt_avg_rates_wrr_PS_rm,
                                         plt_avg_rates_wrr_PS_dr,
                                         plt_avg_rates_wrr_Q_rm,
                                         plt_avg_rates_wrr_Q_dr,
                                         layout_matrix=rbind(c(1,3,5,7),
                                                             c(2,4,6,8)))
    }else{
      plt_avg_invest_wrr <- grid.arrange(plt_avg_invest_wrr_rm,
                                         plt_avg_invest_wrr_dr,
                                         plt_avg_rates_wrr_rm,
                                         plt_avg_rates_wrr_dr,
                                         QC_heat_invest_rel,
                                         QC_heat_invest_rates,
                                         QC_heat_rates_rel,
                                         QC_heat_rates_rates,
                                         layout_matrix=rbind(c(1,3),
                                                             c(2,4)))
    }
    
  }
  
  
  
  return(plt_avg_invest_wrr)
}

#Generate Third Sensitivity Analysis Figure (Institutional Sensitivity)
wrr_fig3 <- fig3(badShocks=TRUE,badCutoff_PS=-0.6, noRel=TRUE, demRates=FALSE,
                 withHeat=TRUE,sepQC=FALSE)

####Sensitivity Analysis Fig 4: Box-Whisker Plot####
#Variables to include in additional sensitivity analysis
vs_addtl <- c(10,4,6,7,13,15,16,17,18,19,20,21,22)

#Generate Fourth Sensitivity Analysis Figure (Box-Whisker Plot of Additional Analysis)
fig4 <- comp_bw(cs, vs_addtl,noavg=FALSE) + ylim(-1.25,3)
fig4
