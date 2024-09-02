setwd("C:/Users/awiec/Dropbox/UrbanWaterTransitions/Analysis/Dynamic Systems Model/Outputs/SensAnalysis/SensAnal_PHX_May2023")

library(tidyverse)
library(gridExtra)
library(cowplot)
library(scales)
library(ggnewscale)

####Source Code####

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

#This script contains the functions used to generate the figures for the RPR paper
source("./analysis_scripts/RPR_Figs_code.R")

#Vector with PMA city indices (1=PHX)
cs <- c(1)

####Fig 1: Default Phoenix CAP Magnitude Sensitivity####
#Generate Magnitude-Based Figure
rpr_fig1_mag <- plt_sens_CAP_default(cs,metrics=c("Rel_Min","Rates_End"),metricSep=FALSE)

#Generate Gradual Heatmaps
rpr_fig1_heat_rel <- plot_Sens_heat(1,2,1,metric="Rel_Min",cmin=-0.17,cmax=0.17)
legend_fig1_heat <- get_legend(rpr_fig1_heat_rel+theme(legend.key.height = unit(1, "cm")))
rpr_fig1_heat_rel <- rpr_fig1_heat_rel + theme(legend.position="none")+
  ggtitle("Minimum Reliability")

#rpr_fig1_heat_infrast <- plot_Sens_heat(1,2,1,metric="Infrast_End",cmin=-0.01,cmax=0.22) +
#  theme(legend.position="none") + ggtitle("Ending Infrastructure")
rpr_fig1_heat_rates <- plot_Sens_heat(1,2,1,metric="Rates_End",cmin=-0.17,cmax=0.17) +
  theme(legend.position="none") + ggtitle("Ending Rates")
#rpr_fig1_heat_d <- plot_Sens_heat(1,2,1,metric="d_End",cmin=-0.01,cmax=0.22) +
#  theme(legend.position="none") + ggtitle("Ending Demand")

rpr_fig1_heat <- plot_grid(rpr_fig1_heat_rel,rpr_fig1_heat_rates,
                           labels=c("B","C"),
                           label_fontfamily = "serif",
                           label_y = 0.125, ncol=2,label_x=0.02,label_size=18)

rpr_fig1_heat <- plot_grid(rpr_fig1_heat, legend_fig1_heat,
                           ncol=2,rel_widths = c(0.8,0.2))

#Compile into Figure
rpr_fig1 <- plot_grid(rpr_fig1_mag,rpr_fig1_heat,
                      labels=c("A",""), 
                      label_fontfamily = "serif",
                      label_y = 0.125,ncol=1, label_x = 0.02,
                      label_size = 18, rel_heights = c(0.5,0.5))

rpr_fig1

####Fig 2: (H1) Institutional Costs & Rates####
#Get flexibility value close to 110
sens3d_test<-sens_Table(cs,9,12,three_d=TRUE,shock_type=2)

#Rate-Making Action Situation (Rates_End and Rel_Min)
fig2_rates <- rpr_fig2(9,12,9,c(4,22,sens3d_test$Param1_Value[106254]),
                       ms=c("Rates_End","Infrast_End","Rel_Avg"))
fig2_rates

#Rate-Making Action Situation (Rates_End and Rel_Min)
fig2_invest <- rpr_fig2(8,11,8,c(4,22,sens3d_test$Param1_Value[106254]),
                        ms=c("Rates_End","Infrast_End","Rel_Min"))
fig2_invest

####Fig 3: (H2) Flex vs. Rigid Heat & Marg Effects####
##Heatmaps
rpr_fig3_heat_rates_rate <- plot_Sens_heat_highmag(1,9,12,metric="Rates_End",v2_def=22,
                                              cmin=-0.13,cmax=0.13) +
  ggtitle("Ending Rates") +
  theme(legend.position="none") + 
  annotate(geom="text",family="serif",angle=90,size=6,x=0.32,y=50,
           label="Over Response") + 
  annotate(geom="text",family="serif",angle=90,size=6,x=0.1,y=55,
           label="Balanced Response") +
  annotate(geom="text",family="serif",angle=0,size=6,x=0.5,y=5.5,
           label="Balanced Response") 
rpr_fig3_heat_rates_rel <- plot_Sens_heat_highmag(1,9,12,metric="Rel_Min",v2_def=22,
                                             cmin=-0.13, cmax = 0.13) +
  ggtitle("Minimum Reliability") +
  theme(legend.position="none") + 
  annotate(geom="text",family="serif",angle=90,size=6,x=0.8,y=50,
           label="Rigidity Trap")
rpr_fig3_legend <- get_legend(rpr_fig3_heat_rates_rate + theme(legend.position="right"))

rpr_fig3_heat_invest_rates <- plot_Sens_heat_highmag(1,8,11,metric="Rates_End",v2_def=22,
                                                   cmin=-0.13,cmax=0.13) +
  ggtitle("Ending Rates") +
  theme(legend.position="none") + 
  annotate(geom="text",family="serif",angle=0,size=6,x=0.1,y=80,
           label="Over \n Response") + 
  annotate(geom="text",family="serif",angle=0,size=6,x=0.4,y=7,
           label="Balanced \n Response") 
rpr_fig3_heat_invest_rel <- plot_Sens_heat_highmag(1,8,11,metric="Rel_Min",v2_def=22,
                                                  cmin=-0.13, cmax = 0.13) +
  ggtitle("Minimum Reliability") +
  theme(legend.position="none") + 
  annotate(geom="text",family="serif",angle=0,size=6,x=0.4,y=120,
           label="Rigidity \n Trap")


rpr_fig3_heat <- plot_grid(rpr_fig3_heat_rates_rel,rpr_fig3_heat_rates_rate,
                                 rpr_fig3_heat_invest_rel,rpr_fig3_heat_invest_rates,
                            labels=c("A","B","C","D"),
                            label_fontfamily = "serif",
                            label_y = 0.075,ncol=2, label_x = 0.02,
                            label_size = 18)
rpr_fig3_heat <- plot_grid(rpr_fig3_heat,rpr_fig3_legend,ncol=2,
                            rel_widths = c(0.85,0.15)) 
rpr_fig3_heat

##Marginal Effect Figures


#### (H3) Flexibility Trap with High Mag####
rpr_fig4_rates_rate <- plot_Sens_heat_1constant(cs,9, 12, 2, "IF_thresh_r", 0, 
                                                metric="Rates_End",#v2_def=22,
                                                cmin=-0.22,cmax=0.22) +
  ggtitle("Ending Rates")
rpr_fig4_invest_rel <- plot_Sens_heat_1constant(cs,7, 10, 2, "IF_thresh_s", 0.5, 
                                                metric="Rel_Avg") +
  ggtitle("Minimum Reliability")
rpr_fig4_invest_rel
rpr_fig4_rates_rate
