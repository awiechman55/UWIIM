setwd("C:/Users/awiec/Dropbox/UrbanWaterTransitions/Analysis/Dynamic Systems Model/Outputs/SensAnalysis/SensAnal_PHX_May2024")

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


####Figure 2: Baseline robustness of reliability and rate pressure given#### 
#institutional costs threshold in investment and rate-setting action situations
wow_fig2_invAS_rateMet <- plot_Sens_heat_1constant(1,8,11,1,"IF_sens_l",22,
                                                   metric="Rates_End",
                                                   cmin=-0.25,cmax=0.25,
                                                   v2_def = -0.285,
                                                   v1_def = 1.2,
                                                   IF_thresh=TRUE,goal=1.2) +
  annotate(geom="text",family="serif",angle=0,size=6,x=0.8,y=-0.325,
                       label="Tier 3 Shortage") +
  annotate(geom="text",family="serif",angle=90,size=6,x=1.22,y=-0.75,
           label="No Inst. Costs") +
  annotate(geom="text",family="serif",angle=90,size=6,x=0.92,y=-0.75,
           label="High Inst. Costs") +
  geom_vline(xintercept=0.9,lty="dashed") +
  xlab("Action Threshold (Goal + Inst. Costs)") +
  #ggtitle("Ending Rates") +
  theme(legend.position="none") 
wow_fig2_invAS_relMet <- plot_Sens_heat_1constant(1,8,11,1,"IF_sens_l",22,
                                                   metric="Shortage",
                                                   cmin=-1,cmax=1,
                                                  v2_def=-0.285,
                                                  v1_def = 1.2,
                                                  IF_thresh=TRUE,goal=1.2) +
  ggtitle("Investment Action Situation") +
  #geom_hline(yintercept=-0.6) +
  guides(fill=guide_legend(title=paste("Shortage","\nSensitivity"))) +
  annotate(geom="text",family="serif",angle=0,size=6,x=0.75,y=-0.4,
           label="Rigidity\nTrap") +
  annotate(geom="text",family="serif",angle=0,size=6,x=0.9,y=-0.9,
           label="Shock >> Capacity") +
  xlab("Action Threshold (Goal + Inst. Costs)") +
  geom_vline(xintercept=0.9,lty="dashed") +
  #xlim(0,0.3) +
  theme(legend.position="none") 
wow_fig2_rateAS_rateMet <- plot_Sens_heat_1constant(1,9,12,1,"IF_sens_r",22,
                                                   metric="Rates_End",
                                                   cmin=-0.25,cmax=0.25,
                                                   v2_def = -0.285,
                                                   v1_def = 2,
                                                   IF_thresh=TRUE,goal=2) +
  guides(fill=guide_legend(title=paste("Rate Burden","\nSensitivity"))) +
  geom_vline(xintercept=1.3,lty="dashed") +
  xlab("Action Threshold (Goal + Inst. Costs)") +
  theme(legend.position="none") 
wow_fig2_rateAS_relMet <- plot_Sens_heat_1constant(1,9,12,1,"IF_sens_r",22,
                                                  metric="Shortage",
                                                  cmin=-1,cmax=1,
                                                  v2_def=-0.285,
                                                  v1_def = 2,
                                                  IF_thresh=TRUE,goal=2) +
  ggtitle("Rate-making Action Situation") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.15,y=-0.125,
           label="Rigidity\nTrap") +
  geom_vline(xintercept=1.3,lty="dashed") +
  xlab("Action Threshold (Goal + Inst. Costs)") +
  theme(legend.position="none") 

fig2_colorbar_rate <- get_legend(wow_fig2_rateAS_rateMet + theme(legend.position="right"))
fig2_colorbar_short <- get_legend(wow_fig2_invAS_relMet + theme(legend.position="right"))

fig2_colorbars <- plot_grid(fig2_colorbar_short,fig2_colorbar_rate,ncol=1)

wow_fig2_figs <- plot_grid(wow_fig2_invAS_relMet,wow_fig2_rateAS_relMet,
                           wow_fig2_invAS_rateMet,wow_fig2_rateAS_rateMet,
                      labels=c("A)", "B)", "C)", "D)"),
                      label_x = 0.02, label_y = 0.05,
                      hjust = -0.5, vjust = -0.5,
                      label_fontfamily = "serif",
                      ncol=2, label_size=20)

wow_fig2 <- plot_grid(wow_fig2_figs,fig2_colorbars,rel_widths=c(0.85,0.15))
wow_fig2

#### Figure 3: Effect of Flexibility####
wow_fig3_noIC_invAS_rateMet <- plot_Sens_heat_1constant(1,8,11,1,
                                                         "IF_thresh_l",0,
                                                         metric="Rates_End",
                                                         cmin=-0.3,cmax=0.3,
                                                         v1_def = 22,
                                                         IF_thresh=FALSE,goal=1.2) +
  #annotate(geom="text",family="serif",angle=0,size=6,x=0.4,y=-0.325,
  #         label="Tier 3 Shortage") +
  theme(legend.position="none") 
wow_fig3_noIC_invAS_relMet <- plot_Sens_heat_1constant(1,8,11,1,
                                                        "IF_thresh_l",0,
                                                        metric="Shortage",
                                                        cmin=-1,cmax=1,
                                                        v1_def=22,
                                                        IF_thresh=FALSE,goal=1.2) +
  theme(legend.position="none") +
  ggtitle("No Institutional Costs")
wow_fig3_noIC_rateAS_rateMet <- plot_Sens_heat_1constant(1,9,12,1,
                                                        "IF_thresh_r",0,
                                                        metric="Rates_End",
                                                        cmin=-0.3,cmax=0.3,
                                                        v1_def = 22,
                                                        IF_thresh=FALSE,goal=1.2) +
  #annotate(geom="text",family="serif",angle=0,size=6,x=0.4,y=-0.325,
  #         label="Tier 3 Shortage") +
  theme(legend.position="none") 
wow_fig3_noIC_rateAS_relMet <- plot_Sens_heat_1constant(1,9,12,1,
                                                       "IF_thresh_r",0,
                                                       metric="Shortage",
                                                       cmin=-1,cmax=1,
                                                       v1_def=22,
                                                       IF_thresh=FALSE,goal=1.2) +
  theme(legend.position="none") +
  ggtitle("No Institutional Costs")
wow_fig3_modIC_invAS_rateMet <- plot_Sens_heat_1constant(1,8,11,1,
                                                         "IF_thresh_l",0.15,
                                                   metric="Rates_End",
                                                   cmin=-0.3,cmax=0.3,
                                                   v1_def = 22,
                                                   IF_thresh=FALSE,goal=1.2) +
  #annotate(geom="text",family="serif",angle=0,size=6,x=0.4,y=-0.325,
  #         label="Tier 3 Shortage") +
  theme(legend.position="none") 
wow_fig3_modIC_invAS_relMet <- plot_Sens_heat_1constant(1,8,11,1,
                                                        "IF_thresh_l",0.15,
                                                  metric="Shortage",
                                                  cmin=-1,cmax=1,
                                                  v1_def=22,
                                                  IF_thresh=FALSE,goal=1.2) +
  theme(legend.position="none") +
  ggtitle("Moderate Institutional Costs")
wow_fig3_modIC_rateAS_rateMet <- plot_Sens_heat_1constant(1,9,12,1,
                                                          "IF_thresh_r",0.3,
                                                    metric="Rates_End",
                                                    cmin=-0.3,cmax=0.3,
                                                    v1_def = 22,
                                                    IF_thresh=FALSE,goal=1.2) +
  theme(legend.position="none") 
wow_fig3_modIC_rateAS_relMet <- plot_Sens_heat_1constant(1,9,12,1,
                                                         "IF_thresh_r",0.3,
                                                   metric="Shortage",
                                                   cmin=-1,cmax=1,
                                                   v1_def=22,
                                                   IF_thresh=FALSE,goal=1.2) +
  theme(legend.position="none") +
  ggtitle("Moderate Institutional Costs")

wow_fig3_highIC_invAS_rateMet <- plot_Sens_heat_1constant(1,8,11,1,
                                                          "IF_thresh_l",0.3,
                                                         metric="Rates_End",
                                                         cmin=-0.3,cmax=0.3,
                                                         v1_def = 22,
                                                         IF_thresh=FALSE,goal=1.2) +
  #annotate(geom="text",family="serif",angle=0,size=6,x=0.4,y=-0.325,
  #         label="Tier 3 Shortage") +
  theme(legend.position="none") 
wow_fig3_highIC_invAS_relMet <- plot_Sens_heat_1constant(1,8,11,1,
                                                         "IF_thresh_l",0.3,
                                                        metric="Shortage",
                                                        cmin=-1,cmax=1,
                                                        v1_def=22,
                                                        IF_thresh=FALSE,goal=1.2) +
  guides(fill=guide_legend(title=paste("Shortage","\nSens"))) +
  theme(legend.position="none") +
  ggtitle("High Institutional Costs")
wow_fig3_highIC_rateAS_rateMet <- plot_Sens_heat_1constant(1,9,12,1,
                                                           "IF_thresh_r",0.7,
                                                          metric="Rates_End",
                                                          cmin=-0.3,cmax=0.3,
                                                          v1_def = 22,
                                                          IF_thresh=FALSE,goal=1.2) +
  guides(fill=guide_legend(title=paste("Rate Burden","\nSens"))) +
  theme(legend.position="none") 
wow_fig3_highIC_rateAS_relMet <- plot_Sens_heat_1constant(1,9,12,1,
                                                          "IF_thresh_r",0.7,
                                                         metric="Shortage",
                                                         cmin=-1,cmax=1,
                                                         v1_def=22,
                                                         IF_thresh=FALSE,goal=1.2) +
  theme(legend.position="none") +
  ggtitle("High Institutional Costs")

fig3_colorbar_rate <- get_legend(wow_fig3_highIC_rateAS_rateMet + theme(legend.position="right"))
fig3_colorbar_short <- get_legend(wow_fig3_highIC_invAS_relMet + theme(legend.position="right"))
fig3_colorbars <- plot_grid(fig3_colorbar_short, fig3_colorbar_rate,ncol=1)

wow_fig3_figs <- plot_grid(wow_fig3_noIC_invAS_relMet,wow_fig3_modIC_invAS_relMet,wow_fig3_highIC_invAS_relMet,
                           wow_fig3_noIC_invAS_relMet,wow_fig3_modIC_invAS_rateMet,wow_fig3_highIC_invAS_rateMet,
                           #wow_fig3_tier3_rateAS_relMet,wow_fig3_highIC_rateAS_relMet,
                           #wow_fig3_tier3_rateAS_rateMet,wow_fig3_highIC_rateAS_rateMet,
                           #labels=c("A)", "B)", "C)", "D)","E)", "F)", "G)", "H)"),
                           labels=c("A)", "B)", "C)", "D)","E)", "F)"),
                           label_x = 0.02, label_y = 0.05,
                           hjust = -0.5, vjust = -0.5,
                           label_fontfamily = "serif",
                           ncol=3, label_size=20)

Met_Rel <- ggdraw() + draw_label("Total\nShortage", 
                            fontfamily="serif",x=0.5, size = 20)
Met_Rates <- ggdraw() + draw_label("Ending\nRates", 
                                   fontfamily="serif",x=0.5, size = 20)

Met_labels <- plot_grid(Met_Rel,Met_Rates,ncol=1)

wow_fig3 <- plot_grid(Met_labels,wow_fig3_figs,fig3_colorbars,rel_widths=c(0.1,0.8,0.1),
                      ncol=3)
wow_fig3

##Rate-Setting
wow_fig4_figs <- plot_grid(wow_fig3_noIC_rateAS_relMet,wow_fig3_modIC_rateAS_relMet,wow_fig3_highIC_rateAS_relMet,
                           wow_fig3_noIC_rateAS_relMet,wow_fig3_modIC_rateAS_rateMet,wow_fig3_highIC_rateAS_rateMet,
                           labels=c("A)", "B)", "C)", "D)","E)", "F)"),
                           label_x = 0.02, label_y = 0.05,
                           hjust = -0.5, vjust = -0.5,
                           label_fontfamily = "serif",
                           ncol=3, label_size=20)

wow_fig4 <- plot_grid(Met_labels,wow_fig4_figs,fig3_colorbars,rel_widths=c(0.1,0.8,0.1),
                      ncol=3)
wow_fig4

#######Figure 5: Rigidity Trap############
wow_fig5_highmag_invAS_rel <- plot_Dev_heat(1,8,11,1,metric="Shortage",
                                            cmin=-1,cmax=1,
                                            IF_thresh=TRUE,goal=1.2,
                                            v1_def = 1.2, v2_def = 22) +
  guides(fill=guide_legend(title=paste("Avg. Dev\n(Over Shocks)\nShortage"))) +
  theme(legend.position="none") +
  geom_vline(xintercept=0.9,lty="dashed") +
  xlab("Action Threshold\n(Goal + Inst. Costs)") +
  ylab("Response Elasticity") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.18,y=8,
           label="No Inst. Costs") +
  annotate(geom="text",family="serif",angle=0,size=6,x=0.88,y=8,
           label="High Inst. Costs") +
  annotate(geom="text",family="serif",angle=0,size=6,x=0.8,y=80,
           label="Rigidity\nTrap") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.05,y=8.5,
           label="High Flexibility") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.05,y=70,
           label="Low Flexibility") +
  ggtitle("Investment Action Situation") + coord_flip()
wow_fig5_highmag_invAS_rate <- plot_Dev_heat(1,8,11,1,metric="Rates_End",
                                            cmin=-0.4,cmax=0.4,
                                            IF_thresh=TRUE,goal=1.2,
                                            v1_def = 1.2, v2_def = 22) +
  guides(fill=guide_legend(title=paste("Avg. Dev\n(Over Shocks)\nRate Burden"))) +
  theme(legend.position="none") +
  geom_vline(xintercept=0.9,lty="dashed") +
  xlab("Action Threshold\n(Goal + Inst. Costs)") +
  ylab("Response Elasticity") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.18,y=8,
           label="No Inst. Costs") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.05,y=8.5,
           label="High Flexibility") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.05,y=70,
           label="Low Flexibility") +
  #ggtitle("Investment Action Situation") + 
  coord_flip()
wow_fig5_highmag_rateAS_rel <- plot_Dev_heat(1,9,12,1,metric="Shortage",
                                             cmin=-1,cmax=1,
                                             IF_thresh=TRUE,goal=2,
                                             v1_def = 2, v2_def=22) +
  guides(fill=guide_legend(title=paste("Rel_Avg","\nAvg. Dev\n(Over Shocks)"))) +
  theme(legend.position="none") +
  xlab("Action Threshold\n(Goal + Inst. Costs)") +
  ylab("Response Elasticity") +
  geom_vline(xintercept=1.3,lty="dashed") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.96,y=8,
           label="No Inst. Costs") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.25,y=8,
           label="High Inst. Costs") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.2,y=100,
           label="Rigidity\nTrap") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.65,y=8.5,
           label="High Flexibility") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.65,y=70,
           label="Low Flexibility") +
  ggtitle("Rate-making Action Situation") + coord_flip()

wow_fig5_highmag_rateAS_rate <- plot_Dev_heat(1,9,12,1,metric="Rates_End",
                                             cmin=-0.4,cmax=0.4,
                                             IF_thresh=TRUE,goal=2,
                                             v1_def = 2, v2_def=22) +
  guides(fill=guide_legend(title=paste("Avg. Dev\n(Over Shocks)\nRate Burden"))) +
  theme(legend.position="none") +
  xlab("Action Threshold\n(Goal + Inst. Costs)") +
  ylab("Response Elasticity") +
  geom_vline(xintercept=1.3,lty="dashed") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.96,y=8,
           label="No Inst. Costs") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.25,y=8,
           label="High Inst. Costs") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.2,y=100,
           label="Rigidity\nTrap") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.65,y=8.5,
           label="High Flexibility") +
  annotate(geom="text",family="serif",angle=0,size=6,x=1.65,y=70,
           label="Low Flexibility") +
  coord_flip()

fig5_heat_legend_short <- get_legend(wow_fig5_highmag_invAS_rel + theme(legend.position="right"))
fig5_heat_legend_rate <- get_legend(wow_fig5_highmag_rateAS_rate + theme(legend.position="right"))

wow_fig5_highmag_invAS_line <- comp_sens_3d_avg(1,8,11,1,avg_v=11,
                                                badCutoff_PS = 0,avg_all=TRUE,
                                                def_val=22,ms=c("Rates_End","Shortage"),
                                                ymin = -0.15, ymax=0.15) +
  annotate(geom="text",family="serif",angle=0,size=6,y=0.15,x=8.5,
           label="High Flexibility") +
  annotate(geom="text",family="serif",angle=0,size=6,y=0.15,x=70,
           label="Low Flexibility") +
  annotate(geom="text",family="serif",angle=0,size=6,y=-0.075,x=100,
           label="*Averaged over all shocks\nand action thresholds") +
  #ylab("Deviation of\nAverage Sensitivity") +
  #theme(legend.position="none") +
  xlab("Response Elasticity") +
  ggtitle("Investment Situation") 
  #ylim(-0.15,0.15)
wow_fig5_highmag_rateAS_line <- comp_sens_3d_avg(1,9,12,1,avg_v=12,
                                                badCutoff_PS = 0,avg_all=TRUE,
                                                def_val=22,ms=c("Rates_End","Shortage"),
                                                ymin= -0.15, ymax=0.15) +
  guides(fill=guide_legend(title=paste("Rate Burden","\nAvg. Dev\n(Over Shocks)"))) +
  annotate(geom="text",family="serif",angle=0,size=6,y=0.15,x=8.5,
           label="High Flexibility") +
  annotate(geom="text",family="serif",angle=0,size=6,y=0.15,x=70,
           label="Low Flexibility") +
  #ylab("Deviation of\nAverage Sensitivity") +
  theme(legend.position="none") +
  ggtitle("Rate-making Action Situation") 
  #ylim(-0.15,0.15)

fig5_line_legend <- get_legend(wow_fig5_highmag_invAS_line + theme(legend.position="right"))
fig5_legends <- plot_grid(fig5_heat_legend_short,fig5_heat_legend_rate,ncol=1)

wow_fig5_figs <- plot_grid(wow_fig5_highmag_invAS_rel,wow_fig5_highmag_rateAS_rel,
                      wow_fig5_highmag_invAS_rate,wow_fig5_highmag_rateAS_rate,
                      labels=c("A)", "B)", "C)", "D)"),
                      label_x = 0.02, label_y = 0.05,
                      hjust = -0.5, vjust = -0.5,
                      label_fontfamily = "serif",
                      ncol=2, label_size=20,
                      rel_heights = c(0.5,0.5))
wow_fig5 <- plot_grid(wow_fig5_figs,fig5_legends,ncol=2,rel_widths=c(0.85,0.15))
wow_fig5

#####Fig 6: Flexibility and Over-Response####
wow_fig6_highIC_invAS_rates <- plot_Sens_heat_1constant(1,8,11,1,
                                                        "IF_thresh_l",0.3,
                                                        metric="Rates_End",
                                                        cmin=-0.25,cmax=0.25,
                                                        v1_def = 22,
                                                        v2_def = -0.285,
                                                        IF_thresh=FALSE,goal=1.2) +
  guides(fill=guide_legend(title=paste("Rates_End","\nSens"))) +
  ylab("Magnitude of\nColorado River Shock") +
  xlab("Response Elasticity") +
  annotate(geom="text",family="serif",angle=0,size=6,x=7,y=-0.22,
           label="Low Rate\nSensitivity") +
  annotate(geom="text",family="serif",angle=0,size=6,x=5.5,y=-0.02,
           label="High Rate\nSensitivity") +
  annotate(geom="text",family="serif",angle=0,size=6,x=60,y=-0.75,
           label="High Investment\nNeeded") +
  ggtitle("High Inst. Costs (Supply_Thresh = 0.9)") +
  theme(legend.position="none")

wow_fig6_highIC_rateAS_rates <- plot_Sens_heat_1constant(1,9,12,1,
                                                        "IF_thresh_r",0.7,
                                                        metric="Rates_End",
                                                        cmin=-0.25,cmax=0.25,
                                                        v1_def = 22,
                                                        v2_def= -0.285,
                                                        IF_thresh=FALSE,goal=1.2) +
  guides(fill=guide_legend(title=paste("Rate Burden","\nSensitivity"))) +
  ylab("Magnitude of\nColorado River Shock") +
  xlab("Response Elasticity") +
  annotate(geom="text",family="serif",angle=0,size=6,y=-0.32,x=150,
           label="Tier 3") +
  #geom_hline(yintercept=-0.525,lty="dashed") + 
  #annotate(geom="text",family="serif",angle=0,size=6, 
  #         y=-0.522, x = 150, label="D-SEIS \n Share") +
  #geom_hline(yintercept=-0.96,lty="dashed") + 
  #annotate(geom="text",family="serif",angle=0,size=6, 
  #         y=-0.957, x = 150, label="D-SEIS \n Priority") +
  #annotate(geom="text",family="serif",angle=0,size=6,y=-0.75,x=6,
  #         label="Flexibility\nTrap") +
  annotate(geom="text",family="serif",angle=0,size=6,y=-0.075,x=70,
           label="Rigidity\nTrap") +
  ggtitle("High Inst. Costs (DSCR_Thresh=1.3)") +
  theme(legend.position="none")

fig6_legend <- get_legend(wow_fig6_highIC_rateAS_rates + theme(legend.position="right"))

fig6_legends <- plot_grid(fig5_line_legend,fig6_legend,ncol=1)

fig6 <- plot_grid(wow_fig5_highmag_invAS_line,wow_fig5_highmag_rateAS_line,
                  wow_fig6_highIC_invAS_rates,wow_fig6_highIC_rateAS_rates,
                  labels=c("A)", "B)","C)","D)"),
                  label_x = 0.02, label_y = 0.05,
                  hjust = -0.5, vjust = -0.5,
                  label_fontfamily = "serif",
                  label_size=20,
                  ncol=2,rel_heights=c(0.4,0.6))
fig6 <- plot_grid(fig6,fig6_legends,ncol=2,rel_widths=c(0.9,0.1))
fig6
