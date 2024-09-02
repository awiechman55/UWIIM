####Old RPR_Figs
rpr_fig2 <- comp_sens_2d_avg(1, 12, badCutoff_PS=0, sameFig=TRUE)

rpr_fig2 <- rpr_fig2 + annotate(geom="text",family="serif",angle=0,size=6, 
                                x=0.125, y = -0.05, 
                                label="Less Rate \n Sensitivity \n (Reduced Response)") +
  annotate(geom="text",family="serif",angle=0,size=6,x=0.35,y=0.06,
           label="More Rate Sensitivity (Over Response)") +
  annotate(geom="text",family="serif",angle=0,size=6,x=0.8,y=-0.04,
           label="No Response (Rigidity Trap)") +
  xlab(expression("Institutional Dependencies/Costs in Rate-Setting (\u03F5"[3]*")"))

rpr_fig2

rpr_fig2_rel <- comp_sens_2d_avg(1, 11, badCutoff_PS=0, sameFig=TRUE)

#rpr_fig2_rel <- rpr_fig2_rel + annotate(geom="text",family="serif",angle=0,size=6, 
#                                x=0.125, y = -0.05, 
#                                label="Less Rate \n Sensitivity \n (Reduced Response)") +
#  annotate(geom="text",family="serif",angle=0,size=6,x=0.35,y=0.06,
##           label="More Rate Sensitivity (Over Response)") +#
#  annotate(geom="text",family="serif",angle=0,size=6,x=0.8,y=-0.04,
#           label="No Response (Rigidity Trap)") +
#  xlab(expression("Institutional Dependencies/Costs in Rate-Setting (\u03F5"[3]*")"))

rpr_fig2_rel

####Fig 3: (H2) Institutional Voids & Rates####
rpr_fig3_invest <- comp_sens_2d_avg(1, 8, def_val = 22, badCutoff_PS=0, sameFig=TRUE) + 
  theme(legend.position = "none") + ylim(-0.1276,0.06) +
  annotate(geom="text",family="serif",angle=0,size=6, x=10, y = -0.04, 
           label="Balanced Response") +
  xlab(expression("Response Elasticity in `Supply Sufficiency`"))
rpr_fig3_rate <- comp_sens_2d_avg(1, 9, def_val = 22, badCutoff_PS=0, sameFig=TRUE) + 
  theme(legend.position = "none") + ylim(-0.1276,0.06) +
  annotate(geom="text",family="serif",angle=0,size=6, x=4.2, y = -0.03, 
           label="Balanced \n Response") +
  annotate(geom="text",family="serif",angle=0,size=6, x=15, y = 0.04, 
           label="Over Response") +
  annotate(geom="text",family="serif",angle=0,size=6, x=100, y = -0.04, 
           label="Balanced Response") +
  xlab(expression("Response Elasticity in Rate-Making"))

rpr_fig3_legend <- get_legend(rpr_fig3_invest + theme(legend.position="top"))

rpr_fig3 <- plot_grid(rpr_fig3_legend, rpr_fig3_invest, rpr_fig3_rate,
                      ncol=1, rel_heights = c(0.1,0.45,0.45))

rpr_fig3

####Fig 4: (H3) Institutional Costs & Ambiguity Deviations####
rpr_fig4_rates_rate <- plot_Sens_heat_highmag(1,9,12,metric="Rates_End",v2_def=22,
                                              cmin=-0.3,cmax=0.3) +
  ggtitle("Ending Rates") +
  theme(legend.position="none") + 
  annotate(geom="text",family="serif",angle=90,size=6,x=0.32,y=50,
           label="Over Response") + 
  annotate(geom="text",family="serif",angle=90,size=6,x=0.1,y=55,
           label="Balanced Response") +
  annotate(geom="text",family="serif",angle=0,size=6,x=0.5,y=5.5,
           label="Balanced Response") 
rpr_fig4_rates_rel <- plot_Sens_heat_highmag(1,9,12,metric="Rel_Min",v2_def=22,
                                             cmin=-0.3, cmax = 0.3) +
  ggtitle("Minimum Reliability") +
  theme(legend.position="none") + 
  annotate(geom="text",family="serif",angle=90,size=6,x=0.8,y=50,
           label="Rigidity Trap")
rpr_rates_legend <- get_legend(rpr_fig4_rates_rate + theme(legend.position="right"))


rpr_fig4_rates <- plot_grid(rpr_fig4_rates_rel,rpr_fig4_rates_rate,
                            labels=c("A","B"),
                            label_fontfamily = "serif",
                            label_y = 0.075,ncol=2, label_x = 0.02,
                            label_size = 18)
rpr_fig4_rates <- plot_grid(rpr_fig4_rates,rpr_rates_legend,ncol=2,
                            rel_widths = c(0.8,0.2)) 

rpr_fig4_rates

####Fig 5: Institutional Costs & Ambiguity with Pace Change####
#No Inst Costs
rpr_fig5_rates_rel_0 <- plot_Sens_heat_1constant(1,9, 12, 3, "IF_thresh_r",0, 
                                                 metric="Rel_Min",v2_def=22,
                                                 cmin=-0.185,cmax=0.185) +
  ggtitle("Minimum Reliability") + theme(legend.position = "none") +
  annotate(geom="text",family="serif",angle=0,size=6,x=15,y=50,
           label="Inst. Costs = 0") 

rpr_fig5_legend <- get_legend(rpr_fig5_rates_rel + theme(legend.position="right"))

rpr_fig5_rates_rate_0 <- plot_Sens_heat_1constant(1,9, 12, 3, "IF_thresh_r",0, 
                                                  metric="Rates_End",v2_def=22, 
                                                  cmin=-0.185,cmax=0.185) +
  ggtitle("Ending Rates") + theme(legend.position = "none") +
  annotate(geom="text",family="serif",angle=0,size=6,x=15,y=10,
           label="Balance Response") +
  annotate(geom="text",family="serif",angle=0,size=6,x=2.5,y=10,
           label="Over \n Response") 

#Moderate Inst Costs (0.34)
rpr_fig5_rates_rel_m <- plot_Sens_heat_1constant(1,9, 12, 3, "IF_thresh_r",0.34, 
                                                 metric="Rel_Min",v2_def=22,
                                                 cmin=-0.185,cmax=0.185) +
  theme(legend.position = "none") +
  annotate(geom="text",family="serif",angle=0,size=6,x=15,y=50,
           label="Inst. Costs = 0.34") 

rpr_fig5_rates_rate_m <- plot_Sens_heat_1constant(1,9, 12, 3, "IF_thresh_r",0.34, 
                                                  metric="Rates_End",v2_def=22, 
                                                  cmin=-0.185,cmax=0.185) +
  theme(legend.position = "none") +
  annotate(geom="text",family="serif",angle=0,size=6,x=2.5,y=70,
           label="Over \n Response") 

#High Inst Costs (0.74)
rpr_fig5_rates_rel_h <- plot_Sens_heat_1constant(1,9, 12, 3, "IF_thresh_r",0.74, 
                                                 metric="Rel_Min",v2_def=22,
                                                 cmin=-0.185,cmax=0.185) +
  theme(legend.position = "none") +
  annotate(geom="text",family="serif",angle=0,size=6,x=15,y=50,
           label="Inst. Costs = 0.74") 

rpr_fig5_rates_rate_h <- plot_Sens_heat_1constant(1,9, 12, 3, "IF_thresh_r",0.74, 
                                                  metric="Rates_End",v2_def=22, 
                                                  cmin=-0.185,cmax=0.185) +
  theme(legend.position = "none")


rpr_fig5_rates <- plot_grid(rpr_fig5_rates_rel_0,rpr_fig5_rates_rate_0,
                            rpr_fig5_rates_rel_m,rpr_fig5_rates_rate_m,
                            rpr_fig5_rates_rel_h,rpr_fig5_rates_rate_h,
                            labels=c("A","B","C","D","E","F"),
                            label_fontfamily = "serif",
                            label_y = 0.15,ncol=2, label_x = 0.02,
                            label_size = 18)
rpr_fig5_rates <- plot_grid(rpr_fig5_rates,rpr_fig5_legend,ncol=2,
                            rel_widths = c(0.8,0.2))

rpr_fig5_rates

#Investment Ambiguity
rpr_fig5_invest_rel <- plot_Sens_heat_1constant(1,8, 12, 3, "IF_thresh_r",0, 
                                                metric="Rates_End",v2_def=22,
                                                cmin=-0.17,cmax=0.17) +
  ggtitle("Minimum Reliability")

rpr_fig5_invest_rel
