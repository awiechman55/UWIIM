####Necessary global variables#######################
p_names <- read.csv("p_names.csv")
cities <- read.csv("cities.csv")
metrics <- read.csv("metrics.csv")
p_labels <- p_names

#Fix p labels
p_labels$p[1]<-"Magnitude of CAP Shock"
p_labels$p[2]<-"Pace of CAP Shock (yrs)"
p_labels$p[4]<-"Fixed Proportion of Rates"
p_labels$p[6]<-"Projection Years"
p_labels$p[7]<-"Elasticity (Curtail)"
p_labels$p[8]<-"Elasticity (Invest)"
p_labels$p[9]<-"Elasticity (Rates)"
p_labels$p[10]<-"Inst. Costs (Curtail)"
p_labels$p[11]<-"Inst. Costs (Invest)"
p_labels$p[12]<-"Inst. Costs (Rates)"
p_labels$p[13]<-"Long-Term Supply Goal"
p_labels$p[15]<-"Max Annual Rate Increse"
p_labels$p[16]<-"Background Demand Decay"
p_labels$p[17]<-"Hard Infrastructure Decay"
p_labels$p[18]<-"Investment Lead Time"
p_labels$p[19]<-"Conservation Sensitivity"
p_labels$p[20]<-"Bond Life"
p_labels$p[21]<-"Goal Debt Service Ratio"
p_labels$p[22]<-"Bond Interest Rate"

#Fix metrics symbols
metrics$Metric[1]<-"Reliab_Total"
metrics$Metric[2]<-"Reliab_Avg"
metrics$Metric[3]<-"Reliab_Min"
metrics$Metric[4]<-"Ending_PC_Demand"
metrics$Metric[5]<-"Ending_PC_Revenue"
metrics$Metric[6]<-"MagCAPShock_1"